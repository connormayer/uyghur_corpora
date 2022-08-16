import apertium
import argparse
import hfst
import io
import pandas as pd
import os
import re
import time

from zipfile import ZipFile

# A bunch of default arguments
METADATA_FILENAME = 'metadata.csv'
DOCS_ZIP = 'pa_documents.zip'
LAT_DOCS_ZIP = 'lat_documents.zip'
DISCARD_FILE = 'discarded.csv'
MULTIPARSE_FILE = 'multiparse_parses.csv'
CONSERVATIVE_FILE = 'conservative_parses.csv'
OUTPUT_DIR = 'output'
FST_DIR = "/mnt/e/apertium-uig"
ORTHO_FST_PATH = "dev/ortho/ara-lat.hfst"
LATIN_FST_PATH = "dev/ortho/lat-ara-test.hfst"
MAIN_FST_PATH = "uig.automorf.hfst"
LAST_PROCESSED = 'last_processed.txt'
PRINT_EVERY = 100

# Some useful constants for text processing
VOWELS = ['a', 'e', 'i', 'o', 'u', 'ö', 'ü', 'é']
TRANSPARENT_VOWELS = ['i', 'é']
BACK_VOWELS = ['a', 'o', 'u']
FRONT_VOWELS = [ 'ö', 'ü', 'e']
CONSONANTS_ONE = [
	'b', 'p', 't', 'j', 'x', 'd', 'r', 'z', 's', 'f', 'q', 'k', 'g', 'l', 'm',
	'n', 'h', 'w', 'y'
]
CONSONANTS_TWO = ['ng', 'sh', 'zh', 'ch']
FRONT_CONSONANTS = ['k', 'g']
BACK_CONSONANTS = ['q', 'gh']

UNIGRAM_CONS_RE = "b|p|t|d|j|x|r|z|s|f|q|k|g|l|m|n|h|w|y"
BIGRAM_CONS_RE = "ng|sh|gh|zh|ch"
VOWELS_RE = "a|e|i|o|u|é|ü|ö"

# Column names for output file
COLNAMES = [
	'root', 'word', 'tags', 'author', 'count', 'back_count', 'front_count', 'raising_candidate', 
	'raised', 'detailed_template', 'template', 'fine_template', 'last_two', 'last_two_distance',
	'root_suffix_distance'
]

def clean_document(doc):
	"""
	Removes punctuation and other characters from corpus file
	"""
	doc = re.sub(r"[,،\/#!$%\^&\*;:{}=_`~()\"<>–\?“—”«»\.]", "", doc)
	doc = re.sub(r"[\u200f\xa0\r\n]", "", doc)
	return doc

def load_transducer(fst_path):
	"""
	Loads the orthographic fst we use to convert Perso-Arabic to Latin
	""" 
	istr = hfst.HfstInputStream(fst_path)
	transducer = istr.read_all()[0]
	transducer.lookup_optimize()
	return transducer

def get_backness_counts(tags):
	"""
	Gets counts of front/back suffix forms
	"""
	back_count = 0
	front_count = 0

	if '-b' in tags:
		back_count += 1
	if '-f' in tags:
		front_count += 1

	return back_count, front_count

def get_templates(root):
	"""
	Gets schematic templates for each root that can be used to group them.

	The detailed template classifies segments as:
		- B: back vowel
		- F: front vowel
		- N: neutral vowel
		- Q: back consonant
		- K: front consonant
		- C: neutral consonant

	The template excludes neutral consonants and the fine template
	excludes neutral vowels.
	"""
	# Make detailed template
	detailed_template = root
	for char in BACK_VOWELS:
		detailed_template = detailed_template.replace(char, 'B')
	for char in FRONT_VOWELS:
		detailed_template = detailed_template.replace(char, 'F')
	for char in BACK_CONSONANTS:
		detailed_template = detailed_template.replace(char, 'Q')
	for char in CONSONANTS_TWO:
		detailed_template = detailed_template.replace(char, 'C')
	for char in FRONT_CONSONANTS:
		detailed_template = detailed_template.replace(char, 'K')
	for char in TRANSPARENT_VOWELS:
		detailed_template = detailed_template.replace(char, 'N') 
	for char in CONSONANTS_ONE:
		detailed_template = detailed_template.replace(char, 'C')
	detailed_template = detailed_template.replace("'", "")
	detailed_template = detailed_template.replace('-', '')

	template = detailed_template.replace('C', '')

	fine_template = template.replace('N', '')

	if not fine_template:
		fine_template = 'N'

	return detailed_template, template, fine_template

def is_raising_candidate(root, word):
	"""
	Checks whether the context is sufficient to allow raising to occur,
	regardless of whether it actually does or not.
	"""
	vowels = list(re.finditer(VOWELS_RE, root))

	if not vowels:
		return False

	last_vowel = vowels[-1]

	if not last_vowel.group() in ['a', 'e']:
		return False

	root_len = len(root)

	if len(word) <= root_len:
		return False

	next_char = word[root_len]
	if root_len == last_vowel.end():
		# Possible raiser is root-final
		next_two_chars = word[root_len: root_len + 2]
		if re.search(BIGRAM_CONS_RE, next_two_chars) and len(word) > root_len + 2:
			# First segment in suffix is a bigram consonant
			following = word[root_len + 2]
			return following in VOWELS
		elif len(word) > root_len + 1:
			following = word[root_len + 1]
			return following in VOWELS
	else:
		# Possible raiser is not word-final
		remaining_root = root[last_vowel.end():]
		if root_len - last_vowel.end() == 1:
			# Root ends in C corresponding to one character
			return next_char in VOWELS 
		elif root_len - last_vowel.end() == 2 and re.search(BIGRAM_CONS_RE, remaining_root):
			# Root ends in C corresponding to two characters
			return next_char in VOWELS

	return False

def convert_latin(word, ortho_transducer):
	"""
	Converts a word from some orthography to Latin
	"""
	try:
		return ortho_transducer.lookup(word)[-1][0]
	except:
		return []

def is_raised(root, word):
	"""
	Checks whether a root occurs in its raised form in a word.
	"""
	vowels = list(re.finditer(VOWELS_RE, root))

	if not vowels:
		return False

	last_vowel = vowels[-1]

	if not last_vowel.group() in ['a', 'e']:
		return False

	i_version = root[:last_vowel.end() - 1] + 'i' + root[last_vowel.end():]
	e_version = root[:last_vowel.end() - 1] + 'é' + root[last_vowel.end():]

	return i_version in word or e_version in word

def get_last_two_harmonizer_distance(detailed_template):
	"""
	Get the distance in segments between the final two harmonizing vowels
	Returns -1 if there's only a single harmonizer
	"""
	harmonizer_idxs = list(re.finditer('F|B', detailed_template))
	if len(harmonizer_idxs) < 2:
		return -1

	return harmonizer_idxs[-1].end() - harmonizer_idxs[-2].end()

def get_root_harm_tag_distance(tags):
	"""
	Gets the distance in tags from the root to the first
	harmonizing tag.

	Returns -1 if there's no harmonizing tag
	"""
	tag_list = tags.split('_')
	# Some roots correspond to multiple tags, we don't want
	# to count these
	if len(tag_list) > 1:
		if tag_list[0] == 'n' and tag_list[1] != 'لىق':
			start = 1
		elif tag_list[0] == 'adj' and tag_list[1] != 'subst':
			start = 1
		elif tag_list[0] == 'adv':
			start = 1
		else:
			start = 2
	else:
		return -1

	try:
		end = list(('-b' in tag or '-f' in tag) for tag in tag_list).index(True)
	except:
		return -1

	return end - start

def write_parses(word, readings, f, ortho_transducer, author):
	"""
	Calculates features and writes parses to file
	"""
	num_parses = len(readings)

	latin_word = convert_latin(word, ortho_transducer)

	if not latin_word:
		return

	for parse in readings:
		parse = parse[0]
		parts = list(filter(None, re.split('<|>', parse)))
		baseform, tags = parts[0], parts[1:]

		# work in Latin
		latin_root = convert_latin(baseform, ortho_transducer)
		if not latin_root:
			return

		tags = '_'.join(tags)
		count = 1 / num_parses

		# Get variables of interest about root/word
		raising_candidate = is_raising_candidate(latin_root, latin_word)
		raised = is_raised(latin_root, latin_word)
		detailed_template, template, fine_template = get_templates(latin_root)
		back_count, front_count = get_backness_counts(tags)
		back_count /= num_parses
		front_count /= num_parses
		last_two = fine_template[-2:]
		last_two_harmonizer_distance = get_last_two_harmonizer_distance(detailed_template)
		intervening_suffix_count = get_root_harm_tag_distance(tags)

		vals = map(str, [
			latin_root, latin_word, tags, author, count, back_count, front_count, 
			raising_candidate, raised, detailed_template, template, fine_template, last_two,
			last_two_harmonizer_distance, intervening_suffix_count
		])

		f.write(','.join(vals) + '\n')

def initialize_files(corpus_dir):
	"""
	Initializes output files by clearing them and writing headers
	"""
	output_path = os.path.join(corpus_dir, OUTPUT_DIR)
	if not os.path.isdir(output_path):
		os.mkdir(output_path)

	open(os.path.join(output_path, DISCARD_FILE), 'w').close()

	with open(os.path.join(output_path, MULTIPARSE_FILE), 'w') as f:
		f.write(','.join(COLNAMES) + '\n')

	with open(os.path.join(output_path, CONSERVATIVE_FILE), 'w') as f:
		f.write(','.join(COLNAMES) + '\n')

def remove_raised_root(readings, ortho_transducer):
	"""
	apertium-uig sometimes lists roots in their raised form as well as their
	unraised form. This function checks whether a list of parses includes raised
	versions, and removes them if so.

	This function erroneously excludes all tokens of 'bir' because the transducer
	erroneously suggests 'ber' as a possible root. This is an error, but it's 
	ok for our purposes here.
	"""
	if len(readings) == 1:
		return readings

	baseforms = [convert_latin(list(filter(None, re.split('<|>', r[0])))[0], ortho_transducer) for r in readings]

	new_readings = []

	for i, reading in enumerate(readings):
		latin_base = baseforms[i]
		if latin_base:
			vowels = list(re.finditer(VOWELS_RE, latin_base))

			if not vowels:
				new_readings.append(reading)
			else:
				last_vowel = vowels[-1]
				last_v_index = last_vowel.end()

				if last_vowel.group() in TRANSPARENT_VOWELS:
					a_version = latin_base[:last_v_index - 1] + 'a' + latin_base[last_v_index:]
					e_version = latin_base[:last_v_index - 1] + 'e' + latin_base[last_v_index:]

					if not (a_version in baseforms or e_version in baseforms):
						new_readings.append(reading)
				else:
					new_readings.append(reading)

	return new_readings

def get_latin_parse(lat_transducer, x):
	"""
	Parses a word using input Latin orthography
	"""
	parse = lat_transducer.lookup(x)
	if parse:
		return parse[0][0]
	else:
		return ''

def parse_corpus(corpus_dir, latin_input, resume, fst_dir, print_every):
	"""
	Top-level logic for parinsg the corpus
	"""
	# This lets us stop part way through and resume if necessary
	if not resume:
		initialize_files(corpus_dir)
	else:
		found_start = False
		with open(os.path.join(corpus_dir, OUTPUT_DIR, LAST_PROCESSED), 'r') as f:
			last_processed = f.read()

	# Load metadata and various transducers
	metadata = pd.read_csv(os.path.join(corpus_dir, METADATA_FILENAME))
	parser = load_transducer(os.path.join(fst_dir, MAIN_FST_PATH))

	if latin_input:
		lat_transducer = load_transducer(
			os.path.join(fst_dir, LATIN_FST_PATH)
		)

	ortho_transducer = load_transducer(
		os.path.join(fst_dir, ORTHO_FST_PATH)
	)

	# Keep track of parse counts
	total_words = 0
	failed_words = 0
	multiparse_words = 0
	single_parse_words = 0

	# Get corpus file
	if not latin_input:
		zip_file = os.path.join(corpus_dir, DOCS_ZIP)
	else:
		zip_file = os.path.join(corpus_dir, LAT_DOCS_ZIP)

	with open(os.path.join(corpus_dir, OUTPUT_DIR, DISCARD_FILE), 'a') as discard_f, \
		 open(os.path.join(corpus_dir, OUTPUT_DIR, MULTIPARSE_FILE), 'a') as mp_f, \
		 open(os.path.join(corpus_dir, OUTPUT_DIR, CONSERVATIVE_FILE), 'a') as conservative_f, \
		 ZipFile(zip_file) as corpus_zip:

		# Iterate through 
		for index, row in metadata.iterrows():
			if resume and not found_start:
				if row['filename'] == last_processed:
					found_start = True
				else:
					print("{} already processed, skipping".format(row['filename']))
					continue

			if index % print_every == 0:
				print("Processing article {} of {}: {}".format(
					index, len(metadata), row['filename'])
				)

			with open(os.path.join(corpus_dir, OUTPUT_DIR, LAST_PROCESSED), 'w') as last_processed_f:
				last_processed_f.write(row['filename'])

			author_latin = 'None_{}'.format(os.path.split(corpus_dir)[-1])

			if not latin_input:
				try:
					author_latin = ' '.join([ortho_transducer.lookup(part)[1][0] for part in row.author.split(' ')])
				except:
					pass
			else:
				if not pd.isna(row.author):
					author_latin = row.author

			# Load data for a single article
			doc_file = corpus_zip.open(os.path.split(zip_file)[1].split('.')[0] + '/' + row['filename'])
			document = io.TextIOWrapper(doc_file, 'utf-8').read()
			
			clean_doc = clean_document(document)
			if latin_input:
				clean_doc = ' '.join(map(lambda x: get_latin_parse(lat_transducer, x), clean_doc.split(' ')))

			word_parses = {x: parser.lookup(x) for x in clean_doc.split(' ')}

			for word, readings in word_parses.items():
				total_words += 1
				if not readings:
					# Can't be parsed
					failed_words += 1
					discard_f.write(word + '\n')
					continue

				readings = remove_raised_root(readings, ortho_transducer)

				if len(readings) > 1:
					# Word can be parsed as more than one root. This is tricky because the transducer
					# frequently lists inflected forms as roots. E.g., tallashqan can be parsed
					# as 'talla + sh + qan' or 'tallash + qan'. We make a maximally conservative 
					# output by discarding EVERY parse with more than one root, and a less
					# conservative one by keeping them but weighting them by their probability.
					write_parses(word, readings, mp_f, ortho_transducer, author_latin)
					multiparse_words += 1
				else:
					# One parse
					write_parses(word, readings, conservative_f, ortho_transducer, author_latin)
					single_parse_words += 1

	print("Total words: {}".format(total_words))
	print("Failed words: {}".format(failed_words))
	print("Multiparse words: {}".format(multiparse_words))
	print("Single parse words: {}".format(single_parse_words))

	with ZipFile(os.path.join(corpus_dir, OUTPUT_DIR, DISCARD_FILE.split('.')[0] + '.zip'), 'w') as my_zip:
		my_zip.write(os.path.join(corpus_dir, OUTPUT_DIR, DISCARD_FILE), arcname=DISCARD_FILE)

	os.remove(os.path.join(corpus_dir, OUTPUT_DIR, DISCARD_FILE))

	with ZipFile(os.path.join(corpus_dir, OUTPUT_DIR, CONSERVATIVE_FILE.split('.')[0] + '.zip'), 'w') as my_zip:
		my_zip.write(os.path.join(corpus_dir, OUTPUT_DIR, CONSERVATIVE_FILE), arcname=CONSERVATIVE_FILE)

	os.remove(os.path.join(corpus_dir, OUTPUT_DIR, CONSERVATIVE_FILE))

	with ZipFile(os.path.join(
			corpus_dir, OUTPUT_DIR, MULTIPARSE_FILE.split('.')[0] + '.zip'), 'w') as my_zip:
		my_zip.write(os.path.join(corpus_dir, OUTPUT_DIR, MULTIPARSE_FILE), arcname=MULTIPARSE_FILE)

	os.remove(os.path.join(corpus_dir, OUTPUT_DIR, MULTIPARSE_FILE))

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="Parse corpus"
	)
	parser.add_argument(
		'corpus_dir', type=str, help='The corpus directory to search.'
	)
	parser.add_argument(
		'--latin_input', action='store_true'
	)
	parser.add_argument(
		'--resume', action='store_true'
	)
	parser.add_argument(
		'--fst_dir', type=str, default=FST_DIR,
		help='The path to the root directory of the FST'
	)
	parser.add_argument(
		'--print_every', type=int, default=PRINT_EVERY,
		help="Print status every n files."
	)
	args = parser.parse_args()
	results = parse_corpus(
		args.corpus_dir, args.latin_input, args.resume, args.fst_dir,
		args.print_every
	)
