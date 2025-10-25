# uyghur_corpora

This contains the code for "A large-scale corpus study of phonological opacity in Uyghur" by Connor Mayer, which will be published in *Phonology*.

* The `corpora` folder contains each of the corpora. The `XX_documents.zip` file in each directory contains the raw text of every article, and the `metadata.csv` file contains the listing of articles and corresponding metadata. The `output` folder in each directory contains the output of the `apertium` transducer run on each corpus. The `raising_candidates.zip` file is a subset of the parses in each corpus that contains all the BB, FB, BF, and BB root tokens that can potentially undergo raising and include at least one harmonizing suffix. The `full_data_maxent`.zip` file is all the BB, FB, BF, and BB root tokens with harmonizing suffixes, regardless of whether they undergo raising. The corpora and parses are stored in zip files for space and efficiency reasons. The scripts that operate on this data automatically zip/unzip them.
  
* The `figures` folder contains several figures used in the paper.
  
* The `maxent_data` folder contains the tableaux used in the MaxEnt analysis, as well as the results of the cross-validation procedure.

* The `parse_corpus.py` script runs the transducer on a corpus, processes the parses as described in the paper, and does some simple text processing to identify different variables of interest. If you want to run this yourself, you need to check out the modified parser from [[https://anonymous.4open.science/r/apertium-uig-EF36](https://anonymous.4open.science/r/apertium-uig-EF36)](https://anonymous.4open.science/r/apertium-uig-C8EB), compile it by running `make` in the root directory, and then change the paths to the transducer in `parse_corpus.py`. The apertium transducer has a number of dependencies that must be installed in order to compile it: see [https://wiki.apertium.org/wiki/Installation](https://wiki.apertium.org/wiki/Installation) for more detail. 

* The `opaque_raising_stats.R` script does some additional data filtering and text processing, carries out the statistical analysis presented in the paper, and generates the figures. You should be able to run this by simply downloading the required R libraries.

* The `maxent_analysis.R` script converts corpus data into tableaux format and then fits MaxEnt models.

The webscrapers used to produce the corpora can be found below:
* [RFA scraper](https://anonymous.4open.science/r/RFA-Scraper-5841)
* [Uyghur Awazi scraper](https://anonymous.4open.science/r/uyghur_tools-E388)
* [Uyghur akadémiyisi scraper](https://anonymous.4open.science/status/UyghurAcademyWebsiteSpider-UG-70D7)
