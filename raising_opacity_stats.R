require(archive)
require(lme4)
require(tidyverse)

# Load in our data. If you want to run this script yourself, you'll need
# to set the working directory to wherever you've checked out the corpora
setwd("E:/git_repos/uyghur_corpora")
# We store the parses in zip files to save space
rfa_file <- archive_read("corpora/rfa/output/conservative_parses.zip")
rfa_data <- read_csv(rfa_file, col_types = cols())
rfa_data$source <- 'rfa'
awazi_file <- archive_read("corpora/awazi/output/conservative_parses.zip")
awazi_data <- read_csv(awazi_file, col_types = cols())
awazi_data$source <- 'awazi'
akademiye_file <- archive_read("corpora/akademiye/output/conservative_parses.zip")
akademiye_data <- read_csv(akademiye_file, col_types = cols())
akademiye_data$source <- 'akademiye'
full_data <- rbind(rfa_data, awazi_data, akademiye_data)

# Total count of tokens
n = sum(full_data$count)

# Get total count, count normalized per million words, and the proportion
# of times a root occurs in a raised form
counts <- full_data %>% 
  group_by(root) %>%
  summarize(total_count = sum(count),
            norm_count = 1000000 * (total_count / n),
            raised_form_prop = sum(raised) / total_count)

# Glom these counts onto our parses so we can use them as fixed effects in
# our model
full_data <- inner_join(full_data, counts, by="root")

# Get only tokens that have the potential to general opaque harmony. These
# are tokens where:
# - the final two harmonizing vowels in the root clash
# - the final vowel COULD raise (that is, it's either 'e' or 'a' and it occurs
#   in an open syllable
# - the token includes at least one harmonizing suffix
possible_opaque_harmonizers <- full_data %>%
  filter(last_two == 'BF' | last_two == 'FB') %>%
  filter(raising_candidate) %>%
  filter(back_count + front_count > 0)

# We remove tokensthat have a number of suffixes that block harmony and impose
# their own, or a number of suffixes that have strange harmonic behavior.
# 'che' and 'ane' routinely harmonize transparently
possible_opaque_harmonizers <- possible_opaque_harmonizers %>% 
  mutate(
    has_blocker = grepl("gpr_rsub3|ger_past2|prog|px2pl", possible_opaque_harmonizers$tags),
    #has_ye = grepl('ye$', possible_opaque_harmonizers$root),
    has_ane = grepl("ane$", possible_opaque_harmonizers$root),
    has_che = grepl("che$", possible_opaque_harmonizers$root),
    has_name = grepl('name$', possible_opaque_harmonizers$root),
  )

# Roots that end in 'che' that aren't the '-che' suffix
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'anche',]$has_che <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'bunche',]$has_che <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'qanche',]$has_che <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'birqanche',]$has_che <- FALSE

# roots that end in ane, but not the -ane suffix
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'bahane',]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'epsane',]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == "i'ane",]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'jerimane',]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'perwane',]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'perghane',]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'péshane',]$has_ane <- FALSE
possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'zamane',]$has_ane <- FALSE

# possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'dékoratsiye',]$has_ye <- FALSE
# possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'diplomatiye',]$has_ye <- FALSE
# possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'énsiklopédiye',]$has_ye <- FALSE
# possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'énsiklopédiye',]$has_ye <- FALSE
# possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'énsiklopédiye',]$has_ye <- FALSE
# possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'énsiklopédiye',]$has_ye <- FALSE
# possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'énsiklopédiye',]$has_ye <- FALSE
# possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'énsiklopédiye',]$has_ye <- FALSE
# possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'nahiye',]$has_ye <- FALSE
# possible_opaque_harmonizers[possible_opaque_harmonizers$root == 'pozitsiye',]$has_ye <- FALSE

# Filter out a few bogus roots
possible_opaque_harmonizers <- possible_opaque_harmonizers %>%
  filter(
     # Not a real word
     root != 'qure' & 
     # Correct Uyghur word is komitét, but seems to often be written
     # as kométit
     root != 'kométet' &
     # Correct root is teqqasla
     root != 'teqqasle' &
     # Correct root is etiwarla
     root != 'etiwarle' &
     # Correct root is mejburla
     root != 'mejburle' &
     # Eli is the name 'Ali', not a raised form of 'ela'
     root != 'ela' &
     # Alem is a real word, but these all look to be cases of alim 'scholar'
     root != 'alem' &
     # Nur'ela is a name, but these stories refer to a man named Nur'eli
     root != "nur'ela" &
     # These correspond to esme 'memory', not esma
     root != 'esma' &
     # Correct root is tuzla
     root != 'tuzle' &
     # Not a real root
     root != 'bade' &
     # This is an interesting case: 'alahide' means 'special', but you can
     # find cases where this is produced as 'alahida'. The parser treats this
     # as 'alahe' + '-da', but I don't think this is the right decomposition
     # of this word. Something interesting with harmony here, but I don't think
     # this can be thought of as suffix alternation, so I exclude it here
     root != 'alahe' &
     # These tokens actually correspond to esir 'century'
     root != 'esra' &
     # This is actually a token of 'resim' spelled incorrectly
     root != 'ressam' &
     # Not a real root
     root != 'peza' &
     # Actual root is xadim
     root != 'xadime' &
     # This is Sewde, one of the wives of Mohammad, not sewda 'passion'
     root != 'sewda' &
     # This is a fixed expression sayisida 'because of'
     root != 'saye' &
     # Not a real root
     root != 'paye' &
     # This is a valid root, but could also correspond to the city Astana
     root != 'astane'
  )

tokens <- possible_opaque_harmonizers %>% 
  filter(!has_blocker) %>%
  select(-has_blocker)

# Does opacity vary as a result of overall frequency
harmonic_raisers <- tokens %>%
  filter(raised) %>%
  mutate(
    opaque = case_when(back_count > 0 & last_two == 'BF' ~ 1,
                       front_count > 0 & last_two == 'FB' ~ 1,
                       TRUE ~ 0),
    log_total_count = log(total_count),
    log_norm_count = log(norm_count)
  )

harmonic_raisers <- raiser_in

root_agg <- harmonic_raisers %>% 
  group_by(root, last_two) %>%
  summarise(n=sum(count),
            percent_back = sum(back_count) / n,
            total_count = mean(total_count),
            raised_form_prop = mean(raised_form_prop))

# Test for significance 
full_model <- glmer(
  opaque ~ log_norm_count + raised_form_prop + last_two + has_name + has_ane + has_che + (1|root) + (1|source), 
  data=harmonic_raisers, 
  family="binomial",
  # The bobyqa optimizer has more luck converging than the default Nelder_Mead
  control=glmerControl(optimizer = 'bobyqa')
)

# Test for significance of source random effects. These both end up being
# significant, so we keep them in the model
no_source <-  glmer(opaque ~ log_norm_count + raised_form_prop + last_two + (1|root), data=harmonic_raisers, family="binomial")
anova(full_model, no_source, test='Chisq')

no_word <- glmer(opaque ~ log_norm_count + raised_form_prop + last_two + (1|source), data=harmonic_raisers, family="binomial")
anova(full_model, no_word, test='Chisq')

# Test for significance of fixed effects
no_last_two <- glmer(opaque ~ log_norm_count + raised_form_prop + (1|root) + (1|source), data=harmonic_raisers, family="binomial")
anova(full_model, no_last_two, test='Chisq')

no_raised_form <- glmer(opaque ~ log_norm_count + last_two + (1|root) + (1|source), data=harmonic_raisers, family="binomial")
anova(full_model, no_raised_form, test='Chisq')

no_freq <- glmer(opaque ~ raised_form_prop + last_two + (1|root) + (1|source), data=harmonic_raisers, family="binomial")
anova(full_model, no_freq, test='Chisq')


