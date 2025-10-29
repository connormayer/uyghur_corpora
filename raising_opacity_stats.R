require(archive)
require(forcats)
require(lme4)
require(tidyverse)
library(brms)
library(bayesplot)
library(DHARMa)
library(viridis)

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
  summarise(total_count = sum(count),
            norm_count = 1000000 * (total_count / n),
            raised_form_prop = sum(raised) / total_count)

# Glom these counts onto our parses so we can use them as fixed effects in
# our model
full_data <- inner_join(full_data, counts, by="root")

# Get only tokens that have the potential to general opaque harmony. These
# are tokens where:
# - the rightmost two harmonizing vowels in the root clash
# - the rightmost vowel COULD raise (that is, it's either 'e' or 'a' and it 
#   occurs in a word-final open syllable
# - the token includes at least one harmonizing suffix

# We flag tokens that have suffixes that block harmony and impose
# their own, and suffixes that have strange harmonic behavior.
# 'che' and 'ane' routinely harmonize transparently
full_data <- full_data %>% 
  mutate(
    has_blocker = grepl("gpr_rsub3|ger_past2|prog|px2pl|loc_subst|iver", full_data$tags),
    has_ane = grepl("ane$", full_data$root),
    has_che = grepl("che$", full_data$root),
    has_name = grepl('name$', full_data$root),
  )

# Roots that end in 'che' that aren't the '-che' suffix
full_data[full_data$root == 'anche',]$has_che <- FALSE
full_data[full_data$root == 'bunche',]$has_che <- FALSE
full_data[full_data$root == 'qanche',]$has_che <- FALSE
full_data[full_data$root == 'birqanche',]$has_che <- FALSE

# roots that end in ane, but not the -ane suffix
full_data[full_data$root == 'bahane',]$has_ane <- FALSE
full_data[full_data$root == 'epsane',]$has_ane <- FALSE
full_data[full_data$root == "i'ane",]$has_ane <- FALSE
full_data[full_data$root == 'jerimane',]$has_ane <- FALSE
full_data[full_data$root == 'perwane',]$has_ane <- FALSE
full_data[full_data$root == 'perghane',]$has_ane <- FALSE
full_data[full_data$root == 'péshane',]$has_ane <- FALSE
full_data[full_data$root == 'zamane',]$has_ane <- FALSE

# Filter out a few bogus roots
full_data <- full_data %>%
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
     root != 'astane' &
     # Should be bolghuchi
     root != 'bolghuche'&
     # muncha is a real root, but munchilik is mistakenly parsed as muncha + liQ
     root != 'muncha' &
     # Should be 'ipalide'
     root != 'ipalida'&
     # This is a typo for üst-i-ge, not underlying usta
     word != 'ustige' &
     # Should be bataréye
     root != 'bataréya'
  )

# Remove words with harmony blockers like '-wat'
full_data <- full_data %>% 
  filter(!has_blocker) %>%
  select(-has_blocker)

# Remove parses with both front and back suffixes. These are often typos
# or transducer errors
full_data <- full_data %>% 
  filter(!(front_count > 0 & back_count > 0))

# For maxent analysis, keep only BF/FB/BB/FF root tokens with harmonizing
# suffixes
maxent_data <- full_data %>%
  filter(last_two %in% c('BF', 'FB', 'BB', 'FF')) %>%
  filter(back_count + front_count > 0) %>%
  mutate(
    log_total_count = log(total_count),
    log_norm_count = log(norm_count)
  )

# write.csv(maxent_data, 'corpora/full_data_maxent.csv')

# For analysis below, keep only tokens that could undergo raising
raising_candidates <- maxent_data %>%
  filter(raising_candidate) 

raisers <- raising_candidates %>%
  filter(raised)

raisers %>% 
  filter(last_two == 'BB') %>% 
  group_by(root) %>% 
  count()

raisers %>% 
  filter(last_two == 'FF') %>% 
  group_by(root) %>% 
  count()

# Find BF and FB forms that undergo raising
opaque_raisers <- raisers %>%
  filter(last_two %in% c('BF', 'FB')) %>%
  mutate(
    opaque = case_when(back_count > 0 & last_two == 'FB' ~ 1,
                       front_count > 0 & last_two == 'BF' ~ 1,
                       TRUE ~ 0)
  )

# Get number of possible raising roots
opaque_raisers %>% group_by(root) %>% count()

# Get number of roots that vacillate
opaque_raisers %>% 
  group_by(root) %>% 
  summarise(front_sum = sum(front_count),
            back_sum = sum(back_count)) %>%
  filter(front_sum > 0 & back_sum > 0)

# Test for significance

# Chi squared tests for whether raised FF/BF and BB/FB roots behave
# differently.
chi_data_b <- full_data %>% 
  filter(last_two %in% c('BB', 'FB')) %>%
  filter(back_count + front_count > 0) %>%
  filter(raised)

chi_data_f <- full_data %>% 
  filter(last_two %in% c('FF', 'BF')) %>%
  filter(back_count + front_count > 0) %>%
  filter(raised)

chisq.test(table(chi_data_b$last_two, chi_data_b$back_count))
chisq.test(table(chi_data_f$last_two, chi_data_f$back_count))

# # Code for a Bayesian analysis
# full_model_brm <- brm(
#   opaque ~ log_norm_count + raised_form_prop + last_two + last_two_distance + root_suffix_distance + has_name + has_ane + has_che + (1|root) + (1|source/author),
#   data=opaque_raisers,
#   family="bernoulli",
#   cores = 4,
#   iter = 6000,
#   control = list(adapt_delta=0.999)
# )
# # Save model so we don't have to train it again :p
# saveRDS(full_model_brm, 'brms_model.rds')

full_model_brm <- readRDS('brms_model.rds')

# Run residual diagnostics using DHARMa
model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(full_model_brm)),
  observedResponse = opaque_raisers$opaque,
  fittedPredictedResponse = apply(t(posterior_epred(full_model_brm)), 1, mean),
  integerResponse = TRUE
)

windowsFonts(Times = windowsFont("Times New Roman"))

plot(model.check)

ggplot(data.frame(residuals = sort(model.check$scaledResiduals)), aes(sample = residuals)) +
  stat_qq(distribution = qunif) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "QQ Plot of Uniform Residuals",
       x = "Expected",
       y = "Observed") +
  theme_minimal()

ggplot(data.frame(predicted = model.check$fittedPredictedResponse, residuals = model.check$scaledResiduals),
       aes(x = predicted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(title = "Residuals vs. Predicted Values",
       x = "Predicted Values",
       y = "Scaled Residuals") +
  theme_minimal()

ggplot(data.frame(predicted = model.check$fittedPredictedResponse, residuals = model.check$scaledResiduals), aes(x = predicted, y = residuals)) +
  geom_point(alpha = 0.5) + # Add points with some transparency
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") + # Expected mean
  geom_smooth(method = "loess", se = FALSE, color = "blue") + # Add a smoother
  labs(
    title = "DHARMa Scaled Residuals vs. Predicted Values",
    x = "Predicted Values",
    y = "Scaled Residuals (0-1)"
  ) +
  theme_minimal()


testOutliers(model.check)
testDispersion(model.check)
testUniformity(model.check)
testZeroInflation(model.check)

# Plot posterior samples
posterior <- as.matrix(full_model_brm)
mcmc_areas(
   posterior,
   pars = c("b_Intercept",
            "b_log_norm_count", 
            "b_raised_form_prop",
            "b_last_twoFB",
            "b_last_two_distance",
            "b_root_suffix_distance",
            "b_has_nameTRUE",
            "b_has_aneTRUE",
            "b_has_cheTRUE"),
   area_method = "equal height",
   prob=0.95,
   point_est = "mean"
  ) +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        legend.title = element_blank(),
        legend.text=element_text(size=25),
        text = element_text(family = "Times")) +
  scale_y_discrete(
    labels = c(
       "b_Intercept" = "Intercept",
       "b_log_norm_count" = "Log frequency", 
       "b_raised_form_prop" = "Proportion raised",
       "b_last_twoFB" = "Last two (ref. BF)",
       "b_last_two_distance" = "Last two distance",
       "b_root_suffix_distance" = "Root-suffix distance",
       "b_has_nameTRUE" = "Has -name",
       "b_has_aneTRUE" = "Has -ane",
       "b_has_cheTRUE" = "Has -che"
    ) 
  ) 
ggsave("figures/bayesian_plot.pdf", height=8, width=12, units="in")

# GRAPHS
graph_data <- raisers %>%
  group_by(last_two, back_count) %>%
  count()

# Break down opacity rates by template
ggplot(data=graph_data, aes(x=fct_relevel(last_two, c('FF', 'BF', 'FB', 'BB')), y=n, fill=factor(back_count))) +
  geom_col(position = 'fill') +
  geom_text(size=9, aes(label=paste("n = ",n)), 
            position=position_fill(vjust=0.5),
            family="Times") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        legend.title = element_blank(),
        legend.text=element_text(size=25),
        text = element_text(family = "Times")) +
  xlab("Underlying root template") +
  ylab("Proportion of tokens")  +
  scale_fill_viridis_d(labels=c("Front suffix", "Back suffix"), begin=0.6)
ggsave("figures/full_harmonic_raisers.pdf", height=8, width=12, units="in")

suffix_agg <- opaque_raisers %>%
  mutate(suffix = case_when(has_ane ~ '-ane',
                            has_che ~ '-che',
                            has_name ~ '-name',
                            TRUE ~ 'other BF')) %>%
  filter(last_two == 'BF') %>%
  group_by(last_two, back_count, suffix) %>%
  count()

ggplot(suffix_agg, aes(x=fct_relevel(suffix, levels=c('-che', '-ane', '-name', 'other BF')), y=n, fill=factor(back_count))) +
  geom_col(position = 'fill') +
  geom_text(size=9, aes(label=paste("n = ",n)), 
            position=position_fill(vjust=0.5),
            family="Times") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        legend.title = element_blank(),
        legend.text=element_text(size=25),
        text = element_text(family = "Times")) +
  xlab("Suffix identity") +
  ylab("Proportion of tokens")  +
  scale_fill_viridis_d(labels=c("Front suffix", "Back suffix"), begin=0.6)
ggsave("figures/suffix_raising.pdf", height=10, width=12, units="in")

# Histograms
# Aggregate rate of opacity by root
root_agg <- opaque_raisers %>%
  group_by(root, last_two) %>%
  summarize(
    opaque_rate = mean(opaque),
    log_norm_count = mean(log_norm_count),
    raised_form_prop=mean(raised_form_prop),
    back_rate = mean(back_count))

root_agg %>%
  ggplot() +
  geom_histogram(aes(back_rate), binwidth=0.025) +
  xlab("Proportion of back responses") +
  ylab("Root type count") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        legend.title = element_blank(),
        legend.text=element_text(size=25),
        text = element_text(family = "Times"),
        strip.text.x = element_text(size = 25)) +
  facet_wrap(~ last_two, scales="free_y")

ggsave('figures/root_histogram.pdf', height=6, width=12, units="in")


