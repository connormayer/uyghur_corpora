require(tidyverse)
require(archive)

# Load in our data. If you want to run this script yourself, you'll need
# to set the working directory to wherever you've checked out the corpora
setwd("E:/git_repos/uyghur_corpora")

# Load data we want to analyze
input_file <- archive_read('corpora/raising_candidates.zip')
raising_candidates <- read_csv(input_file, col_types = cols())

# Create variable corresponding to final vowel
raising_candidates <- raising_candidates %>%
  mutate(last_one = substr(last_two, 2, 2))

# Train simple logistic regression model to calculate P(hc|x)
simple_model <- glm(
  back_count ~ last_one * log_norm_count + last_one * raised_form_prop + has_name + has_ane + has_che,
  data=raising_candidates,
  family="binomial"
)
raising_candidates$predicted_simple <- predict(simple_model)

# Aggregate data for simpler analysis
root_agg <- raising_candidates %>%
  filter(raised) %>%
  group_by(
    root, log_norm_count, raised_form_prop, last_two, last_two_distance, 
    has_name, has_ane, has_che, predicted_simple #, predicted_full
  ) %>%
  summarize(n = sum(back_count + front_count),
            percent_back = sum(back_count) / n)


# Create OTSoft format input for suface-true model
headers <- c('', '', '', 'VAgree', '*Unraised')
n_constraints = 2
output_file <- 'maxent_data/surface_true_output.csv'

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')
write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  no_raise_f <- c(row$root, paste(row$root, '-F', sep=''))
  no_raise_b <- c('', paste(row$root, '-B', sep=''))
  raise_f <- c('', paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c('', paste(row$root, '-RAISE-B', sep=''))
  
  # Add frequencies
  no_raise_f <- c(no_raise_f, 0)
  no_raise_b <- c(no_raise_b, 0)
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))
  
  # Calculate VAgree violations
  if (row$last_two == 'BB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$last_two == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$last_two == 'FB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$last_two == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  
  # Do raising violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_b <- c(raise_b, '')
  
  write.table(matrix(no_raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}

# Create OTSoft format input for opaque model
output_file <- 'maxent_data/opaque_output.csv'
headers <- c('', '', '', '*Unraised', 'HarmonizeBack', 'HarmonizeFront')
n_constraints = 3

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')
write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  no_raise_f <- c(row$root, paste(row$root, '-F', sep=''))
  no_raise_b <- c('', paste(row$root, '-B', sep=''))
  raise_f <- c('', paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c('', paste(row$root, '-RAISE-B', sep=''))
  
  # Add frequencies
  no_raise_f <- c(no_raise_f, 0)
  no_raise_b <- c(no_raise_b, 0)
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))
  
  # Do raising violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeBack violation
  no_raise_f <- c(no_raise_f, ifelse(row$last_two %in% c('FB', 'BB'), 1, ''))
  no_raise_b <- c(no_raise_b, '')
  raise_f <- c(raise_f, ifelse(row$last_two %in% c('FB', 'BB'), 1, ''))
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeFront violation
  no_raise_f <- c(no_raise_f, '')
  no_raise_b <- c(no_raise_b, ifelse(row$last_two %in% c('BF', 'FF'), 1, ''))
  raise_f <- c(raise_f, '')
  raise_b <- c(raise_b, ifelse(row$last_two %in% c('BF', 'FF'), 1, ''))
  
  write.table(matrix(no_raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}

# Create OTSoft format input for lexical model
output_file <- 'maxent_data/indexed_no_phon_output.csv'
headers <- c('', '', '', '*Unraised', 'HarmonizeBack', 'HarmonizeFront')
n_constraints = 3

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')
write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)
#write.table(matrix(c('', '', '', rep(0, n_constraints)), nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  no_raise_f <- c(row$root, paste(row$root, '-F', sep=''))
  no_raise_b <- c('', paste(row$root, '-B', sep=''))
  raise_f <- c('', paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c('', paste(row$root, '-RAISE-B', sep=''))
  
  # Add frequencies
  no_raise_f <- c(no_raise_f, 0)
  no_raise_b <- c(no_raise_b, 0)
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))

  # Do raising violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeBack violation
  inv_logit <- 1 / (1 + exp(-row$predicted_simple))
  h_b <- inv_logit
  
  no_raise_f <- c(no_raise_f, h_b)
  no_raise_b <- c(no_raise_b, '')
  raise_f <- c(raise_f, h_b)
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeFront violation
  h_f <- 1 - inv_logit
  no_raise_f <- c(no_raise_f, '')
  no_raise_b <- c(no_raise_b, h_f)
  raise_f <- c(raise_f, '')
  raise_b <- c(raise_b, h_f)
  
  write.table(matrix(no_raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}

# Create OTSoft format input for opaque-surface model
output_file <- 'maxent_data/opaque_phon_output.csv'
headers <- c('', '', '', 'VAgree', '*Unraised', 'HarmonizeBack', 'HarmonizeFront')
n_constraints = 4

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')
write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  no_raise_f <- c(row$root, paste(row$root, '-F', sep=''))
  no_raise_b <- c('', paste(row$root, '-B', sep=''))
  raise_f <- c('', paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c('', paste(row$root, '-RAISE-B', sep=''))
  
  # Add frequencies
  no_raise_f <- c(no_raise_f, 0)
  no_raise_b <- c(no_raise_b, 0)
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))
  
  # Calculate VAgree violations
  if (row$last_two == 'BB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$last_two == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$last_two == 'FB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$last_two == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  
  # Do raising violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeBack violation
  no_raise_f <- c(no_raise_f, ifelse(row$last_two %in% c('FB', 'BB'), 1, ''))
  no_raise_b <- c(no_raise_b, '')
  raise_f <- c(raise_f, ifelse(row$last_two %in% c('FB', 'BB'), 1, ''))
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeFront violation
  no_raise_f <- c(no_raise_f, '')
  no_raise_b <- c(no_raise_b, ifelse(row$last_two %in% c('BF', 'FF'), 1, ''))
  raise_f <- c(raise_f, '')
  raise_b <- c(raise_b, ifelse(row$last_two %in% c('BF', 'FF'), 1, ''))
  
  write.table(matrix(no_raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}

# Create OTSoft format input for lexical-surface model
output_file <- 'maxent_data/indexed_output.csv'
headers <- c('', '', '', 'VAgree', '*Unraised', 'HarmonizeBack', 'HarmonizeFront')
n_constraints = 4

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')
write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  no_raise_f <- c(row$root, paste(row$root, '-F', sep=''))
  no_raise_b <- c('', paste(row$root, '-B', sep=''))
  raise_f <- c('', paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c('', paste(row$root, '-RAISE-B', sep=''))
  
  # Add frequencies
  no_raise_f <- c(no_raise_f, 0)
  no_raise_b <- c(no_raise_b, 0)
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))
  
  # Calculate VAgree violations
  if (row$last_two == 'BB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$last_two == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$last_two == 'FB') {
    no_raise_f <- c(no_raise_f, 1)
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$last_two == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, 1)
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  
  # Do raising violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeBack violation
  inv_logit <- 1 / (1 + exp(-row$predicted_simple))
  h_b <- inv_logit

  no_raise_f <- c(no_raise_f, h_b)
  no_raise_b <- c(no_raise_b, '')
  raise_f <- c(raise_f, h_b)
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeFront violation
  #h_f <- 1/-log(1-inv_logit)
  
  h_f <- 1 - inv_logit
  no_raise_f <- c(no_raise_f, '')
  no_raise_b <- c(no_raise_b, h_f)
  raise_f <- c(raise_f, '')
  raise_b <- c(raise_b, h_f)
  
  write.table(matrix(no_raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}

# Oracle output
oracle <- root_agg %>% 
  select(n, percent_back) %>%
  mutate(back_count = n * percent_back,
         front_count = n - back_count,
         log_prob = case_when(back_count == 0 ~ front_count * log(1 - percent_back),
                              front_count == 0 ~ back_count * log(percent_back),
                              TRUE ~ back_count * log(percent_back) + front_count * log(1 - percent_back)))

oracle_ll <- sum(oracle$log_prob)
# k is the number of roots * 3 because there are four possible outcomes for each
# root (unraised-b, raised-b, unraised-f, raised-f) and therefore 3 probabilities
# to set for each root.
oracle_k <- nrow(root_agg) * 3
oracle_bic <- log(sum(root_agg$n)) * oracle_k - 2 * oracle_ll

# Fit maxent models
surface_true_m <- optimize_weights('maxent_data/surface_true_output.csv', in_sep=',', upper_bound=50, mu_scalar=0, sigma_scalar=10)
opaque_m <- optimize_weights('maxent_data/opaque_output.csv', in_sep=',', upper_bound = 50)
opaque_phon_m <- optimize_weights('maxent_data/opaque_phon_output.csv', in_sep=',', upper_bound = 50)
indexed_no_phon_m <- optimize_weights('maxent_data/indexed_no_phon_output.csv', in_sep=',', upper_bound = 50)
indexed_m <- optimize_weights('maxent_data/indexed_output.csv', in_sep=',', upper_bound=50)

# Compare models
compare_models(surface_true_m, opaque_m, opaque_phon_m, indexed_no_phon_m, indexed_m, method='bic')

# More conservative estimates for lexical models
df <- simple_model$df.null - simple_model$df.residual
lexical_bic <- log(sum(root_agg$n)) * (df + indexed_no_phon_m$k) - 2 * indexed_no_phon_m$loglik
lexical_surface_bic <- log(sum(root_agg$n)) * (df + indexed_m$k)  - 2 * indexed_m$loglik
