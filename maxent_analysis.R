require(tidyverse)
require(archive)
require(devtools)
install_github("connormayer/maxent.ot")
require(maxent.ot)

# Load in our data. If you want to run this script yourself, you'll need
# to set the working directory to wherever you've checked out the corpora
# setwd("E:/git_repos/uyghur_corpora")
setwd("C:/Users/conno/git_repos/uyghur_corpora")

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
# simple_model_me <- glmer(
#   back_count ~ last_one * log_norm_count + last_one * raised_form_prop + has_name + has_ane + has_che + (1|root) + (1|source),
#   data=raising_candidates,
#   family="binomial"
# )
template_model <- glm(
  back_count ~ template * log_norm_count + template * raised_form_prop + has_name + has_ane + has_che,
  data=raising_candidates,
  family="binomial"
)
# template_model_me <- glmer(
#   back_count ~ template * log_norm_count + template * raised_form_prop + has_name + has_ane + has_che + (1|root) + (1|source),
#   data=raising_candidates,
#   family="binomial"
# )
raising_candidates$predicted_simple <- predict(simple_model, type='response')
#raising_candidates$predicted_simple_me <- predict(simple_model_me, type='response')
raising_candidates$predicted_template <- predict(template_model, type='response')
#raising_candidates$predicted_template_me <- predict(template_model_me, type='response')

# Aggregate data for simpler analysis
root_agg <- raising_candidates %>%
  filter(raised) %>%
  group_by(
    root, log_norm_count, raised_form_prop, last_two, last_two_distance, 
    has_name, has_ane, has_che, predicted_simple, #predicted_simple_me, 
    predicted_template, #predicted_template_me
  ) %>%
  summarize(n = sum(back_count + front_count),
            percent_back = sum(back_count) / n)

# Create OTSoft format input for suface-true model
headers <- c('Input', 'Output', 'Frequency', 'VAgreeSurface')
output_file <- 'maxent_data/surface_true_output.csv'

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  raise_f <- c(row$root, paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c(row$root, paste(row$root, '-RAISE-B', sep=''))
  
  # Add frequencies
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))
  
  # Calculate VAgree violations
  if (row$last_two == 'BB') {
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$last_two == 'FF') {
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$last_two == 'FB') {
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$last_two == 'BF') {
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  
  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}

# Create OTSoft format input for opaque model
output_file <- 'maxent_data/opaque_output.csv'
headers <- c('Input', 'Output', 'Frequency', 'VAgreeUnderlying')

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  raise_f <- c(row$root, paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c(row$root, paste(row$root, '-RAISE-B', sep=''))
  
  # Add frequencies
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))
  
  # Calculate VAgreeUnderlying violation
  raise_f <- c(raise_f, ifelse(row$last_two %in% c('FB', 'BB'), 1, ''))
  raise_b <- c(raise_b, ifelse(row$last_two %in% c('BF', 'FF'), 1, ''))

  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}

# Create OTSoft format input for lexical model
output_file <- 'maxent_data/lexical_output.csv'
headers <- c('Input', 'Output', 'Frequency', 'HarmonicUniformity')

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  raise_f <- c(row$root, paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c(row$root, paste(row$root, '-RAISE-B', sep=''))
  
  # Add frequencies
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))
  
  # Calculate HarmonicUniformity violation
  raise_f <- c(raise_f, row$predicted_template)
  raise_b <- c(raise_b, 1 - row$predicted_template)
  
  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}

# Create OTSoft format input for opaque-surface model
output_file <- 'maxent_data/opaque_surface_output.csv'
headers <- c('Input', 'Output', 'Frequency', 'VAgreeSurface', 'VAgreeUnderlying')

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  raise_f <- c(row$root, paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c(row$root, paste(row$root, '-RAISE-B', sep=''))
  
  # Add frequencies
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))
  
  # Calculate VAgreeSurface violations
  if (row$last_two == 'BB') {
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$last_two == 'FF') {
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$last_two == 'FB') {
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$last_two == 'BF') {
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  
  # Calculate VAgreeUnderlying violation
  raise_f <- c(raise_f, ifelse(row$last_two %in% c('FB', 'BB'), 1, ''))
  raise_b <- c(raise_b, ifelse(row$last_two %in% c('BF', 'FF'), 1, ''))
  
  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}

# Create OTSoft format input for lexical-surface model
output_file <- 'maxent_data/lexical_surface_output.csv'
headers <- c('Input', 'Output', 'Frequency', 'VAgreeSurface', 'HarmonicUniformity')

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  raise_f <- c(row$root, paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c(row$root, paste(row$root, '-RAISE-B', sep=''))
  
  # Add frequencies
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))
  
  # Calculate VAgreeSurface violations
  if (row$last_two == 'BB') {
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  if (row$last_two == 'FF') {
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$last_two == 'FB') {
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, 1)
  }
  if (row$last_two == 'BF') {
    raise_f <- c(raise_f, 1)
    raise_b <- c(raise_b, '')
  }
  
  # Calculate HarmonicUniformity violation
  raise_f <- c(raise_f, row$predicted_template)
  raise_b <- c(raise_b, 1 - row$predicted_template)

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
# k is the number of roots because there are two possible outcomes for each
# root (unraised-b, raised-b, unraised-f, raised-f) and therefore 3 probabilities
# to set for each root.
oracle_k <- nrow(root_agg)
oracle_bic <- log(sum(root_agg$n)) * oracle_k - 2 * oracle_ll

# Fit maxent models

surface_true_output <- read_csv('maxent_data/surface_true_output.csv')
surface_true_m <- optimize_weights(surface_true_output)
surface_true_preds <- predict_probabilities(surface_true_output, surface_true_m$weights)

opaque_output <- read_csv('maxent_data/opaque_output.csv')
opaque_m <- optimize_weights(opaque_output)
opaque_preds <- predict_probabilities(opaque_output, opaque_m$weights)

opaque_surface_output <- read_csv('maxent_data/opaque_surface_output.csv')
opaque_surface_m <- optimize_weights(opaque_surface_output)
opaque_surface_preds <- predict_probabilities(opaque_surface_output, opaque_surface_m$weights)

lexical_output <- read_csv('maxent_data/lexical_output.csv')
lexical_m <- optimize_weights(lexical_output)
lexical_preds <- predict_probabilities(lexical_output, lexical_m$weights)

lexical_surface_output <- read_csv('maxent_data/lexical_surface_output.csv')
lexical_surface_m <- optimize_weights(lexical_surface_output)
lexical_surface_preds <- predict_probabilities(lexical_surface_output, lexical_surface_m$weights)

# Compare models
compare_models(surface_true_m, opaque_m, opaque_phon_m, indexed_no_phon_m, indexed_m, method='bic')

# More conservative estimates for lexical models
df <- simple_model$df.null - simple_model$df.residual
lexical_bic <- log(sum(root_agg$n)) * (df + indexed_no_phon_m$k) - 2 * indexed_no_phon_m$loglik
lexical_surface_bic <- log(sum(root_agg$n)) * (df + indexed_m$k)  - 2 * indexed_m$loglik
