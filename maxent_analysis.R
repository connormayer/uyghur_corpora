require(tidyverse)
require(archive)

# Load in our data. If you want to run this script yourself, you'll need
# to set the working directory to wherever you've checked out the corpora
# setwd("E:/git_repos/uyghur_corpora")
setwd("//wsl.localhost/Ubuntu/home/connor/uyghur_corpora")

input_file <- archive_read('corpora/raising_candidates.zip')
raising_candidates <- read_csv(input_file, col_types = cols())

# Train a logistic regression model to predict backness class
# full_model <- glmer(
#     back_count ~ log_norm_count + raised_form_prop + last_two + last_two_distance + has_name + has_ane + has_che + (1|root),
#     data=raising_candidates,
#     family="binomial"#,
#     # The bobyqa optimizer has more luck converging than the default Nelder_Mead
#     #control=glmerControl(optimizer = 'bobyqa')
#   )

simple_model <- glm(
  back_count ~ last_two * log_norm_count + last_two * raised_form_prop + last_two_distance + has_name + has_ane + has_che,
  data=raising_candidates,
  family="binomial"#,
  # The bobyqa optimizer has more luck converging than the default Nelder_Mead
  #control=glmerControl(optimizer = 'bobyqa')
)

raising_candidates$predicted_simple <- predict(simple_model)
#raising_candidates$predicted_full <- predict(full_model)

root_agg <- raising_candidates %>%
  filter(raised) %>%
  group_by(
    root, log_norm_count, raised_form_prop, last_two, last_two_distance, 
    has_name, has_ane, has_che, predicted_simple #, predicted_full
  ) %>%
  summarize(n = sum(back_count + front_count),
            percent_back = sum(back_count) / n)

headers <- c('', '', '', 'VAgree', '*Unraised')
n_constraints = 2

# Fully surface-true output
output_file <- 'maxent_data/surface_true_output.csv'

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

# Fully opaque output
output_file <- 'maxent_data/opaque_output.csv'
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


# Indexed output
output_file <- 'maxent_data/indexed_output.csv'
headers <- c('', '', '', 'VAgree', '*Unraised', 'HarmonizeBack', 'HarmonizeFront')
n_constraints = 4

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
output_file <- 'maxent_data/fully_indexed_output.csv'
n_general_constraints <- 2
n_total_constraints <- n_general_constraints + nrow(root_agg)
headers <- c('', '', '', 'VAgree', '*Unraised')

for (i in 1:nrow(root_agg)) {
  headers <- c(
    headers, 
    paste('HarmonizeBack-', raising_candidates[i,]$root, sep=''),
    paste('HarmonizeFront-', raising_candidates[i,]$root, sep=''))
}
write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')
write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)
#write.table(matrix(c('', '', '', rep(0, n_total_constraints)), nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)

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
  
  # Add empty spaces for previous constraints
  n_prev_constraints <- 2* (i - 1)
  prev_vals <- rep('', n_prev_constraints)
  no_raise_f <- c(no_raise_f, prev_vals)
  no_raise_b <- c(no_raise_b, prev_vals)
  raise_f <- c(raise_f, prev_vals)
  raise_b <- c(raise_b, prev_vals)
  
  # Calculate HarmonizeBack violation
  no_raise_f <- c(no_raise_f, 1)
  no_raise_b <- c(no_raise_b, '')
  raise_f <- c(raise_f, 1)
  raise_b <- c(raise_b, '')
  
  
  # Calculate HarmonizeFront violation
  no_raise_f <- c(no_raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_f <- c(raise_f, 1)
  raise_b <- c(raise_b, 1)
  
  # Add empty spaces for following constraints
  n_following_constraints <- 2 * (n_total_constraints - i)
  following_vals <- rep('', n_following_constraints)
  no_raise_f <- c(no_raise_f, following_vals)
  no_raise_b <- c(no_raise_b, following_vals)
  raise_f <- c(raise_f, following_vals)
  raise_b <- c(raise_b, following_vals)
  
  write.table(matrix(no_raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}


# Fit maxent models
# These take a while to run
surface_true_m <- optimize_weights('maxent_data/surface_true_output.csv', in_sep=',', upper_bound=50)
opaque_m <- optimize_weights('maxent_data/opaque_output.csv', in_sep=',', upper_bound = 50)
indexed_m <- optimize_weights('maxent_data/indexed_output.csv', in_sep=',', upper_bound=50)

compare_models(surface_true_m, opaque_m, indexed_m, method='bic')
oracle_m <- optimize_weights('maxent_data/fully_indexed_output.csv', in_sep=',', upper_bound=50)
