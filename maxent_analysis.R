require(tidyverse)

# Load in our data. If you want to run this script yourself, you'll need
# to set the working directory to wherever you've checked out the corpora
setwd("E:/git_repos/uyghur_corpora")

input_file <- archive_read('corpora/raising_candidates.zip')
raising_candidates <- read_csv(input_file, col_types = cols())

# Train a logistic regression model to predict backness class
full_model <- glmer(
    back_count ~ log_norm_count + raised_form_prop + last_two + last_two_distance + has_name + has_ane + has_che + (1|root),
    data=raising_candidates,
    family="binomial",
    # The bobyqa optimizer has more luck converging than the default Nelder_Mead
    control=glmerControl(optimizer = 'bobyqa')
  )

root_agg <- raising_candidates %>%
  group_by(
    root, log_norm_count, raised_form_prop, last_two, last_two_distance, has_name, has_ane, has_che
  ) %>%
  summarize(n = sum(back_count + front_count),
            percent_back = sum(back_count) / n)

headers <- c('', '', '', 'VAgree', '*Unraised', 'HarmonizeBack', 'HarmonizeFront')
n_constraints = 4

# Fully surface-true output
output_file <- 'maxent_data/surface_true_output.csv'

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')
write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)
write.table(matrix(c('', '', '', rep(0, n_constraints)), nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  no_raise_f <- c(row$root, paste(row$root, '-F', sep=''))
  no_raise_b <- c(row$root, paste(row$root, '-B', sep=''))
  raise_f <- c(row$root, paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c(row$root, paste(row$root, '-RAISE-B', sep=''))
  
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
  no_raise_f <- c(no_raise_f, '')
  no_raise_b <- c(no_raise_b, '')
  raise_f <- c(raise_f, '')
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeFront violation
  no_raise_f <- c(no_raise_f, '')
  no_raise_b <- c(no_raise_b, '')
  raise_f <- c(raise_f, '')
  raise_b <- c(raise_b, '')
  
  write.table(matrix(no_raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}

# Fully opaque output
output_file <- 'maxent_data/opaque_output.csv'

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')
write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)
write.table(matrix(c('', '', '', rep(0, n_constraints)), nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  no_raise_f <- c(row$root, paste(row$root, '-F', sep=''))
  no_raise_b <- c(row$root, paste(row$root, '-B', sep=''))
  raise_f <- c(row$root, paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c(row$root, paste(row$root, '-RAISE-B', sep=''))
  
  # Add frequencies
  no_raise_f <- c(no_raise_f, 0)
  no_raise_b <- c(no_raise_b, 0)
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))
  
  # Calculate VAgree violations
  if (row$last_two == 'BB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$last_two == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$last_two == 'FB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$last_two == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  
  # Do raising violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeBack violation
  no_raise_f <- c(no_raise_f, 1)
  no_raise_b <- c(no_raise_b, '')
  raise_f <- c(raise_f, 1)
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeFront violation
  no_raise_f <- c(no_raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_f <- c(raise_f, '')
  raise_b <- c(raise_b, 1)
  
  write.table(matrix(no_raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}


# Indexed output
output_file <- 'maxent_data/indexed_output.csv'

write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',')
write.table(matrix(headers, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)
write.table(matrix(c('', '', '', rep(0, n_constraints)), nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, sep=',', append=TRUE)

for (i in 1:nrow(root_agg)) {
  row <- root_agg[i,]
  no_raise_f <- c(row$root, paste(row$root, '-F', sep=''))
  no_raise_b <- c(row$root, paste(row$root, '-B', sep=''))
  raise_f <- c(row$root, paste(row$root, '-RAISE-F', sep=''))
  raise_b <- c(row$root, paste(row$root, '-RAISE-B', sep=''))
  
  # Add frequencies
  no_raise_f <- c(no_raise_f, 0)
  no_raise_b <- c(no_raise_b, 0)
  raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
  raise_b <- c(raise_b, (row$n * row$percent_back))
  
  # Calculate VAgree violations
  if (row$last_two == 'BB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$last_two == 'FF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$last_two == 'FB') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  if (row$last_two == 'BF') {
    no_raise_f <- c(no_raise_f, '')
    no_raise_b <- c(no_raise_b, '')
    raise_f <- c(raise_f, '')
    raise_b <- c(raise_b, '')
  }
  
  # Do raising violations
  no_raise_f <- c(no_raise_f, 1)
  raise_f <- c(raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeBack violation
  no_raise_f <- c(no_raise_f, 1)
  no_raise_b <- c(no_raise_b, '')
  raise_f <- c(raise_f, 1)
  raise_b <- c(raise_b, '')
  
  # Calculate HarmonizeFront violation
  no_raise_f <- c(no_raise_f, '')
  no_raise_b <- c(no_raise_b, 1)
  raise_f <- c(raise_f, '')
  raise_b <- c(raise_b, 1)
  
  write.table(matrix(no_raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(no_raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_f, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
  write.table(matrix(raise_b, nrow=1), file=output_file, row.names=FALSE, col.names = FALSE, append=TRUE, sep=',')
}
