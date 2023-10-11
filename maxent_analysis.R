library(tidyverse)
library(archive)
# library(devtools)
# install_github("connormayer/maxent.ot")
library(maxent.ot)
library(lme4)

# THINGS TO DO:
# - switch logistic regression model to apply to proportions by using weights
# - set up k-fold validation for lexical models

# Load in our data. If you want to run this script yourself, you'll need
# to set the working directory to wherever you've checked out the corpora
setwd("E:/git_repos/uyghur_corpora")
# setwd("C:/Users/conno/git_repos/uyghur_corpora")

############################
# LOAD AND PREPROCESS DATA #
############################

# Load data we want to analyze
input_file <- archive_read('corpora/raising_candidates.zip')
raising_candidates <- read_csv(input_file, col_types = cols())

# Create variable corresponding to final vowel
raising_candidates <- raising_candidates %>%
  mutate(last_one = substr(last_two, 2, 2))

# Train simple logistic regression model to calculate P(hc|x)
lexical_model <- glm(
  back_count ~ last_one * log_norm_count * raised_form_prop + has_name + has_ane + has_che,
  data=raising_candidates,
  family="binomial"
)
# Get number of fitted parameters in lexical model for use in BIC later
num_coefs <- length(coef(lexical_model))

# lexical_model_me <- glmer(
#   back_count ~ last_one * log_norm_count * raised_form_prop + has_name + has_ane + has_che + (1|root),
#   data=raising_candidates,
#   family="binomial"
# )

raising_candidates$p_hc_b <- predict(lexical_model, type='response')
# raising_candidates$p_hc_b_me <- predict(simple_model_me, type='response')

# Aggregate data for simpler analysis
root_agg <- raising_candidates %>%
  filter(raised) %>%
  group_by(
    root, log_norm_count, raised_form_prop, last_two, last_two_distance, 
    has_name, has_ane, has_che, p_hc_b
  ) %>%
  summarize(n = sum(back_count + front_count),
            percent_back = sum(back_count) / n)

###########################
# BUILD AND SAVE TABLEAUX #
###########################

# Function that builds relevant tableaux given a list of constraints
create_tableaux <- function(data, constraints, output_file) {
  headers <- c('Input', 'Output', 'Frequency', constraints)
  write.table(
    matrix(headers, nrow=1), 
    file=output_file, 
    row.names=FALSE, 
    col.names = FALSE, 
    sep=','
  )
  for (i in 1:nrow(data)) {
    row <- data[i,]
    raise_f <- c(row$root, paste(row$root, '-RAISE-F', sep=''))
    raise_b <- c(row$root, paste(row$root, '-RAISE-B', sep=''))
    
    # Add frequencies
    raise_f <- c(raise_f, (row$n * (1 - row$percent_back)))
    raise_b <- c(raise_b, (row$n * row$percent_back))
    
    if ('VAgreeSurface' %in% constraints) {
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
    }
    
    if ('VAgreeUnderlying' %in% constraints) {
      raise_f <- c(raise_f, ifelse(row$last_two %in% c('FB', 'BB'), 1, ''))
      raise_b <- c(raise_b, ifelse(row$last_two %in% c('BF', 'FF'), 1, ''))
    }
    
    if ('HarmonicUniformity' %in% constraints) {
      raise_f <- c(raise_f, row$p_hc_b)
      raise_b <- c(raise_b, 1 - row$p_hc_b)
    }
    
    write.table(
      matrix(raise_f, nrow=1), 
      file=output_file, 
      row.names=FALSE, 
      col.names = FALSE, 
      append=TRUE, 
      sep=','
    )
    write.table(
      matrix(raise_b, nrow=1), 
      file=output_file, 
      row.names=FALSE, 
      col.names = FALSE, 
      append=TRUE, 
      sep=','
    )
  }
}

# Create tableaux for surface-true model
create_tableaux(
  root_agg, 
  c('VAgreeSurface'), 
  'maxent_data/surface_true_output.csv'
)

# Create tableaux for input-sensitive model
create_tableaux(
  root_agg, 
  c('VAgreeUnderlying'), 
  'maxent_data/opaque_output.csv'
)

# Create tableaux for lexical model
create_tableaux(
  root_agg, 
  c('HarmonicUniformity'), 
  'maxent_data/lexical_output.csv'
)

# Create tableaux for input-surface model
create_tableaux(
  root_agg, 
  c('VAgreeSurface', 'VAgreeUnderlying'), 
  'maxent_data/opaque_surface_output.csv'
)

# Create tableaux for lexical-surface model
create_tableaux(
  root_agg, 
  c('VAgreeSurface', 'HarmonicUniformity'), 
  'maxent_data/lexical_surface_output.csv'
)

# Create tableaux for lexical-opaque-surface model
create_tableaux(
  root_agg, 
  c('VAgreeSurface', 'VAgreeUnderlying', 'HarmonicUniformity'), 
  'maxent_data/lexical_surface_opaque_output.csv'
)

##########################
# FIT AND COMPARE MODELS #
##########################

sigmas_to_try <- c(
  5, 4, 3, 2, 1, 0.5, 0.1
)
mus_to_try <- rep(0, length(sigmas_to_try))

# Fit surface-true model
surface_true_output <- read_csv('maxent_data/surface_true_output.csv')
surface_true_cv_m <- cross_validate(
  surface_true_output,
  k = 10,
  mu_values = mus_to_try,
  sigma_values = sigmas_to_try
)

# Fit input-oriented model
opaque_output <- read_csv('maxent_data/opaque_output.csv')
opaque_cv_m <- cross_validate(
  opaque_output,
  k = 10,
  mu_values = mus_to_try,
  sigma_values = sigmas_to_try
)

# Fit input-surface-model
opaque_surface_output <- read_csv('maxent_data/opaque_surface_output.csv')
opaque_surface_cv_m <- cross_validate(
  opaque_surface_output,
  k = 10,
  mu_values = mus_to_try,
  sigma_values = sigmas_to_try
)

# Fit lexical model
lexical_output <- read_csv('maxent_data/lexical_output.csv')
lexical_cv_m <- cross_validate(
  lexical_output,
  k = 10,
  mu_values = mus_to_try,
  sigma_values = sigmas_to_try
)

# Fit lexical-surface model
lexical_surface_output <- read_csv('maxent_data/lexical_surface_output.csv')
lexical_surface_cv_m <- cross_validate(
  lexical_surface_output,
  k = 10,
  mu_values = mus_to_try,
  sigma_values = sigmas_to_try
)

# Fit lexical-input-surface model
lexical_surface_opaque_output <- read_csv(
  'maxent_data/lexical_surface_opaque_output.csv'
)
lexical_surface_opaque_cv_m <- cross_validate(
  lexical_surface_opaque_output,
  k = 10,
  mu_values = mus_to_try,
  sigma_values = sigmas_to_try
)

surface_true_cv_m
opaque_cv_m
opaque_surface_cv_m
lexical_cv_m
lexical_surface_cv_m
lexical_surface_opaque_cv_m
