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
# setwd("E:/git_repos/uyghur_corpora")
setwd("C:/Users/conno/git_repos/uyghur_corpora")

############################
# LOAD AND PREPROCESS DATA #
############################

# Load data we want to analyze
input_file <- archive_read('corpora/full_data_maxent.zip')
maxent_data <- read_csv(input_file, col_types = cols())

# Create variable corresponding to final vowel
maxent_data <- maxent_data %>%
  mutate(last_one = substr(last_two, 2, 2))

# Aggregate data for simpler analysis
root_agg <- maxent_data %>%
  group_by(
    root, log_norm_count, raised_form_prop, last_two, last_two_distance, last_one, 
    has_name, has_ane, has_che, raising_candidate, raised
  ) %>%
  summarize(n = sum(back_count + front_count),
            percent_back = sum(back_count) / n) %>%
  mutate(raised = ifelse(raised, 'raised_n', 'unraised_n')) %>%
  pivot_wider(names_from = raised, values_from = c(n, percent_back), values_fill = 0) 

###########################
# BUILD AND SAVE TABLEAUX #
###########################

# Function that builds relevant tableaux given a list of constraints
create_tableaux <- function(data, constraints, output_file) {
  headers <- c('Input', 'Output', 'Frequency', constraints, '*Unraised', 'ID')
  write.table(
    matrix(headers, nrow=1), 
    file=output_file, 
    row.names=FALSE, 
    col.names = FALSE, 
    sep=','
  )
  for (i in 1:nrow(data)) {
    row <- data[i,]
    # Check whether we expect root to raise. Two conditions: 
    # - root occurs in environment that conditions raising
    # - root is a root that raises. We can't identify this directly, but we
    #   estimate it based on proportion of raised forms in corpus. Non-raisers
    #   should never show any raising.
    should_raise <- row$raising_candidate & row$raised_form_prop > 0
    
    root_name <- str_c(
      row$root, 
      ifelse(should_raise, 'RAISER', 'NON_RAISER'),
      sep='-'
    )
    raise_f <- c(root_name, paste(root_name, '-RAISED-F', sep=''))
    raise_b <- c(root_name, paste(root_name, '-RAISED-B', sep=''))
    unraise_f <- c(root_name, paste(root_name, '-UNRAISED-F', sep=''))
    unraise_b <- c(root_name, paste(root_name, '-UNRAISED-B', sep=''))
    
    # Add frequencies
    raise_f <- c(raise_f, round(row$n_raised_n * (1 - row$percent_back_raised_n)))
    raise_b <- c(raise_b, round(row$n_raised_n * row$percent_back_raised_n))
    unraise_f <- c(unraise_f, round(row$n_unraised_n * (1 - row$percent_back_unraised_n)))
    unraise_b <- c(unraise_b, round(row$n_unraised_n * row$percent_back_unraised_n)) 
    
    if ('VAgreeSurface' %in% constraints) {
      if (row$last_two == 'BB') {
        raise_f <- c(raise_f, 1)
        raise_b <- c(raise_b, '')
        unraise_f <- c(unraise_f, 1)
        unraise_b <- c(unraise_b, '')
      }
      if (row$last_two == 'FF') {
        raise_f <- c(raise_f, '')
        raise_b <- c(raise_b, 1)
        unraise_f <- c(unraise_f, '')
        unraise_b <- c(unraise_b, 1)
      }
      if (row$last_two == 'FB') {
        raise_f <- c(raise_f, '')
        raise_b <- c(raise_b, 1)
        unraise_f <- c(unraise_f, 1)
        unraise_b <- c(unraise_b, '')
      }
      if (row$last_two == 'BF') {
        raise_f <- c(raise_f, 1)
        raise_b <- c(raise_b, '')
        unraise_f <- c(unraise_f, '')
        unraise_b <- c(unraise_b, 1)
      }
    }
    if ('VAgreeUnderlying' %in% constraints) {
      raise_f <- c(raise_f, ifelse(row$last_two %in% c('FB', 'BB'), 1, ''))
      raise_b <- c(raise_b, ifelse(row$last_two %in% c('BF', 'FF'), 1, ''))
      unraise_f <- c(unraise_f, ifelse(row$last_two %in% c('FB', 'BB'), 1, ''))
      unraise_b <- c(unraise_b, ifelse(row$last_two %in% c('BF', 'FF'), 1, ''))
    }
    if ('HarmonicUniformity' %in% constraints) {
      raise_f <- c(raise_f, row$p_hc_b)
      raise_b <- c(raise_b, 1 - row$p_hc_b)
      unraise_f <- c(unraise_f, row$p_hc_b)
      unraise_b <- c(unraise_b, 1 - row$p_hc_b)
    }
      
    # *Unraised and Max violations
    if (should_raise) {
      # Occurs in a context where raising should occur and root is a raiser
      raise_f <- c(raise_f, '', 1)
      raise_b <- c(raise_b, '', 1)
      unraise_f <- c(unraise_f, 1, '')
      unraise_b <- c(unraise_b, 1, '')
    } else {
      # Occurs in a context where raising shouldn't apply or root is a non-raiser
      raise_f <- c(raise_f, '', 1)
      raise_b <- c(raise_b, '', 1)
      unraise_f <- c(unraise_f, '', '')
      unraise_b <- c(unraise_b, '', '')
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
    write.table(
      matrix(unraise_f, nrow=1), 
      file=output_file, 
      row.names=FALSE, 
      col.names = FALSE, 
      append=TRUE, 
      sep=','
    )
    write.table(
      matrix(unraise_b, nrow=1), 
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
  #5, 4, 3, 2, 1, 0.5, 0.1
  2
)
mus_to_try <- rep(0, length(sigmas_to_try))
k <- 5

# Fit surface-true model
surface_true_output <- read_csv('maxent_data/surface_true_output.csv')
surface_true_cv_m <- cross_validate(
  surface_true_output,
  k = k,
  mu_values = mus_to_try,
  sigma_values = sigmas_to_try
)
write_csv(surface_true_cv_m, 'maxent_data/model_results/surface_true.csv')

# Fit input-oriented model
opaque_output <- read_csv('maxent_data/opaque_output.csv')
opaque_cv_m <- cross_validate(
  opaque_output,
  k = k,
  mu_values = mus_to_try,
  sigma_values = sigmas_to_try
)
write_csv(opaque_cv_m, 'maxent_data/model_results/input.csv')

# Fit input-surface-model
opaque_surface_output <- read_csv('maxent_data/opaque_surface_output.csv')
opaque_surface_cv_m <- cross_validate(
  opaque_surface_output,
  k = k,
  mu_values = mus_to_try,
  sigma_values = sigmas_to_try
)
write_csv(opaque_surface_cv_m, 'maxent_data/model_results/input_surface.csv')

###################################
# BESPOKE K-FOLD CROSS-VALIDATION #
###################################

# Set seed for reproducibility
set.seed(123456789)

do_cross_validation <- function(data, k, constraints, model_name, model_folder, mus, sigmas, create_files=TRUE) {
  row_count <- nrow(data)
  randomized_data <- data[sample(row_count),]
  randomized_data$partition <- rep(seq_len(k), ceiling(row_count / k))[1:row_count]
  
  for (held_out in 1:k) {
    training <- randomized_data %>%
      filter(partition != held_out) %>%
      group_by(
        root, log_norm_count, raised_form_prop, last_two, last_two_distance, last_one, 
        has_name, has_ane, has_che, raising_candidate, raised
      ) %>%
      summarize(n = sum(back_count + front_count),
                percent_back = sum(back_count) / n) %>%
      mutate(raised = ifelse(raised, 'raised_n', 'unraised_n')) %>%
      pivot_wider(names_from = raised, values_from = c(n, percent_back), values_fill = 0) 
    
    test <- randomized_data %>%
      filter(partition == held_out) %>%
      group_by(
        root, log_norm_count, raised_form_prop, last_two, last_two_distance, last_one, 
        has_name, has_ane, has_che, raising_candidate, raised
      ) %>%
      summarize(n = sum(back_count + front_count),
                percent_back = sum(back_count) / n) %>%
      mutate(raised = ifelse(raised, 'raised_n', 'unraised_n')) %>%
      pivot_wider(names_from = raised, values_from = c(n, percent_back), values_fill = 0) 
    
    # Train simple logistic regression model to calculate P(hc|x). We're using
    # proportions weighted by count rather than individual tokens to speed things
    # up, because we need to train this model quite a few times.
    # TODO FIX THIS
    lexical_model <- glm(
      percent_back ~ last_one * log_norm_count * raised_form_prop + has_name + has_ane + has_che,
      data=training,
      family="binomial",
      weights=training$n
    )
    # Apply model to predict data
    training$p_hc_b <- predict(lexical_model, type='response')
    test$p_hc_b <- predict(lexical_model, newdata=test, type='response')
    
    # Create training tableau
    training_file <- str_c(
      model_folder, '/', model_name, '_training_', held_out, '.csv'
    )
    if (create_files) {
      create_tableaux(training, constraints, training_file)
    }
    
    # Create test tableau
    test_file <- str_c(
      model_folder, '/', model_name, '_test_', held_out, '.csv'
    )
    if (create_files) {
      create_tableaux(test, constraints, test_file)
    }
    for (i in 1:length(mus)) {
      training_tableaux <- read_csv(training_file)
      test_tableaux <- read_csv(test_file)
      fitted_m <- optimize_weights(training_tableaux, mu_scalar=mus[i], sigma_scalar=sigmas[i])
      predictions <- predict_probabilities(test_tableaux, fitted_m$weights)
    }
  }
}


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
