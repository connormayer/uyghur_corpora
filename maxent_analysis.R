library(tidyverse)
library(archive)
if (!require(maxent.ot)) {
  library(devtools)
  install_github("connormayer/maxent.ot")
  library(maxent.ot)
}
library(lme4)
library(doParallel)
library(foreach)

options(dplyr.summarise.inform = FALSE)

# Set up parallel environment
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type="PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

# If you want to run this script yourself, you'll need
# to set the working directory to wherever you've checked out the corpora
setwd("E:/git_repos/uyghur_corpora")

############################
# LOAD AND PREPROCESS DATA #
############################

# Load data we want to analyze
input_file <- archive_read('corpora/full_data_maxent.zip')
maxent_data <- read_csv(input_file, col_types = cols())

# Create variable corresponding to final vowel
maxent_data <- maxent_data |>
  mutate(last_one = substr(last_two, 2, 2))

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
    
    root_name <- stringr::str_c(
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

##########################
# FIT AND COMPARE MODELS #
##########################

# Can't use maxent.ot cross-validation because we need to fit a logistic
# regression model to each training set. Also, the data sets are so large
# here that this takes ages unless it's parallelized, which this function does.
create_partitions <- function(data, k, model_details) {
  # Define grid of mu/sigmas x held-out partition
  param_grid <- expand.grid(1:length(model_details), 1:k)
  colnames(param_grid) <- c("model_idx", "held_out")
  
  foreach(row_idx=1:nrow(param_grid), .export = 'create_tableaux') %dopar% {
    model_row <- model_details[[param_grid[row_idx,]$model_idx]]
    model_name <- model_row[1]
    constraints <- model_row[2:length(model_row)]
    model_folder <- stringr::str_c('maxent_data/', model_name)
    
    held_out <- param_grid[row_idx,]$held_out
    
    training <- data |>
      dplyr::filter(partition != held_out) |>
      dplyr::group_by(
        root, log_norm_count, raised_form_prop, last_two, last_two_distance, 
        last_one, has_name, has_ane, has_che, raising_candidate, raised
      ) |>
      dplyr::summarize(n = sum(back_count + front_count),
                       percent_back = sum(back_count) / n) |>
      dplyr::mutate(raised = ifelse(raised, 'raised_n', 'unraised_n')) |>
      tidyr::pivot_wider(
        names_from = raised, values_from = c(n, percent_back), values_fill = 0
      )
    
    test <- data |>
      dplyr::filter(partition == held_out) |>
      dplyr::group_by(
        root, log_norm_count, raised_form_prop, last_two, last_two_distance, 
        last_one, has_name, has_ane, has_che, raising_candidate, raised
      ) |>
      dplyr::summarize(n = sum(back_count + front_count),
                       percent_back = sum(back_count) / n) |>
      dplyr::mutate(raised = ifelse(raised, 'raised_n', 'unraised_n')) |>
      tidyr::pivot_wider(
        names_from = raised, values_from = c(n, percent_back), values_fill = 0
      )
    
    if ("HarmonicUniformity" %in% constraints) {
      # Need a separate training data set for logistic regression that combines
      # raised and unraised variants
      training_for_regression <- data |>
        dplyr::filter(partition != held_out) |>
        dplyr::group_by(
          root, log_norm_count, raised_form_prop, last_two, last_two_distance, last_one,
          has_name, has_ane, has_che
        ) |>
        dplyr::summarize(n = sum(back_count + front_count),
                         percent_back = sum(back_count) / n)
      
      # Train a mixed-effects logistic regression model to calculate P(hc=B|x). 
      # We're using proportions weighted by count rather than individual tokens 
      # to speed things up, because we need to train this model quite a few times.
      lexical_model <- lme4::glmer(
        percent_back ~ last_one * log_norm_count * raised_form_prop + has_name + has_ane + has_che + (1|root),
        data=training_for_regression,
        family="binomial",
        weights=training_for_regression$n,
        control=lme4::glmerControl(optimizer = 'bobyqa')
      )
      # Apply model to predict data
      training$p_hc_b <- predict(
        lexical_model, newdata=training, type='response')
      test$p_hc_b <- predict(
        lexical_model, newdata=test, type='response', allow.new.levels = TRUE
      )
    }
    
    # Create training tableau
    training_file <- stringr::str_c(
      model_folder, '/', model_name, '_training_', held_out, '.csv'
    )
    create_tableaux(training, constraints, training_file)
    
    # Create test tableau
    test_file <- stringr::str_c(
      model_folder, '/', model_name, '_test_', held_out, '.csv'
    )
    create_tableaux(test, constraints, test_file)
    return(NA)
  }
  return(NA)
}

do_cross_validation <- function(k, model_name, model_folder, mus, sigmas) {
  # Define grid of mu/sigmas x held-out partition
  param_grid <- expand.grid(1:length(mus), 1:k)
  colnames(param_grid) <- c("param_idx", "held_out")
  
  # Fit models in parallel for efficiency
  result <- foreach(row_idx=1:nrow(param_grid), .combine=rbind,
                    .export = 'create_tableaux') %dopar% {
    param_row <- param_grid[row_idx,]
    mu <- mus[param_row$param_idx]
    sigma <- sigmas[param_row$param_idx]
    held_out <- param_row$held_out
  
    print(stringr::str_c(
      "Running fold ", held_out, " on mu = ", mu, " and sigma = ", sigma)
    )

    training_file <- stringr::str_c(
      model_folder, '/', model_name, '_training_', held_out, '.csv'
    )

    test_file <- stringr::str_c(
      model_folder, '/', model_name, '_test_', held_out, '.csv'
    )

    training_tableaux <- readr::read_csv(training_file, show_col_types=FALSE)
    test_tableaux <- readr::read_csv(test_file, show_col_types=FALSE)
    print("Training model...")
    fitted_m <- maxent.ot::optimize_weights(
      training_tableaux, mu_scalar=mu, sigma_scalar=sigma
    )
    print("Testing on held-out fold")
    predictions <- maxent.ot::predict_probabilities(test_tableaux, fitted_m$weights)
    return(data.frame(model_name = model_name, mu=mu, sigma=sigma, folds=k,
                      held_out=held_out, ll_train=fitted_m$loglik,
                      ll_test=predictions$loglik))
  }
  
  agg_result <- result |>
    group_by(model_name, mu, sigma, folds) |>
    summarize(mean_ll_train=mean(ll_train),
              mean_ll_test=mean(ll_test),
              sd_ll_train=sd(ll_train),
              sd_ll_test=sd(ll_test))
  output_file_name <- str_c('maxent_data/model_results/', model_name, '.csv', sep='')
  write_csv(agg_result, output_file_name)
  return(agg_result)
}

sigmas_to_try <- c(
  1, 10, 20, 50, 100, 200, 500
)
mus_to_try <- rep(0, length(sigmas_to_try))
k <- 10

# Set this to TRUE and change the seed to produce new partitions
CREATE_FILES <- TRUE

if (CREATE_FILES) {
  # Set seed for reproducibility
  set.seed(495869402)
  
  # Define partitions, we'll use the same ones for each model
  row_count <- nrow(maxent_data)
  randomized_data <- maxent_data[sample(row_count),]
  randomized_data$partition <- rep(seq_len(k), ceiling(row_count / k))[1:row_count]
  
  # Create training and test tableaux for each partition of each model
  # This means we use the same partitions for each model
  # We do this in advance to speed up run-time of cross-validation, which lets
  # use easily try a larger range of sigmas
  model_details <- list(
    c("surface", "VAgreeSurface"),
    c("input", "VAgreeUnderlying"),
    c("input_surface", "VAgreeSurface", "VAgreeUnderlying"),
    c("lexical", "HarmonicUniformity"),
    c("lexical_surface", "VAgreeSurface", "HarmonicUniformity"),
    c("lexical_input_surface", "VAgreeSurface", "VAgreeUnderlying",
     "HarmonicUniformity")
  )
  create_partitions(randomized_data, k, model_details)
}

# Fit surface-true model
# surface_m <- do_cross_validation(
#   k, 'surface', 'maxent_data/surface', mus_to_try, sigmas_to_try
# )
# print(surface_m)

# Fit input-oriented model
# input_m <- do_cross_validation(
#   k, 'input', 'maxent_data/input', mus_to_try, sigmas_to_try
# )
# print(input_m)

# Fit input-surface-model
# input_surface_m <- do_cross_validation(
#   k, 'input_surface', 'maxent_data/input_surface', mus_to_try, sigmas_to_try
# )
# print(input_surface_m)

# Fit lexical model
# lexical_m <- do_cross_validation(
#   k, 'lexical', 'maxent_data/lexical', mus_to_try, sigmas_to_try
# )
# print(lexical_m)

# Fit lexical model
lexical_m_me <- do_cross_validation(
  k, 'lexical_me', 'maxent_data/lexical_me', mus_to_try, sigmas_to_try
)
print(lexical_m_me)

# # Fit lexical-surface model
# lexical_surface_m <- do_cross_validation(
#   k, 'lexical_surface', 'maxent_data/lexical_surface', mus_to_try, sigmas_to_try
# )
# print(lexical_surface_m)

# Fit lexical-surface model
lexical_surface_m_me <- do_cross_validation(
  k, 'lexical_surface_me', 'maxent_data/lexical_surface_me', mus_to_try, sigmas_to_try
)
print(lexical_surface_m_me)

# Fit lexical-input-surface model
# lexical_input_surface_m <- do_cross_validation(
#   k, 'lexical_input_surface', 'maxent_data/lexical_input_surface', mus_to_try,
#   sigmas_to_try
# )
# print(lexical_input_surface_m)

lexical_input_surface_m_me <- do_cross_validation(
  k, 'lexical_input_surface_me', 'maxent_data/lexical_input_surface_me', mus_to_try,
  sigmas_to_try
)
print(lexical_input_surface_m_me)

parallel::stopCluster(cl = my.cluster)

# Fit model to full data set to get parameters and generalization cases
root_agg <- maxent_data |>
  group_by(
    root, log_norm_count, raised_form_prop, last_two, last_two_distance, last_one,
    has_name, has_ane, has_che
  ) |>
  summarize(n = sum(back_count + front_count),
                   percent_back = sum(back_count) / n)

lexical_model <- glmer(
  percent_back ~ last_one * log_norm_count * raised_form_prop + has_name + has_ane + has_che + (1|root),
  data=root_agg,
  family="binomial",
  weights=root_agg$n,
  control=glmerControl(optimizer = 'bobyqa')
)

# Generalize to B-final root seen once
wug_words <- tribble(
  ~root, ~last_one, ~log_norm_count, ~raised_form_prop, ~has_name, ~has_ane, ~has_che,
  "front_no_raised", "B", 0, 0, FALSE, FALSE, FALSE
)

predict(lexical_model, wug_words, type="response", allow.new.levels = TRUE)
