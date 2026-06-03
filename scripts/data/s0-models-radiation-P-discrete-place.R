# ---
#' @description
#' Fit discrete radiation models and basic diagnostics for one location.
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse, mixtools
#'
#' @arguments
#'   - param[1] (place): reference location, e.g. "Bologna".
#'   - param[2] (save_data): write generated data ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/data/s0-models-radiation-P-discrete-place.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/data/s0-models-radiation-P-discrete-place.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'
#' @outputs
#'   - `data/models/radiation/P/{place}/discrete_models.RData` (discrete radiation models)
#'   - `data/models/radiation/P/{place}/discrete_models_tests.RData` (diagnostic tests)
#'
#' @depends
#'   - none
#'
#' @tags
#'   - data
#'   - main
#'   - appendix
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
library(mixtools)
print_script_info(file.path("scripts", "data", "s0-models-radiation-P-discrete-place.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Bologna"
# Save the output
save_data <- FALSE
# ***************************** Fixed Arguments ********************************
# Reference name
common_name <- c("discrete_models", "discrete_models_tests")
# Years of training
nyears <- outputs$nyears
# ******************************************************************************
#                             Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly=TRUE)
if (!purrr::is_empty(args)) {
  # Reference location
  place <- ifelse(is.na(args[1]), place, args[1])
  # Save output 
  save_data <- ifelse(is.na(args[2]), save_data, args[2])
  print_script_args(place = place, save_data = save_data)
}
# ******************************************************************************
#                         Fit the discrete models  
# ******************************************************************************
# Default model specification 
spec <- solarModel_spec$new()
spec$set_transform(min_pos = 1, max_pos = 1, delta = 0.05, link = "invgumbel")
spec$set_mean.model(arOrder = 1, maOrder = 0)
spec$set_mixture.model(method = "mixtools")
spec$set_variance.model(garchOrder = 0, archOrder = 0)

# Initialize the list for the models
discrete_models <- list()
discrete_models_tests <- list()
for(nyear in nyears){
  # Test year 
  nyear_test <- nyears + 1
  # Last train year 
  nyear <- as.character(nyear)
  # Max date 
  max_date <- paste0(nyear, "-12-31")
  cli::cli_alert_success(paste0("Fitting model for year: ", nyear, " (Max date: ", max_date, ")"))
  # Model specification
  spec$specification(place, min_date="2005-01-01", from = "2005-01-01", to = max_date)
  # Fit the model
  model <- solarModel$new(spec)
  model$fit()
  # Standard tests
  test_PIT <- solarModel_PIT_test(model, type = "full", nyears = nyear_test)
  test_autocorr <- solarModel_test_autocorr(model, type = "full", method = "bp", lag.max = 10, nyears = nyear_test)
  discrete_models_tests[[nyear]] <- dplyr::tibble(
    Year_train = nyear, 
    Year_test = nyear_test,
    PIT = test_PIT$KS, 
    PIT_p.value = test_PIT$p.value,
    lb =  tail(test_autocorr, 1)$statistic, 
    lb_p.value =  tail(test_autocorr, 1)$p.value
  )
  # Store the model
  discrete_models[[nyear]] <- model$clone(TRUE)
}
# ******************************************************************************
#                                Save the data 
# ******************************************************************************
if (save_data) {
  # Initialize a new folder for the specific output
  dir_output <- file.path(outputs$dir$data$models$radiation$P, place)
  make_new_directory(dir_output)
  # ******************************************
  # Save discrete models 
  save_new_file(dir_output = dir_output, file.name = common_name[1], 
                file.format = "RData", quiet = FALSE, discrete_models)
  # ******************************************
  # Save the tests  
  save_new_file(dir_output = dir_output, file.name = common_name[2], 
                file.format = "RData", quiet = FALSE, discrete_models)
}
