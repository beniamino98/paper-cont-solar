# ---
#' @description
#' Fit monthly DTMC radiation models for one location.
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
#' Rscript scripts/data/s2a-models-radiation-P-DTMC-place.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/data/s2a-models-radiation-P-DTMC-place.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/discrete_models.RData` (discrete radiation models)
#'
#' @outputs
#'   - `data/models/radiation/P/{place}/DTMC_models.RData` (monthly DTMC models)
#'
#' @depends
#'   - `scripts/data/s0-models-radiation-P-discrete-place.R`
#'
#' @tags
#'   - data
#'   - main
#'   - appendix
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "data", "s2a-models-radiation-P-DTMC-place.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Palermo"
# Save the output
save_data <- FALSE
# ***************************** Fixed Arguments ********************************
# Model's directory
dir_models_P <- outputs$dir$data$models$radiation$P
# Reference names
common_name <- "DTMC_models"
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
#                               Load data 
# ******************************************************************************
# Load discrete models
load_data(file.path(dir_models_P, place), "discrete_models")
# ******************************************************************************
#                             Fit HMM models 
# ******************************************************************************
# 2) Fit HMM models 
DTMC_models <- list()
nyear <- "2016"
for(nyear in nyears){
  max_date <- as.Date(paste0(nyear, "-12-31"))
  message("Fitting solar model for training year: ", nyear, " (Max date: ",  as.character(max_date), ")")
  # Reference year 
  nyear <- as.character(nyear)
  # Reference discrete model 
  model <- discrete_models[[nyear]]
  # Continuous model 
  DTMC_model <- radiationModel_CTMC_fit(model, p0 = c(0.5, 0.5), maxit = 1000, tol = 0.0005)
  # Store EM parameters 
  DTMC_model$params$mu_EM <- DTMC_model$params$mu
  DTMC_model$params$sig_EM <- DTMC_model$params$sig
  # Continuous model 
  DTMC_models[[nyear]] <- DTMC_model
}
#nyear <- "2022"
# Reference discrete model 
#model <- discrete_models[[nyear]]
#DTMC_1 <- radiationModel_CTMC_fit_direct(model, p0 = c(0.5, 0.5), maxit = 20, tol = 0.005) 
#DTMC_2 <- radiationModel_CTMC_fit(model, p0 = c(0.5, 0.5), maxit = 1000, tol = 0.005) 
# ******************************************************************************
#                             Save the data 
# ******************************************************************************
if (save_data) {
  # Initialize a new folder for the specific output
  dir_output <- file.path(outputs$dir$data$models$radiation$P, place)
  # ******************************************
  # Save HMM models 
  save_new_file(dir_output = dir_output, file.name = common_name[1], 
                file.format = "RData", quiet = FALSE, DTMC_models)
  
}
