# ---
#' @description
#' Fit CTMC radiation models for one location.
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
#' Rscript scripts/data/s2b-models-radiation-P-CTMC-place.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/data/s2b-models-radiation-P-CTMC-place.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/radiation_models.RData` (IID radiation models)
#'   - `data/models/radiation/P/{place}/DTMC_models.RData` (monthly DTMC models)
#'
#' @outputs
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData` (CTMC radiation models)
#'
#' @depends
#'   - `scripts/data/s1-models-radiation-P-IID-place.R`
#'   - `scripts/data/s2a-models-radiation-P-DTMC-place.R`
#'
#' @tags
#'   - data
#'   - main
#'   - appendix
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "data", "s2b-models-radiation-P-CTMC-place.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Bologna"
# Save the output
save_data <- FALSE
# ***************************** Fixed Arguments ********************************
# Model's directory
dir_models_P <- outputs$dir$data$models$radiation$P
# Reference names
common_name <- "radiation_models_CTMC"
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
# Load IID models
load_data(file.path(dir_models_P, place), "radiation_models")
# Load discrete models
load_data(file.path(dir_models_P, place), "DTMC_models")
# ******************************************************************************
#                          Update radiation models 
# ******************************************************************************
radiation_models_CTMC <- list()
for(nyear in nyears){
  max_date <- as.Date(paste0(nyear, "-12-31"))
  cli::cli_alert_info(paste0("Fitting solar model for training year: ", nyear, " (Max date: ",  as.character(max_date), ") \n"))
  # Reference year 
  nyear <- as.character(nyear)
  # Reference discrete model 
  model_Rt <- radiation_models[[nyear]]$clone(deep = TRUE)
  # DTMC model (A)
  DTMC <- DTMC_models[[nyear]]
  radiation_models_CTMC[[nyear]] <- radiationModel_CTMC$new(model_Rt, DTMC)
}
# ******************************************************************************
#                             Save the data 
# ******************************************************************************
if (save_data) {
  # Initialize a new folder for the specific output
  dir_output <- file.path(outputs$dir$data$models$radiation$P, place)
  # ******************************************
  # Save radiation models 
  save_new_file(dir_output = dir_output, file.name = common_name, 
                file.format = "RData", quiet = FALSE, radiation_models_CTMC)
}
