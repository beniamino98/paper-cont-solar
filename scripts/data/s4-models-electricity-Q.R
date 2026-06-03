# ---
#' @description
#' Fit electricity spot-price models under the risk-neutral measure Q.
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse
#'
#' @arguments
#'   - param[1] (save_data): write generated data ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/data/s4-models-electricity-Q.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/electricity/P/electricity_models.RData` (electricity models under P)
#'
#' @outputs
#'   - `data/models/electricity/Q/electricity_models.RData` (electricity models under Q)
#'
#' @depends
#'   - `scripts/data/s3-models-electricity-P.R`
#'
#' @tags
#'   - data
#'   - main
#'   - appendix
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "data", "s4-models-electricity-Q.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Save the output
save_data <- FALSE
# ***************************** Fixed Arguments ********************************
# Directory for PUN data
dir_GME <- outputs$dir_GME
# Reference name
common_name <- "electricity_models"
# Years of training
nyears <- outputs$nyears
# ******************************************************************************
#                             Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly=TRUE)
if (!purrr::is_empty(args)) {
  # Save output 
  save_data <- ifelse(is.na(args[1]), save_data, args[1])
  print_script_args(save_data = save_data)
}
# ******************************************************************************
#                                 Load data 
# ******************************************************************************
# Load P models 
load_data(outputs$dir$data$models$electricity$P, common_name)
# ******************************************************************************
#                               Generate data 
# ******************************************************************************
# Calibrate Q models
i <- 1
for(i in 1:length(electricity_models)){
  message("Calibrating dQdP model for year: ", names(electricity_models[i]))
  model_Et <- electricity_models[[i]]$clone(deep = TRUE)
  electricity_models[[i]] <- calibrate_dQdP_electricity(tau_max = 365, r = 0, model_Et, tau_hl = 65, quiet = FALSE)
  electricity_models[[i]]$change_measure("Q")
}
# ******************************************************************************
#                              Save the data 
# ******************************************************************************
if (save_data) {
  # Initialize a new folder for the specific output
  dir_output <- outputs$dir$data$models$electricity$Q
  # ******************************************
  # Save radiation models 
  save_new_file(dir_output = dir_output, file.name = common_name, 
                file.format = "RData", quiet = FALSE, electricity_models)
  
}
