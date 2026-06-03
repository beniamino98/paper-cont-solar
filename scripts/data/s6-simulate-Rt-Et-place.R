# ---
#' @description
#' Simulate joint solar-radiation and electricity scenarios for one location.
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse, mixtools
#'
#' @arguments
#'   - param[1] (place): reference location, e.g. "Bologna".
#'   - param[2] (nsim): number of Monte Carlo scenarios.
#'   - param[3] (save_data): write generated data ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/data/s6-simulate-Rt-Et-place.R "Bologna" "5000" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/data/s6-simulate-Rt-Et-place.R "$place" "5000" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData` (CTMC radiation models)
#'   - `data/models/electricity/Q/electricity_models.RData` (electricity models under Q)
#'   - `data/models/rho/{place}/rho_MC.RData` (monthly residual correlations)
#'
#' @outputs
#'   - `data/scenarios/{place}/scenarios.RData` (joint scenarios)
#'
#' @depends
#'   - `scripts/data/s2b-models-radiation-P-CTMC-place.R`
#'   - `scripts/data/s4-models-electricity-Q.R`
#'   - `scripts/data/s5-models-rho-place.R`
#'
#' @tags
#'   - data
#'   - main
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "data", "s6-simulate-Rt-Et-place.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Bologna"
# Number of scenarios 
nsim <- 5000
# Save the output
save_data <- FALSE
# ***************************** Fixed Arguments ********************************
# Initial seed
seed <- 1
# Reference names
common_name <- "scenarios"
# Years of training
nyears <- as.character(outputs$nyears)
# ******************************************************************************
#                             Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly=TRUE)
if (!purrr::is_empty(args)) {
  # Reference location
  place <- ifelse(is.na(args[1]), place, args[1])
  # Number of scenarios 
  nsim <- as.numeric(ifelse(is.na(args[2]), nsim, args[2]))
  # Save output 
  save_data <- ifelse(is.na(args[3]), save_data, args[3])
  print_script_args(place = place, nsim = nsim, save_data = save_data)
}
# ******************************************************************************
#                                  Load data 
# ******************************************************************************
# Load radiation models
load_data(file.path(outputs$dir$data$models$radiation$P, place), "radiation_models_CTMC")
# Load electricity models
load_data(outputs$dir$data$models$electricity$Q, "electricity_models")
# Load Correlation 
load_data(file.path(outputs$dir$data$models$rho, place), "rho_MC")
# ******************************************************************************
#                                Generate data 
# ******************************************************************************
# Generate scenarios
scenarios <- list()
for(nyear in nyears){
  # Max train year
  nyear <- as.character(nyear)
  # Test year
  nyear_test <- as.numeric(nyear)+1
  print(paste0("Simulating year: ", nyear_test, " with a model trained up to ", nyear))
  # Inputs 
  model_Et <- electricity_models[[nyear]]
  model_Rt <- radiation_models_CTMC[[nyear]]
  rho <- flatten(rho_MC[[nyear]]$rho_model)
  # Conditioning date 
  t_now <- as.Date(paste0(as.numeric(nyear)+1, "-01-01")) - 1
  # Horizon date 
  t_hor <- as.Date(paste0(as.numeric(nyear)+1, "-12-31"))
  scenarios[[nyear]] <- scenarios_ER(model_Et, model_Rt, rho, t_now, t_hor, nsim, seed = seed)
  seed <- seed + 1 
}
# ******************************************************************************
#                                  Save data 
# ******************************************************************************
if (save_data) {
  # Initialize output directory 
  dir_output <- file.path(outputs$dir$data$main, "scenarios")
  make_new_directory(dir_output)
  # Initialize output directory 
  dir_output <- file.path(outputs$dir$data$main, "scenarios", place)
  make_new_directory(dir_output)
  # **********************************************************************
  # Save output 
  save_new_file(dir_output, file.name = common_name[1],
                file.format = "RData", quiet = FALSE, scenarios)
}
