# ---
#' @description
#' Generate simulation-based moments for SoRad and SoREd contracts at one location.
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
#' Rscript scripts/data/s7-solarOptions-moments-place.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/data/s7-solarOptions-moments-place.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData` (CTMC radiation models)
#'   - `data/scenarios/{place}/scenarios.RData` (joint scenarios)
#'
#' @outputs
#'   - `data/solarOptions/moments/{place}/moments_index.RData` (index-contract moments)
#'   - `data/solarOptions/moments/{place}/moments_strip.RData` (strip-contract moments)
#'
#' @depends
#'   - `scripts/data/s6-simulate-Rt-Et-place.R`
#'   - `scripts/data/s2b-models-radiation-P-CTMC-place.R`
#'
#' @tags
#'   - data
#'   - main
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "data", "s7-solarOptions-moments-place.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Bologna"
# Save the output
save_data <- FALSE
# ***************************** Fixed Arguments ********************************
# Reference names
common_name <- c("moments_index", "moments_strip")
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
  # Save output 
  save_data <- ifelse(is.na(args[2]), save_data, args[2])
  print_script_args(place = place, save_data = save_data)
}
# ******************************************************************************
#                                  Load data 
# ******************************************************************************
# Load radiation models
load_data(file.path(outputs$dir$data$models$radiation$P, place), "radiation_models_CTMC")
# Load scenarios 
load_data(file.path("data/scenarios", place), "scenarios")
# ******************************************************************************
#                                Generate data 
# ******************************************************************************
# Same strike
K_fun <- function(n) radiation_models_CTMC[[1]]$Rt_bar(n)
moments_index <- list()
moments_strip <- list()
for(nyear in nyears){
  message("Computing moments for year: ", nyear, "\r", appendLF = FALSE)
  nyear <- as.character(nyear)
  preproc <- scenarios[[nyear]]
  moments_index[[nyear]] <- solarOptionIDX_moments_scenario(preproc, K_fun = K_fun)
  moments_strip[[nyear]] <- solarOption_moments_scenario(preproc, K_fun = K_fun)
}
# ******************************************************************************
#                                  Save data 
# ******************************************************************************
if (save_data) {
  # Initialize output directory 
  dir_output <- file.path(outputs$dir$data$main, "solarOptions")
  make_new_directory(dir_output)
  # Initialize output directory 
  dir_output <- file.path(dir_output, "moments")
  make_new_directory(dir_output)
  # Initialize output directory 
  dir_output <- file.path(dir_output, place)
  make_new_directory(dir_output)
  # **********************************************************************
  # Save outputs
  save_new_file(dir_output, file.name = common_name[1],
                file.format = "RData", quiet = FALSE, moments_index)
  save_new_file(dir_output, file.name = common_name[2],
                file.format = "RData", quiet = FALSE, moments_strip)
}
