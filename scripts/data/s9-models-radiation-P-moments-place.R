# ---
#' @description
#' Generate one-day and full-year CTMC radiation moments for one location.
#'
#' @section `supplementary-material`
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr
#'
#' @arguments
#'   - param[1] (place): reference location, e.g. "Bologna".
#'   - param[2] (save_data): write generated data ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/data/s9-models-radiation-P-moments-place.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/data/s9-models-radiation-P-moments-place.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData` (CTMC radiation models)
#'
#' @outputs
#'   - `data/models/radiation/P/{place}/moments/moments_short.RData` (one-day moments)
#'   - `data/models/radiation/P/{place}/moments/moments_long.RData` (full-year moments)
#'
#' @depends
#'   - `scripts/data/s2b-models-radiation-P-CTMC-place.R`
#'
#' @tags
#'   - data
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "data", "s9-models-radiation-P-moments-place.R"))
# C functions to speed up the computations
source("scripts/functions/radiationModel/ctmc-integrals-C-wrappers.R")
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference place
place <- "Bologna"
# Save the output
save_data <- TRUE
# ***************************** Fixed Arguments ********************************
# Model's directory
dir_models_P <- outputs$dir$data$models$radiation$P
# Reference years
nyears <- outputs$nyears
# Name
common_name <- c("moments_short", "moments_long")
# ******************************************************************************
#                             Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly=TRUE)
if (!purrr::is_empty(args)) {
  # Reference place
  place <- ifelse(is.na(args[1]), place, args[1])
  # Save output 
  save_data <- ifelse(is.na(args[2]), save_data, args[2])
  print_script_args(place = place, save_data = save_data)
}
# ******************************************************************************
#                                Load data
# ******************************************************************************
# Load CTMC models
load_data(file.path(dir_models_P, place), "radiation_models_CTMC")
# ******************************************************************************
#                               Generate data  
# ******************************************************************************
moments_short   <- list()
moments_long <- list()
for(nyear in nyears){
  # Reference model 
  model_Rt <- radiation_models_CTMC[[as.character(nyear)]]
  # Sequence of dates
  t_now <- as.Date(paste0(as.numeric(nyear), "-12-31"))
  t_hor <- as.Date(paste0(as.numeric(nyear)+1, "-12-31"))
  t_seq <- seq.Date(t_now+1, t_hor, 1)
  # Compute conditional moments 
  print(paste0("(Short term) Conditional moments for year: ", nyear))
  moments_short[[as.character(nyear)]] <- purrr::map_df(t_seq, ~radiationMoments(t_now = .x-1, t_hor = .x, model_Rt = model_Rt, R0 = NULL))
  # Compute conditional moments 
  print(paste0("(Long term) Conditional moments for year: ", nyear))
  moments_long[[as.character(nyear)]] <- purrr::map_df(t_seq, ~radiationMoments(t_now = t_now, t_hor = .x, model_Rt = model_Rt, R0 = NULL))
}
# ******************************************************************************
#                             Save the data 
# ******************************************************************************
if (save_data) {
  # Output directory 
  dir_output <- file.path(outputs$dir$data$models$radiation$P, place, "moments")
  make_new_directory(dir_output)
  # ******************************************
  # Save moments
  save_new_file(dir_output = dir_output, file.name = common_name[1], 
                file.format = "RData", quiet = FALSE, moments_short)
  # Save moments
  save_new_file(dir_output = dir_output, file.name = common_name[2], 
                file.format = "RData", quiet = FALSE, moments_long)
}




