# ---
#' @description
#' Generate short- and long-horizon out-of-sample radiation forecasts at one location.
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
#' Rscript scripts/data/SM2-models-radiation-P-forecasts-place.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/data/SM2-models-radiation-P-forecasts-place.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData` (CTMC radiation models)
#'   - `data/models/radiation/P/{place}/moments/moments_short.RData` (one-day moments)
#'   - `data/models/radiation/P/{place}/moments/moments_long.RData` (full-year moments)
#'
#' @outputs
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC_forecasts.RData` (forecast data)
#'
#' @depends
#'   - `scripts/data/s2b-models-radiation-P-CTMC-place.R`
#'   - `scripts/data/s9-models-radiation-P-moments-place.R`
#'
#' @tags
#'   - data
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "data", "SM2-models-radiation-P-forecasts-place.R"))
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
common_name <- "radiation_models_CTMC_forecasts"
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
# Load CTMC models
load_data(file.path(dir_models_P, place), "radiation_models_CTMC")
# Short-term moments
load_data(file.path(dir_models_P, place, "moments"), "moments_short")
# Long-term moments
load_data(file.path(dir_models_P, place, "moments"), "moments_long")
# ******************************************************************************
#                                 Helper 
# ******************************************************************************
# Auxiliary function to compute expectation 
fun_e_R <- function(moments_nyear){
  e_Rt <- function(mom){
    M_Y <- c(mom$M_Y1, mom$M_Y0)
    S_Y <- c(mom$S_Y1, mom$S_Y0)
    p_T <- c(mom$p1, mom$p1)
    alpha <- mom$alpha
    beta  <- mom$beta
    C_T   <- mom$Ct
    # Pdf of Yt  
    pdf_Y <- function(x) dmixnorm(x, M_Y, S_Y, p_T)
    # pdf of Rt 
    pdf_R <-  dsolarGHI(x, C_T, alpha, beta, pdf_Y, link = "invgumbel")
    # Compute x*pdf_R(x)
    integrate(function(x) x * pdf_R(x), lower = mom$RT_min, upper = mom$RT_max)$value
  }
  purrr::map_dbl(1:nrow(moments_nyear), ~e_Rt(moments_nyear[.x,]))
}
# ******************************************************************************
#                              Generate data 
# ******************************************************************************
# Initialize a list
radiation_models_CTMC_forecasts <- list(short = list(), long = list())
# ******************************************************************************
moments <- moments_short
for(nyear in names(moments)){
  e_R <- fun_e_R(moments[[nyear]])
  R <- filter(radiation_models_CTMC[[nyear]]$model$data, Year == as.numeric(nyear)+1)$GHI
  radiation_models_CTMC_forecasts$short[[nyear]] <- tibble(Year = nyear, e_R = e_R, R = R)
}
# ******************************************************************************
# Long term-forecast
moments <- moments_long
for(nyear in names(moments)){
  e_R <- fun_e_R(moments[[nyear]])
  R <- filter(radiation_models_CTMC[[nyear]]$model$data, Year == as.numeric(nyear)+1)$GHI
  radiation_models_CTMC_forecasts$long[[nyear]] <- tibble(Year = nyear, e_R = e_R, R = R)
}
# ******************************************************************************
#                             Save the data 
# ******************************************************************************
if (save_data) {
  # Output directory 
  dir_output <- file.path(outputs$dir$data$models$radiation$P, place)
  # ******************************************
  # Save radiation models 
  save_new_file(dir_output = dir_output, file.name = common_name[1], 
                file.format = "RData", quiet = FALSE, radiation_models_CTMC_forecasts)
}

