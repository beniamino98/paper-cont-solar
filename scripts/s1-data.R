# ---
#' @description
#' Run the complete data-generation pipeline for the replication archive.
#'
#' @author Beniamino Sartini
#' @created 2026-05-30
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse
#'
#' @arguments
#'   - param[1] (save_data): write generated data ("TRUE" or "FALSE").
#'   - param[2] (nsim): number of Monte Carlo scenarios, e.g. "5000".
#'   - param[3] (nyear_bounds): training year for two-Gaussian bounds, e.g. "2022".
#'   - param[4] (h_bounds): comma-separated horizons for bounds, e.g. "1,2,3,5,10,15,30".
#'
#' @example
#' Rscript scripts/s1-data.R "TRUE" "5000" "2022" "1,2,3,5,10,15,30"
#'
#' @inputs
#'   - `outputs.RData`
#'   - scripts in `scripts/data`
#'
#' @outputs
#'   - fitted model, scenario, moment, bound, diagnostic, and forecast objects under `data`
#'
#' @depends
#'   - `scripts/s0-load.R`
#'   - `scripts/data/s0-models-radiation-P-discrete-place.R`
#'   - `scripts/data/s1-models-radiation-P-IID-place.R`
#'   - `scripts/data/s2a-models-radiation-P-DTMC-place.R`
#'   - `scripts/data/s2b-models-radiation-P-CTMC-place.R`
#'   - `scripts/data/s3-models-electricity-P.R`
#'   - `scripts/data/s4-models-electricity-Q.R`
#'   - `scripts/data/s5-models-rho-place.R`
#'   - `scripts/data/s6-simulate-Rt-Et-place.R`
#'   - `scripts/data/s7-solarOptions-moments-place.R`
#'   - `scripts/data/s8-solarOptions-mv-place.R`
#'   - `scripts/data/s9-models-radiation-P-moments-place.R`
#'   - `scripts/data/s10-bounds-2gm.R`
#'   - `scripts/data/SM1-models-radiation-P-diagnostic-place.R`
#'   - `scripts/data/SM2-models-radiation-P-forecasts-place.R`
#'
#' @tags
#'   - pipeline
#'   - data
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "s1-data.R"))
# ******************************************************************************
#                                  Inputs
# ******************************************************************************
# Save generated data
save_data <- "TRUE"
# Number of Monte Carlo scenarios
nsim <- "5000"
# Reference year for bounds
nyear_bounds <- "2022"
# Horizons for bounds
h_bounds <- c(1, 2, 3, 5, 10, 15, 30)
# ***************************** Fixed Arguments ********************************
# Reference locations
places <- c("Bologna", "Roma", "Palermo")
# Scripts directory
dir_scripts <- "scripts/data"
# R executable used by the orchestration scripts
r_bin <- Sys.getenv("RSCRIPT_BIN", unset = "")
if (!nzchar(r_bin)) {
  r_bin <- if (nzchar(Sys.which("RScript"))) "RScript" else "Rscript"
}
# ******************************************************************************
#                             Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly = TRUE)
if (!purrr::is_empty(args)) {
  # Save output
  save_data <- ifelse(is.na(args[1]), save_data, args[1])
  # Number of simulations
  nsim <- ifelse(length(args) >= 2 && !is.na(args[2]), args[2], nsim)
  # Reference year for bounds
  nyear_bounds <- ifelse(length(args) >= 3 && !is.na(args[3]), args[3], nyear_bounds)
  # Horizons for bounds
  if (length(args) >= 4 && !is.na(args[4])) {
    h_bounds <- as.numeric(strsplit(args[4], ",", fixed = TRUE)[[1]])
  }
  print_script_args(save_data = save_data, nsim = nsim, nyear_bounds = nyear_bounds,
                    h_bounds = paste(h_bounds, collapse = ","))
}
# ******************************************************************************
#                                 Functions
# ******************************************************************************
# Run an R pipeline script and stop if it fails.
run_pipeline_script <- function(script, ...) {
  script <- file.path(dir_scripts, script)
  script_args <- as.character(c(...))
  if (!file.exists(script)) {
    stop("Missing pipeline script: ", script)
  }
  message("Running: ", paste(c(r_bin, script, script_args), collapse = " "))
  status <- system2(r_bin, args = c(script, script_args))
  if (!identical(as.integer(status), 0L)) {
    stop("Pipeline script failed: ", script)
  }
  invisible(status)
}
# ******************************************************************************
#                            Generate data: s0-s10
# ******************************************************************************
# s0: discrete radiation models
for (place in places) {
  run_pipeline_script("s0-models-radiation-P-discrete-place.R", place, save_data)
}
# s1: IID radiation models
for (place in places) {
  run_pipeline_script("s1-models-radiation-P-IID-place.R", place, save_data)
}
# s2a: DTMC radiation models
for (place in places) {
  run_pipeline_script("s2a-models-radiation-P-DTMC-place.R", place, save_data)
}
# s2b: CTMC radiation models
for (place in places) {
  run_pipeline_script("s2b-models-radiation-P-CTMC-place.R", place, save_data)
}
# s3: electricity model under P
run_pipeline_script("s3-models-electricity-P.R", save_data)
# s4: electricity model under Q
run_pipeline_script("s4-models-electricity-Q.R", save_data)
# s5: residual correlations
for (place in places) {
  run_pipeline_script("s5-models-rho-place.R", place, save_data)
}
# s6: joint scenarios
for (place in places) {
  run_pipeline_script("s6-simulate-Rt-Et-place.R", place, nsim, save_data)
}
# s7: solar-option moments
for (place in places) {
  run_pipeline_script("s7-solarOptions-moments-place.R", place, save_data)
}
# s8: mean-variance and hedging data
for (place in places) {
  run_pipeline_script("s8-solarOptions-mv-place.R", place, save_data)
}
# s9: radiation moment diagnostics
for (place in places) {
  run_pipeline_script("s9-models-radiation-P-moments-place.R", place, save_data)
}
# s10: two-Gaussian bounds
for (h in h_bounds) {
  run_pipeline_script("s10-bounds-2gm.R", nyear_bounds, h, save_data)
}
# ******************************************************************************
#                            Generate data: SM1-SM2
# ******************************************************************************
# SM1: CTMC diagnostics
for (place in places) {
  run_pipeline_script("SM1-models-radiation-P-diagnostic-place.R", place, save_data)
}
# SM2: CTMC forecasts
for (place in places) {
  run_pipeline_script("SM2-models-radiation-P-forecasts-place.R", place, save_data)
}
