# ---
#' @description
#' Run the complete table-generation pipeline for the replication archive.
#'
#' @author Beniamino Sartini
#' @created 2026-05-30
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse
#'
#' @arguments
#'   - param[1] (save_data): write generated tables ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/s2-tables.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - scripts in `scripts/tables`
#'
#' @outputs
#'   - R table objects and LaTeX table strings stored in `outputs.RData`
#'
#' @depends
#'   - `scripts/s0-load.R`
#'   - `scripts/tables/zzz.R`
#'   - `scripts/tables/tbl-model_Et.R`
#'   - `scripts/tables/tbl-model_Et-moments.R`
#'   - `scripts/tables/tbl-bounds-2gm.R`
#'   - `scripts/tables/tbl-model_Rt-place.R`
#'   - `scripts/tables/tbl-correlations.R`
#'   - `scripts/tables/tbl-mv-place.R`
#'   - `scripts/tables/tbl-mv-full.R`
#'   - `scripts/tables/tbl-SM-model_Rt-diagnostic.R`
#'   - `scripts/tables/tbl-SM-model_Rt-forecasts.R`
#'   - `scripts/tables/tbl-SM-bounds-2gm.R`
#'
#' @tags
#'   - pipeline
#'   - tables
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "s2-tables.R"))
# ******************************************************************************
#                                  Inputs
# ******************************************************************************
# Save generated tables
save_data <- "TRUE"
# ***************************** Fixed Arguments ********************************
# Reference locations
places <- c("Bologna", "Roma", "Palermo")
# Scripts directory
dir_scripts <- "scripts/tables"
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
  print_script_args(save_data = save_data)
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
#                              Appendix tables
# ******************************************************************************
# Electricity model tables
run_pipeline_script("tbl-model_Et.R", save_data)
run_pipeline_script("tbl-model_Et-moments.R", save_data)
# Two-Gaussian approximation tables
run_pipeline_script("tbl-bounds-2gm.R", save_data)
# ******************************************************************************
#                         Main and location-level tables
# ******************************************************************************
for (place in places) {
  # Radiation model tables
  run_pipeline_script("tbl-model_Rt-place.R", place, save_data)
  # Correlation tables
  run_pipeline_script("tbl-correlations.R", place, save_data)
  # Mean-variance tables by location
  run_pipeline_script("tbl-mv-place.R", place, save_data)
}
# Aggregate mean-variance tables across locations
run_pipeline_script("tbl-mv-full.R", save_data)
# ******************************************************************************
#                         Supplementary Material tables
# ******************************************************************************
# CTMC diagnostic tables
run_pipeline_script("tbl-SM-model_Rt-diagnostic.R", save_data)
# Forecast tables
run_pipeline_script("tbl-SM-model_Rt-forecasts.R", save_data)
# Multi-horizon two-Gaussian bound tables
run_pipeline_script("tbl-SM-bounds-2gm.R", save_data)
