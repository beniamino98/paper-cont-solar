# ---
#' @description
#' Run the complete figure-generation pipeline for the replication archive.
#'
#' @author Beniamino Sartini
#' @created 2026-05-30
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse
#'
#' @arguments
#'   - param[1] (save_plot): write generated figures ("TRUE" or "FALSE").
#'   - param[2] (nyear_acf): training year for the ACF figure, e.g. "2022".
#'   - param[3] (nyear_mv): reference year for the demand-supply figures, e.g. "2014".
#'   - param[4] (run_extra_figures): run optional SoRadIDX/SoREdIDX demand-supply figures.
#'
#' @example
#' Rscript scripts/s3-figures.R "TRUE" "2022" "2014" "FALSE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - scripts in `scripts/figs`
#'
#' @outputs
#'   - figures saved under `figs`
#'   - example objects saved under `data/exm` by selected figure scripts
#'
#' @depends
#'   - `scripts/s0-load.R`
#'   - `scripts/figs/fig-PUN-price.R`
#'   - `scripts/figs/fig-acf.R`
#'   - `scripts/figs/fig-GHI-sim.R`
#'   - `scripts/figs/fig-mv-sorad.R`
#'   - `scripts/figs/fig-mv-sored.R`
#'   - `scripts/figs/fig-hedged-vs-unhedged.R`
#'   - `scripts/figs/fig-cum-net-ret-sorad.R`
#'   - `scripts/figs/fig-cum-net-ret-sored.R`
#'   - `scripts/figs/fig-mv-soradidx-soredidx.R`
#'
#' @tags
#'   - pipeline
#'   - figures
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "s3-figures.R"))
# ******************************************************************************
#                                  Inputs
# ******************************************************************************
# Save generated figures
save_plot <- "TRUE"
# Training year for the ACF figure
nyear_acf <- "2022"
# Reference year for demand-supply figures
nyear_mv <- "2014"
# Generate extra figures not referenced in the main manuscript
run_extra_figures <- FALSE
# ***************************** Fixed Arguments ********************************
# Reference locations
places <- c("Bologna", "Roma", "Palermo")
# Scripts directory
dir_scripts <- "scripts/figs"
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
  save_plot <- ifelse(is.na(args[1]), save_plot, args[1])
  # ACF year
  nyear_acf <- ifelse(length(args) >= 2 && !is.na(args[2]), args[2], nyear_acf)
  # Demand-supply year
  nyear_mv <- ifelse(length(args) >= 3 && !is.na(args[3]), args[3], nyear_mv)
  # Extra figures
  run_extra_figures <- ifelse(length(args) >= 4 && !is.na(args[4]), as.logical(args[4]), run_extra_figures)
  print_script_args(save_plot = save_plot, nyear_acf = nyear_acf, nyear_mv = nyear_mv,
                    run_extra_figures = run_extra_figures)
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
#                                  Figures
# ******************************************************************************
# Electricity price figure
run_pipeline_script("fig-PUN-price.R", save_plot)
# Radiation diagnostics and simulations
for (place in places) {
  run_pipeline_script("fig-acf.R", place, save_plot, nyear_acf)
  run_pipeline_script("fig-GHI-sim.R", place, save_plot)
}
# Demand-supply figures for SoRad and SoREd
for (place in places) {
  run_pipeline_script("fig-mv-sorad.R", place, nyear_mv, save_plot, TRUE)
  run_pipeline_script("fig-mv-sored.R", place, nyear_mv, save_plot, TRUE)
}
# Hedging figures
run_pipeline_script("fig-hedged-vs-unhedged.R", save_plot)
# Cumulative net-return 
run_pipeline_script("fig-cum-net-ret-sored.R", save_plot)
run_pipeline_script("fig-cum-net-ret-sorad.R", save_plot)
# Extra top-level figures not currently referenced in the main manuscript.
if (isTRUE(run_extra_figures)) {
  for (place in places) {
    run_pipeline_script("fig-mv-soradidx-soredidx.R", place, nyear_mv, save_plot)
  }
}
