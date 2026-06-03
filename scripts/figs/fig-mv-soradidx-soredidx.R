# ---
#' @description
#' Generate demand-supply figures for SoRadIDX and SoREdIDX contracts.
#'
#' @section `main`
#' @label `fig-mv-soradidx`, `fig-mv-soredidx`
#' @name `fig-mv-soradidx`, `fig-mv-soredidx`
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse, mixtools
#'
#' @arguments
#'   - param[1] (place): reference location, e.g. "Bologna".
#'   - param[2] (nyear): reference year, e.g. "2014".
#'   - param[3] (save_plot): write the generated figures ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/figs/fig-mv-soradidx-soredidx.R "Bologna" "2014" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/figs/fig-mv-soradidx-soredidx.R "$place" "2014" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/scenarios/{place}/scenarios.RData` (joint scenarios)
#'
#' @outputs
#'   - `figs/fig-mv/soradidx/{place}/fig-mv-soradidx.pdf`
#'   - `figs/fig-mv/soradidx/{place}/fig-mv-soredidx.pdf`
#'
#' @depends
#'   - `scripts/data/s6-simulate-Rt-Et-place.R`
#'
#' @tags
#'   - figures
#'   - main
# ---
# Load the required functions
suppressMessages(source(file.path("scripts", "s0-load.R")))
source(file.path("scripts", "figs", "functions", "fig_demand_supply.R"))
load("outputs.RData")
print_script_info(file.path("scripts", "figs", "fig-mv-soradidx-soredidx.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference place 
place <- "Bologna"
# Reference year
nyear <- "2017"
# Save plot
save_plot <- FALSE
# ***************************** Fixed Arguments ********************************
# Common name 
common_name <- "fig-mv"
# Risk aversion for SoRadIDX 
nu_sorad <- c(nu_b = 0.01, nu_s = 0.01)
# Risk aversion for SoREdIDX 
nu_sored <- c(nu_b = 0.1, nu_s = 0.1)
# ******************************************************************************
#                             Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly=TRUE)
if (!purrr::is_empty(args)) {
  # Reference location
  place <- ifelse(is.na(args[1]), place, args[1])
  # Reference train year
  nyear <- as.character(ifelse(is.na(args[2]), nyear, args[2]))
  # Save figure
  save_plot <- ifelse(is.na(args[3]), save_plot, args[3])
  print_script_args(place = place, nyear = nyear, save_plot = save_plot)
}
# ******************************************************************************
#                                  Load data 
# ******************************************************************************
# Load scenarios 
load(file.path("data/scenarios", place, "scenarios.RData"))
scenario <- scenarios[[nyear]]
# Compute moments 
moments <- solarOptionIDX_mv_scenario(scenario)
# Moments for contract
sorad <- moments$sorad
sored <- moments$sored
# ******************************************************************************
#                              Generate figure
# ******************************************************************************
# Risk aversion 
nu_s <- nu_sorad[1]
nu_b <- nu_sorad[2]
# Supply and demand for SoRadIDX
supply_demand_sorad <- supplyDemand_mv(sorad$M_Gamma, sorad$v_Gamma, sorad$S_R_Gamma, sorad$v_Gamma, r = 0, tau = 365, w_Gamma = 1)
# Figure 
fig <- fig_demand_supply_sorad(supply_demand_sorad, nu_b, nu_s)
# ******************************************************************************
#                               Save figure
# ******************************************************************************
if (save_plot){
  control <- outputs$dir$figs$control
  dir_output <- file.path(outputs$dir$figs$main, common_name)
  make_new_directory(dir_output)
  dir_output <- file.path(dir_output, "soradidx")
  make_new_directory(dir_output)
  dir_output <- file.path(dir_output, place)
  make_new_directory(dir_output)
  # Save figure
  save_new_fig(dir_output, fig, fig.name = paste0(common_name, "-soradidx"), file.format = stringr::str_remove_all(control$format, "\\."), 
               quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
  if (place == "Bologna"){
    save_new_fig(outputs$dir$figs$main, fig, fig.name = paste0(common_name, "-soradidx"), file.format = stringr::str_remove_all(control$format, "\\."), 
                 quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
  }
}
# ******************************************************************************
#                              Generate figure
# ******************************************************************************
# Risk aversion 
nu_s <- nu_sored[1]
nu_b <- nu_sored[2]
# Supply and demand for SoRadIDX
supply_demand_sored_uh <- supplyDemand_mv(sored$M_E_Gamma, sored$v_E_Gamma, sored$S_ER_EGamma, sored$v_E_Gamma, r = 0, tau = 365, w_Gamma = 1)
supply_demand_sored_h <- supplyDemand_mv(sored$M_E_Gamma, sored$v_E_Gamma, sored$S_ER_E_mid_E_strip, sored$v_E_Gamma_mid_E_strip, r = 0, tau = 365, w_Gamma = 1)
# Figure 
fig <- fig_demand_supply_sored(supply_demand_sored_h, supply_demand_sored_uh, nu_b, nu_s,  r_min_max = c(0.98, 1.02), digits = 1)
# ******************************************************************************
#                               Save figure
# ******************************************************************************
if (save_plot){
  control <- outputs$dir$figs$control
  dir_output <- file.path(outputs$dir$figs$main, common_name)
  make_new_directory(dir_output)
  dir_output <- file.path(dir_output, "soradidx")
  make_new_directory(dir_output)
  dir_output <- file.path(dir_output, place)
  make_new_directory(dir_output)
  # Save figure
  save_new_fig(dir_output, fig, fig.name = paste0(common_name, "-soredidx"), file.format = stringr::str_remove_all(control$format, "\\."), 
               quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
  if (place == "Bologna"){
    # Save figure
    save_new_fig(outputs$dir$figs$main, fig, fig.name = paste0(common_name, "-soredidx"), file.format = stringr::str_remove_all(control$format, "\\."), 
                 quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
  }
  
  
}
