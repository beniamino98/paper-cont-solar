# ---
#' @description
#' 
#' Generate demand-supply figures for SoRad contracts.
#'
#' @section `main`
#' @label `fig-mv-sorad`
#' @name `fig-mv-sorad`
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
#'   - param[4] (save_exm): write the generated example ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/figs/fig-mv-sorad-sored.R "Bologna" "2014" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/figs/fig-mv-sorad-sored.R "$place" "2014" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData` (CTMC radiation models)
#'
#' @outputs
#'   - `figs/fig-mv/sorad/{place}/fig-mv-sorad.pdf`
#'   - `data/exm/exm_sorad.RData`
#'
#' @depends
#'   - `scripts/data/s2b-models-radiation-P-CTMC-place.R`
#'
#' @tags
#'   - figures
#'   - main
# ---
# Load the required functions
suppressMessages(source(file.path("scripts", "s0-load.R")))
source(file.path("scripts", "figs", "functions", "fig_demand_supply.R"))
load("outputs.RData")
print_script_info(file.path("scripts", "figs", "fig-mv-sorad-sored.R"))
# C functions to speed up the computations
source("scripts/functions/radiationModel/ctmc-integrals-C-wrappers.R")
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference place 
place <- "Bologna"
# Reference year
nyear <- "2014"
# Save plot
save_plot <- FALSE
# Save the example
save_exm <- FALSE
# ***************************** Fixed Arguments ********************************
# Common name 
common_name <- "fig-mv"
# Risk aversion 
nu_s <- 0.2
nu_b <- 0.2
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
  # Save the example
  save_exm <- ifelse(is.na(args[4]), save_exm, args[4])
  print_script_args(place = place, nyear = nyear, save_plot = save_plot, save_exm = save_exm)
}
# ******************************************************************************
#                                  Load data 
# ******************************************************************************
# Load radiation models
load_data(file.path(outputs$dir$data$models$radiation$P, place), "radiation_models_CTMC")
ref_nyear <- as.numeric(nyear)
# Reference model
model_Rt <- radiation_models_CTMC[[as.character(ref_nyear)]]$clone(TRUE)
# ******************************************************************************
#                               Generate data 
# ******************************************************************************
# Today date 
t_now <- as.Date(paste0(as.numeric(nyear), "-12-31"))
# Horizon date 
t_hor <- as.Date(paste0(as.numeric(nyear)+1, "-05-07"))
# Compute moments 
sorad <- solarOption_moments_model_sorad(model_Rt, t_now, t_hor, seq_date = FALSE)
# ******************************************************************************
#                              Generate Example 
# ******************************************************************************
# Radiation data 
R0 <- filter(model_Rt$model$data, date == t_now)$GHI
RT <- filter(model_Rt$model$data, date == t_hor)$GHI
print(paste0("Radiation in t_now (", as.character(t_now), "): ", round(R0, 4), " kWh/m2\n"))
print(paste0("Radiation in t_hor (", as.character(t_hor), "): ", round(RT, 4), " kWh/m2\n"))
# Derivatives data 
K <- sorad$K
print(paste0("Strike price for maturity t_hor (", as.character(t_hor), "): ", round(K, 4), " kWh/m2\n"))
print(paste0("Realized SoRad payoff at maturity t_hor (", as.character(t_hor), "): ", round(solarOption_payoff(RT, K), 4), " kWh/m2\n"))
print(paste0("Expected SoRad payoff at maturity t_hor (", as.character(t_hor), "): ", round(sorad$M_Gamma, 4), " kWh/m2\n"))
sorad$M_Gamma
sorad$v_Gamma
sorad$S_R_Gamma
sorad$cr_R_Gamma
# ******************************************************************************
#                          Store Data for Example
# ******************************************************************************
exm_sorad <- list(
  R0 = R0, 
  RT = RT, 
  K  = K, 
  t_now = t_now, 
  t_hor = t_hor, 
  nu_s = nu_s, 
  nu_b = nu_b
)
exm_sorad <- append(exm_sorad, as.list(unlist(sorad)))
# ******************************************************************************
#                               Save the Example
# ******************************************************************************
if (save_exm) {
  dir_output <- file.path(outputs$dir$data$main, "exm")
  make_new_directory(dir_output)
  save_new_file(dir_output, file.name = "exm_sorad", 
                file.format = "RData", quiet = FALSE, exm_sorad)
}
# ******************************************************************************
#                             Generate Figure 
# ******************************************************************************
# Supply and demand
supply_demand_sorad <- supplyDemand_mv(sorad$M_Gamma, sorad$v_Gamma, sorad$S_R_Gamma, sorad$v_Gamma, r = 0, tau = 365, w_Gamma = 1)
# Figure 
fig <- fig_demand_supply_sorad(supply_demand_sorad, nu_b, nu_s)
fig
# ******************************************************************************
#                               Save figure
# ******************************************************************************
# Save the figure
if (save_plot){
  control <- outputs$dir$figs$control
  dir_output <- file.path(outputs$dir$figs$main, common_name)
  make_new_directory(dir_output)
  dir_output <- file.path(dir_output, "sorad")
  make_new_directory(dir_output)
  dir_output <- file.path(dir_output, place)
  make_new_directory(dir_output)
  # Save figure
  save_new_fig(dir_output, fig, fig.name = paste0(common_name, "-sorad"), file.format = stringr::str_remove_all(control$format, "\\."), 
               quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
  if (place == "Bologna"){
    save_new_fig(outputs$dir$figs$main, fig, fig.name = paste0(common_name, "-sorad"), file.format = stringr::str_remove_all(control$format, "\\."), 
                 quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
  }
}