# ---
#' @description
#' Load package dependencies and source project functions used by the replication pipeline.
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies readxl, solarr, purrr, dplyr, Rcpp
#'
#' @inputs
#'   - `scripts/functions`
#'   - `scripts/functions/scripts`
#'
#' @outputs
#'   - packages and project functions loaded into the active R session
#'
#' @depends
#'   - none
#'
#' @tags
#'   - pipeline
#'   - internal
# ---
# Packages
suppressMessages(library(readxl))
suppressMessages(library(solarr))
suppressMessages(library(purrr))
suppressMessages(library(dplyr))
# install.packages(c("Rcpp"))
suppressMessages(library(Rcpp))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference directory for the functions
dir_functions <- file.path("scripts", "functions")
# Reference directory for utility functions
dir_functions_utils <- file.path(dir_functions, "scripts")
# ******************************************************************************
#                              Load functions 
# ******************************************************************************
# Utils functions
source(file.path(dir_functions_utils, "make_new_directory.R"))
source(file.path(dir_functions_utils, "save_new_file.R"))
source(file.path(dir_functions_utils, "save_new_fig.R"))
source(file.path(dir_functions_utils, "print_script_args.R"))
source(file.path(dir_functions_utils, "print_script_info.R"))
source(file.path(dir_functions_utils, "load_data.R"))
# ******************************************************************************
# C functions
# dyn.load(file.path(dir_functions, "C", "bounds_kernels.so"))
# ******************************************************************************
# Electricity Model
source(file.path(dir_functions, "electricityModel", "logOUModel.R"))
source(file.path(dir_functions, "electricityModel", "electricityModel.R"))
source(file.path(dir_functions, "electricityModel", "exp_decay_fp.R"))
# ******************************************************************************
# Solar Model: IID bernoullis
source(file.path(dir_functions, "radiationModel", "radiationModel-internals.R"))
source(file.path(dir_functions, "radiationModel", "radiationModel-R6.R"))
# ******************************************************************************
# Solar Model: CTMC 
source(file.path(dir_functions, "radiationModel", "dtmc.R"))
source(file.path(dir_functions, "radiationModel", "ctmc.R"))
source(file.path(dir_functions, "radiationModel", "radiationModel_CTMC-R6.R"))
# ******************************************************************************
# Moments function 
source(file.path(dir_functions, "radiationModel", "radiationMoments.R"))
# ******************************************************************************
# Integrals function 
source(file.path(dir_functions, "radiationModel", "ctmc-integrals-C.R"))
# ******************************************************************************
# Compute scenarios with solar radiation and electricity models
source(file.path(dir_functions, "scenarios_ER.R"))
# Continuous time 
source(file.path(dir_functions, "radiationModel", "scenarios_radiationModel_CT.R"))
# ******************************************************************************
# MV demand and supply 
source(file.path(dir_functions, "solarOption", "supplyDemand_mv.R"))
# Moments 
source(file.path(dir_functions, "solarOption", "solarOption_mv.R"))
source(file.path(dir_functions, "solarOption", "solarOption_moments.R"))
source(file.path(dir_functions, "solarOption", "solarOptionIDX_moments.R"))
# Hedging 
source(file.path(dir_functions, "solarOption", "solarOption_hedging.R"))
# ******************************************************************************
# Figures 
source(file.path(dir_functions, "figure_theme.R"))
source(file.path(dir_functions, "zzz.R"))
# ******************************************************************************
