# ---
#' @description
#' Fit electricity spot-price models under the historical probability measure P.
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse
#'
#' @arguments
#'   - param[1] (save_data): write generated data ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/data/s3-models-electricity-P.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/GME_daily.csv` (PUN price data)
#'
#' @outputs
#'   - `data/models/electricity/P/electricity_models.RData` (electricity models under P)
#'
#' @depends
#'   - none
#'
#' @tags
#'   - data
#'   - main
#'   - appendix
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "data", "s3-models-electricity-P.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Save the output
save_data <- FALSE
# ***************************** Fixed Arguments ********************************
# Directory for PUN data
dir_GME <- outputs$dir_GME
# Reference name
common_name <- "electricity_models"
# Years of training
nyears <- outputs$nyears
# ******************************************************************************
#                           Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly=TRUE)
if (!purrr::is_empty(args)) {
  # Save output 
  save_data <- ifelse(is.na(args[1]), save_data, args[1])
  print_script_args(save_data = save_data)
}
# ******************************************************************************
#                               Load data 
# ******************************************************************************
# Load GME data 
PUN_data <- data <- readr::read_csv(dir_GME, show_col_types = FALSE)  %>%
  dplyr::select(date, PUN) %>%
  dplyr::mutate(PUN = as.numeric(PUN)/1000) %>%
  dplyr::filter(date <= as.Date("2024-01-01") & date >= as.Date("2005-01-01"))
# ******************************************************************************
#                         Fit the electricity models 
# ******************************************************************************
# Initialize the list for the models
electricity_models <- setNames(vector("list", length(nyears)), as.character(nyears))
# Routine
for(i in 1:length(nyears)){
  # Maximum reference date 
  max_date <- as.Date(paste0(nyears[i], "-12-31")) 
  message("Fitting model for year: ", nyears[i], " (Max date: ", max_date, ")")
  # Fit electricity model 
  electricity_models[[i]] <- electricityModel$new(PUN_data, max_date = max_date)
}
# ******************************************************************************
#                             Save the data 
# ******************************************************************************
if (save_data) {
  # Initialize a new folder for the specific output
  dir_output <- outputs$dir$data$models$electricity$P
  # ******************************************
  # Save radiation models 
  save_new_file(dir_output = dir_output, file.name = common_name, 
                file.format = "RData", quiet = FALSE, electricity_models)
  
}
