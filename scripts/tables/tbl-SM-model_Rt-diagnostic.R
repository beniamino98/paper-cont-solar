# ---
#' @description
#' Generate Supplementary Material diagnostic tables for CTMC radiation models.
#'
#' @section `supplementary-material`
#' @label `tbl-SM-model-Rt-diagnostic`
#' @name `tbl_SM_model_Rt_diagnostic`
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr
#'
#' @arguments
#'   - param[1] (save_data): write generated tables ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/tables/tbl-SM-model_Rt-diagnostic.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/Bologna/radiation_models_CTMC_diagnostic.RData`
#'   - `data/models/radiation/P/Roma/radiation_models_CTMC_diagnostic.RData`
#'   - `data/models/radiation/P/Palermo/radiation_models_CTMC_diagnostic.RData`
#'
#' @outputs
#'   - `outputs$table[["tbl_SM_model_Rt_diagnostic"]][[place]]`
#'   - `outputs$tex[["tbl_SM_model_Rt_diagnostic"]][[place]]`
#'
#' @depends
#'   - `scripts/data/SM1-models-radiation-P-diagnostic-place.R`
#'   - `scripts/tables/zzz.R`
#'
#' @tags
#'   - tables
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
source(file.path("scripts", "tables", "zzz.R"))
print_script_info(file.path("scripts", "tables", "tbl-SM-model_Rt-diagnostic.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Save the output
save_data <- TRUE
# ***************************** Fixed Arguments ********************************
# Model's directory
dir_models_P <- outputs$dir$data$models$radiation$P
# Table labels
tbl_labels <- "tbl-SM-model-Rt-diagnostic"
# Table names
tbl_names <- stringr::str_replace_all(tbl_labels, "-", "_")
# Reference locations
places <- c("Bologna", "Palermo", "Roma")
# ******************************************************************************
#                             Inputs (Command line)
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
# Load diagnostic data
radiation_models_diagnostic <- list()
for(place in places){
  load_data(file.path(dir_models_P, place), "radiation_models_CTMC_diagnostic")
  radiation_models_diagnostic[[place]] <- radiation_models_CTMC_diagnostic
}
# ******************************************************************************
#                         Generating functions
# ******************************************************************************
# Table with KS and Ljung-Box tests
tab_SM_diagnostic <- function(place, radiation_models_diagnostic){
  # Diagnostic 
  data_diagnostic <- dplyr::bind_rows(radiation_models_diagnostic[[place]]) %>%
    select(Train, Test,
           `$\\text{LB}(\\tilde{\\varepsilon}_t)$` = "LB_u",
           `$\\text{LB}(\\tilde{\\varepsilon}^2_t)$` = "LB2_u",
           `$\\text{LB}(\\tilde{z}_t)$` = "LB_r",
           `$\\text{LB}(\\tilde{z}_t^2)$` = "LB2_r",
           `$\\text{KS}(U_{Y_t})$` = "ks_Y", 
           `$\\text{KS}(U_{R_t})$` = "ks_R",
           `$\\ell\\text{(Train)}$` = "loglik.train",
           `$\\ell\\text{(Test)}$` = "loglik.test") 
  data_diagnostic
}
# ******************************************************************************
#                          Generate R Tables 
# ******************************************************************************
tab <- list()
for(place in places){
  # Tables in R 
  tab[[place]] <- tab_SM_diagnostic(place, radiation_models_diagnostic = radiation_models_diagnostic)
}
outputs$table[[tbl_names[1]]] <- tab
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[1], "\033[0m generated!", "\n"))
# ******************************************************************************
#                               Save data
# ******************************************************************************
if (save_data) {
  save(outputs, file = "outputs.RData")
  cli::cli_alert_success(paste0("Tables: ", paste0(purrr::map_chr(tbl_names, ~paste0("\033[1;35m", .x, "\033[0m")), collapse = " - "),  " saved!", "\n"))
}
