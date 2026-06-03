# ---
#' @description
#' Generate Supplementary Material forecast tables for CTMC radiation models.
#'
#' @section `supplementary-material`
#' @label `tbl-SM-model-Rt-forecasts-short`, `tbl-SM-model-Rt-forecasts-long`
#' @name `tbl_SM_model_Rt_forecasts_short`, `tbl_SM_model_Rt_forecasts_long`
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
#' Rscript scripts/tables/tbl-SM-model_Rt-forecasts.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/Bologna/radiation_models_CTMC_forecasts.RData`
#'   - `data/models/radiation/P/Roma/radiation_models_CTMC_forecasts.RData`
#'   - `data/models/radiation/P/Palermo/radiation_models_CTMC_forecasts.RData`
#'
#' @outputs
#'   - `outputs$table[["tbl_SM_model_Rt_forecasts_short"]][[place]]`
#'   - `outputs$table[["tbl_SM_model_Rt_forecasts_long"]][[place]]`
#'   - `outputs$tex[["tbl_SM_model_Rt_forecasts_short"]][[place]]`
#'   - `outputs$tex[["tbl_SM_model_Rt_forecasts_long"]][[place]]`
#'
#' @depends
#'   - `scripts/data/SM2-models-radiation-P-forecasts-place.R`
#'   - `scripts/tables/zzz.R`
#'
#' @tags
#'   - tables
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "tables", "tbl-SM-model_Rt-forecasts.R"))
source(file.path("scripts", "tables", "zzz.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Save the output
save_data <- TRUE
# ***************************** Fixed Arguments ********************************
# Model's directory
dir_models_P <- outputs$dir$data$models$radiation$P
# Table labels
tbl_labels <- c("tbl-SM-model-Rt-forecasts-short", "tbl-SM-model-Rt-forecasts-long")
# Table names
tbl_names <- stringr::str_replace_all(tbl_labels, "-", "_")
# Reference locations
places <- c("Bologna", "Palermo", "Roma")
# ******************************************************************************
#                               Load data 
# ******************************************************************************
radiation_models_forecasts  <- list()
for(place in places){
  # Load forecasts data
  load_data(file.path(dir_models_P, place), "radiation_models_CTMC_forecasts")
  radiation_models_forecasts[[place]]  <- radiation_models_CTMC_forecasts
}
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
#                            Generating functions 
# ******************************************************************************
# Table with Forecasts metric on the out-of-sample year 
tab_SM_forecasts <- function(place, radiation_models_forecasts){
  # Reference years
  nyears <- as.numeric(names(radiation_models_forecasts[[place]]$short))
  # Short term-forecasts
  tab_short <- purrr::map_df(radiation_models_forecasts[[place]]$short, ~forecasts_metrics(.x$R, .x$e_R))
  tab_long  <- purrr::map_df(radiation_models_forecasts[[place]]$long, ~forecasts_metrics(.x$R, .x$e_R))
  
  tab_short <- bind_cols(`Train` = paste0("2005-", nyears),
            `Test` = paste0(nyears+1), tab_short) %>%
    mutate_if(is.numeric, format, digits= 3)
  
  tab_long <- bind_cols(`Train` = paste0("2005-", nyears),
                         `Test` = paste0(nyears+1), tab_long) %>%
    mutate_if(is.numeric, format, digits= 3)
  
  list(long  = tab_long, 
       short = tab_short)
}
# ******************************************************************************
#                             Generate R Tables
# ******************************************************************************
tab1 <- tab2 <- list()
for(place in places){
  # Save R tables
  tab_R <- tab_SM_forecasts(place, radiation_models_forecasts = radiation_models_forecasts)
  tab1[[place]] <- tab_R$short
  tab2[[place]] <- tab_R$long
}
# ******************************************************************************
# Save R tables
outputs$table[[tbl_names[1]]] <- tab1
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[1],  "\033[0m (", place, ")", " generated!", "\n"))
# Save R tables
outputs$table[[tbl_names[2]]] <- tab2
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[2],  "\033[0m (", place, ")", " generated!", "\n"))
# ******************************************************************************
#                                 Save data
# ******************************************************************************
if (save_data) {
  save(outputs, file = "outputs.RData")
  cli::cli_alert_success(paste0("Tables: ", paste0(purrr::map_chr(tbl_names, ~paste0("\033[1;35m", .x, "\033[0m")), collapse = " - "),  " saved!", "\n"))
}
