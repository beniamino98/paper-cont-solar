# ---
#' @description
#' Generate Supplementary Material arrays for multi-horizon two-Gaussian bounds.
#'
#' @section `supplementary-material`
#' @label `tbl-SM-bounds-2gm-cdf`, `tbl-SM-bounds-2gm-moments`
#' @name `tab_SM_bounds_2gm_cdf`, `tab_SM_bounds_2gm_moments`
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
#' Rscript scripts/tables/tbl-SM-bounds-2gm.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/bounds/2022/bounds_2gm-{h}.RData`
#'
#' @outputs
#'   - `outputs$table[["tab_SM_bounds_2gm_cdf"]]`
#'   - `outputs$table[["tab_SM_bounds_2gm_moments"]]`
#'
#' @depends
#'   - `scripts/data/s10-bounds-2gm.R`
#'   - `scripts/tables/zzz.R`
#'
#' @tags
#'   - tables
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
source(file.path("scripts", "tables", "zzz.R"))
print_script_info(file.path("scripts", "tables", "tbl-SM-bounds-2gm.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Save the output
save_data <- TRUE
# ***************************** Fixed Arguments ********************************
# Model's directory
dir_models_P <- outputs$dir$data$models$radiation$P
# Reference year
nyear = "2022"
# Available horizons 
horizon <- c(1, 2, 3, 5, 10, 15, 30)
# Table labels
tbl_labels <- c("tab-SM-bounds-2gm-cdf", "tab-SM-bounds-2gm-moments")
# Table names
tbl_names <- stringr::str_replace_all(tbl_labels, "-", "_")
# Reference locations
places <- setNames(outputs$places, c("Bologna","Rome", "Palermo"))
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
#                           Generating function
# ******************************************************************************
tab_bounds_2gm_cdf <- function(place = "Bologna", nyear = "2022", h = 1){
  load_data(file.path(dir_models_P, place, "bounds", nyear), paste0("bounds_2gm", "-", h))
  data <- head(bounds_2gm$bounds, 365)
  # Median p_diff
  med_p_diff <- quantile(data$p_diff, 0.75)

  kab <- bind_rows(data %>%
                     summarise_if(is.numeric, mean), 
                   data %>%
                     summarise_if(is.numeric, quantile, probs = 0.05),
                   data %>%
                     summarise_if(is.numeric, quantile, probs = 0.95),
                   data %>%
                     mutate(Statistic = paste0("Mean (", ifelse(p_diff < med_p_diff, "<", ">"),  round(med_p_diff*100, 0), "%)")) %>%
                     group_by(Statistic) %>%
                     summarise_if(is.numeric, mean) %>% 
                     select(-Statistic)) %>%
    select(bound, max_abs_err, mean_abs_err, JS_Y, JS_R, integral_part, remainder, M_bound,) %>%
    mutate_all(round, digits = 4) %>%
    mutate(
      bound = format_perc(bound, 2),
      integral_part = format_perc(integral_part, 2),
      remainder = format_perc(remainder, 2),
      max_abs_err = format_perc(max_abs_err, 2),
      mean_abs_err = format_perc(mean_abs_err, 2),
      JS_Y = format_perc(JS_Y/log(2), 2),
      JS_R = format_perc(JS_R/log(2), 2),
    )
  nobs <- c(365, round(365*0.05), round(365*0.05), sum(data$p_diff < med_p_diff), sum(data$p_diff >= med_p_diff))
  kab <- bind_cols(Statistic = c("Mean", "Quantile ($5\\%$)", "Quantile ($95\\%$)", paste0("Mean (", c("<", ">"),  round(med_p_diff*100, 0), "%)")), Obs = nobs,  kab)
  kab
}
# ******************************************************************************
#                             Generate R tables 
# ******************************************************************************
# Extract tables for colnames
tab <- tab_bounds_2gm_cdf(places[1], nyear = "2022", h = 1)
# Initialization 
tab <- array(0,dim = c(5, ncol(tab), length(horizon), length(places)), dimnames = list(NULL, colnames(tab), as.character(horizon), places))
for(place in places){
  for(h in horizon){
    tab[,, as.character(h), as.character(place)] <- as.matrix(tab_bounds_2gm_cdf(place, nyear = "2022", h = h))
  }
}
outputs$table[[tbl_names[1]]] <- tab
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[1], "\033[0m generated!", "\n"))
# ******************************************************************************
#                           Generating function
# ******************************************************************************
tab_bounds_2gm_moments <- function(place = "Bologna", nyear = "2022", h = 1){
  load_data(file.path(dir_models_P, place, "bounds", nyear), paste0("bounds_2gm", "-", h))
  data <- head(bounds_2gm$moments, 365)
  # Standard column names 
  col.names <- c("e_Y", "v_Y", "e_R", "v_R", "e_Gamma", "v_Gamma")
  # True (avg)
  kab_true <- data %>%
    select(e_Y_true, v_Y_true, e_R_true, v_R_true, e_Gamma_true, v_Gamma_true)
  colnames(kab_true) <- col.names
  # 2GM (avg)
  kab_2gm <- data %>%
    select(e_Y_2gm, v_Y_2gm, e_R_2gm, v_R_2gm, e_Gamma_2gm, v_Gamma_2gm)
  colnames(kab_2gm) <- col.names
  # Compute Errors
  errors <- as.matrix(kab_true - kab_2gm)
  rel_errors <- errors / as.matrix(kab_true)
  # RMSE
  RMSE <- apply(errors, 2, function(x) sqrt(mean(x^2)))
  # MAE
  MAE <- apply(errors, 2, function(x) mean(abs(x)))
  # MAPE
  MAPE <- apply(rel_errors, 2, function(x) mean(abs(x)*100))
  # Combine data 
  kab <- bind_rows(summarise_all(kab_true, mean),
                   summarise_all(kab_2gm, mean),
                   RMSE, MAE, MAPE)
  # Add statistic name 
  kab <- bind_cols(Statistic = c("Numerical", "2GM", "RMSE", "MAE", "MAPE ($\\%$)"), kab) %>%
    mutate_if(is.numeric, round, digits = 3) %>%
    mutate_if(is.numeric, format, digits = 3)
  # Column names 
  colnames(kab) <- c(" ", "$\\bar{M}_Y$", "$\\bar{S}^2_Y$", "$\\bar{M}_R$", "$\\bar{S}^2_R$", "$\\bar{M}_{\\Gamma}$", "$\\bar{S}^2_{\\Gamma}$")
  kab
}
# ******************************************************************************
#                           Generate  R Tables 
# ******************************************************************************
# Extract tables for colnames
tab <- tab_bounds_2gm_moments(places[1], nyear = "2022", h = 1)
# Initialization 
tab <- array(0,dim = c(5, ncol(tab), length(horizon), length(places)), 
             dimnames = list(NULL, colnames(tab), as.character(horizon), places))
for(place in places){
  for(h in horizon){
    tab[,, as.character(h), as.character(place)] <- as.matrix(tab_bounds_2gm_moments(place, nyear = "2022", h = h))
  }
}
# ******************************************************************************
outputs$table[[tbl_names[2]]] <- tab
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[2], "\033[0m generated!", "\n"))
# ******************************************************************************
#                                Save data
# ******************************************************************************
if (save_data) {
  save(outputs, file = "outputs.RData")
  cli::cli_alert_success(paste0("Tables: ", paste0(purrr::map_chr(tbl_names, ~paste0("\033[1;35m", .x, "\033[0m")), collapse = " - "),  " saved!", "\n"))
}
