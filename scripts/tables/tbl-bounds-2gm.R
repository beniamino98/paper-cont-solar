# ---
#' @description
#' Generate Appendix tables validating the two-Gaussian approximation.
#'
#' @section `appendix`
#' @label `tbl-bounds-2gm-cdf`, `tbl-bounds-2gm-moments`
#' @name `tbl_bounds_2gm_cdf`, `tbl_bounds_2gm_moments`
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
#' Rscript scripts/tables/tbl-bounds-2gm.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/bounds/2022/bounds_2gm-1.RData`
#'
#' @outputs
#'   - `outputs$table[["tbl_bounds_2gm_cdf"]]`
#'   - `outputs$table[["tbl_bounds_2gm_moments"]]`
#'   - `outputs$tex[["tbl_bounds_2gm_cdf"]]`
#'   - `outputs$tex[["tbl_bounds_2gm_moments"]]`
#'
#' @depends
#'   - `scripts/data/s10-bounds-2gm.R`
#'   - `scripts/tables/zzz.R`
#'
#' @tags
#'   - tables
#'   - appendix
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
source(file.path("scripts", "tables", "zzz.R"))
print_script_info(file.path("scripts", "tables", "tbl-bounds-2gm.R"))
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
# Table labels
tbl_labels <- c("tbl-bounds-2gm-cdf", "tbl-bounds-2gm-moments")
# Table names
tbl_names <- stringr::str_replace_all(tbl_labels, "-", "_")
# Reference locations
places <- setNames(outputs$places, c("Bologna","Rome", "Palermo"))
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
bounds_2gm_list <- list()
for(place in places){
  load_data(file.path(dir_models_P, place, "bounds", nyear), "bounds_2gm-1")
  bounds_2gm_list[[place]] <- bounds_2gm
}
# ******************************************************************************
#                           Generating function
# ******************************************************************************
# CDF approximation
tab_bounds_2gm_cdf <- function(place = "Bologna", bounds_2gm_list){
  
  data <- head(bounds_2gm_list[[place]]$bounds, 365)
  kab <- bind_rows(data %>%
                     summarise_if(is.numeric, mean), 
                   data %>%
                     summarise_if(is.numeric, quantile, probs = 0.05),
                   data %>%
                     summarise_if(is.numeric, quantile, probs = 0.95),
                   data %>%
                     mutate(Statistic = ifelse(p_diff < 0.7, "Mean (<80%)", "Mean (>80%)"))%>%
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
  nobs <- c(365, round(365*0.05), round(365*0.05), sum(data$p_diff < 0.7), sum(data$p_diff >= 0.7))
  kab <- bind_cols(Statistic = c("Mean", "Quantile ($5\\%$)", "Quantile ($95\\%$)", "Mean given state probability ($<70\\%$)", "Mean given state probability ($\\ge70\\%$)"), Obs = nobs,  kab)
  colnames(kab) <- c("$\\text{ }$", "Obs", "$d_K$", "$d_{\\max}$", "$d_{\\text{avg}}$", "JS ($Y$)", "JS ($R$)", 
                     "Integral", "Remainder", "$M_{f}^{\\small \\text{2GM}}$")
  kab
}
# Moments approximation
tab_bounds_2gm_moments <- function(place = "Bologna", bounds_2gm_list){
  data <- head(bounds_2gm_list[[place]]$moments, 365)
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
  kab <- bind_cols(Statistic = c("Numerical", "2GM", "RMSE", "MAE", "MAPE ($\\%$)"), kab)
  
  kab %>%
    mutate_if(is.numeric, round, digits = 3)%>%
    mutate_if(is.numeric, format, digits = 3)
  
  # Column names 
  colnames(kab) <- c(" ", "$\\bar{M}_Y$", "$\\bar{S}^2_Y$", "$\\bar{M}_R$", "$\\bar{S}^2_R$", "$\\bar{M}_{\\Gamma}$", "$\\bar{S}^2_{\\Gamma}$")
  kab
}
# ******************************************************************************
#                     Generate table: tbl-bounds-2gm-cdf
# ******************************************************************************
# R table 
data_kable <- purrr::map_df(places, ~tab_bounds_2gm_cdf(.x, bounds_2gm_list))
outputs$table[[tbl_names[1]]] <- data_kable
# ******************************************************************************
# TeX Table 
outputs$tex[[tbl_names[1]]]<- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 2) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::pack_rows(names(places[1]), 1, 5, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>% 
  kableExtra::pack_rows(names(places[2]), 6, 10, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>% 
  kableExtra::pack_rows(names(places[3]), 11, 15, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>% 
  kableExtra::column_spec(2, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(7, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::row_spec(5, hline_after = TRUE) %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  #kableExtra::add_header_above(c(" ", "SoREdIDX" = 4, "Electricity futures" = 3, "P&L" = 3, "Hedging Effect" = 4), escape = TRUE, bold = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 7)
outputs$tex[[tbl_names[1]]][1] <- stringr::str_replace(outputs$tex[[tbl_names[1]]][1], "\\[!h\\]", "")
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[1], "\033[0m generated!", "\n"))
# ******************************************************************************
#                   Generate table: tbl-bounds-2gm-moments
# ******************************************************************************
data_kable <- purrr::map_df(places, ~tab_bounds_2gm_moments(.x, bounds_2gm_list))
outputs$table[[tbl_names[2]]] <- data_kable
# ******************************************************************************
# TeX Table 
outputs$tex[[tbl_names[2]]]<- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 2) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::pack_rows(names(places[1]), 1, 5, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>% 
  kableExtra::pack_rows(names(places[2]), 6, 10, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>% 
  kableExtra::pack_rows(names(places[3]), 11, 15, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>% 
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(3, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::row_spec(5, hline_after = TRUE) %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::add_header_above(c(" ", "$Y$" = 2, "$R$" = 2, "$\\\\Gamma$" = 2), escape = FALSE, bold = FALSE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 7)
# ******************************************************************************
outputs$tex[[tbl_names[2]]][1] <- stringr::str_replace(outputs$tex[[tbl_names[2]]][1], "\\[!h\\]", "")
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[2], "\033[0m generated!", "\n"))
# ******************************************************************************
#                               Save data
# ******************************************************************************
if (save_data) {
  save(outputs, file = "outputs.RData")
  cli::cli_alert_success(paste0("Tables: ", paste0(purrr::map_chr(tbl_names, ~paste0("\033[1;35m", .x, "\033[0m")), collapse = " - "),  " saved!", "\n"))
}
