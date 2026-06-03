# ---
#' @description
#' Aggregate location-level empirical application tables into full-sample tables.
#'
#' @section `main`
#' @label `tbl-avg-moments-sorad`, `tbl-hedging-sorad`, `tbl-hedging-soradidx`,
#'        `tbl-avg-moments-sored`, `tbl-hedging-sored-uh`, `tbl-hedging-sored-h`, `tbl-hedging-soredidx-strip`
#' @name `tbl_avg_moments_sorad`, `tbl_hedging_sorad`, `tbl_hedging_soradidx`,
#'       `tbl_avg_moments_sored`, `tbl_hedging_sored_uh`, `tbl_hedging_sored_h`, `tbl_hedging_soredidx_strip`
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
#' Rscript scripts/tables/tbl-mv-full.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - `outputs$table[["tbl_avg_moments_sorad"]][[place]]`
#'   - `outputs$table[["tbl_hedging_sorad"]][[place]]`
#'   - `outputs$table[["tbl_hedging_soradidx"]][[place]]`
#'   - `outputs$table[["tbl_avg_moments_sored"]][[place]]`
#'   - `outputs$table[["tbl_hedging_sored_uh"]][[place]]`
#'   - `outputs$table[["tbl_hedging_sored_h"]][[place]]`
#'   - `outputs$table[["tbl_hedging_soredidx_strip"]][[place]]`
#'
#' @outputs
#'   - `outputs$table[["tbl_avg_moments_sorad"]][["full"]]`
#'   - `outputs$table[["tbl_hedging_sorad"]][["full"]]`
#'   - `outputs$table[["tbl_hedging_soradidx"]][["full"]]`
#'   - `outputs$table[["tbl_avg_moments_sored"]][["full"]]`
#'   - `outputs$table[["tbl_hedging_sored_uh"]][["full"]]`
#'   - `outputs$table[["tbl_hedging_sored_h"]][["full"]]`
#'   - `outputs$table[["tbl_hedging_soredidx_strip"]][["full"]]`
#'
#' @depends
#'   - `scripts/tables/tbl-mv-place.R`
#'   - `scripts/tables/zzz.R`
#'
#' @tags
#'   - tables
#'   - main
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
source(file.path("scripts", "tables", "zzz.R"))
print_script_info(file.path("scripts", "tables", "tbl-mv-full.R"))
# ******************************************************************************
#                               Inputs 
# ******************************************************************************
place <- "full"
# Reference places 
places <- c(Bologna = "Bologna", Rome = "Roma", Palermo = "Palermo")
# Save the output
save_data <- FALSE
# ***************************** Fixed Arguments ********************************
# Table labels
tbl_labels <- c(
  "tbl_avg_moments_sorad",
  "tbl_hedging_sorad",
  "tbl_hedging_soradidx",
  "tbl_avg_moments_sored",
  "tbl_hedging_sored_uh",
  "tbl_hedging_sored_h",
  "tbl_hedging_soredidx_strip"
)
# Table names
tbl_names <- stringr::str_replace_all(tbl_labels, "-", "_")
# ******************************************************************************
#                        Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly=TRUE)
if (!purrr::is_empty(args)) {
  # Save output 
  save_data <- ifelse(is.na(args[1]), save_data, args[1])
  print_script_args(save_data = save_data)
}
# ******************************************************************************
#                   Table: Average daily moments (SoRad)
# ******************************************************************************
# Table name
tbl_name <- "tbl_avg_moments_sorad"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# R Table 
data_kable <- purrr::map_df(places, ~outputs$table[[tbl_name]][[.x]] )
outputs$table[[tbl_name]][[place]] <- data_kable
# ******************************************************************************
# TeX Table 
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::pack_rows(names(places[1]), 1, 11, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[2]), 12, 22, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[3]), 23, 33, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::row_spec(11, hline_after = TRUE) %>%
  kableExtra::row_spec(21, hline_after = TRUE) %>%
  kableExtra::row_spec(22, hline_after = TRUE) %>%
  kableExtra::row_spec(32, hline_after = TRUE) %>%
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(9, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 6)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
# ******************************************************************************
#                          Table: SoRad Supply and Demand 
# ******************************************************************************
# Table name
tbl_name <- "tbl_hedging_sorad"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# R Table 
data_kable <- purrr::map_df(places, ~outputs$table[[tbl_name]][[.x]] )
outputs$table[[tbl_name]][[place]] <- data_kable
# ******************************************************************************
# TeX Table
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::pack_rows(names(places[1]), 1, 11, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[2]), 12, 22, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[3]), 23, 33, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::row_spec(11, hline_after = TRUE) %>%
  kableExtra::row_spec(21, hline_after = TRUE) %>%
  kableExtra::row_spec(22, hline_after = TRUE) %>%
  kableExtra::row_spec(32, hline_after = TRUE) %>%
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(8, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(10, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::add_header_above(c(" ", "SoRad" = 4, "Profit and Loss" = 3, "Hedging effect buyer" = 4), escape = TRUE, bold = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 6)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
# ******************************************************************************
#                         Table: SoRadIDX Supply and Demand 
# ******************************************************************************
# Table name
tbl_name <- "tbl_hedging_soradidx"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# R Table 
data_kable <- purrr::map_df(places, ~outputs$table[[tbl_name]][[.x]] )
outputs$table[[tbl_name]][[place]] <- data_kable
# ******************************************************************************
# TeX Table 
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::pack_rows(names(places[1]), 1, 11, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[2]), 12, 22, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[3]), 23, 33, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::row_spec(11, hline_after = TRUE) %>%
  kableExtra::row_spec(21, hline_after = TRUE) %>%
  kableExtra::row_spec(22, hline_after = TRUE) %>%
  kableExtra::row_spec(32, hline_after = TRUE) %>%
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(8, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(10, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::add_header_above(c(" ", "SoRadIDX" = 4, "Profit and Loss" = 3, "Hedging effect buyer" = 4), escape = TRUE, bold = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 6)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
# ******************************************************************************
#                                Table: SoREd Moments 
# ******************************************************************************
# Table name
tbl_name <- "tbl_avg_moments_sored"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# R Table 
data_kable <- purrr::map_df(places, ~outputs$table[[tbl_name]][[.x]] )
outputs$table[[tbl_name]][[place]] <- data_kable
# ******************************************************************************
# TeX Table 
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::pack_rows(names(places[1]), 1, 11, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[2]), 12, 22, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[3]), 23, 33, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::row_spec(11, hline_after = TRUE) %>%
  kableExtra::row_spec(21, hline_after = TRUE) %>%
  kableExtra::row_spec(22, hline_after = TRUE) %>%
  kableExtra::row_spec(32, hline_after = TRUE) %>%
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(9, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 6)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
# ******************************************************************************
#               Table: SoREd Supply and Demand (unhedged)
# ******************************************************************************
# Table name
tbl_name <- "tbl_hedging_sored_uh"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# R Table 
data_kable <- purrr::map_df(places, ~outputs$table[[tbl_name]][[.x]] )
outputs$table[[tbl_name]][[place]] <- data_kable
# ******************************************************************************
# TeX Table 
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::pack_rows(names(places[1]), 1, 11, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[2]), 12, 22, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[3]), 23, 33, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::row_spec(11, hline_after = TRUE) %>%
  kableExtra::row_spec(21, hline_after = TRUE) %>%
  kableExtra::row_spec(22, hline_after = TRUE) %>%
  kableExtra::row_spec(32, hline_after = TRUE) %>%
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(8, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::add_header_above(c(" ", "$\\\\textbf{SoREd}$" = 4, "$\\\\textbf{Profit and Loss}$" = 3, "$\\\\textbf{Hedging effect buyer}$ ($\\\\times 10^{-2}$)" = 4), escape = FALSE, bold = FALSE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 6)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
# ******************************************************************************
#                 Table: SoREd Supply and Demand (hedged)
# ******************************************************************************
# Table name
tbl_name <- "tbl_hedging_sored_h"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# R Table 
data_kable <- purrr::map_df(places, ~outputs$table[[tbl_name]][[.x]] )
outputs$table[[tbl_name]][[place]] <- data_kable
# ******************************************************************************
# TeX Table 
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::pack_rows(names(places[1]), 1, 11, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[2]), 12, 22, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[3]), 23, 33, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::row_spec(11, hline_after = TRUE) %>%
  kableExtra::row_spec(21, hline_after = TRUE) %>%
  kableExtra::row_spec(22, hline_after = TRUE) %>%
  kableExtra::row_spec(32, hline_after = TRUE) %>%
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(8, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(11, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::add_header_above(c(" ", "$\\\\textbf{SoREd}$" = 4, "$\\\\textbf{Electricity futures}$" = 3, "$\\\\textbf{Profit and Loss}$" = 3, "$\\\\textbf{Hedging effect buyer}$ ($\\\\times 10^{-2}$)" = 4), escape = FALSE, bold = FALSE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 5)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
# ******************************************************************************
#                   Table: SoREdIDX Supply and Demand 
# ******************************************************************************
# Table name
tbl_name <- "tbl_hedging_soredidx_strip"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# R Table 
data_kable <- purrr::map_df(places, ~outputs$table[[tbl_name]][[.x]] )
outputs$table[[tbl_name]][[place]] <- data_kable
# ******************************************************************************
# TeX Table 
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::pack_rows(names(places[1]), 1, 11, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[2]), 12, 22, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::pack_rows(names(places[3]), 23, 33, escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "l", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>%
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(8, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(11, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::row_spec(11, hline_after = TRUE) %>%
  kableExtra::row_spec(21, hline_after = TRUE) %>%
  kableExtra::row_spec(22, hline_after = TRUE) %>%
  kableExtra::row_spec(32, hline_after = TRUE) %>%
  kableExtra::add_header_above(c(" ", "$\\\\textbf{SoREdIDX}$" = 4, "$\\\\textbf{Electricity futures}$" = 3, "$\\\\textbf{Profit and Loss}$" = 3, "$\\\\textbf{Hedging effect buyer}$ ($\\\\times 10^{-2}$)" = 4), escape = FALSE, bold = FALSE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 5)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
#cli::cli_alert_success(paste0("Tables: ", paste0(purrr::map_chr(tbl_names, ~paste0("\033[1;35m", .x, "\033[0m")), collapse = " - "),  " saved!", "\n"))
