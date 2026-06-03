# ---
#' @description
#' Generate R and LaTeX tables for estimated electricity/radiation residual correlations.
#'
#' @section `appendix` and `supplementary-material`
#' @label `tbl-correlations`
#' @name `tbl_correlations`
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr
#'
#' @arguments
#'   - param[1] (place): reference location, e.g. "Bologna".
#'   - param[2] (save_data): write generated tables ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/tables/tbl-correlations.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/tables/tbl-correlations.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/rho/{place}/rho_MC.RData`
#'
#' @outputs
#'   - `outputs$table[["tbl_correlations"]][[place]]`
#'   - `outputs$tex[["tbl_correlations"]][[place]]`
#'
#' @depends
#'   - `scripts/data/s5-models-rho-place.R`
#'   - `scripts/tables/zzz.R`
#'
#' @tags
#'   - tables
#'   - appendix
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
source(file.path("scripts", "tables", "zzz.R"))
print_script_info(file.path("scripts", "tables", "tbl-correlations.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Bologna"
# Save the output
save_data <- TRUE
# ***************************** Fixed Arguments ********************************
# Model's directory
dir_models_P <- outputs$dir$data$models$rho
# Table labels
tbl_labels <- "tbl-correlations"
# Table names
tbl_names <- stringr::str_replace_all(tbl_labels, "-", "_")
# ******************************************************************************
#                             Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly=TRUE)
if (!purrr::is_empty(args)) {
  # Reference location
  place <- ifelse(is.na(args[1]), place, args[1])
  # Save output 
  save_data <- ifelse(is.na(args[2]), save_data, args[2])
  print_script_args(place = place, save_data = save_data)
}
# ******************************************************************************
#                               Load data 
# ******************************************************************************
# Load correlations 
load_data(file.path(dir_models_P, place), "rho_MC")
rho_models <- rho_MC
# ******************************************************************************
#                           Generating function
# ******************************************************************************
#' Estimated correlations parameters 
#' @param rho_models List of fitted monthly correlation objects.
#' @param digits Integer. Number of significant digits.
#' @return A tibble with state-dependent and unconditional correlations.
#' @examples
#' digits = 3
#' tab_correlations(rho_models, 3)
tab_correlations <- function(rho_models, digits = 3){
  # References years 
  nyears <- names(rho_models)
  # References months
  nmonths <- lubridate::month(1:12, label = TRUE)
  df_rho_1 <- list()
  df_rho_2 <- list()
  df_rho_tot <- list()
  for(nyear in nyears){
    nyear <- as.character(nyear)
    rho_year <-  rho_models[[nyear]]
    # Dataset with correlations
    df_rho <- purrr::map_df(rho_year$rho_model, ~bind_rows(.x))
    # Detect significance 
    rho1 <- paste0(round(df_rho$rho1, digits = digits), purrr::map_chr(rho_year$tests$rho_1.p.value, detect_signif))
    names(rho1) <- nmonths
    rho2 <- paste0(round(df_rho$rho2, digits = digits), purrr::map_chr(rho_year$tests$rho_2.p.value, detect_signif))
    names(rho2) <- nmonths
    rho_tot <- paste0(round(rho_year$tests$rho.emp, digits = digits), purrr::map_chr(rho_year$tests$rho.emp_p.value, detect_signif))
    names(rho_tot) <- nmonths
    
    # Complete dataset
    df_rho_1[[nyear]] <- dplyr::bind_cols(Year = nyear, dplyr::bind_rows(rho1))
    df_rho_2[[nyear]] <- dplyr::bind_cols(Year = nyear, dplyr::bind_rows(rho2))
    df_rho_tot[[nyear]] <- dplyr::bind_cols(Year = nyear, dplyr::bind_rows(rho_tot))
  }
  structure(
    list(
      rho1 = dplyr::bind_rows(df_rho_1),
      rho2 = dplyr::bind_rows(df_rho_2),
      rho_tot = dplyr::bind_rows(df_rho_tot)
    )
  )}
# ******************************************************************************
#                           Generate R Table  
# ******************************************************************************
cat(paste0("Generating table: ", "\033[1;35m", tbl_labels, "\033[0m\n"))
# Construct the R object
tbl_corr <- tab_correlations(rho_models)
# Store the table for a specific location 
outputs$table[[tbl_names[1]]][[place]] <- dplyr::bind_rows(tbl_corr$rho1, tbl_corr$rho2, tbl_corr$rho_tot)
# ******************************************************************************
#                           Generate LaTeX Table  
# ******************************************************************************
# Construct the LaTeX object
outputs$tex[[tbl_names[1]]][[place]] <- outputs$table[[tbl_names[1]]][[place]]  %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               linesep = "",
               format = "latex",
               caption = NA,
               escape = FALSE) %>%
  kableExtra::pack_rows("$\\rho_1 \\\\approx \\\\mathbb{C}r\\\\{dW_t, dW_{1,t}\\\\}$", 1, 10, escape = FALSE, bold = FALSE, latex_gap_space = "0.3em", latex_align = "c", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>% 
  #kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::pack_rows("$\\rho_0 \\\\approx \\\\mathbb{C}r\\\\{dW_t, dW_{0,t}\\\\}$", 11, 20, escape = FALSE, bold = FALSE, latex_gap_space = "0.3em", latex_align = "c", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>% 
  #kableExtra::row_spec(20, hline_after = TRUE) %>%
  kableExtra::pack_rows("$\\rho^{\\\\mathrm{avg}} \\\\approx \\\\mathbb{C}r\\\\{dW_t, dM_{t}\\\\}$", 21, nrow(outputs$table[[tbl_names[1]]][[place]]), escape = FALSE, bold = TRUE, 
                        latex_gap_space = "0.3em", latex_align = "c", hline_after = TRUE, hline_before = TRUE, background = "black!7") %>% 
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 6)
# ******************************************************************************
outputs$tex[[tbl_names[1]]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_names[1]]][[place]][1], "\\[!h\\]", "")
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[1], "\033[0m generated!", "\n"))
# ******************************************************************************
#                             Save the data 
# ******************************************************************************
if (save_data) {
  save(outputs, file = "outputs.RData")
  cli::cli_alert_success(paste0("Tables: ", paste0(purrr::map_chr(tbl_names, ~paste0("\033[1;35m", .x, "\033[0m")), collapse = " - "),  " saved!", "\n"))
}
