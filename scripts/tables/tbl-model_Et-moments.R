# ---
#' @description
#' Generate historical and model-implied expectation/variance tables for electricity models.
#'
#' @section `appendix`
#' @label `tbl-model_Et-moments`
#' @name `tbl_model_Et_moments`
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
#' Rscript scripts/tables/tbl-model_Et-moments.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/electricity/P/electricity_models.RData`
#'   - `data/models/electricity/Q/electricity_models.RData`
#'
#' @outputs
#'   - `outputs$table[["tbl_model_Et_moments"]]`
#'   - `outputs$tex[["tbl_model_Et_moments"]]`
#'
#' @depends
#'   - `scripts/data/s3-models-electricity-P.R`
#'   - `scripts/data/s4-models-electricity-Q.R`
#'   - `scripts/tables/zzz.R`
#'
#' @tags
#'   - tables
#'   - appendix
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
source(file.path("scripts", "tables", "zzz.R"))
print_script_info(file.path("scripts", "tables", "tbl-model_Et-moments.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Save the output
save_data <- FALSE
# ***************************** Fixed Arguments ********************************
# Table labels
tbl_labels <- "tbl-model_Et-moments"
# Table names
tbl_names <- stringr::str_replace_all(tbl_labels, "-", "_")
# ******************************************************************************
#                            Inputs (Command line)
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
# Load electricity models under P
load_data(outputs$dir$data$models$electricity$P, "electricity_models")
electricity_models_P <- purrr::map(electricity_models, ~.x$clone(TRUE))
# Load electricity models under lambda^E
load_data(outputs$dir$data$models$electricity$Q, "electricity_models")
electricity_models_lambda_E <- purrr::map(electricity_models, ~.x$clone(TRUE))
# ******************************************************************************
#                           Generating functions  
# ******************************************************************************
#' Build long-horizon electricity-price moment summaries.
#'
#' @param model_Et_P Electricity model under the physical measure.
#' @param model_Et_lambda_E Electricity model with calibrated pricing premium.
#' @param digits Integer. Number of significant digits.
#'
#' @return A tibble with historical and model-implied long-horizon moments.
#' @examples
#' model_Et_P <- electricity_models_P[[1]]
#' model_Et_lambda_E <- electricity_models_lambda_E[[1]]
#' digits <- 3
#' tbl_model_Et_moments(model_Et_P, electricity_models_lambda_E, digits)
tab_model_Et_moments <- function(model_Et_P, model_Et_lambda_E, digits = 3){
  # Training year
  nyear_train <- lubridate::year(max(model_Et_P$date_train))
  # Long term moments under P 
  e_Et_P <- model_Et_P$F_E("2010-01-01", "2014-01-01", 0.05)
  v_Et_P <- model_Et_P$S_E("2010-01-01", "2014-01-01", 0.05)
  # Long term moments (lambda^E)
  e_Et_lambda_E <- model_Et_lambda_E$F_E("2010-01-01", "2014-01-01", 0.05)
  v_Et_lambda_E <- model_Et_lambda_E$S_E("2010-01-01", "2014-01-01", 0.05)
  # Historical average moments 
  data <- dplyr::filter(model_Et_P$data, isTrain)
  e_Et_train <- mean(data$PUN)
  v_Et_train <- var(data$PUN)
  # Moments last train year
  data <- dplyr::filter(model_Et_P$data, Year == nyear_train)
  e_Et_last <- mean(data$PUN)
  v_Et_last <- var(data$PUN)
  # Historical moments test year
  data <- dplyr::filter(model_Et_P$data, Year == nyear_train+1)
  e_Et_test <- mean(data$PUN)
  v_Et_test <- var(data$PUN)

  mom <- dplyr::tibble(
    Year = paste0("2005-", nyear_train),
    Test = as.character(nyear_train + 1),
    e_Et_hist = e_Et_train, 
    sd_Et_hist = sqrt(v_Et_train),
    e_Et_P = e_Et_P, 
    sd_Et_P = sqrt(v_Et_P),
    lambda_E = model_Et_lambda_E$model$lambda,
    e_Et_lambda_E = e_Et_lambda_E,
    sd_Et_lambda_E = sqrt(v_Et_lambda_E),
    e_Et_year = e_Et_last, 
    sd_Et_year = sqrt(v_Et_last),
    e_Et_test =  e_Et_test,
    v_Et_test = sqrt(v_Et_test)
  )
  names(mom) <-  c("Train", "Test", 
                   "$\\mathbb{E}\\{E_t\\}$", "$\\sqrt{\\mathbb{V}\\{E_t\\}}$", 
                   "$\\mathbb{E}\\{E_t\\}$ ", "$\\sqrt{\\mathbb{V}\\{E_t\\}}$ ", 
                   "$\\lambda^E$", "$\\mathbb{E}\\{E_t\\}$  ", "$\\sqrt{\\mathbb{V}\\{E_t\\}}$  ", 
                   "$\\mathbb{E}\\{E_t\\}$   ", "$\\sqrt{\\mathbb{V}\\{E_t\\}}$   ", 
                   "$\\mathbb{E}\\{E_t\\}$    ", "$\\sqrt{\\mathbb{V}\\{E_t\\}}$    ")
  return(mom)
}
# ******************************************************************************
#                           Generate R Table  
# ******************************************************************************
# Construct the R object
data_kable <- purrr::map2_df(electricity_models_P, electricity_models_lambda_E, tab_model_Et_moments)
# Add average between the period 
avg <- dplyr::bind_cols(Train = "Average", Test = "$\\text{}$", dplyr::summarise_if(data_kable, is.numeric, mean))
# Save the final table 
outputs$table[[tbl_names]] <- dplyr::bind_rows(data_kable, avg)
# ******************************************************************************
#                         Generate LaTeX Table  
# ******************************************************************************
# Construct the LaTeX object
outputs$tex[[tbl_names]] <- outputs$table[[tbl_names]]  %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>%
  kableExtra::add_header_above(c(" ", " ", "Train data (historical)" = 2, "Model ($\\\\lambda^{\\\\small \\\\text{E}} = 0$)" = 2, "Model ($\\\\lambda^{\\\\small \\\\text{E}} \\\\neq 0$)" = 3, "Last train year" = 2, "Test year" = 2), 
                               escape = FALSE, bold = TRUE) %>%
  kableExtra::column_spec(2, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(4, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(6, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(9, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(11, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 6)
# ******************************************************************************
# Extra // generated after row_spec removed 
outputs$tex[[tbl_names]][1] <- stringr::str_replace(outputs$tex[[tbl_names]][1], stringr::fixed("\\midrule\\\\\n"), "\\midrule\n")
# Remove "\\[!h\\]" that creates conflicts in quarto
outputs$tex[[tbl_names]][1] <- stringr::str_replace(outputs$tex[[tbl_names]][1], "\\[!h\\]", "")
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[1], "\033[0m generated!", "\n"))
# ******************************************************************************
#                             Save the data 
# ******************************************************************************
if (save_data) {
  save(outputs, file = "outputs.RData")
  cli::cli_alert_success(paste0("Tables: ", paste0(purrr::map_chr(tbl_names, ~paste0("\033[1;35m", .x, "\033[0m")), collapse = " - "),  " saved!", "\n"))
}
