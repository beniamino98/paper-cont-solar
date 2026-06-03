# ---
#' @description
#' Generate estimated-parameter tables for electricity models.
#'
#' @section `appendix`
#' @label `tbl-model_Et`
#' @name `tbl_model_Et`
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
#' Rscript scripts/tables/tbl-model_Et.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/electricity/P/electricity_models.RData`
#'
#' @outputs
#'   - `outputs$table[["tbl_model_Et"]]`
#'   - `outputs$tex[["tbl_model_Et"]]`
#'
#' @depends
#'   - `scripts/data/s3-models-electricity-P.R`
#'   - `scripts/tables/zzz.R`
#'
#' @tags
#'   - tables
#'   - appendix
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
source(file.path("scripts", "tables", "zzz.R"))
print_script_info(file.path("scripts", "tables", "tbl-model_Et.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Save the output
save_data <- TRUE
# ***************************** Fixed Arguments ********************************
# Table labels
tbl_labels <- "tbl-model_Et"
# Table names
tbl_names <- stringr::str_replace_all(tbl_labels, "-", "_")
# ******************************************************************************
#                         Inputs (Command line)
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
# Load Electricity models
load_data(outputs$dir$data$models$electricity$P, "electricity_models")
# ******************************************************************************
#                            Generating function 
# ******************************************************************************
#' Parameters of the electricity models
#' @param model_Et A fitted `electricityModel` object.
#' @param digits Integer. Number of significant digits.
#' @return A tibble with sample size, AR(1), and continuous-time OU parameters.
#' @examples
#' model_Et <- electricity_models[[1]]
#' digits = 3
#' tbl_model_Et(model_Et, digits)
tbl_model_Et <- function(model_Et, digits = 3){
  # Train Years 
  nyear_train <- lubridate::year(max(model_Et$date_train))
  nyears_train <- paste0("2005-", nyear_train)
  # Number of observations
  nobs <- nrow(model_Et$model$data)
  # log-OU model
  OU <- model_Et$model$model
  # AR parameters 
  params.AR <- OU$params.AR
  std.errors.AR <- OU$std.errors.AR
  params.AR  <- format_param_error(params.AR, std.errors.AR, digits = 3)
  names(params.AR) <-  c("$\\phi_0$", "$\\phi_1$","$\\hat{\\sigma_u}$")
  # OU parameters 
  params.OU <- OU$params.OU
  std.errors.OU <- OU$std.errors.OU
  params.OU  <- format_param_error(params.OU, std.errors.OU, digits = 3)
  names(params.OU) <- c("$\\mu_X$", "$\\kappa$", "$\\sigma_X$")
  # Dataset 
  tbl_params <-  dplyr::bind_cols(`Train` = nyears_train,
                                  `Obs.` = nobs, 
                                  dplyr::bind_rows(params.AR), 
                                  dplyr::bind_rows(params.OU))
                           
  return(tbl_params)}
# ******************************************************************************
#                            Generate R Table  
# ******************************************************************************
# Construct the R object
outputs$table[[tbl_names[1]]] <- purrr::map_df(electricity_models, tbl_model_Et)
# ******************************************************************************
#                           Generate LaTeX Table  
# ******************************************************************************
# Construct the LaTeX object
outputs$tex[[tbl_names[1]]] <- outputs$table[[tbl_names[1]]]  %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>%
  kableExtra::add_header_above(c(" ", " ", "AR parameters" = 3, "OU parameters" = 3), escape = FALSE, bold = TRUE) %>%
  kableExtra::column_spec (2, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec (5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 6)
# ******************************************************************************
# Remove "\\[!h\\]" that creates conflicts in quarto
outputs$tex[[tbl_names[1]]][1] <- stringr::str_replace(outputs$tex[[tbl_names[1]]][1], "\\[!h\\]", "")
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[1], "\033[0m generated!", "\n"))
# ******************************************************************************
#                               Save the data 
# ******************************************************************************
if (save_data) {
  save(outputs, file = "outputs.RData")
  cli::cli_alert_success(paste0("Tables: ", paste0(purrr::map_chr(tbl_names, ~paste0("\033[1;35m", .x, "\033[0m")), collapse = " - "),  " saved!", "\n"))
}
