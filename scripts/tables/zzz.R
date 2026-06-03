# ---
#' @description
#' Define formatting helpers shared by the table-generation scripts.
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse, kableExtra
#'
#' @inputs
#'   - none
#'
#' @outputs
#'   - helper functions used by `scripts/tables/*.R`
#'
#' @depends
#'   - none
#'
#' @tags
#'   - tables
#'   - internal
# ---
format_param_error <- function(value, error = NULL, digits = 3, scientific = FALSE){
  x_value <- format(value, digits = digits, scientific = scientific)
  if (is.null(error)) {
    out_value <- x_value 
  } else {
    x_error <- format(error, digits = digits, scientific = scientific)
    out_value <- paste0("$\\underset{", "(", x_error, ")", "}{", x_value, "}$")
  }
  return(out_value)
}

detect_signif <- function(x){
  if (is.nan(x)) {
    return("")
  }
  if(x < 0.01){
    "***"
  } else if(x >= 0.01 & x <= 0.05){
    "**"
  } else if(x > 0.05){
    ""
  }
}

format_perc <- function(x, digits = 1){
  paste0(round(x*100, digits = digits), "\\%")
}

# Convert R tables in LaTeX 
tab_from_R_to_TeX <- function(data, digits = 3) {
  kab <- data %>%
    dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
    knitr::kable(booktabs = TRUE,
                 longtable = FALSE,
                 label = NA,
                 format = "latex",
                 linesep = "",
                 caption = NA,
                 escape = FALSE) %>% 
    kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 6)
  kab[1] <- stringr::str_replace(kab[1], "\\[!h\\]", "")
  kab}

# Helper forecasts metrics
forecasts_metrics <- function(x_true, x_pred){
  # Errors 
  errors <- na.omit(x_true - x_pred)
  abs.errors <- abs(errors)
  rel.errors <- errors / x_true
  tibble(
    SSE = sum(errors^2),
    MSE = mean(errors^2),
    RMSE = sqrt(MSE),
    MAE = mean(abs.errors),
    MAPE = mean(abs(rel.errors))
  )
}
