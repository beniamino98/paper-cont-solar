# ---
#' @description
#' Generate estimated-parameter tables for CTMC radiation models at one location.
#'
#' @section `main` and `supplementary-material`
#' @label `tbl-model_Rt-mean`, `tbl-model_Rt-variance`, `tbl-model_Rt-mixture`
#' @name `tbl_model_Rt_mean`, `tbl_model_Rt_variance`, `tbl_model_Rt_mixture`
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
#' Rscript scripts/tables/tbl-model_Rt-place.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/tables/tbl-model_Rt-place.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData`
#'
#' @outputs
#'   - `outputs$table[["tbl_model_Rt_mean"]][[place]]`
#'   - `outputs$table[["tbl_model_Rt_variance"]][[place]]`
#'   - `outputs$table[["tbl_model_Rt_mixture"]][[place]]`
#'   - `outputs$tex[["tbl_model_Rt_mean"]][[place]]`
#'   - `outputs$tex[["tbl_model_Rt_variance"]][[place]]`
#'   - `outputs$tex[["tbl_model_Rt_mixture"]][[place]]`
#'
#' @depends
#'   - `scripts/data/s2b-models-radiation-P-CTMC-place.R`
#'   - `scripts/tables/zzz.R`
#'
#' @tags
#'   - tables
#'   - main
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
source(file.path("scripts", "tables", "zzz.R"))
print_script_info(file.path("scripts", "tables", "tbl-model_Rt-place.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Bologna"
# Save the output
save_data <- TRUE
# ***************************** Fixed Arguments ********************************
# Model's directory
dir_models_P <- outputs$dir$data$models$radiation$P
# Table labels
tbl_labels <- c("tbl-model_Rt-mean", "tbl-model_Rt-variance", "tbl-model_Rt-mixture")
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
# Load HMM models
load_data(file.path(dir_models_P, place), "radiation_models_CTMC")
radiation_models <- radiation_models_CTMC
# ******************************************************************************
#                           Generating functions  
# ******************************************************************************
#' Solar model parameters (Step 1.)
#' Clear-sky, alpha and beta, seasonal mean of Y
#' @param model_Rt A fitted radiation model.
#' @param digits Integer. Number of significant digits.
#' @return A tibble with clear-sky, transformation, and seasonal-mean parameters.
#' @examples
#' model_Rt <- radiation_models[[1]]
#' digits = 3
tbl_model_Rt_mean <- function(model_Rt, digits = 3){
  # Model's specification 
  spec <- model_Rt$model$spec
  # Reference year 
  nyear_train <- paste0("2005-",lubridate::year(max(spec$dates$train$to)))
  # Number of observations
  nobs <- spec$dates$train$nobs
  # Transform parameters
  alpha <- spec$transform$alpha
  beta  <- spec$transform$beta
  # Seasonal model clear-sky
  Ct_tidy <- spec$seasonal_model_Ct$tidy
  Ct_tidy <- purrr::map2_chr(Ct_tidy$estimate, Ct_tidy$std.error, ~format_param_error(.x, .y, digits = 3))
  Ct_tidy <- dplyr::bind_rows(setNames(Ct_tidy, paste0("$\\delta_", 0:3, "$")))
  # Seasonal model clear-sky
  Yt_bar_tidy <- spec$seasonal.mean$tidy
  Yt_bar_tidy <- purrr::map2_chr(Yt_bar_tidy$estimate, Yt_bar_tidy$std.error, ~format_param_error(.x, .y, digits = 3))
  Yt_bar_tidy <- dplyr::bind_rows(setNames(Yt_bar_tidy, paste0("$a_", 0:2, "$")))
  
  df_params <- dplyr::bind_cols(
    `Train years` = nyear_train,
    `Obs.` = nobs,
    Ct_tidy, 
    `$\\alpha$` = format_param_error(alpha, error = NULL, digits = digits),
    `$\\beta$` = format_param_error(beta, error = NULL, digits = digits),
    Yt_bar_tidy
  )
  return(df_params)}
# ******************************************************************************
#                             Generate R Table 
# ******************************************************************************
# Construct the R object
outputs$table[[tbl_names[1]]][[place]] <- purrr::map_df(radiation_models, ~tbl_model_Rt_mean(.x))
# ******************************************************************************
#                           Generate LaTeX Table  
# ******************************************************************************
# Construct the LaTeX object
outputs$tex[[tbl_names[1]]][[place]] <- outputs$table[[tbl_names[1]]][[place]]  %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 2) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>%
  # Cols 1,2, (Train Years, Number of obs)
  kableExtra::column_spec(2, border_left = FALSE, border_right = TRUE) %>%
  # Cols 3,4,5,6 (Clear-sky)
  kableExtra::column_spec(6, border_left = FALSE, border_right = TRUE) %>%
  # Cols 7,8 (alpha, beta)
  kableExtra::column_spec(8, border_left = FALSE, border_right = TRUE) %>%
  # Cols 9,10,11 (sesonal mean)
  kableExtra::add_header_above(c(" ", " ", "Clear-sky ($C_t$)" = 4, "Bounds" = 2, "Seasonal mean ($\\\\bar{Y}_t$)" = 3), escape = FALSE, bold = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 7.5)
# ******************************************************************************
# Remove "\\[!h\\]" that creates conflicts in quarto
outputs$tex[[tbl_names[1]]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_names[1]]][[place]][1], "\\[!h\\]", "")
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[1], "\033[0m generated!", "\n"))
# ******************************************************************************
#                           Generating functions  
# ******************************************************************************
#' Solar model mean reversion and variance parameters
#' @param model_Rt A fitted radiation model.
#' @param digits Integer. Number of significant digits.
#' @return A tibble with mean-reversion and seasonal-variance parameters.
#' @examples
#' model_Rt <- radiation_models[[1]]
#' digits = 3
tbl_model_Rt_variance <- function(model_Rt, digits = 3){
  # Model's specification 
  spec <- model_Rt$model$spec
  # Reference year 
  nyear_train <- paste0("2005-",lubridate::year(max(spec$dates$train$to)))
  # AR parameters 
  phi <- spec$mean.model$phi
  theta <- -log(phi)
  std.error.phi <- spec$mean.model$std.errors[2]
  std.error.theta <- abs(1 / phi) * std.error.phi
  # **********************************************************************
  # Seasonal variance reparametrization 
  reparam <- model_Rt$seasonal_variance$extra_params$reparam
  # Original parameters 
  b_tidy <- spec$seasonal.variance$tidy
  b_tidy <- purrr::map2_chr(b_tidy$estimate, b_tidy$std.error, ~format_param_error(.x, .y, digits = 3))
  b_tidy <- dplyr::bind_rows(setNames(b_tidy, paste0("$b_", 0:2, "$")))
  # Reparametrized parameters
  gamma_tidy <- reparam$gamma
  # Numeric Jacobian for std. errors 
  params <- c(spec$seasonal.variance$tidy$estimate, theta)
  J_gamma <- numDeriv::jacobian(function(params) reparam_seasonal_function(params[1:3], params[4], omega = 2*base::pi/365)$gamma, params)
  gamma.std.errors <- sqrt(diag(J_gamma %*% diag(c(spec$seasonal.variance$tidy$std.error, std.error.theta)^2) %*% t(J_gamma)))
  gamma_tidy <- purrr::map2_chr(gamma_tidy, gamma.std.errors, ~format_param_error(.x, .y, digits = 3))
  gamma_tidy <- dplyr::bind_rows(setNames(gamma_tidy, paste0("$\\gamma_", 0:2, "$")))
  # Continuous time parameters 
  c_tidy <- reparam$c_
  params <- c(spec$seasonal.variance$tidy$estimate, theta)
  J_c <- numDeriv::jacobian(function(params) reparam_seasonal_function(params[1:3], params[4], omega = 2*base::pi/365)$c_, params)
  c.std.errors <- sqrt(diag(J_c %*%  diag(c(spec$seasonal.variance$tidy$std.error, std.error.theta)^2)%*% t(J_c)))
  c_tidy <- purrr::map2_chr(c_tidy, c.std.errors, ~format_param_error(.x, .y, digits = 3))
  c_tidy <- dplyr::bind_rows(setNames(c_tidy, paste0("$c_", 0:2, "$")))
  # Psi params
  psi_1 <- 1 - exp(-2*theta*1) * cos(spec$seasonal.variance$omega * 1)
  psi_2 <- exp(-2*theta*1) * sin(spec$seasonal.variance$omega * 1)
  # Delta method for std. errors 
  dpsi1_dtheta <-  2 * 1 * exp(-2 * theta * 1) * cos(spec$seasonal.variance$omega * 1)
  dpsi2_dtheta <- -2 * 1 * exp(-2 * theta * 1) * sin(spec$seasonal.variance$omega * 1)
  # Std. errors 
  std.error.psi1 <- abs(dpsi1_dtheta) * std.error.theta
  std.error.psi2 <- abs(dpsi2_dtheta) * std.error.theta
  
  df_params <- dplyr::bind_cols(
    `Train years` = nyear_train,
    `$\\theta$` = format_param_error(theta, std.error.theta, digits = 3),
    `$\\psi_1$` = format_param_error(psi_1, std.error.psi1, digits = 3),
    `$\\psi_2$` = format_param_error(psi_2, std.error.psi2, digits = 3),
    b_tidy,
    c_tidy,
  )
  return(df_params)}
# ******************************************************************************
#                             Generate R Table 
# ******************************************************************************
# Construct the R object
outputs$table[[tbl_names[2]]][[place]] <- purrr::map_df(radiation_models, ~tbl_model_Rt_variance(.x))
# ******************************************************************************
#                           Generate LaTeX Table  
# ******************************************************************************
# Construct the LaTeX object
outputs$tex[[tbl_names[2]]][[place]] <- outputs$table[[tbl_names[2]]][[place]]  %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 2) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>%
  # Col 1, Train Years
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  # Col 2 Theta
  kableExtra::column_spec(2, border_left = FALSE, border_right = TRUE) %>%
  # Col 3,4 psi
  kableExtra::column_spec(4, border_left = FALSE, border_right = TRUE) %>%
  # Col 5,6,7 b_tidy
  kableExtra::column_spec(7, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 7.5)
# ******************************************************************************
# Remove "\\[!h\\]" that creates conflicts in quarto
outputs$tex[[tbl_names[2]]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_names[2]]][[place]][1], "\\[!h\\]", "")
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[2], "\033[0m generated!", "\n"))
# ******************************************************************************
#                           Generating functions  
# ******************************************************************************
#' CTMC parameters 
#' @param model_Rt A fitted `radiationModel_CTMC` object.
#' @param digits Integer. Number of significant digits.
#' @return A tibble with monthly CTMC emission and transition parameters.
#' @examples
#' model_Rt <- radiation_models[[1]]
#' digits = 3
tbl_model_Rt_mixture <- function(model_Rt, digits = 3){
  # Mixture models 
  NM_model <- model_Rt$model$spec$mixture.model
  # Number of observations 
  nobs <- purrr::map_dbl(NM_model$model, ~nrow(.x$responsabilities))
  
  if (class(model_Rt)[1] == "radiationModel_CTMC"){
    # Mean parameters
    mu <- bind_rows(model_Rt$CTMC$params$mu_EM)
    sd <- bind_rows(model_Rt$CTMC$params$sig_EM)
    # Prob: P(1 | 1)
    p11 <- sapply(model_Rt$CTMC$params$Pm, function(x)x[1,1])
    p01 <- sapply(model_Rt$CTMC$params$Pm, function(x)x[1,2])
    p10 <- sapply(model_Rt$CTMC$params$Pm, function(x)x[2,1])
    p00 <- sapply(model_Rt$CTMC$params$Pm, function(x)x[2,2])
    probs <- tibble(p11 = p11, p01 = p01, p10 = p10, p00 = p00)
    df_display <- bind_cols(Month = lubridate::month(1:12, label = TRUE), mu, sd, probs, n = nobs)%>%
      dplyr::select(Month, n, p11, p01, p10, p00, mu1, mu2, sd1, sd2)
    #dplyr::select(Month, n, p11, p01, p10, p00, mu1, sd1, mu2, sd2)
    colnames(df_display) <- c("Month", "Obs.", "$p_{1 \\mid 1}$", "$p_{0 \\mid 1}$",
                              "$p_{1 \\mid 0}$", "$p_{0 \\mid 0}$",
                              "$\\mu_{1}$", "$\\mu_{0}$", 
                              "$\\sigma_{1}$", "$\\sigma_{0}$")
    #"$\\mu_{1}$", "$\\sigma_{1}$", 
    #"$\\mu_{0}$", "$\\sigma_{0}$")
    
  } else {
    df_display <- NM_model$coefficients %>%
      dplyr::left_join(dplyr::select(NM_model$moments, Month, mean, variance, skewness, kurtosis), by = "Month") %>%
      dplyr::mutate(Month = lubridate::month(Month, label = TRUE), n = nobs) %>%
      dplyr::select(Month, n, p1, mu1, sd1, mu2, sd2, mean, variance, skewness, kurtosis)
    colnames(df_display) <- c("Month", "N", "$\\mathbb{P}(B = 1)$",
                              "$\\mu_{1}$", "$\\sigma_{1}$", 
                              "$\\mu_{0}$", "$\\sigma_{0}$",
                              "Mean", "Variance", "Skewness", "Kurtosis")
  }
  
  df_display <- dplyr::mutate_if(df_display, is.numeric, round, digits)
  return(df_display)}
# ******************************************************************************
#                             Generate R Table 
# ******************************************************************************
# Construct the R object
outputs$table[[tbl_names[3]]][[place]] <- list()
for(nyear in names(radiation_models)){
  outputs$table[[tbl_names[3]]][[place]][[nyear]] <- tbl_model_Rt_mixture(radiation_models[[nyear]])
}
# ******************************************************************************
#                           Generate LaTeX Table  
# ******************************************************************************
# Construct the LaTeX object
outputs$tex[[tbl_names[3]]][[place]] <- outputs$table[[tbl_names[3]]][[place]][["2022"]]  %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 2) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>%
  kableExtra::add_header_above(c(" ", " ", "Transition probabilities" = 4, "$\\\\mu_{b}$" = 2, "$\\\\sigma_{b}$" = 2), escape = FALSE, bold = TRUE) %>%
  kableExtra::column_spec(2, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(6, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(8, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 7.5)
# ******************************************************************************
# Remove "\\[!h\\]" that creates conflicts in quarto
outputs$tex[[tbl_names[3]]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_names[3]]][[place]][1], "\\[!h\\]", "")
cli::cli_alert_info(paste0("Table: ", "\033[1;35m", tbl_names[3], "\033[0m generated!", "\n"))
# ******************************************************************************
#                             Save the data 
# ******************************************************************************
if (save_data) {
  save(outputs, file = "outputs.RData")
  cli::cli_alert_success(paste0("Tables: ", paste0(purrr::map_chr(tbl_names, ~paste0("\033[1;35m", .x, "\033[0m")), collapse = " - "),  " saved!", "\n"))
}
