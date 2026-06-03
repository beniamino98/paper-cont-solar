# ---
#' @description
#' Compute out-of-sample diagnostic tests for CTMC radiation models at one location.
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse, mixtools
#'
#' @arguments
#'   - param[1] (place): reference location, e.g. "Bologna".
#'   - param[2] (save_data): write generated data ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/data/SM1-models-radiation-P-diagnostic-place.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/data/SM1-models-radiation-P-diagnostic-place.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData` (CTMC radiation models)
#'
#' @outputs
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC_diagnostic.RData` (diagnostic results)
#'
#' @depends
#'   - `scripts/data/s2b-models-radiation-P-CTMC-place.R`
#'
#' @tags
#'   - data
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "data", "SM1-models-radiation-P-diagnostic-place.R"))
source("scripts/functions/radiationModel/radiationModel-CTMC-density-C-wrappers.R")
source("scripts/functions/radiationModel/ctmc-integrals-C-wrappers.R")
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Palermo"
# Save the output
save_data <- FALSE
# ***************************** Fixed Arguments ********************************
# Model's directory
dir_models_P <- outputs$dir$data$models$radiation$P
# Reference names
common_name <- "radiation_models_CTMC_diagnostic"
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
source("scripts/functions/radiationModel/radiationModel-CTMC-density-C-wrappers.R")
# ******************************************************************************
#                               Functions
# ******************************************************************************
#' Compute one-step diagnostic statistics for a fitted CTMC radiation model.
#'
#' @param model_Rt A fitted `radiationModel_CTMC` object.
#'
#' @return A one-row tibble with PIT, Ljung-Box, and KS diagnostics.
#' @examples
#' model_Rt <- radiation_models_CTMC[[5]]
#' radiation_model_CTMC_diagnostic(model_Rt)
radiation_model_CTMC_diagnostic <- function(model_Rt){
  # Reference test year
  nyear_test <- lubridate::year(model_Rt$model$spec$dates$train$to)+1
  # Index train/test 
  isTrain <- which(model_Rt$model$data$Year < nyear_test)
  isTest  <- which(model_Rt$model$data$Year == nyear_test)
  isFull  <- c(isTrain, isTest)
  # Monthly index 
  t_month <- model_Rt$model$data[c(isTrain, isTest),]$Month
  # HMM model 
  CTMC <- model_Rt$CTMC
  # Log-likelihood 
  loglik.train <- dtmc_monthly_EM(CTMC$data[c(isTrain),], CTMC$params, update_emissions = FALSE, update_transitions = FALSE)$loglik
  loglik.test  <- dtmc_monthly_EM(CTMC$data[c(isTest),], CTMC$params, update_emissions = FALSE, update_transitions = FALSE)$loglik
  # Extract data (full)
  data <- model_Rt$model$data[c(isFull),]
  # Clear-sky
  Ct <- data$Ct
  # Solar transform parameters
  alpha <- model_Rt$model$spec$transform$alpha
  beta  <- model_Rt$model$spec$transform$beta
  # Realized radiation
  Rt <- data$GHI
  # Realized transformed variable 
  Yt <- data$Yt
  # Seasonal mean Y
  Yt_bar <- data$Yt_bar
  # Mean-reversion parameter
  theta <- model_Rt$theta
  # Seasonal std. deviation 
  J <- data$sigma_bar
  # Standardized residuals 
  u_tilde <- data$u_tilde
  # Reference dates
  ref_dates <- data$date
  # Predictive probabilities
  alpha_ctmc <- CTMC$alpha[CTMC$data$date %in% ref_dates,]
  # Initialization 
  z_tilde <- c()
  grades_U <- c()
  grades_R <- c()
  i <- 1
  for(i in isTest){
    message("Train: 2005-", nyear_test-1, " - Test (", nyear_test, ") t: ", i, "/", max(isTest), "\r", appendLF = FALSE)
    # Monthly parameters
    mu_m <- CTMC$params$mu[[t_month[i]]]
    sd_m <- CTMC$params$sig[[t_month[i]]]
    pT_m <- drop(alpha_ctmc[i-1,] %*% CTMC$params$Pm[[t_month[i]]])
    # Moments 
    mom <- GM_moments(mu_m, sd_m, pT_m)
    # Standardized residuals 
    z_tilde[i] <- (u_tilde[i] - mom$mean) / sqrt(mom$variance) 
    # Grades of eps
    grades_U[i] <- pmixnorm(u_tilde[i], mu_m, sd_m, pT_m)
    # State-dependent expectation of Y
    M_Y <- Yt_bar[i] + (Yt[i-1] - Yt_bar[i-1]) * exp(-theta) + mu_m * J[i]
    # State-dependent std.deviation of Y
    S_Y <- J[i] * sd_m
    # Grades of eps
    grades_U[i] <- pmixnorm(Yt[i], M_Y, S_Y, pT_m)
    # Grades of R 
    grades_R[i] <- psolarGHI(Rt[i], Ct[i], alpha, beta, function(x) pmixnorm(x, M_Y, S_Y, pT_m))
  }
  # KS-test grades
  ks_Y   <- ks_test(grades_U[isTest], punif)
  ks_R <- ks_test(grades_R[isTest], punif)
  # LB tests 
  LB_u  <- Box.test(na.omit(u_tilde[isTest]), lag = 30)
  LB2_u <- Box.test(na.omit(u_tilde[isTest])^2, lag = 30)
  # LB tests 
  LB_r  <- Box.test(na.omit(z_tilde[isTest]), lag = 30)
  LB2_r <- Box.test(na.omit(z_tilde[isTest])^2, lag = 30)
  
  format_bp_test <- function(bp_test){
    paste0(format(bp_test$statistic, digits = 2, scientific = FALSE), format_pval(bp_test$p.value)) 
  }
  format_ks_test <- function(ks_test){
    paste0(format(ks_test$KS, digits = 2, scientific = FALSE), format_pval(ks_test$p.value)) 
  }
  
  dplyr::tibble(
    place = model_Rt$model$spec$place, 
    Train = paste0("2005-", nyear_test-1),
    Test = nyear_test,
    LB_r = LB_r %>% format_bp_test,
    LB2_r = LB2_r %>% format_bp_test,
    LB_u = LB_u %>% format_bp_test,
    LB2_u = LB2_u %>% format_bp_test,
    # KS-test grades Y
    ks_Y = ks_Y %>% format_ks_test,
    # KS-test grades R
    ks_R = ks_R %>% format_ks_test,
    # Log-likelihoods
    loglik.train = loglik.train, 
    loglik.test = loglik.test)
}
# ******************************************************************************
#                               Generate data
# ******************************************************************************
# Diagnostic tests
radiation_models_CTMC_diagnostic <- purrr::map(radiation_models_CTMC, radiation_model_CTMC_diagnostic)
# ******************************************************************************
#                             Save the data 
# ******************************************************************************
if (save_data) {
  # Output directory 
  dir_output <- file.path(outputs$dir$data$models$radiation$P, place)
  # ******************************************
  # Save radiation models 
  save_new_file(dir_output = dir_output, file.name = common_name[1], 
                file.format = "RData", quiet = FALSE, radiation_models_CTMC_diagnostic)
}
