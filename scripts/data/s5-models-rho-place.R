# ---
#' @description
#' Estimate monthly correlations between electricity and solar-radiation residuals for one location.
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse
#'
#' @arguments
#'   - param[1] (place): reference location, e.g. "Bologna".
#'   - param[2] (save_data): write generated data ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/data/s5-models-rho-place.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/data/s5-models-rho-place.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/electricity/P/electricity_models.RData` (electricity models under P)
#'   - `data/models/radiation/P/{place}/radiation_models.RData` (IID radiation models)
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData` (CTMC radiation models)
#'
#' @outputs
#'   - `data/models/rho/{place}/rho_GM.RData` (IID-regime correlations)
#'   - `data/models/rho/{place}/rho_MC.RData` (CTMC-regime correlations)
#'
#' @depends
#'   - `scripts/data/s3-models-electricity-P.R`
#'   - `scripts/data/s2b-models-radiation-P-CTMC-place.R`
#'
#' @tags
#'   - data
#'   - main
#'   - appendix
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "data", "s5-models-rho-place.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Bologna"
# Save the output
save_data <- FALSE
# ***************************** Fixed Arguments ********************************
# Reference names
common_name <- c("rho_GM", "rho_MC")
# Years of training
nyears <- outputs$nyears
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
# Load radiation models (IID)
load_data(file.path(outputs$dir$data$models$radiation$P, place), "radiation_models")
# Load radiation models (MC)
load_data(file.path(outputs$dir$data$models$radiation$P, place), "radiation_models_CTMC")
# Load electricity models
load_data(outputs$dir$data$models$electricity$P, "electricity_models")
# ******************************************************************************
#                                 Function 
# ******************************************************************************
# Fit monthly correlations
fit_monthly_correlations <- function(model_Et, model_Rt){
  # Electricity model data
  data_Et <- dplyr::select(model_Et$data, date, Et = "PUN", weights, z_E, z_tilde_E)[-1,]
  # Radiation model data
  data_Rt <- dplyr::select(model_Rt$model$data, Rt = "GHI", date, Year, Month, u_tilde, isTrain)[-1,]
  # Joint dataset with for correlation
  data_joint <- dplyr::inner_join(data_Et, data_Rt, by = "date")
  data_joint <- dplyr::filter(data_joint, isTrain & date %in% model_Et$date_train & isTrain)
  # Extract maximum year for train data 
  max_year <- max(data_joint$Year)
  # Compute the monthly correlation with flexible probabilities
  rho_model <- list()
  moments   <- list()
  tests     <- list()
  
  nmonth <- 1
  for(nmonth in 1:12){
    # Filter till the specific month
    data_m <- dplyr::filter(data_joint, Month == nmonth)
    # Target residuals 
    eps_R <-  data_m$u_tilde
    eps_E <- data_m$z_tilde_E
    # Responsabilities 
    if (class(model_Rt)[1] == "radiationModel") {
      gamma_t <- model_Rt$model$spec$mixture.model$model[[nmonth]]$E_step(data_m$u_tilde)
    } else if (class(model_Rt)[1] == "radiationModel_CTMC"){
      gamma_t <- model_Rt$CTMC$alpha[which(model_Rt$CTMC$data$date %in% data_m$date),]
    }
    
    mu_E_k  <- c(mu_E_1 = 0, mu_E_0 = 0)
    mu_R_k  <- c(mu_R_1 = 0, mu_R_0 = 0)
    v_E_k   <- c(v_E_1 = 1, v_E_0 = 1)
    v_R_k   <- c(v_R_1 = 1, v_R_0 = 1)
    cv_RE_k <- c(cv_RE_1 = 0, cv_RE_0 = 0)
    cr_RE_k <- c(cr_RE_1 = 0, cr_RE_0 = 0)
    k <- 1
    for(k in 1:2){
      # Weights
      w_k <- unlist(gamma_t[,k])
      # Normalizing factor
      W_k <- sum(w_k)
      # State-conditional means
      mu_E_k[k] <- sum(eps_E * w_k) / W_k
      mu_R_k[k] <- sum(eps_R * w_k) / W_k
      # State-conditional variances
      v_E_k[k] <- sum((eps_E - mu_E_k[k])^2 * w_k) / W_k
      v_R_k[k] <- sum((eps_R - mu_R_k[k])^2 * w_k) / W_k
      # State-conditional covariance 
      cv_RE_k[k] <- sum((eps_E - mu_E_k[k]) * (eps_R - mu_R_k[k]) * w_k) / W_k
      # State-conditional correlation
      cr_RE_k[k] <- cv_RE_k[k] / sqrt(v_E_k[k] * v_R_k[k])
    }
    
    T_k <- apply(gamma_t, 2, sum)
    # Unconditional probabilities 
    pi_k <- T_k / nrow(gamma_t)
    # Unconditional means 
    mu_E <- sum(pi_k * mu_E_k)
    mu_R <- sum(pi_k * mu_R_k)
    # Unconditional variances
    v_E <- sum(pi_k * (v_E_k + mu_E_k^2)) - mu_E^2
    v_R <- sum(pi_k * (v_R_k + mu_R_k^2)) - mu_R^2
    # Unconditional covariance
    cv_RE <- sum(pi_k  * (cv_RE_k + mu_E_k * mu_R_k)) - mu_E * mu_R
    # Unconditional correlation (Fitted rho)
    cr_RE <- cv_RE / sqrt(v_E * v_R)
    # Empirical rho
    cr.test <- cor.test(eps_E, eps_R)
    
    cr.test_k <- list(statistic = 0, p.value = 0, T_k = T_k)
    cr.test_k$statistic <- cr_RE_k * sqrt((T_k - 2) / (1 - cr_RE_k^2))
    cr.test_k$p.value <- (1 - pt(abs(cr.test_k$statistic), df = T_k - 2)) / 2
    # Store moments 
    moments[[nmonth]] <- dplyr::bind_cols(Month = nmonth, 
                                          dplyr::bind_rows(mu_E_k), mu_E = mu_E, dplyr::bind_rows(v_E_k), v_E = v_E,
                                          dplyr::bind_rows(mu_R_k), mu_R = mu_R, dplyr::bind_rows(v_R_k), v_R = v_R,
                                          dplyr::bind_rows(cv_RE_k), cv_RE = cv_RE, dplyr::bind_rows(cr_RE_k), cr_RE = cr_RE)
    # Store estimates 
    tests[[nmonth]] <- dplyr::bind_cols(Month = nmonth, 
                                        rho_hat = cr_RE,
                                        rho.emp = cr.test$estimate, 
                                        rho.stat = cr.test$statistic,
                                        rho.emp_p.value = cr.test$p.value,
                                        rho_1 = cr_RE_k[1],
                                        N1_eff = cr.test_k$T_k[1],
                                        rho_1.stat = cr.test_k$statistic[1],
                                        rho_1.p.value = cr.test_k$p.value[1], 
                                        rho_2 = cr_RE_k[2],
                                        N2_eff = cr.test_k$T_k[2],
                                        rho_2.stat = cr.test_k$statistic[2],
                                        rho_2.p.value = cr.test_k$p.value[2])
    
    # Correlations
    rho_model[[nmonth]] <- list(c(rho1 = cr_RE_k[[1]], rho2 = cr_RE_k[[2]]))
  }
  structure(
    list(
      rho_model = rho_model,
      moments = dplyr::bind_rows(moments), 
      tests = dplyr::bind_rows(tests)
    )
  )
}
# ******************************************************************************
#                               Generate data 
# ******************************************************************************
# Initialize the list for the models
rho_GM  <- list()
rho_MC <- list()
# Compute the correlations 
i <- 1
for(nyear in nyears){
  nyear <- as.character(nyear)
  max_date <- paste0(nyears[i], "-12-31")
  cli::cli_alert_info(paste0("Fitting correlation for training year: ", nyear, " (Max date: ", max_date, ") \n"))
  # Electricity model 
  model_Et <- electricity_models[[nyear]]
  # Radiation model (IID Bernoullis)
  model_Rt <- radiation_models[[nyear]]
  rho_GM[[nyear]] <- fit_monthly_correlations(model_Et, model_Rt)
  # Radiation model (Markov Chain)
  model_Rt <- radiation_models_CTMC[[nyear]]
  rho_MC[[nyear]] <- fit_monthly_correlations(model_Et, model_Rt)
}
# ******************************************************************************
#                             Save the data 
# ******************************************************************************
if (save_data) {
  # Initialize a new folder for the specific output
  dir_output <- file.path(outputs$dir$data$models$rho, place)
  make_new_directory(dir_output)
  # ******************************************
  # Save correlations 
  save_new_file(dir_output = dir_output, file.name = common_name[1], 
                file.format = "RData", quiet = FALSE, rho_GM)
  # ******************************************
  # Save correlations
  save_new_file(dir_output = dir_output, file.name = common_name[2], 
                file.format = "RData", quiet = FALSE, rho_MC)
  
}
