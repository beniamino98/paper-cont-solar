#' Scenario-Based Mean-Variance Moments for SoRad and SoREd
#'
#' Compute Monte Carlo moments, covariances, hedge ratios, and realized cashflow
#' inputs from joint electricity-radiation scenarios.
#'
#' @param preproc Scenario object returned by `scenarios_ER()`.
#' @param K_fun Optional strike function of day-of-year. If missing, the
#'   radiation seasonal mean is used.
#'
#' @return List containing updated scenarios and SoRad/SoREd daily and annual
#'   moment tables.
#' @export
solarOption_mv_scenario <- function(preproc, K_fun){
  # Matrices with simulated joint scenarios (N x scenarios)
  E_mat <- preproc$scenarios$ET_Q[-c(1),]
  F_mat <- preproc$scenarios$F_net_Q[-c(1),]
  R_mat <- preproc$scenarios$RT[-c(1),]
  ER_mat <- E_mat * R_mat
  # Simulated payoffs
  data_sim <- preproc$data_sim
  # Reference dates
  dates <- data_sim$date[-1]
  # Reference strike 
  if (!missing(K_fun)){
    data_sim$GHI_bar <- K_fun(data_sim$n) 
  } else {
    K_fun <- function(n) preproc$model_Rt$Rt_bar(n)
  }
  K_mat <- data_sim$GHI_bar[-c(1)] * matrix(1, nrow = nrow(R_mat), ncol = ncol(R_mat))
  # Realized payoff
  Gamma_mat <- solarOption_payoff(R_mat, K_mat) 
  E_Gamma_mat <- E_mat * Gamma_mat
  # Future prices for Hedging PnL
  E0 <- data_sim$PUN[1]
  FE_tT <- preproc$model_Et$F_E(data_sim$date[1], data_sim$date[-c(1)], E0)
  PnL <- data_sim$PUN[-1] - FE_tT
  # ***************************************************************
  #                 Total averages in the period 
  # ***************************************************************                   
  # mean(E[R])
  e_M_R <- mean(apply(R_mat, 2, mean))
  # mean(V[R])
  e_v_R <- mean(apply(R_mat, 2, var))
  # mean(E[Gamma])
  e_M_Gamma <- mean(apply(Gamma_mat, 2, mean))
  # mean(V[Gamma])
  e_v_Gamma <- mean(apply(Gamma_mat, 2, var))
  # mean(Cv[Gamma, R])
  e_cv_R_Gamma <- mean(apply(R_mat * Gamma_mat, 2, mean) - apply(R_mat, 2, mean) * apply(Gamma_mat, 2, mean))
  # mean(Cr[Gamma, R])
  e_cr_R_Gamma <- mean((apply(R_mat * Gamma_mat, 2, mean) - apply(R_mat, 2, mean) * apply(Gamma_mat, 2, mean)) / sqrt(apply(R_mat, 2, var) * apply(Gamma_mat, 2, var)))
  # ***************************************************************
  # mean(E[E])
  e_M_E <- mean(apply(E_mat, 2, mean))
  # mean(v[E])
  e_v_E <- mean(apply(E_mat, 2, function(x) mean(x^2))) -  mean(apply(E_mat, 2, mean))^2
  # mean(E[E R])
  e_M_ER <- mean(apply(ER_mat, 2, mean))
  # mean(V[E R])
  e_v_ER <- mean(apply(ER_mat, 2, function(x) mean(x^2))) -  mean(apply(ER_mat, 2, mean))^2
  # mean(E[E Gamma])
  e_M_E_Gamma <- mean(apply(E_Gamma_mat, 2, mean))
  # mean(V[E Gamma])
  e_v_E_Gamma <- mean(apply(E_Gamma_mat, 2, var))
  # mean(Cv[E Gamma, E R])
  e_cv_ER_EGamma <- mean(apply(R_mat * E_Gamma_mat * E_mat, 2, mean) - apply(ER_mat, 2, mean) * apply(E_Gamma_mat, 2, mean))
  # mean(Cr[E Gamma, E R])
  e_cr_ER_EGamma <- mean((apply(R_mat * E_Gamma_mat * E_mat, 2, mean) - apply(ER_mat, 2, mean) * apply(E_Gamma_mat, 2, mean)) / sqrt(apply(ER_mat, 2, var) * apply(E_Gamma_mat, 2, var)))
  # mean(Cv[E Gamma, E ])
  e_cv_EGamma_E <- mean(apply(E_Gamma_mat * E_mat, 2, mean) - apply(E_mat, 2, mean) * apply(E_Gamma_mat, 2, mean))
  # mean(Cr[E Gamma, E])
  e_cr_EGamma_E <- mean((apply(E_Gamma_mat * E_mat, 2, mean) - apply(E_mat, 2, mean) * apply(E_Gamma_mat, 2, mean)) / sqrt(apply(E_mat, 2, var) * apply(E_Gamma_mat, 2, var)))
  # mean(Cr[E R, E])
  e_cv_ER_E <- mean((apply(R_mat * E_mat^2, 2, mean) - apply(E_mat * R_mat, 2, mean) * apply(E_mat, 2, mean))) 
  e_cr_ER_E <- mean((apply(R_mat * E_mat^2, 2, mean) - apply(E_mat * R_mat, 2, mean) * apply(E_mat, 2, mean)) / sqrt(apply(E_mat * R_mat, 2, var) * apply(E_mat, 2, var)))
  # ***************************************************************
  #                 Total averages for each t 
  # ***************************************************************        
  # E[R]
  M_R <- apply(R_mat, 1, mean)
  # V[R]
  v_R <- apply(R_mat, 1, var)
  # E[Gamma]
  M_Gamma <- apply(Gamma_mat, 1, mean)
  # V[Gamma]
  v_Gamma <- apply(Gamma_mat, 1, var)
  # Cv[Gamma, R]
  S_R_Gamma <- diag(cov(t(R_mat), t(Gamma_mat)))
  # ***************************************************************
  # E[E]
  M_E <- apply(E_mat, 1, mean)
  # v[E]
  v_E <- apply(E_mat, 1, var)
  # E[E R]
  M_ER <- apply(ER_mat, 1, mean)
  # V[E R]
  v_ER <- apply(ER_mat, 1, var)
  # E[E R]
  M_E_Gamma <- apply(E_Gamma_mat, 1, mean)
  # V[E R]
  v_E_Gamma <- apply(E_Gamma_mat, 1, var)
  # Cv[E Gamma, E R]
  S_ER_EGamma <- diag(cov(t(E_mat * Gamma_mat), t(E_mat * R_mat)))
  # Cv[E R, E]
  S_ER_E <- diag(cov(t(E_mat * R_mat), t(E_mat)))
  # Cv[E Gamma, E]
  S_EGamma_E <- diag(cov(t(E_mat * Gamma_mat), t(E_mat)))
  # Cv[Gamma, E ]
  S_Gamma_E <- diag(cov(t(Gamma_mat), t(E_mat)))
  # ***************************************************************
  #                   Strip Hedging (SoREd)
  # ***************************************************************
  # Strip by strip hedging 
  beta_sored <- S_EGamma_E / v_E
  # Residual variance
  v_E_Gamma_mid_E <- v_E_Gamma - S_EGamma_E^2 / v_E
  # ***************************************************************
  #                   Strip Hedging (SoRad)
  # ***************************************************************
  # Strip by strip hedging 
  beta_sorad <- S_Gamma_E / v_E
  # Residual variance
  v_Gamma_mid_E <- v_Gamma - S_Gamma_E^2 / v_E
  # ***************************************************************
  # Realized radiation 
  Rt <- data_sim$GHI[-1]
  Et <- data_sim$PUN[-1]
  Gamma <- solarOption_payoff(data_sim$GHI, data_sim$GHI_bar)[-1]
  E_Gamma <- Gamma * Et
  # Structure data for SoRad
  moments_sorad <- tibble(date = dates,
                          Year = lubridate::year(dates),
                          Month = lubridate::month(dates),
                          Day = lubridate::day(dates),
                          Rt = Rt,
                          Gamma = Gamma,
                          M_R = M_R, 
                          v_R = v_R, 
                          M_Gamma = M_Gamma, 
                          v_Gamma = v_Gamma, 
                          S_R_Gamma = S_R_Gamma, 
                          cr_R_Gamma = S_R_Gamma / sqrt(v_R * v_Gamma),
                          v_Gamma_mid_E = v_Gamma_mid_E, 
                          beta = beta_sorad,
                          PnL = beta_sorad * PnL)
  
  # Structure data for SoRad
  moments_sored <- tibble(date = dates,
                          Year = lubridate::year(dates),
                          Month = lubridate::month(dates),
                          Day = lubridate::day(dates),
                          Rt = Rt,
                          Et = Et,
                          F_E = FE_tT,
                          Gamma = E_Gamma, 
                          # Marginal 
                          M_E = M_E,
                          v_E = v_E,
                          M_R = M_R,
                          v_R = v_R,
                          # Joint 
                          M_ER = M_ER, 
                          v_ER = v_ER, 
                          M_E_Gamma = M_E_Gamma, 
                          v_E_Gamma = v_E_Gamma, 
                          # Covariances 
                          S_ER_EGamma = S_ER_EGamma, 
                          cr_ER_EGamma = S_ER_EGamma / sqrt(v_ER * v_E_Gamma),
                          S_EGamma_E = S_EGamma_E, 
                          cr_EGamma_E = S_EGamma_E / sqrt(v_E_Gamma * v_E), 
                          S_ER_E = S_ER_E, 
                          cr_ER_E = S_ER_E / sqrt(v_ER * v_E), 
                          v_E_Gamma_mid_E = v_E_Gamma_mid_E, 
                          beta = beta_sored,
                          PnL = beta_sored * PnL)
  # Add row PnL
  moments_sored$net_PnL <- PnL
  
  moments_sorad_day <- tibble(
    Year = lubridate::year(dates[1]),
    M_R = e_M_R,
    v_R = e_v_R, 
    M_R_emp = mean(Rt), 
    v_R_emp = var(Rt),
    M_Gamma = e_M_Gamma,
    v_Gamma = e_v_Gamma,
    M_Gamma_emp = mean(Gamma),
    v_Gamma_emp = var(Gamma),
    cv_R_Gamma = e_cv_R_Gamma,
    cv_R_Gamma_emp = cov(Rt, Gamma),
    cr_R_Gamma = e_cr_R_Gamma,
    cr_R_Gamma_emp = cor(Rt, Gamma)
  )
  
  moments_sored_day <- tibble(
    Year = lubridate::year(dates[1]),
    M_ER = e_M_ER,
    v_ER = e_v_ER, 
    M_ER_emp = mean(Et*Rt), 
    v_ER_emp = var(Et*Rt),
    M_E_Gamma = e_M_E_Gamma,
    v_E_Gamma = e_v_E_Gamma,
    M_E_Gamma_emp = mean(E_Gamma),
    v_E_Gamma_emp = var(E_Gamma),
    cv_ER_EGamma = e_cv_ER_EGamma,
    cv_ER_EGamma_emp = cov(Et*Rt, E_Gamma),
    cr_ER_EGamma = e_cr_ER_EGamma,
    cr_ER_EGamma_emp = cor(Et*Rt, Et*Gamma),
    cr_EGamma_E = e_cr_EGamma_E,
    cr_EGamma_E_emp = cor(E_Gamma, Et),
    cr_ER_E = e_cr_ER_E,
    cr_ER_E_emp = cor(Et*Rt, Et)
  )
  
  preproc$scenarios$Gamma <- Gamma_mat 
  preproc$scenarios$E_Gamma <- E_Gamma_mat 
  
  structure(
    list(
      preproc = preproc,
      sorad = moments_sorad, 
      sored = moments_sored,
      sorad_day = moments_sorad_day, 
      sored_day = moments_sored_day,
      K_fun = K_fun
    )
  )
}

#' Scenario-Based Mean-Variance Moments for Indexed Contracts
#'
#' Aggregate simulated daily payoffs into index-style SoRadIDX and SoREdIDX
#' moments and hedge inputs.
#'
#' @param preproc Scenario object returned by `scenarios_ER()`.
#' @param K_fun Optional strike function of day-of-year. If missing, the
#'   radiation seasonal mean is used.
#'
#' @return List containing updated scenarios and indexed SoRad/SoREd moment
#'   tables.
#' @export
solarOptionIDX_mv_scenario <- function(preproc, K_fun){
  # Matrices with simulated joint scenarios 
  E_mat <- preproc$scenarios$ET_Q[-c(1),]
  F_mat <- preproc$scenarios$F_net_Q[-c(1),]
  R_mat <- preproc$scenarios$RT[-c(1),]
  ER_mat <- E_mat * R_mat
  # Simulated payoffs
  data_sim <- preproc$data_sim
  # Reference dates
  dates <- data_sim$date[-1]
  # Reference strike 
  if (!missing(K_fun)){
    data_sim$GHI_bar <- K_fun(data_sim$n) 
  } else {
    K_fun <- function(n) preproc$model_Rt$Rt_bar(n)
  }
  K_mat <- data_sim$GHI_bar[-c(1)] * matrix(1, nrow = nrow(R_mat), ncol = ncol(R_mat))
  # Realized payoffs 
  Gamma_mat <- solarOption_payoff(R_mat, K_mat) 
  E_Gamma_mat <- E_mat * Gamma_mat
  # Future prices for Hedging PnL
  E0 <- data_sim$PUN[1]
  FE_tT <- preproc$model_Et$F_E(data_sim$date[1], data_sim$date[-c(1)], E0)
  PnL <- data_sim$PUN[-1] - FE_tT
  # ***************************************************************
  # Scenarios of Gamma
  soradidx_scenarios <- apply(Gamma_mat, 2, sum)
  # E[sum R]
  M_R <- mean(apply(R_mat, 2, sum))
  # V[sum R]
  v_R <- var(apply(R_mat, 2, sum))
  # E[sum Gamma]
  M_Gamma <- mean(soradidx_scenarios)
  # V[sum Gamma]
  v_Gamma <- var(soradidx_scenarios)
  # Cv[sum Gamma, sum RT]
  S_R_Gamma <- sum(cov(t(R_mat), t(Gamma_mat))) 
  # ***************************************************************
  # Scenarios of Gamma
  soredidx_scenarios <- apply(E_Gamma_mat, 2, sum)
  ER_scenarios <- apply(ER_mat, 2, sum) 
  # E[sum R-]
  M_E <- mean(apply(E_mat, 2, sum))
  # v[sum R]
  v_E <- var(apply(E_mat, 2, sum))
  # E[sum ET RT]
  M_ER <- mean(apply(ER_mat, 2, sum))
  # V[sum ET RT]
  v_ER <- var(apply(ER_mat, 2, sum))
  # E[sum ET RT]
  M_E_Gamma <- mean(soredidx_scenarios)
  # V[sum ET RT]
  v_E_Gamma <- var(soredidx_scenarios)
  # Cv[sum ET Gamma, sum ET RT]
  S_ER_EGamma <- sum(cov(t(E_Gamma_mat), t(ER_mat)))
  # Cv[sum ET Gamma, sum RT]
  S_ER_E <- cov(t(ER_mat), t(E_mat))
  # Cv[sum ET Gamma, sum ET]
  S_EGamma_E <- cov(t(E_mat), t(E_Gamma_mat)) 
  # Cv[sum ET Gamma, sum ET]
  S_E_Gamma <- cov(t(E_mat), t(Gamma_mat)) 
  # ***************************************************************
  #                   Global Hedging (SoREdIDX)
  # ***************************************************************
  # Covariance between futures strip 
  S_FF <- cov(t(F_mat))
  S_FI <- cov(t(F_mat), cbind(soredidx_scenarios))
  S_FER <- cov(t(F_mat), cbind(ER_scenarios))
  # Weights in futures 
  beta_total <- c(solve(S_FF) %*% S_FI)
  alpha_total <- c(solve(S_FF) %*% S_FER)
  # Residuals variance
  v_E_Gamma_mid_E_total <- drop(v_E_Gamma - 2 * t(beta_total) %*% S_FI + t(beta_total) %*% S_FF %*% beta_total)
  # Equivalently: valid only in this optimal case
  # v_E_Gamma - beta_total %*% S_FI
  # Residual covariance
  S_ER_E_mid_E_total <- drop(S_ER_EGamma
                             - t(alpha_total) %*% S_FI
                             - t(beta_total)  %*% S_FER
                             + t(alpha_total) %*% S_FF %*% beta_total)
  #S_ER_E_mid_E_total <- drop(S_ER_EGamma - t(alpha_total) %*% S_FI)
  # ***************************************************************
  #                   Strip Hedging (SoREdIDX)
  # ***************************************************************
  beta_strip <- diag(S_EGamma_E) / diag(S_FF)
  alpha_strip <- diag(S_ER_E) / diag(S_FF)
  # Residuals variance
  v_E_Gamma_mid_E_strip <- drop(v_E_Gamma - 2 * t(beta_strip) %*% S_FI + t(beta_strip) %*% S_FF %*% beta_strip)
  # Residual covariance
  S_ER_E_mid_E_strip <- drop(S_ER_EGamma
                             - t(alpha_strip) %*% S_FI
                             - t(beta_strip)  %*% S_FER
                             + t(alpha_strip) %*% S_FF %*% beta_strip)
  
  # ***************************************************************
  # Realized radiation 
  Rt <- data_sim$GHI[-1]
  Et <- data_sim$PUN[-1]
  Gamma <- solarOption_payoff(data_sim$GHI, data_sim$GHI_bar)[-1]
  E_Gamma <- Gamma * Et
  # Structure data for SoRadIDX
  moments_sorad <- tibble(t_hor = max(data_sim$date), 
                          t_now = data_sim$date[1],
                          tau = as.numeric(difftime(t_hor, t_now, units = "days")),
                          Et_Rt = sum(Et*Rt),
                          Rt = sum(Rt),
                          Et = sum(Et),
                          Gamma = sum(Gamma),
                          M_R = M_R, 
                          v_R = v_R, 
                          M_Gamma = M_Gamma, 
                          v_Gamma = v_Gamma, 
                          S_R_Gamma = S_R_Gamma, 
                          cr_R_Gamma = S_R_Gamma / sqrt(v_R * v_Gamma))
  
  # Structure data for SoREdIDX
  moments_sored <- tibble(t_hor = max(data_sim$date), 
                          t_now = data_sim$date[1],
                          tau = as.numeric(difftime(t_hor, t_now, units = "days")),
                          Et_Rt = sum(Et*Rt),
                          Rt = sum(Rt),
                          Et = sum(Et),
                          Gamma = sum(E_Gamma),
                          M_ER = M_ER, 
                          v_ER = v_ER, 
                          M_E = M_E,
                          M_R = M_R,
                          v_E = v_E, 
                          v_R = v_R,
                          M_E_Gamma = M_E_Gamma, 
                          v_E_Gamma = v_E_Gamma, 
                          S_ER_EGamma = S_ER_EGamma, 
                          cr_ER_EGamma = S_ER_EGamma / sqrt(v_ER * v_E_Gamma),
                          S_EGamma_E = sum(S_EGamma_E), 
                          S_ER_E = sum(S_ER_E),
                          v_E_Gamma_mid_E_total = v_E_Gamma_mid_E_total, 
                          v_E_Gamma_mid_E_strip = v_E_Gamma_mid_E_strip, 
                          S_ER_E_mid_E_total = S_ER_E_mid_E_total,
                          S_ER_E_mid_E_strip = S_ER_E_mid_E_strip, 
                          pnl = list(PnL),
                          alpha_total = list(alpha_total),
                          beta_total = list(beta_total),
                          alpha_strip = list(alpha_strip),
                          beta_strip = list(beta_strip))
  
  preproc$scenarios$Gamma <- Gamma_mat 
  preproc$scenarios$E_Gamma <- E_Gamma_mat 
  
  structure(
    list(
      preproc = preproc,
      sorad = moments_sorad, 
      sored = moments_sored,
      K_fun = K_fun
    )
  )
}

