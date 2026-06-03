#' Model-Based Mean-Variance Moments for SoRad
#'
#' Compute analytical moment tables for SoRad contracts using the radiation model 
#'
#' @param model_Rt Radiation HMM model object.
#' @param t_now Date or character scalar. Conditioning date.
#' @param t_hor Date or character scalar. Horizon date.
#' @param seq_date Logical. If `TRUE`, compute moments for every daily maturity
#'   between `t_now + 1` and `t_hor`; otherwise compute only `t_hor`.
#' @param K_fun Strike function
#' @return A table with the moments.
#' @export
solarOption_moments_model_sorad <- function(model_Rt, t_now, t_hor, seq_date = TRUE, K_fun = NULL){
  # Transform parameters
  alpha <- model_Rt$model$spec$transform$alpha
  beta  <- model_Rt$model$spec$transform$beta
  # Generate the approximated density of Y
  generate_pdf_Y <- function(mom){
    M_Y <- c(mom$M_Y1, mom$M_Y0)
    S_Y <- c(mom$S_Y1, mom$S_Y0)
    p_T <- c(mom$p1, 1-mom$p1)
    force(M_Y); force(S_Y); force(p_T)
    # Pdf of Yt  
    function(x){
      dmixnorm(x, M_Y, S_Y, p_T)
    }
  }
  # Compute sorad moments all in once 
  sorad_moments <- function(x, pdf_Y, RT, Gamma){
    # Pre-compute quantites
    pdf_x <- pdf_Y(x)
    RT_x <- RT(x)
    Gamma_x <- Gamma(x)
    # *************************
    # E[RT]
    M_R <- RT_x * pdf_x
    # E[RT2]
    M2_R <- RT_x^2 * pdf_x
    # E[Gamma]
    M_Gamma <- Gamma_x * pdf_x
    # E[Gamma^2]
    M2_Gamma <- Gamma_x^2 * pdf_x
    # E[RT Gamma] 
    M_R_Gamma <- RT_x * Gamma_x * pdf_x
    # Output 
    rbind(M_R, M2_R, M_Gamma, M2_Gamma, M_R_Gamma)
  }
  # Generate moments computator
  generate_eval_moment <- function(maxEval = 100000, tol = 0.00001, vectorInterface = TRUE, lowerLimit = c(-Inf), upperLimit = c(Inf)){
    force(maxEval); force(tol); force(vectorInterface);
    force(lowerLimit); force(upperLimit);
    function(.f, fDim = 5) {
      intg <- cubature::hcubature(
        function(x) .f(x), 
        lowerLimit = lowerLimit,
        upperLimit = upperLimit,
        fDim = fDim,
        maxEval = maxEval,
        tol = tol, 
        vectorInterface = vectorInterface
      )
      out <- intg$integral
      attr(out, "error") <- intg$error
      return(out)
    }
  }
  # Utility functions
  generate_RT <- function(Ct,  model_Rt){
    force(Ct); force(model_Rt);
    function(x){
      model_Rt$model$spec$transform$iRY(x, Ct)
    }
  }
  generate_Gamma <- function(GHI_bar, Ct, model_Rt){
    force(Ct); force(model_Rt);
    function(x){
      RT_x <- model_Rt$model$spec$transform$iRY(x, Ct)
      solarOption_payoff(RT_x, GHI_bar)
    }
  }
  # Strike price function 
  if (is.null(K_fun)) {
    K_fun <- function(t) model_Rt$Rt_bar(t)
  }
  # ***************************************************************************
  #                             Moments 
  # ***************************************************************************  
  t_now <- as.Date(t_now)
  t_hor <- t_seq <- as.Date(t_hor)
  # Sequence of dates for each maturity 
  if (seq_date) {
    t_seq <- seq.Date(t_now + 1, t_hor, 1)
  }
  # Times to maturity
  tau <- as.numeric(t_seq - t_now)
  # Compute marginal moments for Y 
  moments_Y <- radiationMoments(t_now, t_seq, model_Rt, R0 = NULL)
  moments_Y$K <- K_fun(moments_Y$date)
  # ***************************************************************************  
  # Transform function for R
  RT <- purrr::map(1:nrow(moments_Y), ~generate_RT(moments_Y$Ct[.x], model_Rt))
  # Transform function for Gamma
  Gamma <- purrr::map(1:nrow(moments_Y), ~generate_Gamma(moments_Y$K[.x], moments_Y$Ct[.x], model_Rt))
  # Approximated densities
  pdfs_Y <- purrr::map(1:length(t_seq), ~generate_pdf_Y(moments_Y[.x,]))
  # Moments 
  moments_sorad <- purrr::map(1:length(t_seq), ~generate_eval_moment(lowerLimit = c(-Inf), upperLimit = c(Inf))(function(x) sorad_moments(x, pdf_Y = pdfs_Y[[.x]], RT = RT[[.x]], Gamma = Gamma[[.x]]), fDim = 5))
  # ***************************************************************
  # E[RT]
  M_R <- purrr::map_dbl(moments_sorad, ~.x[[1]])
  # V[RT]
  v_R <- purrr::map_dbl(moments_sorad, ~.x[[2]]) - M_R^2
  # E[Gamma]
  M_Gamma <- purrr::map_dbl(moments_sorad, ~.x[[3]])
  # V[Gamma]
  v_Gamma <- purrr::map_dbl(moments_sorad, ~.x[[4]]) - M_Gamma^2
  # cov[ET Gamma] = E[RT Gamma] - E[RT] * E[Gamma]
  S_R_Gamma <- purrr::map_dbl(moments_sorad, ~.x[[5]]) - M_Gamma * M_R
  # ***************************************************************
  # Structure data for SoRad
  R_tT <- dplyr::filter(model_Rt$model$data, date %in% moments_Y$date)$GHI

  # Output data
  dplyr::tibble(date = moments_Y$date,
                M_R = M_R, 
                v_R = v_R, 
                K = moments_Y$K,
                M_Gamma = M_Gamma, 
                v_Gamma = v_Gamma, 
                S_R_Gamma = S_R_Gamma, 
                cr_R_Gamma = S_R_Gamma / sqrt(v_R * v_Gamma),
                Gamma = solarOption_payoff(R_tT, moments_Y$K))
}

#' Model-Based Mean-Variance Moments for SoREd
#'
#' Compute analytical moment tables for SoREd contracts using the
#' radiation and electricity marginal models plus CTMC regime correlations.
#'
#' @param model_Et Electricity model object.
#' @param model_Rt Radiation HMM model object.
#' @param rho List with 12 monthly state-correlation vectors.
#' @param t_now Date or character scalar. Conditioning date.
#' @param t_hor Date or character scalar. Horizon date.
#' @param seq_date Logical. If `TRUE`, compute moments for every daily maturity
#'   between `t_now + 1` and `t_hor`; otherwise compute only `t_hor`.
#' @param K_fun Strike function
#' @return A table with the moments.
#' @export
solarOption_moments_model_sored <- function(model_Et, model_Rt, rho, t_now, t_hor, seq_date = TRUE, K_fun = NULL){
  # ***************************************************************************
  # Generate Joint pdf 
  generate_joint_pdf <- function(mom_Y, mom_X, p_T, S_XY){
    # Marginal means 
    M_XY_1 <- c(mom_X$M_X, mom_Y$M_Y1)
    M_XY_0 <- c(mom_X$M_X, mom_Y$M_Y0)
    # Covariance matrices 
    S_XY_1 <- matrix(c(mom_X$S2_X, S_XY[1], S_XY[1], mom_Y$S_Y1^2), nrow = 2, ncol = 2, byrow = TRUE)
    S_XY_0 <- matrix(c(mom_X$S2_X, S_XY[2], S_XY[2], mom_Y$S_Y0^2), nrow = 2, ncol = 2, byrow = TRUE)
    p_1 <- p_T[1]
    p_0 <- p_T[2]
    function(x){
      p_1 * mvtnorm::dmvnorm(x, mean = M_XY_1, sigma = S_XY_1) +  p_0 * mvtnorm::dmvnorm(x, mean =  M_XY_0, sigma = S_XY_0)
    }
  }
  # Generate the approximated density of Y
  generate_pdf_Y <- function(mom){
    M_Y <- c(mom$M_Y1, mom$M_Y0)
    S_Y <- c(mom$S_Y1, mom$S_Y0)
    p_T <- c(mom$p1, 1-mom$p1)
    force(M_Y); force(S_Y); force(p_T)
    # Pdf of Yt  
    function(x){
      dmixnorm(x, M_Y, S_Y, p_T)
    }
  }
  # Compute sored moments all in once 
  sored_moments <- function(x, pdf_joint, RT, Gamma){
    # Pre-compute quantites
    pdf_x <- pdf_joint(t(x))
    RT_x <- RT(x[2,])
    ET_x <- exp(x[1,])
    Gamma_x <- Gamma(x[2,])
    # *************************
    # E[ET]
    M_E <- ET_x * pdf_x
    # E[RT2]
    M2_E <- ET_x^2 * pdf_x
    # Joint moment
    # E[ET RT]
    M_ER <- ET_x * RT_x * pdf_x
    # E[RT^2 Gamma^2]
    M2_ER <- ET_x^2 * RT_x^2 * pdf_x
    # E[ET Gamma]
    M_E_Gamma <- ET_x * Gamma_x * pdf_x
    # Variances 
    # E[ET^2 Gamma^2]
    M2_E_Gamma <- ET_x^2 * Gamma_x^2 * pdf_x
    # E[ET RT, ET Gamma]
    M_ER_E_Gamma <- ET_x^2 * RT_x  * Gamma_x * pdf_x
    # E[ET Gamma, ET]
    M_EGamma_E <- ET_x^2 * Gamma_x * pdf_x
    # E[E^2 R]
    M_ER_E <- ET_x^2 * RT_x * pdf_x
    # E[Gamma E]
    M_Gamma_E <- Gamma_x * ET_x * pdf_x
    # E[Gamma]
    M_Gamma <- Gamma_x * pdf_x
    # E[Gamma^2]
    M2_Gamma <- Gamma_x^2 * pdf_x
    # E[Gamma]
    M_R <- RT_x * pdf_x
    # E[Gamma^2]
    M2_R <- RT_x^2 * pdf_x
    rbind(M_E, M2_E, M_ER, M2_ER, M_E_Gamma, M2_E_Gamma, M_ER_E_Gamma, M_EGamma_E, M_ER_E, M_Gamma_E,
          M_Gamma, M2_Gamma, M_R, M2_R)
  }
  # Generate moments computator
  generate_eval_moment <- function(maxEval = 100000, tol = 0.00001, vectorInterface = TRUE, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf)){
    force(maxEval); force(tol); force(vectorInterface);
    force(lowerLimit); force(upperLimit);
    function(.f, fDim = 1) {
      intg <- cubature::hcubature(
        function(x) .f(x), 
        lowerLimit = lowerLimit,
        upperLimit = upperLimit,
        fDim = fDim,
        maxEval = maxEval,
        tol = tol, 
        vectorInterface = vectorInterface
      )
      out <- intg$integral
      attr(out, "error") <- intg$error
      return(out)
    }
  }
  # Utility functions
  generate_RT <- function(Ct,  model_Rt){
    force(Ct); force(model_Rt);
    function(x){
      model_Rt$model$spec$transform$iRY(x, Ct)
    }
  }
  generate_Gamma <- function(GHI_bar, Ct, model_Rt){
    force(Ct); force(model_Rt);
    function(x){
      RT_x <- model_Rt$model$spec$transform$iRY(x, Ct)
      solarOption_payoff(RT_x, GHI_bar)
    }
  }
  # Strike price function 
  if (is.null(K_fun)) {
    K_fun <- function(t) model_Rt$Rt_bar(t)
  }
  # ***************************************************************************
  #                             Moments 
  # ***************************************************************************  
  t_now <- as.Date(t_now)
  t_hor <- t_seq <- as.Date(t_hor)
  # Sequence of dates for each maturity 
  if (seq_date) {
    t_seq <- seq.Date(t_now + 1, t_hor, 1)
  }
  # Times to maturity
  tau <- as.numeric(t_seq - t_now)
  # Solar radiation at time t_now 
  R0 <- filter(model_Rt$model$data, date == t_now)$GHI
  # Compute marginal moments for Y 
  moments_Y <- radiationMoments(t_now, t_seq, model_Rt, R0 = R0)
  moments_Y$K <- K_fun(moments_Y$date)
  RT <- purrr::map(1:nrow(moments_Y), ~generate_RT(moments_Y$Ct[.x], model_Rt))
  Gamma <- purrr::map(1:nrow(moments_Y), ~generate_Gamma(moments_Y$K[.x], moments_Y$Ct[.x], model_Rt))
  # Future prices for Hedging PnL
  E0 <- filter(model_Et$data, date == t_now)$PUN
  FE_tT <- model_Et$F_E(t_now, t_seq, E0)
  # Marginal moments for X
  mom_X <- tibble(M_X = model_Et$model$M_X(tau, E0), S2_X = model_Et$model$V_X(tau))
  # Inputs 
  CTMC <- model_Rt$CTMC
  if (!is.null(CTMC$params$Qm)){
    Qm <- CTMC$params$Qm
  } else {
    Qm <- transition_list_to_generator_2state(CTMC$params$Pm)
  }
  sd <- CTMC$params$sig
  theta <- model_Rt$theta
  kappa <- model_Et$model$kappa
  sigma_bar <- model_Rt$seasonal_variance$extra_params$seasonal_function
  sigma_X <- model_Et$model$sigma
  # probability at time t_now 
  p0 <- CTMC$alpha[CTMC$data$date == t_now, ]
  # Generate bounds for integration 
  bounds <- purrr::map(t_seq, ~create_bounds(t_now, .x, Qm))
  # Probability at maturity 
  p_T <- purrr::map(1:length(bounds), ~p0 %*% Phi_C(bounds[[.x]]$n[1], bounds[[.x]]$tau[1], bounds[[.x]])[[1]])
  # Covariance between X and Y
  S_XY <- purrr::map(1:length(bounds), ~integral_E_Yt_Xt(bounds[[.x]], p0, sd, theta, kappa, sigma_bar, sigma_X, rho, maxEval = 10000, tol = 0.00001))
  # Joint density 
  pdf_joint <- purrr::map(1:length(bounds), ~generate_joint_pdf(moments_Y[.x,], mom_X[.x,], p_T[[.x]], S_XY[[.x]]))
  # ***************************************************************
  # Moments 
  sored_mom <- purrr::map(1:length(bounds), ~generate_eval_moment()(function(x) sored_moments(x, pdf_joint = pdf_joint[[.x]], RT = RT[[.x]], Gamma = Gamma[[.x]]), fDim = 14))
  # E[ET]
  M_E <- purrr::map_dbl(sored_mom, ~.x[[1]])
  # V[ET]
  v_E <- purrr::map_dbl(sored_mom, ~.x[[2]]) - M_E^2
  # E[ET RT]
  M_ER <- purrr::map_dbl(sored_mom, ~.x[[3]])
  # V[ET RT]
  v_ER <- purrr::map_dbl(sored_mom, ~.x[[4]]) - M_ER^2
  # E[ET Gamma]
  M_E_Gamma <- purrr::map_dbl(sored_mom, ~.x[[5]])
  # V[ET Gamma]
  v_E_Gamma <- purrr::map_dbl(sored_mom, ~.x[[6]])  - M_E_Gamma^2 
  # Cv[ET RT, ET Gamma]
  S_ER_EGamma <- purrr::map_dbl(sored_mom, ~.x[[7]]) - M_ER * M_E_Gamma
  # Cv[ET Gamma, ET]
  S_EGamma_E <- purrr::map_dbl(sored_mom, ~.x[[8]]) - M_E_Gamma * M_E
  # Cv[E R, E]
  S_ER_E <- purrr::map_dbl(sored_mom, ~.x[[9]]) - M_E * M_ER
  # E[Gamma]
  M_Gamma <- purrr::map_dbl(sored_mom, ~.x[[11]])
  # Cv[Gamma, E ]
  S_Gamma_E <- purrr::map_dbl(sored_mom, ~.x[[10]]) - M_E * M_Gamma
  # V[Gamma]
  S_Gamma <- purrr::map_dbl(sored_mom, ~.x[[12]]) - M_Gamma^2
  # E[Gamma]
  M_R <- purrr::map_dbl(sored_mom, ~.x[[13]])
  # V[Gamma]
  S_R <- purrr::map_dbl(sored_mom, ~.x[[14]]) - M_R^2
  # Strip by strip hedging
  beta_strip <- S_EGamma_E / v_E
  # Residual variance 
  v_E_Gamma_mid_E <- v_E_Gamma - S_EGamma_E^2 / v_E
  # ***************************************************************
  # Realized radiation 
  R_tT <- dplyr::filter(model_Rt$model$data, date %in% moments_Y$date)$GHI
  # Realized payoff (SoRad)
  Gamma <- solarOption_payoff(R_tT, moments_Y$K)
  # Realized electricity 
  E_tT <- dplyr::filter(model_Et$data, date %in% moments_Y$date)$PUN
  # Realized payoff (SoREd)
  Gamma_E <- E_tT * Gamma
  # Futures PnL 
  PnL <- E_tT - FE_tT
  dplyr::tibble(date = t_seq,
                M_ER = M_ER, 
                v_ER = v_ER,
                M_E = M_E, 
                v_E = v_E, 
                M_E_Gamma = M_E_Gamma, 
                v_E_Gamma = v_E_Gamma, 
                S_ER_EGamma = S_ER_EGamma, 
                cr_ER_EGamma = S_ER_EGamma / sqrt(v_ER * v_E_Gamma),
                S_EGamma_E = S_EGamma_E, 
                S_ER_E = S_ER_E,
                S_Gamma_E = S_Gamma_E,
                v_E_Gamma_mid_E = v_E_Gamma_mid_E, 
                beta = beta_strip,
                PnL = beta_strip * PnL, 
                Rt = R_tT, 
                Et = E_tT, 
                FE = FE_tT,
                Gamma = Gamma_E)
}

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
solarOption_moments_scenario <- function(preproc, K_fun){
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
