#' Preprocess Radiation Data to Fit a Monthly DTMC
#'
#' Build standardized residuals, seasonal quantities, weights, and initial DTMC
#' parameter containers for the monthly two-state radiation DTMC.
#'
#' @param Y Numeric vector. Transformed radiation observations.
#' @param dates Date vector aligned with `Y`.
#' @param weights Optional numeric/logical training weights.
#' @param params List containing initial DTMC parameters (`pi`, `Pm`, `a`,
#'   `theta`, `b`, and optionally `mu`, `sig`).
#' @param eps Numeric tolerance used to stabilize probabilities.
#' @param update_emissions Logical. If `TRUE`, update starting emission
#'   parameters from residuals.
#'
#' @return A list with preprocessed `data` and possibly completed `params`.
#' @keywords internal
dtmc_monthly_preprocess <- function(Y, dates, weights, params, eps = 1e-8, update_emissions = FALSE) {
  #------------------------------------------------------------
  # 0) Basic checks and unpack
  #------------------------------------------------------------
  stopifnot(length(Y) == length(dates))
  Tn <- length(Y)
  if (Tn < 2) stop("Need at least 2 observations for transitions.")
  
  if (is.null(params$pi))    stop("params$pi missing")
  if (is.null(params$Pm))    stop("params$Pm (list of 12 2x2) missing")
  if (is.null(params$a))     stop("params$a (seasonal mean coeffs) missing")
  if (is.null(params$theta)) stop("params$theta (OU parameter) missing")
  if (is.null(params$b))     stop("params$b (seasonal variance coeffs) missing")
  
  # small numeric epsilon for probs, rename to avoid clash
  eps_small <- eps
  
  # Month index 1..12
  month_idx <- as.integer(format(as.Date(dates), "%m"))
  # Day index 1...365
  doy <- as.integer(format(as.Date(dates), "%j"))
  # Custom weights 
  if (!missing(weights)){
    stopifnot(Tn == length(weights))
    weights <- ifelse(weights == 0, 0, 1)
  } else {
    weights <- rep(1, Tn)
  }
  #------------------------------------------------------------
  # Seasonal functions 
  seasonal_function <- function(doy, par) {
    omega <- 2 * base::pi / 365
    par[1] + par[2] * sin(omega * doy) + par[3] * cos(omega * doy)
  }
  #------------------------------------------------------------
  # 1) Seasonal mean + OU innovations (a, theta are kept fixed)
  #------------------------------------------------------------
  # Seasonal mean 
  Y_bar <- seasonal_function(doy, params$a)
  # Deseasonalized time series 
  Y_tilde <- Y - Y_bar
  # AR Residuals 
  eps_raw <- c(0, Y_tilde[2:Tn] - exp(-params$theta) * Y_tilde[1:(Tn - 1)])
  #------------------------------------------------------------
  # 2) Seasonal variance: log-variance parametrization
  #    log sigma_bar^2 = X b  -> sigma_bar = exp(0.5 * X b)
  #------------------------------------------------------------
  # Seasonal variance
  b <- seasonalModel_params_to_phi(params$b)
  J_t <- seasonal_function(doy, b)
  # Reparametrization
  reparam <- reparam_seasonal_function(b, params$theta, 2*base::pi/365)
  integral_expectation <- integral_sigma_numeric(params$theta, reparam$c_)
  # Integral for expected value 
  M_t <- purrr::map_dbl(doy, ~integral_expectation(.x-1, .x, .x)) 
  # Ratio 
  kappa_t <- 1 / (M_t / sqrt(J_t))
  # standardized residuals (for mixture updates)
  r_t <- eps_raw / sqrt(J_t)
  # Initialize parameters
  #------------------------------------------------------------
  # 4) Emission parameters (mu, sig) per-month, per-state
  #    Model: eps_t | state=i,month=m ~ N(mu_{i,m} * sigma_bar_t,
  #                                       (sig_{i,m} * sigma_bar_t)^2)
  #    or equivalently r_t ~ N(mu_{i,m}, sig_{i,m}^2)
  #------------------------------------------------------------
  have_emissions <- !is.null(params$mu) && !is.null(params$sig)
  if (!have_emissions && !update_emissions) {
    stop("Provide params$mu and params$sig or set update_emissions=TRUE to estimate them.")
  }
  
  if (!have_emissions) {
    mu  <- vector("list", 12)
    sig <- vector("list", 12)
    for (m in 1:12) {
      idx <- which(month_idx == m)
      if (length(idx) == 0) {
        mu[[m]]  <- c(mean(r_t), mean(r_t))
        sdm      <- sd(r_t)
        if (!is.finite(sdm) || sdm <= 0) sdm <- 1
        sig[[m]] <- c(sdm, sdm)
        next
      }
      mu_m  <- mean(r_t[idx], na.rm = TRUE)
      sig_m <- sd(r_t[idx],   na.rm = TRUE)
      if (!is.finite(sig_m) || sig_m <= 0) sig_m <- 1
      mu[[m]]  <- c(mu_m, mu_m)
      sig[[m]] <- c(sig_m, sig_m)
    }
    params$mu <- mu
    params$sig <- sig 
  } 
  
  structure(
    list(
      data = dplyr::tibble(date = dates, Month = month_idx, weights = weights, Y = Y, Y_bar = Y_bar, Y_tilde = Y_tilde, eps_raw = eps_raw, r_t = r_t, J_t = J_t, M_t = M_t, kappa_t = kappa_t),
      params = params
    )
  )
}

#' Run One Monthly DTMC EM Update
#'
#' Perform one forward-backward E-step and one optional M-step for monthly
#' transition matrices and state-dependent emissions.
#'
#' @param data Preprocessed data returned by `preprocess_HMM_monthly()`.
#' @param params List of DTMC parameters.
#' @param update_emissions Logical. If `TRUE`, update monthly emission means
#'   and standard deviations.
#' @param update_transitions Logical. If `TRUE`, update monthly transition
#'   matrices.
#' @param eps Numeric tolerance used to stabilize probabilities.
#'
#' @return A list containing updated parameters, filtered/smoothed
#'   probabilities, transition posteriors, and log-likelihood.
#' @keywords internal
dtmc_monthly_EM <- function(data, params,
                           update_emissions = TRUE,
                           update_transitions  = TRUE,
                           eps = 1e-8) {
  Tn <- nrow(data)
  month_idx <- data$Month
  dates     <- data$date
  weights   <- data$weights
  eps_raw   <- data$eps_raw
  J_t       <- data$J_t
  M_t       <- data$M_t
  r_t       <- data$r_t
  kappa_t   <- data$kappa_t
  # small numeric epsilon for probs, rename to avoid clash
  eps_small <- eps
  # Helper 
  logSumExp <- function(x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
  }
  #------------------------------------------------------------
  # 3) DTMC structure: pi and Pm
  #------------------------------------------------------------
  # Initial probabilities
  pi <- as.numeric(params$pi)
  if (length(pi) != 2) stop("params$pi must have length 2")
  pi <- pmax(pi, eps_small)
  pi <- pi / sum(pi)
  # Transition matrices (list of 12 2x2)
  Pm <- params$Pm
  if (length(Pm) != 12) stop("params$Pm must be a list of length 12")
  for (m in 1:12) {
    if (!is.matrix(Pm[[m]]) || any(dim(Pm[[m]]) != c(2, 2))) {
      stop("Each Pm[[m]] must be a 2x2 matrix")
    }
    row_sums <- rowSums(Pm[[m]])
    if (any(row_sums <= 0)) stop("Invalid row in Pm[[m]]")
    Pm[[m]] <- pmax(Pm[[m]], eps_small)
    Pm[[m]] <- Pm[[m]] / row_sums
  }
  #------------------------------------------------------------
  # 5) Emission log-densities log f_t(i) on eps_raw
  #    f(eps) = (1/sigma_bar_t) * phi(r_t; mu, sig)
  #------------------------------------------------------------
  mu <- params$mu
  # Extract monthly means 
  mu_1 <- purrr::map_dbl(1:Tn, ~mu[[month_idx[.x]]][1]) # * M_t
  mu_2 <- purrr::map_dbl(1:Tn, ~mu[[month_idx[.x]]][2]) #* M_t
  sig <- params$sig
  # Extract monthly std. deviations 
  sd_1 <- purrr::map_dbl(1:Tn, ~sig[[month_idx[.x]]][1]) #* sqrt(J_t)
  sd_2 <- purrr::map_dbl(1:Tn, ~sig[[month_idx[.x]]][2])# * sqrt(J_t)
  # Emission log-densities
  logf <- cbind(log(dnorm((r_t - mu_1) / sd_1) / sd_1), 
                log(dnorm((r_t - mu_2) / sd_2) / sd_2))
  #------------------------------------------------------------
  # 6) E-step: forward-backward
  #------------------------------------------------------------
  logalpha <- matrix(NA_real_, Tn, 2)
  logbeta  <- matrix(NA_real_, Tn, 2)
  # init forward
  logalpha[1, ] <- log(pi + 1e-12) + logf[1, ]
  # forward recursion
  for (t in 2:Tn) {
    m_prev <- month_idx[t - 1]
    for (j in 1:2) {
      trans_log <- logalpha[t - 1, ] + log(pmax(Pm[[m_prev]][, j], eps_small))
      logalpha[t, j] <- logf[t, j] + logSumExp(trans_log)
    }
  }
  # init backward
  logbeta[Tn, ] <- 0
  # backward recursion
  for (t in (Tn - 1):1) {
    m_curr <- month_idx[t]
    for (i in 1:2) {
      trans_log <- log(pmax(Pm[[m_curr]][i, ], eps_small)) + logf[t + 1, ] + logbeta[t + 1, ]
      logbeta[t, i] <- logSumExp(trans_log)
    }
  }
  # Total log-likelihood 
  loglik <- logSumExp(logalpha[Tn, ])
  
  # Posterior probabilities: gamma_t(i)
  gamma <- matrix(NA_real_, Tn, 2)
  for (t in 1:Tn) {
    z <- logalpha[t, ] + logbeta[t, ]
    gamma[t, ] <- exp(z - logSumExp(z))
  }
  # Joint probabilities in row-stochastic code order: xi_t(previous, next).
  xi <- array(NA_real_, dim = c(Tn - 1, 2, 2))
  for (t in 1:(Tn - 1)) {
    m_curr <- month_idx[t]
    Z <- matrix(NA_real_, 2, 2)
    for (i in 1:2) {
      for (j in 1:2) {
        Z[i, j] <- logalpha[t, i] + log(pmax(Pm[[m_curr]][i, j], eps_small)) + logf[t + 1, j] + logbeta[t + 1, j]
      }
    }
    z_norm <- logSumExp(as.vector(Z))
    xi[t, , ] <- exp(Z - z_norm)
  }
  #------------------------------------------------------------
  # 7) M-step
  #------------------------------------------------------------
  # (1) Initial probabilities
  pi_new <- gamma[1, ]
  pi_new <- pmax(pi_new, eps_small)
  pi_new <- pi_new / sum(pi_new)
  # (2) Transitions
  Pm_new <- Pm
  if (update_transitions) {
    for (m in 1:12) {
      idx_t <- which(month_idx[1:(Tn - 1)] == m)
      if (length(idx_t) == 0) next
      numer <- matrix(0, 2, 2)
      denom <- rep(0, 2)
      for (t in idx_t) {
        numer <- numer + xi[t, , ] * weights[t]
        denom <- denom + gamma[t, ] * weights[t]
      }
      for (i in 1:2) {
        if (denom[i] <= 0) {
          Pm_new[[m]][i, ] <- Pm[[m]][i, ]
        } else {
          row <- numer[i, ] / denom[i]
          row <- pmax(row, eps_small)
          Pm_new[[m]][i, ] <- row / sum(row)
        }
      }
      # Enforce the two-state CTMC embeddability condition.
      s_off <- Pm_new[[m]][1, 2] + Pm_new[[m]][2, 1]
      if (s_off >= 1 - eps_small) {
        scale <- (1 - eps_small) / s_off
        Pm_new[[m]][1, 2] <- Pm_new[[m]][1, 2] * scale
        Pm_new[[m]][2, 1] <- Pm_new[[m]][2, 1] * scale
        Pm_new[[m]][1, 1] <- 1 - Pm_new[[m]][1, 2]
        Pm_new[[m]][2, 2] <- 1 - Pm_new[[m]][2, 1]
      }
    }
  }
  
  # (3) Emissions: update mu, sig on r_t
  mu_new  <- mu
  sig_new <- sig
  if (update_emissions) {
    for (m in 1:12) {
      idx_m <- which(month_idx == m)
      if (length(idx_m) == 0) next
      for (i in 1:2) {
        w    <- gamma[idx_m, i]
        wsum <- sum(weights[idx_m] * w)
        if (wsum <= 0) next
        mu_i  <- sum(w * weights[idx_m] * r_t[idx_m]) / sum(weights[idx_m] * w)
        var_i <- sum(w * weights[idx_m] * (r_t[idx_m] - mu_i)^2) / wsum 
        mu_new[[m]][i]  <- mu_i
        sig_new[[m]][i] <- sqrt(max(var_i, eps_small))
      }
    }
  }
  #------------------------------------------------------------
  # 8) Filtered probs and predictive weights (as before)
  #------------------------------------------------------------
  alpha <- matrix(NA_real_, nrow(gamma), 2)
  lSE <- function(z) { m <- max(z); m + log(sum(exp(z - m))) }
  for (t in 1:nrow(alpha)) {
    z <- logalpha[t, ]
    alpha[t, ] <- exp(z - lSE(z))
  }
  w_pred <- matrix(NA_real_, nrow(alpha), 2)
  w_pred[1, ] <- pi_new / sum(pi_new)
  for (t in 2:nrow(alpha)) {
    m_prev <- as.integer(format(as.Date(dates[t - 1]), "%m"))
    w_pred[t, ] <- as.numeric(alpha[t - 1, ] %*% Pm_new[[m_prev]])
    w_pred[t, ] <- pmax(w_pred[t, ], eps_small)
    w_pred[t, ] <- w_pred[t, ] / sum(w_pred[t, ])
  }
  # Assign standard names 
  colnames(alpha) <- paste0("alpha", 1:ncol(alpha))
  colnames(gamma)  <- paste0("gamma", 1:ncol(gamma))
  colnames(w_pred) <- paste0("w_pred", 1:ncol(w_pred))

  mu_new <- purrr::map(mu_new,  ~setNames(.x, c("mu1", "mu2")))
  sig_new <- purrr::map(sig_new, ~setNames(.x, c("sd1", "sd2")))
  #------------------------------------------------------------
  # 9) Return
  # Y = standardized innovations under final b
  #------------------------------------------------------------
  list(
    data = data, 
    params = list(
      pi    = pi_new,
      Pm    = Pm_new,
      mu    = mu_new,
      sig   = sig_new,
      a     = params$a,
      theta = params$theta,
      b     = params$b
    ),
    alpha = alpha, 
    gamma  = gamma,
    w_pred = w_pred,
    xi     = xi,
    loglik = loglik
  )
}

#' Fit a Monthly Two-State DTMC
#'
#' Fit the monthly DTMC from transformed radiation observations using the
#' discrete radiation mixture model as starting values.
#'
#' @param Y Numeric vector. Transformed radiation observations.
#' @param dates Date vector aligned with `Y`.
#' @param weights Numeric vector of observation weights.
#' @param model Fitted discrete radiation model used for starting values.
#' @param p0 Numeric length-two initial state probability.
#' @param maxit Integer. Maximum EM iterations.
#' @param tol Numeric. Log-likelihood convergence tolerance.
#'
#' @return Fitted DTMC list returned by `dtmc_monthly_EM()`.
#' @keywords internal
dtmc_monthly_fit <- function(Y, dates, weights, model, p0 = c(0.5, 0.5), maxit = 1000, tol = 0.01){
  NM_model <- model$spec$mixture.model
  # Initialize the list of parameters
  params <- list()
  params[["mu"]] <- purrr::map(1:12, ~unlist(NM_model$means[.x,]))
  params[["sig"]] <- purrr::map(1:12, ~unlist(NM_model$sd[.x,]))
  # Initialization of Pm
  build_P_from_probs <- function(pi_vec, lambda0 = 0.3) {
    p1 <- pi_vec[1]
    p2 <- pi_vec[2]
    a  <- lambda0 * p2  # P12
    b  <- lambda0 * p1  # P21
    matrix(c(1 - a, a,
             b, 1 - b),
           nrow = 2, byrow = TRUE)
  }
  probs <- purrr::map(1:12, ~unlist(NM_model$p[.x,]))
  params[["Pm"]] <- purrr::map(1:12, ~build_P_from_probs(probs[[.x]], lambda0 = 0.3))
  
  # Initial guess 
  params[["pi"]] <- p0
  params$a <- model$spec$seasonal.mean$coefficients
  params$theta <- -log(model$spec$mean.model$phi)
  params$b <- seasonalModel_params_to_zeta(model$spec$seasonal.variance$coefficients)
  # Prepare data
  data <- dtmc_monthly_preprocess(Y, dates, weights, params, eps = 1e-8)
  
  # EM-recursion 
  loop_condition <- TRUE
  iter <- 0
  loglik <- -1000000
  while(loop_condition){
    # EM-step 
    em <- dtmc_monthly_EM(data$data, params, 
                          update_emissions = TRUE, 
                          update_transitions = TRUE,
                          eps = 1e-8)
    print(em$loglik) 
    # Update iteration 
    iter <- iter + 1
    # Check break condition 
    loop_condition <- (iter <= maxit) & abs(loglik - em$loglik) > tol
    # Update parameters 
    if (!loop_condition) break
    # Update parameters
    params <- em$params
    # Update log-likelihod 
    loglik <- em$loglik
  }
  return(em)
}
