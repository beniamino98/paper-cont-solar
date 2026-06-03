# install.packages(c("Rcpp"))
# library(Rcpp)
# Compile integrand (only when update the function)
#old <- getwd()
#setwd("functions/C")
#system("R CMD SHLIB bounds_kernels.c")
#setwd(old)
# dyn.load("functions/C/bounds_kernels.so")
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' Matrix Exponential for a Two-State CTMC Generator
#'
#' Compute `expm(k * Q)` after validating that `Q` is a two-state CTMC
#' generator.
#'
#' @param k Numeric scalar. Non-negative time horizon.
#' @param Q Numeric `2 x 2` CTMC generator.
#'
#' @return Numeric `2 x 2` transition matrix.
#' @export
matrix_exponential <- function(k, Q) {
  k <- as.numeric(k)
  Q <- as.matrix(Q)

  if (length(k) != 1 || !is.finite(k)) {
    stop("k must be a finite scalar horizon.")
  }
  if (k < -1e-12) {
    stop("k must be non-negative.")
  }
  if (!ctmc_is_generator_2state(Q)) {
    stop("Q must be a valid 2x2 CTMC generator.")
  }

  k <- max(k, 0)
  unname(as.matrix(expm::expm(k * Q)))
}

#' Identify the Active Monthly CTMC Interval
#'
#' Return the month associated with a CTMC time point using half-open interval
#' bounds.
#'
#' @param tau Numeric vector of time points in day-of-year coordinates.
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#'
#' @return Numeric vector of month indices.
#' @keywords internal
get_month_index_C <- function(tau, bounds) {
  vapply(as.numeric(tau), function(x) {
    idx <- which(bounds$n <= x & x < bounds$N)
    if (length(idx) == 0) {
      if (x >= max(bounds$N)) {
        return(bounds$Month[length(bounds$Month)])
      }
      if (x <= min(bounds$n)) {
        return(bounds$Month[1])
      }
      stop("tau is outside the supplied CTMC bounds.")
    }
    bounds$Month[idx[length(idx)]]
  }, numeric(1))
}

#' Transition Matrix Over One CTMC Interval
#'
#' Compute the ordered product of monthly CTMC transitions over `[a, b)`.
#'
#' @param a Numeric scalar. Start time in day-of-year coordinates.
#' @param b Numeric scalar. End time in day-of-year coordinates.
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#'
#' @return Numeric `2 x 2` transition matrix.
#' @keywords internal
ctmc_phi_one <- function(a, b, bounds) {
  if (b < a) {
    stop("Phi_C expects a <= b.")
  }
  if (isTRUE(all.equal(a, b, tolerance = 1e-14))) {
    return(diag(1, 2, 2))
  }

  idx <- which(bounds$n < b & bounds$N > a)
  if (length(idx) == 0) {
    stop("No CTMC interval overlaps [a, b).")
  }

  P <- diag(1, 2, 2)
  for (j in idx) {
    left <- max(a, bounds$n[j])
    right <- min(b, bounds$N[j])
    if (right - left <= 1e-14) {
      next
    }
    P <- P %*% matrix_exponential(right - left, bounds$Q[[j]])
  }
  unname(P)
}

#' CTMC Transition Product Over Possibly Vectorized Intervals
#'
#' Compute transition matrices over one or more half-open intervals `[a, b)`.
#'
#' @param a Numeric scalar or vector of start times.
#' @param b Numeric scalar or vector of end times.
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#'
#' @return List of `2 x 2` transition matrices.
#' @export
Phi_C <- function(a, b, bounds) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  n_out <- max(length(a), length(b))

  if (length(a) == 1) {
    a <- rep(a, n_out)
  }
  if (length(b) == 1) {
    b <- rep(b, n_out)
  }
  if (length(a) != n_out || length(b) != n_out) {
    stop("a and b must have the same length, unless one is scalar.")
  }

  lapply(seq_len(n_out), function(i) ctmc_phi_one(a[i], b[i], bounds))
}

#' Build Half-Open Monthly CTMC Bounds
#'
#' Partition the interval `[t_now, t_hor)` into monthly subintervals and attach
#' the corresponding CTMC generator to each row.
#'
#' @param t_now Date or character scalar. Start date.
#' @param t_hor Date or character scalar. Horizon date.
#' @param Q List of monthly CTMC generators or monthly transition matrices.
#'
#' @return Tibble with interval starts `n`, ends `N`, horizon `tau`, months, and
#'   monthly generators/products.
#' @export
create_bounds <- function(t_now, t_hor, Q) {
  t_now <- as.Date(t_now)
  t_hor <- as.Date(t_hor)
  total_days <- as.numeric(difftime(t_hor, t_now, units = "days"))

  if (length(total_days) != 1 || !is.finite(total_days) || total_days < 0) {
    stop("t_hor must be a date on or after t_now.")
  }

  Q <- ctmc_as_generator_list(Q)
  start_n <- as.numeric(number_of_day(t_now))
  tau <- start_n + total_days

  if (total_days == 0) {
    bounds <- data.frame(
      Year = lubridate::year(t_now),
      Month = lubridate::month(t_now),
      n = start_n,
      N = start_n,
      tau = tau
    )
  } else {
    dates <- seq.Date(t_now, t_hor - 1, by = "day")
    month <- lubridate::month(dates)
    year <- lubridate::year(dates)
    run_id <- cumsum(c(TRUE, diff(month) != 0 | diff(year) != 0))
    offsets <- seq_along(dates) - 1

    bounds <- do.call(rbind, lapply(split(offsets, run_id), function(offset) {
      first <- offset[1]
      last <- offset[length(offset)]
      active_date <- t_now + first
      data.frame(
        Year = lubridate::year(active_date),
        Month = lubridate::month(active_date),
        n = start_n + first,
        N = start_n + last + 1,
        tau = tau
      )
    }))
  }

  bounds <- dplyr::as_tibble(bounds)
  bounds$n_idx <- bounds$n
  bounds$Q <- lapply(bounds$Month, function(m) Q[[m]])
  bounds$Q_prod <- lapply(seq_len(nrow(bounds)), function(j) {
    matrix_exponential(bounds$N[j] - bounds$n[j], bounds$Q[[j]])
  })
  bounds
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' Expected Mixture Drift at a CTMC Time
#'
#' @param s Numeric vector of time points.
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#' @param p0 Numeric length-two initial probability.
#' @param mu List of monthly state means.
#'
#' @return Numeric vector of unconditional expected state means.
#' @keywords internal
e_mu_B_fast <- function(s, bounds, p0, mu){
  # Ensure numeric
  s <- as.numeric(s)
  # Initial time 
  t_init <- rep(bounds$n[1], length(s))
  # Monthly indexes 
  tm_s <- get_month_index_C(s, bounds)
  # Monthly means 
  muq <- sapply(tm_s, function(x) mu[[x]])
  # Monthly transitions 
  Qq <- Phi_C(t_init, s, bounds)
  # Vectorized product 
  prod_p0_Phi <- function(s, Qq) p0 %*% Qq
  p0_Phi <- t(mapply(prod_p0_Phi, s = s, Qq = Qq))
  diag((p0_Phi %*% muq))
}
#' Expected Mixture Diffusion Variance at a CTMC Time
#'
#' @param s Numeric vector of time points.
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#' @param p0 Numeric length-two initial probability.
#' @param sd2 List of monthly state variances.
#'
#' @return Numeric vector of unconditional expected state variances.
#' @keywords internal
e_sigma2_B_fast <- function(s, bounds, p0, sd2){
  # Ensure numeric
  s <- as.numeric(s)
  # Initial time 
  t_init <- rep(bounds$n[1], length(s))
  # Monthly indexes 
  tm_s <- get_month_index_C(s, bounds)
  # Monthly means 
  
  sd2q <- sapply(tm_s, function(x) sd2[[x]])
  # Monthly transitions 
  Qq <- Phi_C(t_init, s, bounds)
  # Vectorized product 
  prod_p0_Phi <- function(s, Qq) p0 %*% Qq
  p0_Phi <- t(mapply(prod_p0_Phi, s = s, Qq = Qq))
  diag((p0_Phi %*% sd2q))
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' Integrated Expected Drift Contribution
#'
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#' @param p0 Numeric length-two initial probability.
#' @param mu List of monthly state means.
#' @param theta Numeric scalar. Mean reversion speed.
#' @param sigma_bar Function returning seasonal volatility.
#' @param t0 Numeric offset for the lower integration bound.
#' @param T0 Numeric offset for the upper integration bound.
#'
#' @return Numeric scalar with absolute integration error stored as an
#'   attribute.
#' @keywords internal
integral_E_mu_tT_fast <- function(bounds, p0, mu, theta, sigma_bar, t0 = 0, T0 = 0){
  # Integration bounds 
  t_init <- bounds$n[1]
  t_end  <- bounds$tau[1]
  # Integrand of E[mu(t,T)]
  integrand_E_mu_tT <- function(s, t_end, bounds, p0, mu, theta, sigma_bar){
    exp(-theta * (t_end - s)) * sigma_bar(s) * e_mu_B_fast(s, bounds, p0, mu)
  }
  intg <-  integrate(integrand_E_mu_tT, 
                     lower = t_init + t0, upper = t_end + T0, stop.on.error = FALSE, 
                     t_end = t_end, bounds = bounds, p0 = p0, mu = mu,
                     theta = theta, sigma_bar = sigma_bar)
  out <- intg$value
  attr(out, "error") <- intg$abs.error
  return(out)
}
#' Integrated Expected Diffusion Variance Contribution
#'
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#' @param p0 Numeric length-two initial probability.
#' @param sd List of monthly state standard deviations.
#' @param theta Numeric scalar. Mean reversion speed.
#' @param sigma_bar Function returning seasonal volatility.
#' @param t0 Numeric offset for the lower integration bound.
#' @param T0 Numeric offset for the upper integration bound.
#'
#' @return Numeric scalar with absolute integration error stored as an
#'   attribute.
#' @keywords internal
integral_E_sigma_tT_fast <- function(bounds, p0, sd, theta, sigma_bar, t0 = 0, T0 = 0){
  # Integration bounds 
  t_init <- bounds$n[1]
  t_end  <- bounds$tau[1]
  # Compute variances 
  sd2 <- purrr::map(sd, ~.x^2)
  # Integrand of E[mu(t,T)]
  integrand_E_sigma_tT <- function(s, t_end, bounds, p0, sd2, theta, sigma_bar){
    exp(-2*theta * (t_end - s)) * sigma_bar(s)^2 * e_sigma2_B_fast(s, bounds, p0, sd2)
  }
  intg <- integrate(integrand_E_sigma_tT, 
                    lower = t_init + t0, upper = t_end + T0, stop.on.error = FALSE, 
                    t_end = t_end, bounds = bounds, p0 = p0, sd2 = sd2,
                    theta = theta, sigma_bar = sigma_bar)
  
  out <- intg$value
  attr(out, "error") <- intg$abs.error
  return(out)
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' Conditional Expected Mixture Drift at a CTMC Time
#'
#' @param s Numeric vector of time points.
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#' @param p0 Numeric length-two initial probability.
#' @param mu List of monthly state means.
#' @param ei Numeric length-two terminal-state indicator.
#' @param p_T Numeric scalar. Probability of the terminal conditioning event.
#'
#' @return Numeric vector of conditional expected state means.
#' @keywords internal
e_mu_B_cond_fast <- function(s, bounds, p0, mu, ei = c(1,0), p_T){
  # Ensure numeric
  s <- as.numeric(s)
  # Initial time 
  t_init <- rep(bounds$n[1], length(s))
  # End time
  t_end <- rep(bounds$tau[1], length(s))
  # Monthly indexes 
  tm_s <- get_month_index_C(s, bounds)
  # Monthly means 
  mu_m <- sapply(tm_s, function(x) mu[[x]])
  # Monthly transitions t_init -> s 
  probs_t_s <- Phi_C(t_init, s, bounds)
  # Monthly transitions s -> t_end
  probs_s_T <- Phi_C(s, t_end, bounds)
  # Vectorized product 
  prod_probs_s <- function(s, probs_t_s, probs_s_T) p0 %*% probs_t_s %*% diag((probs_s_T %*% ei)[,1])
  prod_probs <- t(mapply(prod_probs_s, s = s, probs_t_s = probs_t_s, probs_s_T = probs_s_T)) 
  # Output
  diag((prod_probs %*% mu_m)) / p_T
}
#' Conditional Expected Mixture Diffusion Variance at a CTMC Time
#'
#' @param s Numeric vector of time points.
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#' @param p0 Numeric length-two initial probability.
#' @param sd2 List of monthly state variances.
#' @param ei Numeric length-two terminal-state indicator.
#' @param p_T Numeric scalar. Probability of the terminal conditioning event.
#'
#' @return Numeric vector of conditional expected state variances.
#' @keywords internal
e_sigma2_B_cond_fast <- function(s, bounds, p0, sd2, ei = c(1,0), p_T){
  # Ensure numeric
  s <- as.numeric(s)
  # Initial time 
  t_init <- rep(bounds$n[1], length(s))
  # End time
  t_end <- rep(bounds$tau[1], length(s))
  # Monthly indexes 
  tm_s <- get_month_index_C(s, bounds)
  # Monthly means 
  sd2_m <- sapply(tm_s, function(x) sd2[[x]])
  # Monthly transitions t_init -> s 
  probs_t_s <- Phi_C(t_init, s, bounds)
  # Monthly transitions s -> t_end
  probs_s_T <- Phi_C(s, t_end, bounds)
  # Vectorized product 
  prod_probs_s <- function(s, probs_t_s, probs_s_T) p0 %*% probs_t_s %*% diag((probs_s_T %*% ei)[,1])
  prod_probs <- t(mapply(prod_probs_s, s = s, probs_t_s = probs_t_s, probs_s_T = probs_s_T)) 
  # Output
  diag((prod_probs %*% sd2_m)) / p_T
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' Integrated Drift Conditional on Terminal Regime
#'
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#' @param p0 Numeric length-two initial probability.
#' @param mu List of monthly state means.
#' @param theta Numeric scalar. Mean reversion speed.
#' @param sigma_bar Function returning seasonal volatility.
#' @param t0 Numeric offset for the lower integration bound.
#' @param T0 Numeric offset for the upper integration bound.
#' @param maxEval Integer. Maximum cubature evaluations.
#' @param tol Numeric cubature tolerance.
#'
#' @return Numeric length-two vector with integration errors stored as an
#'   attribute.
#' @keywords internal
integral_E_mu_tT_cond_fast <- function(bounds, p0, mu, theta, sigma_bar, t0 = 0, T0 = 0, maxEval = 100000, tol = 0.00001){
  # Initial time 
  t_init <- bounds$n[1]
  # End time
  t_end <- bounds$tau[1]
  # Probability at maturity 
  p_T <- p0 %*% Phi_C(t_init, t_end, bounds)[[1]]
  
  # Integrand of E[mu(t,T) | B_T = i]
  integrand_E_mu_tT_cond <- function(s, t0, T_, bounds, p0, mu, theta, sigma_bar, p_T){
    # Ensure numeric
    s <- as.numeric(s)
    # Initial time 
    t_init <- rep(bounds$n[1], length(s))
    # End time
    t_end <- rep(bounds$tau[1], length(s))
    # Monthly indexes 
    tm_s <- get_month_index_C(s, bounds)
    # Monthly means 
    mu_m <- sapply(tm_s, function(x) mu[[x]])
    # Monthly transitions t_init -> s 
    probs_t_s <- Phi_C(t_init, s, bounds)
    # Monthly transitions s -> t_end
    probs_s_T <- Phi_C(s, t_end, bounds)
    # Vectorized product 
    prod_probs_s <- function(s, probs_t_s, probs_s_T, ei) p0 %*% probs_t_s %*% diag((probs_s_T %*% ei)[,1])
    prod_probs_1 <- t(mapply(prod_probs_s, s = s, probs_t_s = probs_t_s, probs_s_T = probs_s_T, MoreArgs = list(ei = c(1, 0)))) 
    prod_probs_0 <- t(mapply(prod_probs_s, s = s, probs_t_s = probs_t_s, probs_s_T = probs_s_T, MoreArgs = list(ei = c(0, 1)))) 
    # Seasonal weights
    w_s <- exp(-theta * (t_end - s)) * sigma_bar(s)
    # Final product 
    e_mu_tT_cond_1 <- w_s * diag((prod_probs_1 %*% mu_m)) / p_T[1]
    e_mu_tT_cond_0 <- w_s * diag((prod_probs_0 %*% mu_m)) / p_T[2]
    
    rbind(e_mu_tT_cond_1, e_mu_tT_cond_0)
  }
  
  intg <- cubature::hcubature(
    integrand_E_mu_tT_cond, 
    lowerLimit = t_init + t0,
    upperLimit = t_end + T0,
    bounds = bounds, 
    p0 = p0, 
    mu = mu, 
    theta = theta, 
    sigma_bar = sigma_bar,
    p_T = c(p_T),
    fDim = 2,
    tol = tol,
    norm = "INDIVIDUAL",
    maxEval = maxEval,
    vectorInterface = TRUE)
  out <- intg$integral
  attr(out, "error") <- intg$error
  return(out)
}
#' Integrated Diffusion Variance Conditional on Terminal Regime
#'
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#' @param p0 Numeric length-two initial probability.
#' @param sd List of monthly state standard deviations.
#' @param theta Numeric scalar. Mean reversion speed.
#' @param sigma_bar Function returning seasonal volatility.
#' @param t0 Numeric offset for the lower integration bound.
#' @param T0 Numeric offset for the upper integration bound.
#' @param maxEval Integer. Maximum cubature evaluations.
#' @param tol Numeric cubature tolerance.
#'
#' @return Numeric length-two vector with integration errors stored as an
#'   attribute.
#' @keywords internal
integral_E_sigma_tT_cond_fast <- function(bounds, p0, sd, theta, sigma_bar, t0 = 0, T0 = 0, maxEval = 100000, tol = 0.00001){
  # Initial time 
  t_init <- bounds$n[1]
  # End time
  t_end <- bounds$tau[1]
  # Probability at maturity 
  p_T <- p0 %*% Phi_C(t_init, t_end, bounds)[[1]]
  # Compute variances 
  sd2 <- purrr::map(sd, ~.x^2)
  # Integrand of E[mu(t,T) | B_T = i]
  integrand_E_sigma_tT_cond <- function(s, t0, T_, bounds, p0, sd2, theta, sigma_bar, p_T){
    # Ensure numeric
    s <- as.numeric(s)
    # Initial time 
    t_init <- rep(bounds$n[1], length(s))
    # End time
    t_end <- rep(bounds$tau[1], length(s))
    # Monthly indexes 
    tm_s <- get_month_index_C(s, bounds)
    # Monthly means 
    sd2_m <- sapply(tm_s, function(x) sd2[[x]])
    # Monthly transitions t_init -> s 
    probs_t_s <- Phi_C(t_init, s, bounds)
    # Monthly transitions s -> t_end
    probs_s_T <- Phi_C(s, t_end, bounds)
    # Vectorized product 
    prod_probs_s <- function(s, probs_t_s, probs_s_T, ei) p0 %*% probs_t_s %*% diag((probs_s_T %*% ei)[,1])
    prod_probs_1 <- t(mapply(prod_probs_s, s = s, probs_t_s = probs_t_s, probs_s_T = probs_s_T, MoreArgs = list(ei = c(1, 0)))) 
    prod_probs_0 <- t(mapply(prod_probs_s, s = s, probs_t_s = probs_t_s, probs_s_T = probs_s_T, MoreArgs = list(ei = c(0, 1)))) 
    # Seasonal weights
    w_s <- exp(-2*theta * (t_end - s)) * sigma_bar(s)^2
    # Final product 
    e_sigma_tT_cond_1 <- w_s * diag((prod_probs_1 %*% sd2_m)) / p_T[1]
    e_sigma_tT_cond_0 <- w_s * diag((prod_probs_0 %*% sd2_m)) / p_T[2]
    
    rbind(e_sigma_tT_cond_1, e_sigma_tT_cond_0)
  }
  
  intg <- cubature::hcubature(
    integrand_E_sigma_tT_cond, 
    lowerLimit = t_init + t0,
    upperLimit = t_end + T0,
    bounds = bounds, 
    p0 = p0, 
    sd2 = sd2, 
    theta = theta, 
    sigma_bar = sigma_bar,
    p_T = c(p_T),
    fDim = 2,
    tol = tol,
    norm = "INDIVIDUAL",
    maxEval = maxEval,
    vectorInterface = TRUE)
  out <- intg$integral
  attr(out, "error") <- intg$error
  return(out)
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' Integrated Conditional Moment Components
#'
#' Compute conditional drift, diffusion variance, and diffusion-mean
#' contribution terms used by `radiationModelHMM_moments()`.
#'
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#' @param p0 Numeric length-two initial probability.
#' @param mu List of monthly state means.
#' @param sd List of monthly state standard deviations.
#' @param theta Numeric scalar. Mean reversion speed.
#' @param sigma_bar Function returning seasonal volatility.
#' @param p_T Optional terminal state probabilities.
#' @param t0 Numeric offset for the lower integration bound.
#' @param T0 Numeric offset for the upper integration bound.
#' @param maxEval Integer. Maximum cubature evaluations.
#' @param tol Numeric cubature tolerance.
#'
#' @return Named numeric vector of six conditional moment components with
#'   integration errors stored as an attribute.
#' @keywords internal
integral_E_cond_fast <- function(bounds, p0, mu, sd, theta, sigma_bar, p_T, t0 = 0, T0 = 0, maxEval = 100000, tol = 0.00001){
  # Initial time 
  t_init <- bounds$n[1]
  # End time
  t_end <- bounds$tau[1]
  # Compute variances 
  sd2 <- purrr::map(sd, ~.x^2)
  # Probability at maturity 
  if (missing(p_T)) {
    p_T <- p0 %*% Phi_C(t_init, t_end, bounds)[[1]]
  }
  # Integrand of E[mu(t,T) | B_T = i]
  integrand_E_cond <- function(s, bounds, p0, mu, sd, sd2, theta, sigma_bar, p_T){
    # Ensure numeric
    s <- as.numeric(s)
    # Initial time 
    t_init <- rep(bounds$n[1], length(s))
    # End time
    t_end <- rep(bounds$tau[1], length(s))
    # Monthly indexes 
    tm_s <- get_month_index_C(s, bounds)
    # Monthly parameters
    mu_m <- sapply(tm_s, function(x) mu[[x]])
    sd_m <- sapply(tm_s, function(x) sd[[x]])
    sd2_m <- sapply(tm_s, function(x) sd2[[x]])
    # Monthly transitions t_init -> s 
    probs_t_s <- Phi_C(t_init, s, bounds)
    # Monthly transitions s -> t_end
    probs_s_T <- Phi_C(s, t_end, bounds)
    # Vectorized product 
    prod_probs_s <- function(s, probs_t_s, probs_s_T, ei) p0 %*% probs_t_s %*% diag((probs_s_T %*% ei)[,1])
    prod_probs_1 <- t(mapply(prod_probs_s, s = s, probs_t_s = probs_t_s, probs_s_T = probs_s_T, MoreArgs = list(ei = c(1, 0)))) 
    prod_probs_0 <- t(mapply(prod_probs_s, s = s, probs_t_s = probs_t_s, probs_s_T = probs_s_T, MoreArgs = list(ei = c(0, 1)))) 
    # Seasonal weights
    w_s <- exp(-theta * (t_end - s)) * sigma_bar(s)
    
    # V[sigma(t, T) | B_T = i]
    v_sigma_tT_cond_1 <- w_s^2 * diag((prod_probs_1 %*% sd2_m)) / p_T[1]
    v_sigma_tT_cond_0 <- w_s^2 * diag((prod_probs_0 %*% sd2_m)) / p_T[2]
    # E[sigma(t, T) | B_T = i]
    e_sigma_tT_cond_1 <- w_s * diag((prod_probs_1 %*% sd_m)) / p_T[1]
    e_sigma_tT_cond_0 <- w_s * diag((prod_probs_0 %*% sd_m)) / p_T[2]
    # E[mu(t, T) | B_T = i]
    e_mu_tT_cond_1 <- w_s * diag((prod_probs_1 %*% mu_m)) / p_T[1]
    e_mu_tT_cond_0 <- w_s * diag((prod_probs_0 %*% mu_m)) / p_T[2]
    
    rbind(e_mu_tT_cond_1, e_mu_tT_cond_0, 
          v_sigma_tT_cond_1, v_sigma_tT_cond_0, 
          e_sigma_tT_cond_1, e_sigma_tT_cond_0)
  }
  
  intg <- cubature::hcubature(
    integrand_E_cond, 
    lowerLimit = t_init + t0,
    upperLimit = t_end + T0,
    bounds = bounds, 
    p0 = p0, 
    mu = mu,
    sd = sd,
    sd2 = sd2, 
    theta = theta, 
    sigma_bar = sigma_bar,
    p_T = c(p_T),
    fDim = 6,
    tol = tol,
    norm = "INDIVIDUAL",
    maxEval = maxEval,
    vectorInterface = TRUE)
  out <- intg$integral
  names(out) <- c("E_mu_tT_1", "E_mu_tT_0", 
                  "V_sigma_tT_1", "V_sigma_tT_0",
                  "E_sigma_tT_1", "E_sigma_tT_0")
  attr(out, "error") <- intg$error
  return(out)
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' Integrated Unconditional Drift Cross-Moment
#'
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#' @param p0 Numeric length-two initial probability.
#' @param mu List of monthly state means.
#' @param theta Numeric scalar. Mean reversion speed.
#' @param sigma_bar Function returning seasonal volatility.
#' @param t0 Numeric length-two lower-bound offsets.
#' @param T0 Numeric length-two upper-bound offsets.
#' @param maxEval Integer. Maximum cubature evaluations.
#' @param tol Numeric cubature tolerance.
#'
#' @return Numeric scalar with integration error stored as an attribute.
#' @keywords internal
integral_E_mu_su_fast <- function(bounds, p0, mu, theta, sigma_bar, t0 = c(0, 0), T0 = c(0, 0), maxEval = 100000, tol = 0.00001){
  # Initial time 
  t_init <- bounds$n[1]
  # End time
  t_end <- bounds$tau[1]
  # Integrand for the drift cross-moment.
  integrand_E_mu_su_vec <- function(su, bounds, p0, mu, theta, sigma_bar, mat = TRUE){
    # Extract variables
    if (mat){
      s_ <- su[1,] 
      u_ <- su[2,]
    } else {
      s_ <- su[1] 
      u_ <- su[2]
    }
    # Switch bounds
    idx_switch <- s_ > u_
    tot_switch <- sum(idx_switch)
    s <- s_; u <- u_
    if (tot_switch != 0){
      s[idx_switch] <- u_[idx_switch]
      u[idx_switch] <- s_[idx_switch]
    }
    # Ensure numeric
    s <- as.numeric(s)
    u <- as.numeric(u)
    # Integration bounds
    t_init <- rep(bounds$n[1], length(s))
    t_end <- rep(bounds$tau[1], length(s))
    # Monthly indexes 
    tm_s <- get_month_index_C(s, bounds)
    tm_u <- get_month_index_C(u, bounds)
    # Monthly means 
    muq_u <- sapply(tm_u, function(x) mu[[x]])
    # Monthly transitions 
    probs_t_s <- Phi_C(t_init, s, bounds)
    probs_s_u <- Phi_C(s, u, bounds)
    # Vectorized product
    prod_probs <- function(tm_s, probs_t_s, probs_s_u) p0 %*% probs_t_s %*% diag(mu[[tm_s]]) %*% probs_s_u 
    prod_probs_su <- t(mapply(prod_probs, tm_s = tm_s, probs_t_s = probs_t_s, probs_s_u = probs_s_u))
    # Conditional cross moment 
    C_su <- diag(prod_probs_su %*% muq_u) 
    # Weights
    w_su <- exp(-theta*(2*t_end - s - u)) * sigma_bar(s) * sigma_bar(u) 
    matrix(w_su * C_su, nrow = 1)
  }
  
  intg <- cubature::hcubature(
    integrand_E_mu_su_vec, 
    lowerLimit = c(t_init + t0[1], t_init + t0[2]), 
    upperLimit = c(t_end + T0[1], t_end + T0[2]),
    bounds = bounds, 
    p0 = p0, 
    mu = mu, 
    theta = theta, 
    sigma_bar = sigma_bar,
    mat = TRUE,
    norm = "INDIVIDUAL",
    maxEval = maxEval,
    tol = tol,
    vectorInterface = TRUE)
  out <- intg$integral
  attr(out, "error") <- intg$error
  return(out)
}
#' Integrated Drift Cross-Moment Conditional on Terminal Regime
#'
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#' @param p0 Numeric length-two initial probability.
#' @param mu List of monthly state means.
#' @param theta Numeric scalar. Mean reversion speed.
#' @param sigma_bar Function returning seasonal volatility.
#' @param p_T Optional terminal state probabilities.
#' @param t0 Numeric length-two lower-bound offsets.
#' @param T0 Numeric length-two upper-bound offsets.
#' @param maxEval Integer. Maximum cubature evaluations.
#' @param tol Numeric cubature tolerance.
#'
#' @return Numeric length-two vector with integration errors stored as an
#'   attribute.
#' @keywords internal
integral_E_mu_su_cond_fast <- function(bounds, p0, mu, theta, sigma_bar, p_T, t0 = c(0, 0), T0 = c(0, 0), maxEval = 100000, tol = 0.00001){
  # Initial time 
  t_init <- bounds$n[1]
  # End time
  t_end <- bounds$tau[1]
  # Probability at maturity 
  if (missing(p_T)) {
    p_T <- p0 %*% Phi_C(t_init, t_end, bounds)[[1]]
  }
  
  # Integral E[mu(t,T)]# It_hor = ntegral E[mu(t,T)]
  integrand_E_mu_su_cond_vec <- function(su, bounds, p0, mu, theta, sigma_bar, p_T, mat = TRUE){
    # Extract variables
    if (mat){
      s_ <- su[1,] 
      u_ <- su[2,]
    } else {
      s_ <- su[1] 
      u_ <- su[2]
    }
    # Switch bounds
    idx_switch <- s_ > u_
    tot_switch <- sum(idx_switch)
    s <- s_; u <- u_
    if (tot_switch != 0){
      s[idx_switch] <- u_[idx_switch]
      u[idx_switch] <- s_[idx_switch]
    }
    # Ensure numeric
    s <- as.numeric(s)
    u <- as.numeric(u)
    # Integration bounds
    t_init <- rep(bounds$n[1], length(s))
    t_end <- rep(bounds$tau[1], length(s))
    # Monthly indexes 
    tm_s <- get_month_index_C(s, bounds)
    tm_u <- get_month_index_C(u, bounds)
    # Monthly means 
    mu_m_u <- sapply(tm_u, function(x) mu[[x]])
    # Monthly transitions 
    probs_t_s <- Phi_C(t_init, s, bounds)
    probs_s_u <- Phi_C(s, u, bounds)
    probs_u_T <- Phi_C(u, t_end, bounds)
    # Vectorized product
    prod_probs <- function(tm_s, probs_t_s, probs_s_u, probs_u_T, ei) p0 %*% probs_t_s %*% diag(mu[[tm_s]]) %*% probs_s_u %*% diag((probs_u_T %*% ei)[,1])
    prod_probs_1 <- t(mapply(prod_probs, tm_s = tm_s, probs_t_s = probs_t_s, probs_s_u = probs_s_u, probs_u_T = probs_u_T, MoreArgs = list(ei = c(1, 0))))
    prod_probs_0 <- t(mapply(prod_probs, tm_s = tm_s, probs_t_s = probs_t_s, probs_s_u = probs_s_u, probs_u_T = probs_u_T, MoreArgs = list(ei = c(0, 1))))
    # Conditional cross moment 
    C_su_1 <- diag(prod_probs_1 %*% mu_m_u) / p_T[1]
    C_su_0 <- diag(prod_probs_0 %*% mu_m_u) / p_T[2]
    # Seasonal weights
    w_su <- exp(-theta*(2*t_end - s - u)) * sigma_bar(s) * sigma_bar(u)
    rbind(w_su * C_su_1, w_su * C_su_0)
  }
  
  intg <- cubature::hcubature(
    integrand_E_mu_su_cond_vec, 
    lowerLimit = c(t_init + t0[1], t_init + t0[2]), 
    upperLimit = c(t_end + T0[1], t_end + T0[2]),
    bounds = bounds, 
    p0 = p0, 
    mu = mu, 
    theta = theta, 
    sigma_bar = sigma_bar,
    p_T = c(p_T),
    fDim = 2,
    mat = TRUE,
    maxEval = maxEval,
    tol = tol,
    norm = "INDIVIDUAL",
    vectorInterface = TRUE)
  
  out <- intg$integral
  attr(out, "error") <- intg$error
  return(out)
}

#' Integrated Conditional Cross-Covariance Between Radiation and Electricity
#'
#' Compute the CTMC-regime weighted covariance contribution between transformed
#' radiation and log-electricity shocks.
#'
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#' @param p0 Numeric length-two initial probability.
#' @param sd List of monthly radiation state standard deviations.
#' @param theta Numeric scalar. Radiation mean reversion speed.
#' @param kappa Numeric scalar. Electricity mean reversion speed.
#' @param sigma_bar Function returning radiation seasonal volatility.
#' @param sigma_X Numeric scalar. Electricity log-volatility.
#' @param rho List of monthly state correlations.
#' @param maxEval Integer. Maximum cubature evaluations.
#' @param tol Numeric cubature tolerance.
#'
#' @return Numeric length-two vector with integration errors stored as an
#'   attribute.
#' @keywords internal
integral_E_Yt_Xt <- function(bounds, p0, sd, theta, kappa, sigma_bar, sigma_X, rho,  maxEval = 100000, tol = 0.00001){
  
  # Initial time 
  t_init <- bounds$n[1]
  # End time
  t_end <- bounds$tau[1]
  # Probability at maturity 
  p_T <- p0 %*% Phi_C(t_init, t_end, bounds)[[1]]
 
  integrand_fast <- function(s, bounds, p0, sd, theta, kappa, sigma_bar, sigma_X, rho, p_T){
    # Ensure numeric
    s <- as.numeric(s)
    # Initial time 
    t_init <- rep(bounds$n[1], length(s))
    # End time
    t_end <- rep(bounds$tau[1], length(s))
    # Monthly indexes 
    tm_s <- get_month_index_C(s, bounds)
    # Monthly product 
    rho_sigma_B <- purrr::map(1:12, ~sd[[.x]] * rho[[.x]])
    # Monthly means 
    rho_sigma_B_s <- sapply(tm_s, function(x) rho_sigma_B[[x]])
    # Monthly transitions t_init -> s 
    probs_t_s <- Phi_C(t_init, s, bounds)
    # Monthly transitions s -> t_end
    probs_s_T <- Phi_C(s, t_end, bounds)
    # Vectorized product 
    prod_probs_s <- function(s, probs_t_s, probs_s_T, ei) p0 %*% probs_t_s %*% diag((probs_s_T %*% ei)[,1])
    prod_probs_1 <- t(mapply(prod_probs_s, s = s, probs_t_s = probs_t_s, probs_s_T = probs_s_T, MoreArgs = list(ei = c(1, 0)))) 
    prod_probs_0 <- t(mapply(prod_probs_s, s = s, probs_t_s = probs_t_s, probs_s_T = probs_s_T, MoreArgs = list(ei = c(0, 1)))) 
    # Seasonal weights
    w_s <- exp(-theta * (t_end- s)) * exp(-kappa * (t_end- s)) * sigma_bar(s) * sigma_X 
    # Final product 
    S_XY_cond_1 <- w_s * diag((prod_probs_1 %*% rho_sigma_B_s)) / p_T[1]
    S_XY_cond_0 <- w_s * diag((prod_probs_0 %*% rho_sigma_B_s)) / p_T[2]
    
    rbind(S_XY_cond_1, S_XY_cond_0)
  }
  
  intg <- cubature::hcubature(
    integrand_fast, 
    lowerLimit = c(t_init), 
    upperLimit = c(t_end),
    bounds = bounds, 
    p0 = p0, 
    p_T = p_T,
    sd = sd, 
    theta = theta, 
    kappa = kappa, 
    sigma_bar = sigma_bar,
    sigma_X = sigma_X,
    rho = rho,
    fDim = 2,
    norm = "INDIVIDUAL",
    maxEval = maxEval,
    tol = tol,
    vectorInterface = TRUE)
  intg
  out <- intg$integral
  attr(out, "error") <- intg$error
  return(out)
}
