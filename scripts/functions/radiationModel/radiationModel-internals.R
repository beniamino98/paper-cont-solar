#' Create a monthly sequence of dates
#'
#' @param t_now Character, today date.
#' @param t_hor Character, horizon date.
#' @param last_day Logical. When `TRUE` the last day will be treated as conditional variance otherwise not.
#' @examples
#' t_now <- "2022-01-01"
#' t_hor <- "2022-03-24"
#' create_monthly_sequence("2022-01-01", "2022-03-24")
#' create_monthly_sequence("2022-01-01", "2022-03-24", TRUE)
#'
#' @return A tibble describing the monthly sub-intervals between `t_now` and
#'   `t_hor`.
#' @keywords internal integral radiationModel
#' @export
#' @noRd
create_monthly_sequence <- function(t_now, t_hor, last_day = FALSE){
  # Convert in dates
  t_now <- as.Date(t_now)
  t_hor <- as.Date(t_hor)
  # Compute year and months for reference dates
  # Standard dates
  t_start <- as.Date(paste0(lubridate::year(t_now), "-", lubridate::month(t_now), "-01"))
  t_end <- as.Date(paste0(lubridate::year(t_hor), "-", lubridate::month(t_hor), "-01"))
  # Sequence of dates
  dates_hor <- dates_now <- seq.Date(t_start, t_end, by = "1 month")
  n_months <- length(dates_now)
  dates_now <- c(t_now, dates_now[-1] - lubridate::day(dates_now[-1]))
  dates_hor <- c(dates_hor[-n_months] - lubridate::day(dates_hor[-n_months]) + lubridate::days(lubridate::days_in_month(dates_hor[-n_months])), t_hor)
  n_of_day <- as.numeric(difftime(dates_hor, dates_now, units = "days"))
  df_dates <- data.frame(Year = lubridate::year(dates_now), Month = lubridate::month(dates_hor), n = n_of_day)
  df_dates$N <- df_dates$n + number_of_day(t_now) + c(0, cumsum(lag(df_dates$n, 1)[-1]))
  df_dates$n <- df_dates$N - df_dates$n
  df_dates$tau <- as.numeric(difftime(t_hor, t_now, units = "days")) + df_dates$n[1]
  if (last_day) {
    df_last <- dplyr::mutate(tail(df_dates, 1), n = tau - 1)
    df_dates[nrow(df_dates),]$N <- df_dates[nrow(df_dates),]$N - 1
    df_dates <- dplyr::bind_rows(df_dates, df_last)
  }
  attr(df_dates, "last_day") <- last_day
  return(df_dates)
}
#' Estimate Mean Reversion by a Martingale Moment Condition
#'
#' Estimate the continuous-time mean reversion parameter from the transformed
#' radiation process and its seasonal mean.
#'
#' @param Yt Numeric vector. Transformed radiation observations.
#' @param Yt_bar Numeric vector. Seasonal mean of `Yt`.
#' @param e_mu Numeric vector or scalar. Optional drift correction.
#'
#' @return List with `theta`, `phi`, and approximate standard errors.
#' @keywords internal
martingale_method_seasonal <- function(Yt, Yt_bar, e_mu = 0){
  # Dataset
  data <- dplyr::tibble(Yt = Yt, Yt_bar = Yt_bar, e_mu = e_mu)
  data$n <- 1:nrow(data)
  # Quadratic variation
  data$dYt2 <- (lag(data$Yt,1) - lag(data$Yt, 2))^2
  # Seeasonal variance from quadratic variation
  seasonal_variance <- seasonalModel$new(formula = "dYt2 ~ 1")
  seasonal_variance$fit(data = data)
  # Estimated seasonal variance
  data$sigma2_bar <- seasonal_variance$predict(1:nrow(data)-1)
  # Martingale estimation
  data$Y_est_L1 <- (lag(data$Yt_bar, 1) - lag(data$Yt, 1)) / data$sigma2_bar
  # Differences from seasonal mean
  data$dYt <-  data$Yt - data$Yt_bar
  data$dYt_L1 <- lag(data$Yt, 1) - lag(data$Yt_bar, 1) - lag(data$e_mu, 1)
  data <- na.omit(data)
  a <- data$Y_est_L1 * data$dYt
  b <- data$Y_est_L1 * data$dYt_L1
  A <- sum(a) 
  B <- sum(b)
  # Bartlett HAC for vector time series X (T x k)
  hac_cov <- function(X, L = NULL){
    X <- as.matrix(X)
    Tn <- nrow(X); 
    k <- ncol(X)
    if(is.null(L)) L <- floor(4*(Tn/100)^(2/9))  # common automatic choice
    Xc <- sweep(X, 2, colMeans(X), "-")
    Gamma0 <- crossprod(Xc)/Tn
    S <- Gamma0
    for(l in 1:L){
      w <- 1 - l/(L+1)  # Bartlett
      Gl <- crossprod(Xc[(l+1):Tn, , drop=FALSE], Xc[1:(Tn-l), , drop=FALSE]) / Tn
      S <- S + w*(Gl + t(Gl))
    }
    # This is HAC for the mean: Var(sqrt(T) * mean) ~ S
    # Hence Var(sum) = Var(T * mean) ~ T * S
    list(S_mean = S, S_sum = Tn * S, L = L)
  }
  
  V_AB <- hac_cov(cbind(a, b))$S_sum
  g <- c(-1/A, 1/B)                            # gradient wrt (A,B)
  var_theta <- drop(t(g) %*% V_AB %*% g)
  std.error.theta  <- sqrt(var_theta)
  theta_N <- -log(A/B)
  phi <- exp(-theta_N)
  std.error.phi <- phi * std.error.theta
  
  list(
    theta = theta_N,
    phi = phi, 
    std.error.theta = std.error.theta,
    std.error.phi = std.error.phi
  )
}

#' Reparameterize the Seasonal Variance Function
#'
#' Convert discrete seasonal variance coefficients into continuous-time
#' coefficients used by the integrated radiation model.
#'
#' @param par Numeric vector of seasonal variance coefficients.
#' @param theta Numeric scalar. Mean reversion speed.
#' @param omega Numeric scalar. Seasonal angular frequency.
#'
#' @return List with original, continuous-time, and integral coefficients.
#' @keywords internal
reparam_seasonal_function <- function(par, theta, omega){
  # Original parameters
  a0 <- par[1]
  a1 <- par[2]
  a2 <- par[3]
  # Correction for long term variance
  c0_long <- a0 * 2 * theta
  c1_long <- a1 * 2 * theta - omega * a2
  c2_long <- a2 * 2 * theta + omega * a1
  # Integral parameters
  gamma0_long <-  c0_long / (2 * theta)
  gamma1_long <- (c1_long * 2 * theta + c2_long * omega) / (4 * theta^2 + omega^2)
  gamma2_long <- (c2_long * 2 * theta - c1_long * omega) / (4 * theta^2 + omega^2)
  # Correction for short term variance
  alpha  <- 1 - exp(-2 * theta) * cos(omega)
  beta  <- exp(-2 * theta) * sin(omega)
  detM  <- alpha^2 + beta^2
  # Dynamic parameters 
  c0 <- (2 * theta * a0) / (1 - exp(-2 * theta))
  c1 <- ((2 * theta * alpha + omega * beta) * a1 + (2 * theta * beta - omega * alpha) * a2) / detM
  c2 <- ((omega * alpha - 2 * theta * beta) * a1 + (omega * beta + 2 * theta * alpha) * a2) / detM
  # Integral parameters
  gamma0 <-  c0 / (2 * theta)
  gamma1 <- (c1 * 2 * theta + c2 * omega) / (4 * theta^2 + omega^2)
  gamma2 <- (c2 * 2 * theta - c1 * omega) / (4 * theta^2 + omega^2)
  # Output 
  structure(
    list(
      alpha = alpha,
      beta = beta,
      detM = detM,
      a_ = c(a0 = a0, a1 = a1, a2 = a2),
      c_ = c(c0 = c0, c1 = c1, c2 = c2),
      c_long = c(c0 = c0_long, c1 = c1_long, c2 = c2_long),
      gamma = c(gamma0 = gamma0, gamma1 = gamma1, gamma2 = gamma2),
      gamma_long = c(gamma0 = gamma0_long, gamma1 = gamma1_long, gamma2 = gamma2_long)
    )
  )
}
#' Build a Numerical Seasonal Volatility Integral
#'
#' Return a vectorized function that computes
#' \eqn{\int_t^s \bar\sigma_u \exp[-\theta(T-u)]du}.
#'
#' @param theta Numeric scalar. Mean reversion speed.
#' @param par Numeric vector. Seasonal variance coefficients.
#' @param omega Numeric scalar. Seasonal angular frequency.
#'
#' @return A function of `(t, s, T_)`.
#' @keywords internal
integral_sigma_numeric <- function(theta, par, omega = 2*base::pi/365){
  seasonal_function <- function(t) par[1] + par[2] * sin(omega * t) + par[3] * cos(omega * t)
  integrand <- function(tau, T_) sqrt(seasonal_function(tau) * exp(-2*theta*(T_-tau)))
  function(t, s, T_){
    result <- c()
    for(i in 1:length(T_)){
      result[i] <- integrate(integrand, lower = t[i], upper = s[i], T_ = T_[i])$value
    }
    return(result)
  }
}
#' Build a Closed-Form Seasonal Variance Integral
#'
#' Return a vectorized function that computes the analytical integral of the
#' squared seasonal volatility term.
#'
#' @param theta Numeric scalar. Mean reversion speed.
#' @param par Numeric vector. Integral seasonal variance coefficients.
#' @param omega Numeric scalar. Seasonal angular frequency.
#'
#' @return A function of `(t, s, T_)`.
#' @keywords internal
integral_sigma2_formula <- function(theta, par, omega = 2*base::pi/365){
  # Functions
  f0 <- function(t, s, T_) exp(-2 * theta * (T_ - s)) - exp(-2 * theta * (T_ - t))
  f1 <- function(t, s, T_) exp(-2 * theta * (T_ - s)) * sin(omega * s) - exp(-2 * theta * (T_ - t)) * sin(omega * t)
  f2 <- function(t, s, T_) exp(-2 * theta * (T_ - s)) * cos(omega * s) - exp(-2 * theta * (T_ - t)) * cos(omega * t)
  
  function(t, s, T_){
    result <- par[1] * f0(t, s, T_) + par[2] * f1(t, s, T_) + par[3] * f2(t, s, T_)
    return(result)
  }
}
#' Build a Numerical Double Integral for Seasonal Volatility Products
#'
#' Return a function that integrates products of seasonal volatility terms over
#' a square time domain.
#'
#' @param theta Numeric scalar. Mean reversion speed.
#' @param par Numeric vector. Seasonal variance coefficients.
#' @param omega Numeric scalar. Seasonal angular frequency.
#'
#' @return A function of `(t, s, T_)`.
#' @keywords internal
integral_sigma_ij_numeric <- function(theta, par, omega = 2*base::pi/365){
  seasonal_function <- function(t) par[1] + par[2] * sin(omega * t) + par[3] * cos(omega * t)
  integrand <- function(x, T_) matrix(seasonal_function(x[1,]) * exp(-theta*(T_-x[1,])) * seasonal_function(x[2,]) * exp(-theta*(T_-x[2,])), nrow = 1)
  function(t, s, T_){
    result <- c()
    for(i in 1:length(T_)){
      intg <- cubature::hcubature(
        integrand, 
        lowerLimit = c(t[i], t[i]), 
        upperLimit = c(s[i], s[i]),
        T_ = T_[i],
        norm = "INDIVIDUAL",
        maxEval = maxEval,
        vectorInterface = TRUE)
    result[i] <- intg$integral
    attr(result[i], "error") <- intg$error
    }
    return(result)
  }
}
#' @description
#' Integral mixture drift of both component of \eqn{Y_t}.
#' @param t_now Character, today date.
#' @param t_hor Character, horizon date.
#' @param model_Rt A `radiationModel` object.
#' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
#' @return Mixture expected value for both component of \eqn{Y_t}.
#' @keywords internal
integral_drift <- function(t_now, t_hor, model_Rt, df_date){
  
  if (missing(df_date)){
    # Create a sequence of dates
    df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
  }
  # Adjust parameters for Q-measure
  df_params <- model_Rt$model$spec$mixture.model$coefficients
  # Expected drift
  df_params$e_mu_B <- df_params$mu1 * df_params$p1 + df_params$mu2 * df_params$p2
  # Expected diffusion drift
  df_params$e_sigma_B <- df_params$sd1 * df_params$p1 + df_params$sd2 * df_params$p2
  # Combine the datasets
  df <- merge(df_date, df_params, by = "Month", all.x = TRUE)
  # Keep monthly intervals in chronological order after the merge.
  df <- df[order(df$n),]
  # Compute the integral
  df$int <- model_Rt$integral_expectation(df$n, df$N, df_date)$int_sigma
  # Index for last day
  nrows <- nrow(df)
  # Total drift from time t up to T-1 for expectation under P
  mix_drift_mu <- df[-nrows,]$int * df[-nrows,]$e_mu_B
  # Conditional drift from time T-1 up to T for expectation under P
  mix_drift_mu_1 <- sum(mix_drift_mu) + df[nrows,]$int * df[nrows,]$mu1
  mix_drift_mu_2 <- sum(mix_drift_mu) + df[nrows,]$int * df[nrows,]$mu2
  # Total drift from time t up to T-1 for expectation under Q
  mix_drift_sigma_B <- df[-nrows,]$int * df[-nrows,]$e_sigma_B
  # Conditional drift from time T-1 up to T for expectation under Q
  mix_drift_sd_1 <- sum(mix_drift_sigma_B) + df[nrows,]$int * df[nrows,]$sd1
  mix_drift_sd_2 <- sum(mix_drift_sigma_B) + df[nrows,]$int * df[nrows,]$sd2
  # Drift depending on lambda
  mix_drift_1 <- mix_drift_mu_1 + mix_drift_sd_1 * model_Rt$lambda
  mix_drift_2 <- mix_drift_mu_2 + mix_drift_sd_2 * model_Rt$lambda
  
  list(e_drift_1 = mix_drift_1, 
       e_drift_2 = mix_drift_2,
       e_mu1 = mix_drift_mu_1, 
       e_mu2 = mix_drift_mu_2,
       e_sd1 = mix_drift_sd_1, 
       e_sd2 = mix_drift_sd_2)
}

#' @description
#' Integral mixture diffusion of both component of \eqn{Y_t}.
#' @param t_now Character, today date.
#' @param t_hor Character, horizon date.
#' @param model_Rt A `radiationModel` object.
#' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
#' @return Mixture expected value for both component of \eqn{Y_t}.
#' @keywords internal
integral_diffusion = function(t_now, t_hor, model_Rt, df_date){
  if (missing(df_date)){
    # Create a sequence of dates
    df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
  }
  # Mixture parameters
  NM_model <- model_Rt$model$spec$mixture.model
  # Adjust parameters for Q-measure
  df_params <- NM_model$coefficients %>%
    dplyr::mutate(mu_diff = mu1 - mu2,
                  sigma_diff = sd1 - sd2,
                  sigma2_1 = sd1^2,
                  sigma2_2 = sd2^2,
                  e_sigma2 = sigma2_1 * p1 + sigma2_2 * p2,
                  v_mu = mu_diff^2 * p1 * p2,
                  v_sigma = sigma_diff^2 * p1 * p2)
  # Combine the datasets
  df <- merge(df_date, df_params, by = "Month", all.x = TRUE)
  df <- df[order(df$n),]
  # Compute the integral
  df$int_sigma2 <- model_Rt$integral_variance(df$n, df$N, df_date)$int_sigma2
  # Index for last day
  nrows <- nrow(df)
  df_t <- df[-nrows,]
  df_T <- df[nrows,]
  variance_drift_P <- variance_drift_Q <- variance_diffusion <- 0
  if (nrow(df_t) > 1){
    # Drift variance P in t, T-1
    variance_drift_P <- sum(df_t$int_sigma2 * df_t$v_mu)
    # Diffusion variance Q in t, T-1
    variance_drift_Q <- sum(df_t$int_sigma2 * df_t$v_sigma)
    # Diffusion second moment in t, T-1
    variance_diffusion <- sum(df_t$int_sigma2 * df_t$e_sigma2)
  }
  # Total variance in t, T-1
  common_variance <- variance_drift_P + variance_diffusion + variance_drift_Q * model_Rt$lambda^2
  # Realized variance for each component between T-1 and T
  variance_1 <- common_variance + df_T$int_sigma2 * df_T$sigma2_1
  variance_2 <- common_variance + df_T$int_sigma2 * df_T$sigma2_2
  
  list(variance_1 = variance_1, 
       variance_2 = variance_2,
       variance_drift_P = variance_drift_P,
       variance_drift_Q = variance_drift_Q,
       common_variance = common_variance + df_T$int_sigma2 * (df_T$v_mu + df_T$e_sigma2),
       variance_diffusion = variance_diffusion,
       last_1 = df_T$int_sigma2 * df_T$sigma2_1, last_2 = df_T$int_sigma2 * df_T$sigma2_2)
  
}

