#' Radiation model
#'
#' @rdname radiationModel
#' @name radiationModel
#' @note Version 1.0.0
#' @export
radiationModel <- R6::R6Class("radiationModel",
                              # ====================================================================================================== #
                              #                                             Public slots
                              # ====================================================================================================== #
                              public = list(
                                #' @field theta Numeric, mean reversion parameter.
                                theta = NA,
                                #' @field lambda_S Numeric, market risk premium Q-measure.
                                lambda_S = 0,
                                #' @description
                                #' Initialize a `radiationModel` object
                                #' @param model A `solarModel` object. See \code{\link{solarModel}}.
                                #' @param martingale.method Logical. If `TRUE`, estimate the OU parameter with the martingale method.
                                #' @param means.correction Logical. If `TRUE`, correct mixture means for the variance scaling.
                                #' @return Initializes the object in place.
                                initialize = function(model, martingale.method = TRUE, means.correction = FALSE){
                                  # Fit continuous model 
                                  model.cont <- radiationModel_fit(model, martingale.method = martingale.method, means.correction = means.correction)
                                  # Store the model
                                  private$..model <- model.cont$model$clone(TRUE)
                                  # Store mean reversion parameter
                                  self$theta <- model.cont$theta
                                  # Store seasonal function  
                                  private$..seasonal_variance <- model.cont$sigma2_bar$clone(TRUE)
                                  # Integral functions
                                  reparam <- model.cont$sigma2_bar$extra_params$reparam
                                  private$..integral_variance <- integral_sigma2_formula(self$theta, reparam$gamma, omega = model.cont$sigma2_bar$omega)
                                  private$..integral_expectation <- integral_sigma_numeric(self$theta, reparam$c_, omega = model.cont$sigma2_bar$omega)
                                },
                                #' @description
                                #' Change the reference probability measure
                                #' @param measure Character, probability measure. Can be `P` or `Q`.
                                change_measure = function(measure){
                                  measure <- match.arg(measure, choices = c("P", "Q"))
                                  private$..measure <- measure
                                  if (measure == "Q"){
                                    private$..lambda <- self$lambda_S
                                  } else {
                                    private$..lambda <- 0
                                  }
                                },
                                #' @description
                                #' Forecast the probability of the first mixture state.
                                #' @param t_hor Date or character scalar. Horizon date.
                                #' @param ... Additional arguments passed through for compatibility.
                                #' @return Numeric probability of state 1 at `t_hor`.
                                prob = function(t_hor, ...){
                                  self$model$spec$mixture.model$prob$predict(t_hor)
                                },
                                #' @description
                                #' Clear sky radiation for a day of the year.
                                #' @param t_now Character, today date.
                                #' @return Clear sky radiation at time t_now.
                                Ct = function(t_now){
                                  t <- number_of_day(t_now)
                                  self$model$spec$seasonal_model_Ct$predict(number_of_day(t_now))
                                },
                                #' @description
                                #' Differential of the Clear-sky radiation for a day of the year.
                                #' @param t_now Character, today date.
                                #' @return Differential Clear sky radiation at time t_now.
                                dCt = function(t_now){
                                  t <- number_of_day(t_now)
                                  self$model$spec$seasonal_model_Ct$differential(t)
                                },
                                #' @description
                                #' Seasonal mean for the transformed variable \eqn{Y_t} for a given day of the year.
                                #' @param t_now Character, today date.
                                #' @return Seasonal mean for \eqn{Y_t} at time t_now.
                                Yt_bar = function(t_now){
                                  t <- number_of_day(t_now)
                                  self$model$spec$seasonal.mean$predict(t)
                                },
                                #' @description
                                #' Differential of the seasonal mean for a day of the year.
                                #' @param t_now Character, today date.
                                #' @return Differential Clear sky radiation at time t_now.
                                dYt_bar = function(t_now){
                                  t <- number_of_day(t_now)
                                  self$model$spec$seasonal.mean$differential(t)
                                },
                                #' @description
                                #' Seasonal mean for the solar radiation for a given day of the year.
                                #' @param t_now Character, today date.
                                #' @return Seasonal mean for Rt.
                                Rt_bar = function(t_now){
                                  t <- number_of_day(t_now)
                                  self$model$spec$transform$iRY(self$Yt_bar(t), self$Ct(t))
                                },
                                #' @description
                                #' Transformed variable instantaneous seasonal std. deviation \eqn{\bar{\sigma_{t}}}.
                                #' @param t_now Character, today date.
                                #' @return Seasonal std. deviation for Yt on date t_now.
                                sigma_bar = function(t_now){
                                  t <- number_of_day(t_now)
                                  #sqrt(private$..seasonal_variance$predict(t))
                                  private$..seasonal_variance$extra_params$seasonal_function(t)
                                },
                                #' @description
                                #' Transformed variable mixture mean drift \eqn{\mu_(B)}.
                                #' @param t_now Character, today date.
                                #' @param B Integer. If `B = 1` it is returned the mean of the first component,
                                #' otherwise if `B = 0` the second.
                                #' @return Mixture seasonal drift for \eqn{Y_t} at time t_now.
                                mu_B = function(t_now, B = 1){
                                  NM_model <- self$model$spec$mixture.model
                                  result <- ifelse(B == 1, NM_model$mu1$predict(t_now), NM_model$mu2$predict(t_now))
                                  return(result)
                                },
                                #' @description
                                #' Transformed variable mixture diffusion drift \eqn{\sigma_(B)}
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @return Mixture seasonal diffusion for \eqn{Y_t}.
                                sigma_B = function(t_now, B = 1){
                                  NM_model <- self$model$spec$mixture.model
                                  result <- ifelse(B == 1, NM_model$sd1$predict(t_now), NM_model$sd2$predict(t_now))
                                  return(result)
                                },
                                #' @description
                                #' Transformed variable drift \eqn{\mu_(Y)}.
                                #' @param Yt Numeric, transformed solar radiation.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @param dt Numeric, time step.
                                #' @return Mixture drift for \eqn{Y_t}.
                                mu_Y = function(Yt, t_now, B = 1, dt = 1){
                                  # Number of the day
                                  t <- number_of_day(t_now)
                                  # Seasonal variance
                                  sigma_bar <- self$sigma_bar(t)
                                  # Drift for Yt
                                  drift_Y <- self$dYt_bar(t) - self$theta * (Yt - self$Yt_bar(t)) + sigma_bar * self$mu_B(t_now, B)
                                  # Add the market price of risk under the Q measure.
                                  drift_Y + sigma_bar * self$sigma_B(t_now, B) * self$lambda * ifelse(self$measure == "Q", 1, 0)
                                },
                                #' @description
                                #' Transformed variable diffusion \eqn{\sigma_(Y)}.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @return Diffusion for \eqn{Y_t}.
                                sigma_Y = function(t_now, B = 1){
                                  # Seasonal variance
                                  sigma_bar <- self$sigma_bar(t_now)
                                  # Diffusion for Rt process
                                  sigma_bar * self$sigma_B(t_now, B)
                                },
                                #' @description
                                #' Solar radiation drift \eqn{\mu_(R)}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @param dt Numeric, time step.
                                #' @return Drift for \eqn{R_t}.
                                mu_R = function(Rt, t_now, B = 1, dt = 1){
                                  # Number of the day
                                  n <- number_of_day(t_now)
                                  # Clear-sky at time t
                                  Ct <- self$Ct(n)
                                  # Convert Rt to Yt
                                  Yt <- self$model$spec$transform$RY(Rt, Ct)
                                  # Drift for Ct
                                  dCt_dt <- self$dCt(n) / dt
                                  # Transform parameters
                                  alpha <- self$model$spec$transform$alpha
                                  beta <- self$model$spec$transform$beta
                                  # Clearness index
                                  Kt <- 1 - alpha - beta * exp(-exp(Yt))
                                  Kt * dCt_dt + Ct * beta * exp(Yt - exp(Yt)) * (self$mu_Y(Yt, t_now, B, dt) + 0.5 * (1 - exp(Yt)) * self$sigma_Y(t_now, B)^2)
                                },
                                #' @description
                                #' Solar radiation diffusion \eqn{\sigma_(R)}.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @return Diffusion for \eqn{R_t}.
                                sigma_R = function(Rt, t_now, B = 1){
                                  # Clear-sky at time t
                                  Ct <- self$Ct(t_now)
                                  # Convert Rt to Yt
                                  Yt <- self$model$spec$transform$RY(Rt, Ct)
                                  # Diffusion for Rt process
                                  Ct * self$model$spec$transform$beta * exp(Yt - exp(Yt)) * self$sigma_Y(t_now, B)
                                },
                                #' @description
                                #' Compute the integral for expectation \eqn{\mu_(t,T)}.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @param last_day Logical. When `TRUE` the last day will be treated as conditional variance otherwise not.
                                integral_expectation = function(t_now, t_hor, df_date, last_day = TRUE){
                                  # Create a sequence of dates for constant monthly parameters
                                  if (missing(df_date)) {
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = last_day)
                                  }
                                  # Compute the integral
                                  df_date$int_sigma <- private$..integral_expectation(df_date$n, df_date$N, df_date$tau)
                                  return(df_date)
                                },
                                #' @description
                                #' Compute the integral for variance \eqn{\sigma^2_(t,T)}.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @param last_day Logical. When `TRUE` the last day will be treated as conditional variance otherwise not.
                                integral_variance = function(t_now, t_hor, df_date, last_day = TRUE){
                                  # Create a sequence of dates for constant monthly parameters
                                  if (missing(df_date)) {
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = last_day)
                                  }
                                  # Compute the integral
                                  df_date$int_sigma2 <- private$..integral_variance(df_date$n, df_date$N, df_date$tau)
                                  return(df_date)
                                },
                                #' @description
                                #' Method print for `radiationModel` object.
                                print = function(){
                                  self$model$print()
                                  cat("MRP: ", self$lambda, "\n")
                                  cat("Measure: ", self$measure, "\n")
                                  cat("Version: ", private$version, "\n")
                                }
                              ),
                              # ====================================================================================================== #
                              #                                             Private slots
                              # ====================================================================================================== #
                              private = list(
                                version = "1.0.0",
                                ..model = NA,
                                ..measure = "P",
                                ..lambda = 0,
                                ..integral_variance = NA,
                                ..integral_expectation = NA,
                                ..seasonal_variance = NA
                              ),
                              # ====================================================================================================== #
                              #                                             Active slots
                              # ====================================================================================================== #
                              active = list(
                                #' @field model An object of the class `solarModel`.
                                model = function(){
                                  private$..model
                                },
                                #' @field measure Character, reference probability measure used.
                                measure = function(){
                                  private$..measure
                                },
                                #' @field lambda Numeric, market risk premium used.
                                lambda = function(){
                                  private$..lambda
                                },
                                #' @field lambda Numeric, market risk premium used.
                                seasonal_variance = function(){
                                  private$..seasonal_variance
                                }
                              )
)


#' Fit the Continuous-Time Radiation Model
#'
#' Clone a fitted discrete radiation model, estimate/reuse the mean-reversion
#' parameter, reparameterize the seasonal variance, and build the integral
#' functions used by `radiationModel`.
#'
#' @param model A fitted discrete solar radiation model.
#' @param martingale.method Logical. If `TRUE`, estimate mean reversion with the
#'   martingale moment condition.
#' @param means.correction Logical. If `TRUE`, rescale mixture means to account
#'   for continuous-time seasonal volatility integration.
#'
#' @return A list containing the fitted cloned model and continuous-time
#'   auxiliary objects.
#' @keywords internal
radiationModel_fit <- function(model, martingale.method = TRUE, means.correction = TRUE){
  
  model.cont <- model$clone(TRUE)
  # 1) Estimate mean reversion parameter
  if (martingale.method) {
    df <- dplyr::filter(model$data, isTrain & weights != 0)
    theta_N <- martingale_method_seasonal(df$Yt, df$Yt_bar)
    theta <- theta_N$theta
    # Update AR model 
    model.cont$spec$mean.model$update(c(phi_1 = theta_N$phi))
    model.cont$spec$mean.model$update_std.errors(c(phi_1 = theta_N$std.error.phi))
    model.cont$filter()
    model.cont$fit_seasonal_variance()
    model.cont$fit_mixture_model()
    model.cont$update_classification()
    model.cont$update_moments()
    model.cont$update_logLik()
  } else {
    theta <- -log(model.cont$spec$mean.model$phi)
  }
  
  # 2) Reparametrize seasonal function with continuous time parameters
  sigma2_bar <- model.cont$spec$seasonal.variance$clone(TRUE)
  # Reparametrize seasonal function with continuous time parameters
  reparam <- reparam_seasonal_function(sigma2_bar$coefficients, theta, omega = sigma2_bar$omega)
  # Update coefficients to match c parameters 
  names(reparam$c_) <- names(sigma2_bar$coefficients)
  sigma2_bar$update(reparam$c_)
  # Store reparametrization 
  sigma2_bar$extra_params$reparam <- reparam
  # Seasonal function for integration 
  seasonal_function <- function(c_, omega){
    force(c_); force(omega);
    function(t){
      sqrt(c_[1] + c_[2]*sin(omega*t) + c_[3]*cos(omega*t))
    }
  }
  sigma2_bar$extra_params$seasonal_function <- seasonal_function(reparam$c_, 2*base::pi/365)
  # Integrals functions
  integral_variance <- integral_sigma2_formula(theta, reparam$gamma, omega = sigma2_bar$omega)
  integral_expectation <- integral_sigma_numeric(theta, reparam$c_, omega = sigma2_bar$omega)
  
  # Correct NM parameters
  k1 <- 1
  if (means.correction) {
    # Adjustment for the mean that multiplies I
    t <- 1:365
    J_t <- integral_variance(t-1, t, t)
    I_t <- integral_expectation(t-1, t, t)
    k1 <- mean(I_t/ sqrt(J_t))
        # Adjusted means
    means <- model.cont$spec$mixture.model$means * k1
    # Update mixture parameters
    model.cont$spec$mixture.model$update(means = means)
    model.cont$update_moments()
    model.cont$update_logLik()
  }
  
  structure(
    list(
      model = model.cont, 
      martingale.method = martingale.method, 
      means.correction = means.correction, 
      theta = theta,
      sigma2_bar = sigma2_bar, 
      k1 = k1, 
      integral_variance = integral_variance, 
      integral_expectation = integral_expectation
    )
  )
}

#' Moment-Matched Radiation Distribution for the IID Mixture Model
#'
#' Compute the two-component Gaussian moment approximation for transformed
#' radiation and expose the corresponding density/CDF functions for both
#' transformed radiation and physical GHI.
#'
#' @param t_now Date or character scalar. Conditioning date.
#' @param t_hor Date or character vector. Horizon date(s).
#' @param model_Rt A `radiationModel` object.
#' @param R0 Optional numeric initial GHI at `t_now`. If missing, it is read
#'   from `model_Rt$model$data`.
#'
#' @return A tibble with horizon dates, component means/standard deviations,
#'   mixture weights, transformed and GHI density/CDF closures, and diagnostic
#'   moment components.
#' @export
radiationModel_moments <- function(t_now, t_hor, model_Rt, R0 = NULL){
  
  # Number of horizons 
  n.ahead <- length(t_hor)
  if (n.ahead > 1 && length(t_now) != 1) {
    warning("The length of `t_now` should be 1! \n Only first element is used!")
    t_now <- t_now[1]
  }
  # Convert in dates
  t_now <- as.Date(rep(t_now, n.ahead))
  t_hor <- as.Date(t_hor) 
  # Time to maturity in days
  tau <- as.numeric(difftime(t_hor, t_now, units = "days"))
  # ************************************************************
  # Radiation at time t_now 
  if (is.null(R0)) {
    R0 <- dplyr::filter(model_Rt$model$data, date %in% t_now)$GHI  
  }
  # Transform bounds 
  alpha = model_Rt$model$spec$transform$alpha
  beta = model_Rt$model$spec$transform$beta
  # Mean reversion parameter
  theta <- model_Rt$theta
  # ************************************************************
  # Clear-sky at time t_now 
  C0 <- model_Rt$Ct(t_now)
  # Transformed variable at time t_now 
  Y0 <- model_Rt$model$spec$transform$RY(R0, C0)
  # Seasonal mean of Y at time t_now and t_hor 
  Y0_bar <- model_Rt$Yt_bar(t_now)
  YT_bar <- model_Rt$Yt_bar(t_hor)
  # Seasonal + AR expectation 
  e_Yt <- YT_bar + (Y0 - Y0_bar) * exp(-theta * tau)
  # ************************************************************
  # 2) Expectation of mixture drifts
  e_mix_drift <- purrr::map2(t_now, t_hor, ~unlist(integral_drift(.x, .y, model_Rt)))
  e_mu_tT_1 <- purrr::map_dbl(e_mix_drift, ~.x["e_mu1"])
  e_mu_tT_0 <- purrr::map_dbl(e_mix_drift, ~.x["e_mu2"])
  e_sigma_tT_1 <- purrr::map_dbl(e_mix_drift, ~.x["e_sd1"])
  e_sigma_tT_0 <- purrr::map_dbl(e_mix_drift, ~.x["e_sd2"])
  # 3) Diffusion 
  e_mix_diffusion <- purrr::map2(t_now, t_hor, ~unlist(integral_diffusion(.x, .y, model_Rt)))
  v_mu_tT_P <- purrr::map_dbl(e_mix_diffusion, ~.x["variance_drift_P"])
  v_mu_tT_Q <- purrr::map_dbl(e_mix_diffusion, ~.x["variance_drift_Q"])
  v_sigma_tT <- purrr::map_dbl(e_mix_diffusion, ~.x["variance_diffusion"])
  last_1 <- purrr::map_dbl(e_mix_diffusion, ~.x["last_1"])
  last_0 <- purrr::map_dbl(e_mix_diffusion, ~.x["last_2"])
  # ************************************************************
  lambda_R <- model_Rt$lambda
  # Total mixture expectations 
  M_Y1 <- e_Yt + e_mu_tT_1 + e_sigma_tT_1 * lambda_R
  M_Y0 <- e_Yt + e_mu_tT_0 + e_sigma_tT_0 * lambda_R
  # Total variance in t, T-1
  common_variance <- v_mu_tT_P + v_sigma_tT + v_mu_tT_Q * lambda_R^2
  # Total mixture variances
  S_Y1 <- sqrt(common_variance + last_1)
  S_Y0 <- sqrt(common_variance + last_0)
  # Mixuture probability 
  p1 <- model_Rt$prob(t_hor = t_hor)
  p_T <- cbind(p1, 1-p1)
  # Realized variance for each component between T-1 and T
  sd_Yt <- sqrt(purrr::map_dbl(e_mix_diffusion, ~.x["common_variance"]))
  # ************************************************************
  # Pdf of Rt and Yt 
  pdf_Yt <- function(M_Y, S_Y, p_T){
    force(M_Y); force(S_Y); force(p_T)
    function(x){
      dmixnorm(x, M_Y, S_Y, p_T)
    }
  } 
  pdf_Rt <- function(Ct, alpha, beta, M_Y, S_Y, p_T){
    pdf_Y <- function(x){dmixnorm(x, M_Y, S_Y, p_T)}
    force(Ct); force(alpha); force(beta)
    function(x){
      dsolarGHI(x, Ct, alpha, beta, pdf_Y, link = "invgumbel")
    }
  } 
  # Cdf of Rt and Yt 
  cdf_Yt <- function(M_Y, S_Y, p_T){
    force(M_Y); force(S_Y); force(p_T)
    function(x){
      pmixnorm(x, M_Y, S_Y, p_T)
    }
  } 
  cdf_Rt <- function(Ct, alpha, beta, M_Y, S_Y, p_T){
    cdf_Y <- function(x){pmixnorm(x, M_Y, S_Y, p_T)}
    force(Ct); force(alpha); force(beta)
    function(x){
      psolarGHI(x, Ct, alpha, beta, cdf_Y, link = "invgumbel")
    }
  } 
  # ************************************************************
  # Clear-sky at time t_hor 
  C_T <- model_Rt$Ct(t_hor)
  # Seasonal mean of R at t_hor 
  GHI_bar <- model_Rt$model$spec$transform$iRY(YT_bar, C_T)
  # Generate pdf and cdf 
  pdf_R <- purrr::map(1:length(t_hor), ~pdf_Rt(C_T[.x], alpha, beta, c(M_Y1[.x], M_Y0[.x]), c(S_Y1[.x], S_Y0[.x]), unlist(p_T[.x, ])))
  cdf_R <- purrr::map(1:length(t_hor), ~cdf_Rt(C_T[.x], alpha, beta, c(M_Y1[.x], M_Y0[.x]), c(S_Y1[.x], S_Y0[.x]), unlist(p_T[.x, ])))
  # Pdf of Yt 
  pdf_Y <- purrr::map(1:length(t_hor), ~pdf_Yt(c(M_Y1[.x], M_Y0[.x]), c(S_Y1[.x], S_Y0[.x]), unlist(p_T[.x, ])))
  cdf_Y <- purrr::map(1:length(t_hor), ~cdf_Yt(c(M_Y1[.x], M_Y0[.x]), c(S_Y1[.x], S_Y0[.x]), unlist(p_T[.x, ])))
  
  # Full dataset 
  dplyr::tibble(
    date = t_hor,
    Month = lubridate::month(date),
    e_Yt = e_Yt, 
    sd_Yt = sd_Yt,
    M_Y1 = M_Y1,
    M_Y0 = M_Y0,
    S_Y1 = S_Y1,
    S_Y0 = S_Y0,
    p1 = p1,
    Ct = C_T,
    GHI_bar = GHI_bar,
    alpha = alpha, 
    beta = beta, 
    RT_min = Ct*(1-alpha-beta),
    RT_max = Ct*(1-alpha),
    pdf_R = pdf_R,
    cdf_R = cdf_R,
    pdf_Y = pdf_Y,
    cdf_Y = cdf_Y,
    e_mu_tT_1 = e_mu_tT_1, 
    e_mu_tT_0 = e_mu_tT_0, 
    e_sigma_tT_1 = e_sigma_tT_1, 
    e_sigma_tT_0 = e_sigma_tT_0, 
    v_mu_tT_P = v_mu_tT_P,
    v_sigma_tT = v_sigma_tT,
    v_mu_tT_Q = v_mu_tT_Q,
    last_1 = last_1, 
    last_0 = last_0
  )
}
