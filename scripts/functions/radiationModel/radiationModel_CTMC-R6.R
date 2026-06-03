#' Radiation Model CTMC
#'
#' @rdname radiationModel_CTMC
#' @name radiationModel_CTMC
#' @description
#' R6 wrapper extending `radiationModel` with a fitted monthly two-state CTMC
#' for radiation residual regimes.
#'
#' @note Version 1.0.0
#' @export
radiationModel_CTMC <- R6::R6Class("radiationModel_CTMC",
                                 inherit = radiationModel,
                                 # ====================================================================================================== #
                                 #                                             Public slots
                                 # ====================================================================================================== #
                                 public = list(
                                   #' @description
                                   #' Initialize a `radiationModel_CTMC` object
                                   #' @param model_Rt A fitted `radiationModel` object.
                                   #' @param DTMC Fitted monthly two-state DTMC used to construct CTMC generators.
                                   #' @return Initializes the object in place.
                                   initialize = function(model_Rt, DTMC){
                                     # Store the discrete time model 
                                     private$..model <- model_Rt$model$clone(TRUE)
                                     # Store mean reversion parameter
                                     self$theta <- model_Rt$theta
                                     # Store seasonal function  
                                     private$..seasonal_variance <- model_Rt$seasonal_variance$clone(TRUE)
                                     # Integral functions
                                     reparam <- private$..seasonal_variance$extra_params$reparam
                                     private$..integral_variance <- integral_sigma2_formula(self$theta, reparam$gamma, omega = self$seasonal_variance$omega)
                                     private$..integral_expectation <- integral_sigma_numeric(self$theta, reparam$c_, omega = self$seasonal_variance$omega)
                                     # DTMC Model in CTMC 
                                     CTMC <- ctmc_from_dtcm(DTMC, delta = 1)
                                     private[["..CTMC"]] <- CTMC
                                     # New means 
                                     mix.means <- as.matrix(dplyr::bind_rows(CTMC$params$mu))
                                     # New std deviations
                                     mix.sd <- as.matrix(dplyr::bind_rows(CTMC$params$sig))
                                     private$..model$spec$mixture.model$update(means = mix.means, sd = mix.sd)
                                   },
                                   #' @description
                                   #' Forecast the probability of the first DTMC state.
                                   #'
                                   #' @param t_hor Date or character scalar. Horizon date.
                                   #' @param t_now Date or character scalar. Conditioning date.
                                   #' @param p0 Optional numeric length-two initial probability.
                                   #'
                                   #' @return Numeric scalar. Forecast probability of state 1 at `t_hor`.
                                   prob = function(t_hor, t_now, p0){
                                     t_now <- as.Date(t_now)
                                     t_hor <- as.Date(t_hor)
                                     if (missing(p0)) {
                                       alpha0 <- self$DTMC$alpha
                                       dates <- self$DTMC$data$date
                                       p0 <- unlist(alpha0[dates == t_now])
                                     }
                                     Qm <- self$CTMC$params$Qm 
                                     bounds <- create_bounds(t_now, t_hor, Qm)
                                     drop(p0 %*% Phi_C(bounds$n[1], bounds$tau[1], bounds)[[1]])[1]
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
                                   ..CTMC = NA
                                 ),
                                 # ====================================================================================================== #
                                 #                                             Active slots
                                 # ====================================================================================================== #
                                 active = list(
                                   #' @field DTMC Fitted monthly two-state DTMC object.
                                   CTMC = function(){
                                     private$..CTMC
                                   }
                                 )
)


#' Fit the Continuous-Time Radiation CTMC 
#'
#' Fit the monthly two-state CTMC on the training sample and then recompute
#' filtered probabilities on the full dataset without updating parameters.
#'
#' @param model Fitted discrete radiation model.
#' @param p0 Numeric length-two initial state probability.
#' @param maxit Integer. Maximum EM iterations.
#' @param tol Numeric. Log-likelihood convergence tolerance.
#'
#' @return Fitted CTMC list with parameters, filtered probabilities, and data.
#' @keywords internal
radiationModel_CTMC_fit <- function(model, p0 = c(0.5, 0.5), maxit = 1000, tol = 0.01){
  # Extract mixture model for starting parameters 
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
  params[["a"]] <- model$spec$seasonal.mean$coefficients # (not optimized)
  params[["theta"]] <- -log(model$spec$mean.model$phi) # (not optimized)
  # Reparametrization to ensure positivity (not optimized)
  params[["b"]] <- seasonalModel_params_to_zeta(model$spec$seasonal.variance$coefficients)
  # ***************************************************
  data_train <- dplyr::filter(model$data, isTrain)
  Y <- data_train$Yt
  dates <- data_train$date
  weights = data_train$weights
  # Prepare data
  data <- dtmc_monthly_preprocess(Y, dates, weights, params, eps = 1e-8)
  # EM-recursion 
  loop_condition <- TRUE
  iter <- 0
  loglik <- -1000000
  while(loop_condition){
    # EM-step 
    em <- dtmc_monthly_EM(data = data$data, 
                          params = params, 
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
  # Full not update on whole data
  data <- model$data
  Y <- data$Yt
  dates <- data$date
  weights = data$weights
  # Prepare data
  preproc <- dtmc_monthly_preprocess(Y, dates, params = params, eps = 1e-8)
  # EM-step 
  em <- dtmc_monthly_EM(preproc$data, preproc$params, 
                        update_emissions = FALSE, 
                        update_transitions = FALSE,
                        eps = 1e-8)
  em
}

#' Approximate CTMC Density for the Radiation DTMC
#'
#' Numerically propagate the transformed radiation density over a grid using the
#' recovered CTMC generators and a Gaussian transition kernel.
#'
#' @param t_now Date or character scalar. Conditioning date.
#' @param t_hor Date or character scalar. Horizon date.
#' @param model_Rt A `radiationModel_CTMC` object.
#' @param R0 Optional numeric initial GHI at `t_now`.
#' @param y_grid Optional numeric grid for transformed radiation.
#' @param n_grid Integer. Grid size used when `y_grid` is not supplied.
#' @param dt Numeric. CTMC/kernel time step in days.
#' @param y_sd_mult Numeric. Width multiplier for automatic grid construction.
#' @param normalize Logical. If `TRUE`, renormalize density mass after each
#'   time step.
#'
#' @return List containing state densities, total density/CDF closures, grid,
#'   initial condition, and diagnostics.
#' @export
radiationModel_CTMC_density <- function(
    t_now,
    t_hor,
    model_Rt,
    R0=NULL,
    y_grid = NULL,
    n_grid = 1201,
    dt = 1/30,
    y_sd_mult = 8,
    normalize = TRUE
) {
  stopifnot(any(class(model_Rt) %in% "radiationModel_CTMC"))
  # Totay date
  t_now <- as.Date(t_now)
  # Horizon date
  t_hor <- as.Date(t_hor)
  # Days to maturity 
  tau <- as.numeric(difftime(t_hor, t_now, units = "days"))
  if (tau <= 0) stop("t_hor must be after t_now.")
  # Mean-reversion 
  theta <- model_Rt$theta
  # CTMC model 
  CTMC <- model_Rt$CTMC
  mu <-  CTMC$params$mu
  sig <- CTMC$params$sig
  if (!is.null(CTMC$params$Qm)){
    Q <- CTMC$params$Qm
  } else {
    Q <- transition_list_to_generator_2state(CTMC$params$Pm)
  }
  # Initial value of radiation 
  if (is.null(R0)) {
    R0 <- dplyr::filter(model_Rt$model$data, date == t_now)$GHI
  }
  # Initial value of Y 
  C0 <- model_Rt$Ct(t_now)
  Y0 <- model_Rt$model$spec$transform$RY(R0, C0)
  # Initial probability 
  p0 <- as.numeric(CTMC$alpha[CTMC$data$date == t_now, ])
  p0 <- p0 / sum(p0)
  # Time grid 
  s_grid <- seq(0, tau, by = dt)
  if (tail(s_grid, 1) < tau) s_grid <- c(s_grid, tau)
  
  n0 <- number_of_day(t_now)
  # Seasonal functions
  Ybar_fun <- function(s) model_Rt$Yt_bar(n0 + s)
  sigbar_fun <- function(s) model_Rt$sigma_bar(n0 + s)
  # Derivative of seasonal mean
  dYbar_fun <- function(s) model_Rt$dYt_bar(n0 + s)
  # dYbar_fun <- function(s, eps = 1e-4) (Ybar_fun(s + eps) - Ybar_fun(s - eps)) / (2 * eps)
  month_fun <- function(s) {
    lubridate::month(t_now + floor(s + 1e-12))
  }
  if (is.null(y_grid)) {
    s_tmp <- seq(0, tau, length.out = 200)
    Ybar_vals <- Ybar_fun(s_tmp)
    sigbar_vals <- sigbar_fun(s_tmp)
    max_sig <- max(unlist(sig), na.rm = TRUE)
    approx_sd <- sqrt(sum((sigbar_vals * max_sig)^2) * tau / length(s_tmp))
    y_grid <- seq(
      min(Y0, Ybar_vals, na.rm = TRUE) - y_sd_mult * approx_sd,
      max(Y0, Ybar_vals, na.rm = TRUE) + y_sd_mult * approx_sd,
      length.out = n_grid
    )
  }
  #y_grid  <- seq(-20, 20, length.out = 2000)
  dy <- y_grid[2] - y_grid[1]
  N <- length(y_grid)
  
  f <- matrix(0, nrow = N, ncol = 2)
  idx0 <- which.min(abs(y_grid - Y0))
  f[idx0, 1] <- p0[1] / dy
  f[idx0, 2] <- p0[2] / dy
  
  for (k in seq_len(length(s_grid) - 1)) {
    s0 <- s_grid[k]
    s1 <- s_grid[k + 1]
    ds <- s1 - s0
    m <- month_fun(s0)
    P <- matrix_exponential(ds, Q[[m]])
    Ybar0 <- Ybar_fun(s0)
    Ybar1 <- Ybar_fun(s1)
    sigbar1 <- sigbar_fun(s1)
    f_new <- matrix(0, nrow = N, ncol = 2)
    for (i in 1:2) {
      
      g_prev <- P[1, i] * f[, 1] + P[2, i] * f[, 2]
      
      mean_next <- Ybar1 +
        exp(-theta * ds) * (y_grid - Ybar0) +
        sigbar1 * mu[[m]][i] * ds
      
      sd_next <- sigbar1 * sig[[m]][i] * sqrt(ds)
      
      # Kernel matrix: rows y_next, columns y_prev
      K <- outer(
        y_grid,
        mean_next,
        function(y_next, mean_col) {
          dnorm(y_next, mean = mean_col, sd = sd_next)
        }
      )
      f_new[, i] <- as.numeric(K %*% (g_prev * dy))
    }
    
    if (normalize) {
      mass <- sum(f_new) * dy
      if (is.finite(mass) && mass > 0) f_new <- f_new / mass
    }
    
    f <- f_new
  }
  
  fY <- rowSums(f)
  cdfY <- cumsum(fY) * dy
  cdfY <- pmin(pmax(cdfY, 0), 1)
  
  list(
    y_grid = y_grid,
    f_state = f,
    f_Y = fY,
    pdf_Y = approxfun(y_grid, fY, yleft = 0, yright = 0, rule = 2),
    cdf_Y = approxfun(y_grid, cdfY, yleft = 0, yright = 1, rule = 2),
    mass = sum(fY) * dy,
    Y0 = Y0,
    p0 = p0,
    t_now = t_now,
    t_hor = t_hor,
    dt = dt
  )
}


#' Moment-Matched Radiation Distribution for the CTMC Model
#'
#' Compute conditional two-Gaussian moment approximations for transformed
#' radiation under the monthly CTMC dynamics and return density/CDF closures
#' for transformed radiation and physical GHI.
#'
#' @param t_now Date or character scalar. Conditioning date.
#' @param t_hor Date or character vector. Horizon date(s).
#' @param model_Rt A `radiationModel_CTMC` object.
#' @param R0 Optional numeric initial GHI at `t_now`. If missing, it is read
#'   from `model_Rt_DTMC$model$data`.
#' @param maxEval Integer. Maximum number of numerical integrand evaluations.
#' @param tol Numeric. Numerical integration tolerance.
#'
#' @return A tibble with horizon dates, component means/standard deviations,
#'   state probabilities, transformed and GHI density/CDF closures, and
#'   diagnostic moment components.
#' @export
radiationModel_CTMC_moments <- function(t_now, t_hor, model_Rt, R0 = NULL, maxEval = 50000, tol = 0.0001){
 
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
  # Model's parameters 
  # Radiation at time t_now 
  if (is.null(R0)) {
    R0 <- dplyr::filter(model_Rt$model$data, date == t_now[1])$GHI  
  }
  # Transform bounds 
  alpha = model_Rt$model$spec$transform$alpha
  beta = model_Rt$model$spec$transform$beta
  # Mean reversion parameter
  theta <- model_Rt$theta
  # Seasonal variance 
  sigma_bar <- model_Rt$seasonal_variance$extra_params$seasonal_function
  # Extract DTMC params
  CTMC <- model_Rt$CTMC
  # Mean parameters
  mu <- CTMC$params$mu
  # Mean parameters
  sd <- CTMC$params$sig
  # Transition matrix. v2 intentionally recovers the CTMC generator from the
  # fitted one-step transition matrices instead of trusting legacy Qm.
  if (!is.null(CTMC$params$Qm)){
    Q <- CTMC$params$Qm
  } else {
    Q <- transition_list_to_generator_2state(CTMC$params$Pm)
  }
  # Extract initial probabilities
  p0 <- purrr::map_df(t_now, ~setNames(CTMC$alpha[CTMC$data$date == .x,], c("alpha1", "alpha2")))
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
  # Pre compute the list of bounds
  bounds_list <- purrr::map2(t_now, t_hor, ~create_bounds(.x, .y, Q))
  # Forecast probability at maturity
  p_T <- purrr::map2(bounds_list, 1:nrow(p0), ~c(unlist(p0[.y,]) %*% Phi_C(.x$n[1], .x$tau[1], .x)[[1]]))
  # Conditional moments 
  E_cond <- purrr::map(1:nrow(p0), ~integral_E_cond_fast(bounds_list[[.x]], unlist(p0[.x,]), mu, sd, theta, sigma_bar, p_T = p_T[[.x]], t0 = c(0), T0 = c(0), maxEval = maxEval, tol = tol))
  E_mu_su_cond <- purrr::map(1:nrow(p0), ~integral_E_mu_su_cond_fast(bounds_list[[.x]], unlist(p0[.x,]), mu, theta, sigma_bar, p_T[[.x]], c(0,0), c(0, 0), maxEval, tol))
  # Probability at maturity B_T = 1
  p1 <- purrr::map_dbl(p_T, ~.x[1])
  # ************************************************************
  # Components for P-moments 
  # Conditional expected values drift mu(t,T)
  M_mu_1 <- purrr::map_dbl(E_cond, ~.x[1])
  M_mu_0 <- purrr::map_dbl(E_cond, ~.x[2])
  M_mu   <- M_mu_1 * p1 + M_mu_0 * (1 - p1)
  # Conditional second moment diffusion sigma(t,T)
  S2_sigma_1 <- purrr::map_dbl(E_cond, ~.x[3]) 
  S2_sigma_0 <- purrr::map_dbl(E_cond, ~.x[4]) 
  S2_sigma   <- S2_sigma_1 * p1 + S2_sigma_0 * (1 - p1) 
  # Cross integral mu(t,T)
  S2_mu_1 <- (purrr::map_dbl(E_mu_su_cond, ~.x[1]) - M_mu_1^2)
  S2_mu_0 <- (purrr::map_dbl(E_mu_su_cond, ~.x[2]) - M_mu_0^2)
  # Total variance mu(t, T)
  S2_mu  <- (p1 * S2_mu_1 + (1 - p1) * S2_mu_0) + (M_mu_1 - M_mu_0)^2 * p1 * (1 - p1)
  # Components for Q-moments 
  # Conditional expected values drift gamma(t,T)
  M_gamma_1 <- purrr::map_dbl(E_cond, ~.x[5])
  M_gamma_0 <- purrr::map_dbl(E_cond, ~.x[6])
  M_gamma <- M_gamma_1 * p1 + M_gamma_0 * (1 - p1) 
  # ************************************************************
  # Q-measure variance adjustment for lambda_R != 0 is not implemented here.
  # Moments of YT
  lambda_R <- model_Rt$lambda
  # Expectations
  M_Y1 <- e_Yt + M_mu_1 + M_gamma_1 * lambda_R
  M_Y0 <- e_Yt + M_mu_0 + M_gamma_0 * lambda_R
  # Std. deviations
  S2_Y <- S2_sigma + S2_mu
  S_Y1 <- sqrt(S2_mu_1 + S2_sigma_1)
  S_Y0 <- sqrt(S2_mu_0 + S2_sigma_0)
  # ************************************************************
  # Clear-sky at time t_hor 
  C_T <- model_Rt$Ct(t_hor)
  # Seasonal mean of R at t_hor 
  GHI_bar <- model_Rt$model$spec$transform$iRY(YT_bar, C_T)
  # ************************************************************
  # Full dataset 
  dplyr::tibble(
    date = t_hor,
    Month = lubridate::month(t_hor),
    e_Yt = e_Yt, 
    sd_Yt = sqrt(S2_Y),
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
    M_mu_1 = M_mu_1, 
    M_mu_0 = M_mu_0, 
    S2_sigma_1 = S2_sigma_1, 
    S2_sigma_0 = S2_sigma_0, 
    M_gamma_1 = M_gamma_1, 
    M_gamma_0 = M_gamma_0, 
    S2_mu_1 = S2_mu_1, 
    S2_mu_0 = S2_mu_0
  )
}

#' Bound Two-Gaussian Approximation Error for CTMC Radiation Density
#'
#' Compare the two-Gaussian approximation with a numerical CTMC density over one
#' horizon and compute Fourier, CDF, and moment discrepancies.
#'
#' @param t_now Date or character scalar. Conditioning date.
#' @param h Integer or numeric scalar. Horizon in days.
#' @param model_Rt A `radiationModel_CTMC` object.
#' @param n_grid Integer. Number of grid points for the transformed radiation
#'   density.
#' @param df Numeric. Grid spacing used for the numerical density benchmark.
#'
#' @return A tibble with approximation bounds, direct density diagnostics, and
#'   moment discrepancies.
#' @keywords internal
radiationModel_CTMC_density_bound <- function(t_now, h = 1, model_Rt, n_grid = 2001, df = 0.1){

    # Trapezoidal integration weights.
  trapz_weights <- function(x) {
    if (length(x) < 2) {
      stop("Need at least two quadrature nodes.")
    }
    dx <- diff(x)
    if (any(dx <= 0)) {
      stop("Quadrature grid must be strictly increasing.")
    }
    w <- numeric(length(x))
    w[1] <- dx[1] / 2
    w[length(w)] <- dx[length(dx)] / 2
    if (length(x) > 2) {
      w[2:(length(x) - 1)] <- (head(dx, -1) + tail(dx, -1)) / 2
    }
    w
  }
  # Build the positive frequency grid for Esseen-type integration.
  make_u_grid <- function(U, du, lower = lower_u) {
    if (!is.finite(U) || U <= 0) {
      stop("U must be positive.")
    }
    if (!is.finite(du) || du <= 0) {
      stop("du must be positive.")
    }
    if (!is.finite(lower) || lower <= 0 || lower >= U) {
      stop("lower_u must be in (0, U).")
    }
    u <- c(lower, seq(du, U, by = du))
    if (tail(u, 1) < U) {
      u <- c(u, U)
    }
    sort(unique(u))
  }
  # Two-Gaussian moment-matched characteristic function.
  two_gm_phi <- function(u, row) {
    p1 <- as.numeric(row$p1)
    M_Y1 <- as.numeric(row$M_Y1)
    M_Y0 <- as.numeric(row$M_Y0)
    S_Y1 <- as.numeric(row$S_Y1)
    S_Y0 <- as.numeric(row$S_Y0)
    p1 * exp(1i * u * M_Y1 - 0.5 * u^2 * S_Y1^2) + (1 - p1) * exp(1i * u * M_Y0 - 0.5 * u^2 * row$S_Y0^2)
  }
  # Build an automatic y-grid from the two-Gaussian approximation.
  default_y_grid <- function(row, n_grid) {
    p1 <- as.numeric(row$p1)
    M <- c(as.numeric(row$M_Y1), as.numeric(row$M_Y0))
    S <- c(as.numeric(row$S_Y1), as.numeric(row$S_Y0))
    pi_vec <- c(p1, 1 - p1)
    mix_mean <- sum(pi_vec * M)
    mix_var <- sum(pi_vec * (S^2 + M^2)) - mix_mean^2
    mix_sd <- sqrt(max(mix_var, 1e-12))
    lo <- min(M - 8 * S, mix_mean - 8 * mix_sd, na.rm = TRUE)
    hi <- max(M + 8 * S, mix_mean + 8 * mix_sd, na.rm = TRUE)
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
      lo <- mix_mean - 10
      hi <- mix_mean + 10
    }
    seq(lo, hi, length.out = n_grid)
  }
  kl_density <- function(y_grid, f_true, f_approx, eps = 1e-300) {
      w <- trapz_weights(y_grid)
      f_true <- pmax(as.numeric(f_true), eps)
      f_approx <- pmax(as.numeric(f_approx), eps)
      # Optional renormalization on the grid
      f_true <- f_true / sum(w * f_true)
      f_approx <- f_approx / sum(w * f_approx)
      sum(w * f_true * log(f_true / f_approx))
  }
  js_density <- function(y_grid, f, g, eps = 1e-300) {
      w <- trapz_weights(y_grid)

      f <- pmax(as.numeric(f), eps)
      g <- pmax(as.numeric(g), eps)

      f <- f / sum(w * f)
      g <- g / sum(w * g)

      m <- 0.5 * (f + g)

      0.5 * sum(w * f * log(f / m)) +
        0.5 * sum(w * g * log(g / m))
  }

  # Horizon date
  t_hor <- as.Date(t_now) + h
  # Radiation at time t_now 
  R0 <- lookup_R0(model_Rt, t_now)
  # CTMC moments and 2GM approximation 
  row <- radiationMoments(
    t_now = t_now,
    t_hor = t_hor,
    model_Rt = model_Rt,
    R0 = R0
  )[1,]
  y_grid <- default_y_grid(row, n_grid)
  dens <- radiationModel_CTMC_density(
    t_now = t_now,
    t_hor = t_hor,
    model_Rt = model_Rt,
    R0 = R0,
    y_grid = y_grid,
    n_grid = n_grid,
    dt = df,
    normalize = TRUE
  )
  # Integration weights
  u_bound <- make_u_grid(U_bound, du_bound)
  w_bound <- trapz_weights(u_bound)
  # Distributions of Y
  F_true <- pmin(pmax(dens$cdf_Y(y_grid), 0), 1)
  F_2gm <- pmixnorm(y_grid, c(row$M_Y1, row$M_Y0), c(row$S_Y1, row$S_Y0), c(row$p1, 1-row$p1))
  # Characteristic functions 
  phi_true_bound <- phi_from_density(u_bound, dens$y_grid, dens$f_Y)
  phi_2gm_bound <- two_gm_phi(u_bound, row)
  # Integral bound 
  integral_part <- (2 / base::pi) * sum(w_bound * abs((phi_true_bound - phi_2gm_bound) / u_bound)) 
  p1 <- as.numeric(row$p1)
  # Remainder 
  pi_vec <- c(p1, 1 - p1)
  V <- c(as.numeric(row$S_Y1)^2, as.numeric(row$S_Y0)^2)
  keep <- pi_vec > 1e-12 & is.finite(V) & V > 0
  M_bound <- sum(pi_vec[keep] / sqrt(2 * base::pi * V[keep]))
  remainder <- 24 * M_bound / (pi * U_bound)
  # Total bound 
  bound <- integral_part + remainder
  # Distributon of R
  r_grid <- seq(row$RT_min, row$RT_max, length.out = n_grid)
  cdf_R_true <- psolarGHI(r_grid, row$Ct, row$alpha, row$beta, dens$cdf_Y) 
  cdf_R_2gm  <- psolarGHI(r_grid, row$Ct, row$alpha, row$beta, row$cdf_Y[[1]])
  # Empirical error on Y 
  max_abs_err_y <- max(abs(F_true - F_2gm))
  max_abs_err_r <- max(abs(cdf_R_true - cdf_R_2gm))
  mean_abs_err_y <- mean(abs(F_true - F_2gm))
  mean_abs_err_r <- mean(abs(cdf_R_true - cdf_R_2gm))

  # Densities of Y 
  f_true <- dens$f_Y
  f_2gm <- dmixnorm(y_grid, c(row$M_Y1, row$M_Y0), c(row$S_Y1, row$S_Y0), c(row$p1, 1-row$p1))
  # Density of R
  pdf_R_true <- dsolarGHI(r_grid, row$Ct, row$alpha, row$beta, dens$pdf_Y) 
  pdf_R_2gm  <- dsolarGHI(r_grid, row$Ct, row$alpha, row$beta, row$pdf_Y[[1]])
  # JS distances (Y,R)
  js_dist_y <- js_density(y_grid, f_true, f_2gm)
  js_dist_r <- js_density(r_grid, pdf_R_true, pdf_R_2gm)
  # KL distances (Y)
  KL_true_to_2gm_y <- kl_density(y_grid, f_true, f_2gm)
  KL_2gm_to_true_y <- kl_density(y_grid, f_2gm, f_true)
  # KL distances (R)
  KL_true_to_2gm_r <- kl_density(r_grid, pdf_R_true, pdf_R_2gm)
  KL_2gm_to_true_r <- kl_density(r_grid, pdf_R_2gm, pdf_R_true)
  
  # Data for the CDF 
  tibble(
    t_now = t_now, 
    t_hor = t_hor, 
    h = h,
    integral_part = integral_part, 
    remainder = remainder, 
    M_bound = M_bound,
    bound = bound, 
    max_abs_err = max_abs_err_y,
    mean_abs_err = mean_abs_err_y,
    JS_Y = js_dist_y, 
    JS_R = js_dist_r,
    KL_true_y = KL_true_to_2gm_y,
    KL_GM_y = KL_2gm_to_true_y,
    KL_true_r = KL_true_to_2gm_r,
    KL_GM_r = KL_2gm_to_true_r,
    p_diff = max(c(row$p1, 1-row$p1)),
    M_diff = row$M_Y1 - row$M_Y0
  ) 
}
