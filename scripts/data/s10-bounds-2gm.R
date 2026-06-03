# ---
#' @description
#' Generate two-Gaussian CDF, moment, and error-bound diagnostics for CTMC radiation models.
#'
#' @section `appendix` and `supplementary-material`
#' @label `bounds-2gm`
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse
#'
#' @arguments
#'   - param[1] (nyear): training year, e.g. "2022".
#'   - param[2] (h): forecast horizon in days.
#'   - param[3] (save_data): write generated data ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/data/s10-bounds-2gm.R "2022" "1" "TRUE"
#' for h in 1 2 3 5 10 15 30; do Rscript scripts/data/s10-bounds-2gm.R "2022" "$h" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData` for each configured location
#'
#' @outputs
#'   - `data/models/radiation/P/{place}/bounds/{nyear}/bounds_2gm-{h}.RData` (CDF, moment, and error diagnostics)
#'
#' @depends
#'   - `scripts/data/s2b-models-radiation-P-CTMC-place.R`
#'
#' @tags
#'   - data
#'   - appendix
#'   - supplementary-material
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "data", "s10-bounds-2gm.R"))
source("scripts/functions/radiationModel/radiationModel-CTMC-density-C-wrappers.R")
source("scripts/functions/radiationModel/ctmc-integrals-C-wrappers.R")
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference year
nyear <- "2022"
# Horizon
h <- 1
# Save the output
save_data <- TRUE
# ***************************** Fixed Arguments ********************************
# Reference places
places <- c("Bologna", "Palermo", "Roma")
# Frequency-domain truncation for the CDF bound.
U_bound <- 2000
# Frequency-domain truncation kept in the output for comparison with older runs.
U_inv <- 2000
# Number of points in the y-grid used to evaluate the true density/CDF.
n_grid <- 2001
# Quadrature grid refinements in the frequency domain.
du_bound <- 0.5
du_inv <- 0.5
# Lower integration bound used to avoid division by zero at u = 0.
lower_u <- 1e-6
# Numerical settings for production radiation moments.
moments_maxEval <- 50000
moments_tol <- 1e-6
# Resume controls.
max_rows <- Inf
resume <- TRUE
# Model's directory
dir_models_P <- outputs$dir$data$models$radiation$P
# Table labels
tbl_labels <- "tbl-model_Et"
# Table names
tbl_names <- stringr::str_replace_all(tbl_labels, "-", "_")
# ******************************************************************************
#                             Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly=TRUE)
if (!purrr::is_empty(args)) {
  # Reference year
  nyear <- ifelse(is.na(args[1]), nyear, args[1])
  # Reference horizon
  h <- ifelse(is.na(args[2]), h, args[2])
  # Save output 
  save_data <- ifelse(is.na(args[3]), save_data, args[3])
  print_script_args(nyear = nyear, h = h, save_data = save_data)
}
# ******************************************************************************
#                           Generating functions  
# ******************************************************************************
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
# Two-Gaussian moment-matched CDF.
two_gm_cdf <- function(y, row) {
  p1 <- as.numeric(row$p1)
  M_Y1 <- as.numeric(row$M_Y1)
  M_Y0 <- as.numeric(row$M_Y0)
  S_Y1 <- as.numeric(row$S_Y1)
  S_Y0 <- as.numeric(row$S_Y0)
  p1 * stats::pnorm(y, M_Y1, S_Y1) + (1 - p1) * stats::pnorm(y, M_Y0, S_Y0)
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
# Extract the observed initial radiation value, when available.
lookup_R0 <- function(model, t_now) {
  idx <- which(as.Date(model$model$data$date) == as.Date(t_now))
  if (length(idx) == 0) {
    return(NULL)
  }
  model$model$data$GHI[idx[1]]
}
# Numerical characteristic function obtained from the production density grid.
phi_from_density <- function(u, y, f) {
  wy <- trapz_weights(y)
  vapply(u, function(ui) sum(wy * exp(1i * ui * y) * f), complex(1))
}
#t_now  <- "2022-01-08"
#h <- 1
#df <- 0.05
# Density bounds
radiationModel_CTMC_density_bound <- function(t_now, h = 1, model_Rt, n_grid = 2001, df = 0.1){
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
  r_grid <- seq(row$RT_min, row$RT_max, length.out = n_grid)
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
  # Distributions 
  F_true <- pmin(pmax(dens$cdf_Y(y_grid), 0), 1)
  F_2gm <- two_gm_cdf(y_grid, row)
  # Distributon of solar radiation
  cdf_R_true <- psolarGHI(r_grid, row$Ct, row$alpha, row$beta, dens$cdf_Y) 
  cdf_R_2gm  <- psolarGHI(r_grid, row$Ct, row$alpha, row$beta, row$cdf_Y[[1]])
 
  # Theoric bound between the CDF
  u_bound <- make_u_grid(U_bound, du_bound)
  w_bound <- trapz_weights(u_bound)
  # Characteristic functions 
  phi_true_bound <- phi_from_density(u_bound, dens$y_grid, dens$f_Y)
  phi_2gm_bound <- two_gm_phi(u_bound, row)
  integral_part <- (2 / pi) *
    sum(w_bound * abs((phi_true_bound - phi_2gm_bound) / u_bound))
  p1 <- as.numeric(row$p1)
  pi_vec <- c(p1, 1 - p1)
  V <- c(as.numeric(row$S_Y1)^2, as.numeric(row$S_Y0)^2)
  keep <- pi_vec > 1e-12 & is.finite(V) & V > 0
  M_bound <- sum(pi_vec[keep] / sqrt(2 * base::pi * V[keep]))
  remainder <- 24 * M_bound / (pi * U_bound)
  # Final bound 
  bound <- integral_part + remainder
  # Empirical error 
  abs_error <- abs(F_true - F_2gm)
  max_idx <- which.max(abs_error)
  F_true_y_star <- F_true[max_idx]
  F_2gm_y_star <- F_2gm[max_idx]
  abs_err_y_star <- abs_error[max_idx]
  abs_err_r_star <- max(abs(cdf_R_true - cdf_R_2gm))

  # Densities 
  f_true <- dens$f_Y
  f_2gm <- dmixnorm(y_grid, c(row$M_Y1, row$M_Y0), c(row$S_Y1, row$S_Y0), c(row$p1, 1-row$p1))
  # Density of solar radiation
  pdf_R_true <- dsolarGHI(r_grid, row$Ct, row$alpha, row$beta, dens$pdf_Y) 
  pdf_R_2gm  <- dsolarGHI(r_grid, row$Ct, row$alpha, row$beta, row$pdf_Y[[1]])
  # JS distances
  js_dist_y <- js_density(y_grid, f_true, f_2gm)
  js_dist_r <- js_density(r_grid, pdf_R_true, pdf_R_2gm)
  # KL distances
  KL_true_to_2gm_y <- kl_density(y_grid, f_true, f_2gm)
  KL_2gm_to_true_y <- kl_density(y_grid, f_2gm, f_true)
  KL_true_to_2gm_r <- kl_density(r_grid, pdf_R_true, pdf_R_2gm)
  KL_2gm_to_true_r <- kl_density(r_grid, pdf_R_2gm, pdf_R_true)
  # Helper 
  fun_e_Y <- function(g_x = function(x) x, pdf_Y) integrate(function(x) g_x(x)*pdf_Y(x), lower = -Inf, upper = Inf, stop.on.error = FALSE)$value
  fun_e_R <- function(g_x = function(x) x, pdf_Y) integrate(function(x) g_x(x)*dsolarGHI(x, row$Ct, row$alpha, row$beta, pdf_Y), lower = row$RT_min, upper = row$RT_max, stop.on.error = FALSE)$value
  # Expectations of Y
  e_Y_true <- fun_e_Y(g_x = function(x) x, dens$pdf_Y)
  e_Y_2gm  <- fun_e_Y(g_x = function(x) x, row$pdf_Y[[1]])
  # Variance of Y
  v_Y_true <- fun_e_Y(g_x = function(x) x^2, dens$pdf_Y) - e_Y_true^2
  v_Y_2gm <- fun_e_Y(g_x = function(x) x^2, row$pdf_Y[[1]]) - e_Y_2gm^2
  # Expectations of R
  e_R_true <- fun_e_R(g_x = function(x) x, dens$pdf_Y)
  e_R_2gm  <- fun_e_R(g_x = function(x) x, row$pdf_Y[[1]])
  # Variance of R
  v_R_true <- fun_e_R(g_x = function(x) x^2, dens$pdf_Y)  - e_R_true^2
  v_R_2gm  <- fun_e_R(g_x = function(x) x^2, row$pdf_Y[[1]])  - e_R_2gm^2
  # Expectations of Gamma
  e_Gamma_true <- fun_e_R(g_x = function(x) solarOption_payoff(x, row$GHI_bar), dens$pdf_Y)
  e_Gamma_2gm  <- fun_e_R(g_x = function(x) solarOption_payoff(x, row$GHI_bar), row$pdf_Y[[1]])
  # Variance of Gamma
  v_Gamma_true <- fun_e_R(g_x = function(x) solarOption_payoff(x, row$GHI_bar)^2, dens$pdf_Y) - e_Gamma_true^2
  v_Gamma_2gm  <- fun_e_R(g_x = function(x) solarOption_payoff(x, row$GHI_bar)^2, row$pdf_Y[[1]])  - e_Gamma_2gm^2
  # Covariance Gamma, R
  C_R_Gamma_true <- fun_e_R(g_x = function(x) solarOption_payoff(x, row$GHI_bar)*x, dens$pdf_Y) - e_Gamma_true*e_R_true
  C_R_Gamma_2gm <- fun_e_R(g_x = function(x) solarOption_payoff(x, row$GHI_bar)*x, row$pdf_Y[[1]])  - e_Gamma_2gm*e_R_2gm
  
  # Data for the CDF 
  data_cdf <- tibble(
    t_now = t_now, 
    t_hor = t_hor, 
    h = h,
    integral_part = integral_part, 
    remainder = remainder, 
    bound = bound, 
    max_abs_err = abs_err_y_star,
    mean_abs_err = mean(abs_error),
    JS_Y = js_dist_y, 
    JS_R = js_dist_r,
    KL_true_y = KL_true_to_2gm_y,
    KL_GM_y = KL_2gm_to_true_y,
    KL_true_r = KL_true_to_2gm_r,
    KL_GM_r = KL_2gm_to_true_r,
    p_T = row$p1,
    M_diff = row$M_Y1 - row$M_Y0,
    F_true = F_true_y_star, 
    F_2gm = F_2gm_y_star) 
  # Data for the moments  
  data_mom <- tibble(
    t_now = t_now, 
    t_hor = t_hor, 
    h = h,
    e_Y_true = e_Y_true, 
    e_Y_2gm = e_Y_2gm,
    v_Y_true = v_Y_true, 
    v_Y_2gm = v_Y_2gm,
    e_R_true = e_R_true, 
    e_R_2gm = e_R_2gm,
    v_R_true = v_R_true, 
    v_R_2gm = v_R_2gm,
    e_Gamma_true = e_Gamma_true, 
    e_Gamma_2gm = e_Gamma_2gm,
    v_Gamma_true = v_Gamma_true, 
    v_Gamma_2gm = v_Gamma_2gm,
    C_R_Gamma_true = C_R_Gamma_true, 
    C_R_Gamma_2gm = C_R_Gamma_2gm,
  )
  # Data for the errors  
  data_err <- tibble(
    t_now = t_now, 
    t_hor = t_hor, 
    h = h,
    e_Y = abs((e_Y_true - e_Y_2gm)/e_Y_true)*100,
    v_Y = abs((v_Y_true - v_Y_2gm)/v_Y_true)*100,
    e_R = abs((e_R_true - e_R_2gm)/e_R_true)*100,
    v_R = abs((v_R_true - v_R_2gm)/v_R_true)*100,
    e_Gamma = abs((e_Gamma_true - e_Gamma_2gm)/e_Gamma_true)*100,
    v_Gamma = abs((v_Gamma_true - v_Gamma_2gm)/v_Gamma_true)*100,
    C_R_Gamma = abs((C_R_Gamma_true - C_R_Gamma_2gm)/C_R_Gamma_true)*100,
  )

  list(
    cdf = data_cdf, 
    mom = data_mom,
    err = data_err
  )
}

radiationModel_CTMC_test_GM_approx <- function(t_now, h = 1, model_Rt, n_grid = 2001, df = 0.1){

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

  # Helper 
  fun_e_Y <- function(g_x = function(x) x, pdf_Y) integrate(function(x) g_x(x)*pdf_Y(x), lower = -Inf, upper = Inf, stop.on.error = FALSE)$value
  fun_e_R <- function(g_x = function(x) x, pdf_Y) integrate(function(x) g_x(x)*dsolarGHI(x, row$Ct, row$alpha, row$beta, pdf_Y), lower = row$RT_min, upper = row$RT_max, stop.on.error = FALSE)$value
  # Expectations of Y
  e_Y_true <- fun_e_Y(g_x = function(x) x, dens$pdf_Y)
  e_Y_2gm  <- fun_e_Y(g_x = function(x) x, row$pdf_Y[[1]])
  # Variance of Y
  v_Y_true <- fun_e_Y(g_x = function(x) x^2, dens$pdf_Y) - e_Y_true^2
  v_Y_2gm <- fun_e_Y(g_x = function(x) x^2, row$pdf_Y[[1]]) - e_Y_2gm^2
  # Expectations of R
  e_R_true <- fun_e_R(g_x = function(x) x, dens$pdf_Y)
  e_R_2gm  <- fun_e_R(g_x = function(x) x, row$pdf_Y[[1]])
  # Variance of R
  v_R_true <- fun_e_R(g_x = function(x) x^2, dens$pdf_Y)  - e_R_true^2
  v_R_2gm  <- fun_e_R(g_x = function(x) x^2, row$pdf_Y[[1]])  - e_R_2gm^2
  # Expectations of Gamma
  e_Gamma_true <- fun_e_R(g_x = function(x) solarOption_payoff(x, row$GHI_bar), dens$pdf_Y)
  e_Gamma_2gm  <- fun_e_R(g_x = function(x) solarOption_payoff(x, row$GHI_bar), row$pdf_Y[[1]])
  # Variance of Gamma
  v_Gamma_true <- fun_e_R(g_x = function(x) solarOption_payoff(x, row$GHI_bar)^2, dens$pdf_Y) - e_Gamma_true^2
  v_Gamma_2gm  <- fun_e_R(g_x = function(x) solarOption_payoff(x, row$GHI_bar)^2, row$pdf_Y[[1]])  - e_Gamma_2gm^2
  # Covariance Gamma, R
  C_R_Gamma_true <- fun_e_R(g_x = function(x) solarOption_payoff(x, row$GHI_bar)*x, dens$pdf_Y) - e_Gamma_true*e_R_true
  C_R_Gamma_2gm <- fun_e_R(g_x = function(x) solarOption_payoff(x, row$GHI_bar)*x, row$pdf_Y[[1]])  - e_Gamma_2gm*e_R_2gm
  
  # Data for the CDF 
  data_cdf <- tibble(
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
  # Data for the moments  
  data_mom <- tibble(
    t_now = t_now, 
    t_hor = t_hor, 
    h = h,
    e_Y_true = e_Y_true, 
    e_Y_2gm = e_Y_2gm,
    v_Y_true = v_Y_true, 
    v_Y_2gm = v_Y_2gm,
    e_R_true = e_R_true, 
    e_R_2gm = e_R_2gm,
    v_R_true = v_R_true, 
    v_R_2gm = v_R_2gm,
    e_Gamma_true = e_Gamma_true, 
    e_Gamma_2gm = e_Gamma_2gm,
    v_Gamma_true = v_Gamma_true, 
    v_Gamma_2gm = v_Gamma_2gm,
    C_R_Gamma_true = C_R_Gamma_true, 
    C_R_Gamma_2gm = C_R_Gamma_2gm,
  )
  # Data for the errors  
  data_err <- tibble(
    t_now = t_now, 
    t_hor = t_hor, 
    h = h,
    e_Y = abs((e_Y_true - e_Y_2gm)/e_Y_true)*100,
    v_Y = abs((v_Y_true - v_Y_2gm)/v_Y_true)*100,
    e_R = abs((e_R_true - e_R_2gm)/e_R_true)*100,
    v_R = abs((v_R_true - v_R_2gm)/v_R_true)*100,
    e_Gamma = abs((e_Gamma_true - e_Gamma_2gm)/e_Gamma_true)*100,
    v_Gamma = abs((v_Gamma_true - v_Gamma_2gm)/v_Gamma_true)*100,
    C_R_Gamma = abs((C_R_Gamma_true - C_R_Gamma_2gm)/C_R_Gamma_true)*100,
  )
  list(
    cdf = data_cdf, 
    mom = data_mom,
    err = data_err
  )
}

# Generate data for  model 
generate_bounds_test <- function(place = "Bologna", nyear = "2022", nyear_test = "2022", h = 1, df = 0.05, n_grid = 2001) {
  # Extract the model
  # Load HMM models
  load_data(file.path(dir_models_P, place), "radiation_models_CTMC")
  model_Rt <- radiation_models_CTMC[[nyear]]
  # Today date
  t_now <- as.Date(paste0(as.numeric(nyear_test)-1, "-12-31"))
  # Horizon date 
  t_hor <- as.Date(paste0(nyear_test, "-12-31"))
  # Sequence of dates 
  date_seq <- seq.Date(t_now, t_hor, 1)
  # Initialization 
  bounds_list <- moments_list <- errors_list <- list()
  err_e_R <- err_e_Gamma <- c()
  for(i in 1:length(date_seq)){
    bounds_new <- radiationModel_CTMC_test_GM_approx(date_seq[i], h =  h, model_Rt, n_grid = n_grid, df = df)
    bounds_list <- append(bounds_list, list(bounds_new$cdf))
    moments_list <- append(moments_list, list(bounds_new$mom))
    errors_list <- append(errors_list, list(bounds_new$err))
    err_e_R[i] <- bounds_new$err$e_R
    err_e_Gamma[i] <- bounds_new$err$e_Gamma
    print(paste0("(", i, "/", length(date_seq),"): ",  "Date: ", date_seq[i], " (h: ", h, ") ",
    " Error R (avg): ", round(mean(err_e_R), 3), "% ", " (", round(bounds_new$err$e_R, 3), "%)", " - ",
     "Error Gamma (avg): ", round(mean(err_e_Gamma), 3), "%", " (", round(bounds_new$err$e_Gamma, 3), "%)"))
  } 
  data_bounds <- bind_rows(bounds_list)
  data_moments <- bind_rows(moments_list)
  data_errors <- bind_rows(errors_list)
  list(
    place = place, 
    nyear = nyear,
    nyear_test = nyear_test, 
    h = h, 
    df = df, 
    n_grid = n_grid,
    bounds = data_bounds, 
    moments = data_moments,
    errors = data_errors
  )
}
# ******************************************************************************
#                               Generate data  
# ******************************************************************************
# Test year
nyear_test <- as.character(as.numeric(nyear)+1)
for(place in places){
  bounds_2gm <- generate_bounds_test(place = place, nyear = nyear, nyear_test = nyear_test, h = as.numeric(h), df = 0.05, n_grid = 2001) 
  if (save_data){
    dir_output <- file.path(dir_models_P, place, "bounds")
    make_new_directory(dir_output)
    dir_output <- file.path(dir_output, nyear)
    make_new_directory(dir_output)
    save_new_file(dir_output, file.name = paste0("bounds_2gm-", h), 
                  file.format = "RData", quiet = FALSE, bounds_2gm)
  }
  
  
}
