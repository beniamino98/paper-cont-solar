# ---
#' @description
#' Shared fixtures and independent reference calculations for CTMC/radiation
#' tests.
#'
#' @details
#' This helper is intentionally independent from the production implementation
#' where possible. It provides:
#' - small two-state CTMC generators and monthly parameter lists;
#' - mock `radiationModel_CTMC` objects for moment, density, and scenario tests;
#' - explicit transition-product and integral references used by function-level
#'   test files.
# ---
suppressPackageStartupMessages({
  library(testthat)
  library(solarr)
  library(dplyr)
  library(purrr)
  library(lubridate)
})

ctmc_find_repo_root <- function() {
  here <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "scripts", "functions", "radiationModel", "ctmc-integrals-C.R"))) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) {
      stop("Cannot find repository root containing scripts/functions/radiationModel/ctmc-integrals-C.R")
    }
    here <- parent
  }
}

ctmc_root <- ctmc_find_repo_root()
ctmc_functions_dir <- file.path(ctmc_root, "scripts", "functions")
ctmc_radiation_dir <- file.path(ctmc_functions_dir, "radiationModel")
source(file.path(ctmc_radiation_dir, "radiationModel-internals.R"))
source(file.path(ctmc_radiation_dir, "ctmc.R"))
source(file.path(ctmc_radiation_dir, "ctmc-integrals-C.R"))
source(file.path(ctmc_radiation_dir, "radiationModel-R6.R"))
source(file.path(ctmc_radiation_dir, "radiationModel_CTMC-R6.R"))
source(file.path(ctmc_radiation_dir, "radiationMoments.R"))
source(file.path(ctmc_radiation_dir, "scenarios_radiationModel_CT.R"))
source(file.path(ctmc_functions_dir, "scenarios_ER.R"))

ctmc_generator <- function(q12 = 0.08, q21 = 0.03) {
  matrix(c(-q12, q12,
           q21, -q21),
         nrow = 2, byrow = TRUE)
}

ctmc_month_list <- function(default, replacements = list()) {
  out <- replicate(12, default, simplify = FALSE)
  for (nm in names(replacements)) {
    out[[as.integer(nm)]] <- replacements[[nm]]
  }
  out
}

ctmc_mu_list <- function(default = c(mu1 = 0.35, mu0 = -0.15),
                         replacements = list()) {
  ctmc_month_list(default, replacements)
}

ctmc_sd_list <- function(default = c(sd1 = 0.7, sd0 = 1.1),
                         replacements = list()) {
  ctmc_month_list(default, replacements)
}

ctmc_sigma_bar_one <- function(s) {
  rep(1, length(s))
}

ctmc_sigma_bar_seasonal <- function(s) {
  1 + 0.1 * sin(2 * pi * s / 365)
}

ctmc_fixture <- function(t_now = "2022-01-15",
                         t_hor = "2022-01-20",
                         Q = ctmc_month_list(ctmc_generator()),
                         mu = ctmc_mu_list(),
                         sd = ctmc_sd_list(),
                         p0 = c(0.4, 0.6),
                         sigma_bar = ctmc_sigma_bar_one) {
  bounds <- create_bounds(t_now, t_hor, Q)
  list(
    t_now = as.Date(t_now),
    t_hor = as.Date(t_hor),
    bounds = bounds,
    Q = Q,
    mu = mu,
    sd = sd,
    p0 = p0,
    sigma_bar = sigma_bar,
    t_init = bounds$n[1],
    t_end = bounds$tau[1]
  )
}

ctmc_month_boundary_fixture <- function(t_now, t_hor) {
  Q <- ctmc_month_list(
    ctmc_generator(0.04, 0.02),
    replacements = list(
      `1` = ctmc_generator(0.18, 0.05),
      `2` = ctmc_generator(0.03, 0.11),
      `12` = ctmc_generator(0.07, 0.13)
    )
  )
  mu <- ctmc_mu_list(
    c(mu1 = 0.35, mu0 = -0.15),
    replacements = list(
      `1` = c(mu1 = 0.5, mu0 = -0.2),
      `2` = c(mu1 = 0.1, mu0 = 0.25),
      `12` = c(mu1 = -0.4, mu0 = 0.6)
    )
  )
  sd <- ctmc_sd_list(
    c(sd1 = 0.7, sd0 = 1.1),
    replacements = list(
      `1` = c(sd1 = 0.8, sd0 = 1.0),
      `2` = c(sd1 = 1.2, sd0 = 0.9),
      `12` = c(sd1 = 0.6, sd0 = 1.4)
    )
  )
  ctmc_fixture(t_now = t_now, t_hor = t_hor, Q = Q, mu = mu, sd = sd)
}

ctmc_zero_transition_fixture <- function(t_now = "2022-01-15",
                                         t_hor = "2022-01-20",
                                         p0 = c(0.4, 0.6)) {
  ctmc_fixture(
    t_now = t_now,
    t_hor = t_hor,
    Q = ctmc_month_list(matrix(0, nrow = 2, ncol = 2)),
    mu = ctmc_mu_list(c(mu1 = 0.25, mu0 = -0.10)),
    sd = ctmc_sd_list(c(sd1 = 0.5, sd0 = 0.9)),
    p0 = p0
  )
}

ctmc_reference_phi <- function(a, b, fix) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (length(a) != 1 || length(b) != 1) {
    stop("ctmc_reference_phi expects scalar a and b")
  }
  if (b < a) {
    stop("ctmc_reference_phi expects a <= b")
  }
  if (isTRUE(all.equal(a, b, tolerance = 1e-14))) {
    return(diag(1, 2, 2))
  }

  start_n <- solarr::number_of_day(fix$t_now)
  offset_a <- a - start_n
  offset_b <- b - start_n
  interior <- if (ceiling(offset_a) <= floor(offset_b)) {
    seq(ceiling(offset_a), floor(offset_b), by = 1)
  } else {
    numeric(0)
  }
  cuts <- sort(unique(c(offset_a, offset_b, interior)))
  cuts <- cuts[cuts >= offset_a - 1e-12 & cuts <= offset_b + 1e-12]
  cuts[1] <- offset_a
  cuts[length(cuts)] <- offset_b

  P <- diag(1, 2, 2)
  for (j in seq_len(length(cuts) - 1)) {
    left <- cuts[j]
    right <- cuts[j + 1]
    if (right - left <= 1e-12) {
      next
    }
    active_date <- fix$t_now + floor(left + 1e-12)
    active_month <- lubridate::month(active_date)
    P <- P %*% expm::expm((right - left) * fix$Q[[active_month]])
  }
  P
}

ctmc_reference_month <- function(s, fix) {
  start_n <- solarr::number_of_day(fix$t_now)
  offset <- as.numeric(s) - start_n
  lubridate::month(fix$t_now + floor(offset + 1e-12))
}

ctmc_reference_pT <- function(fix, p0 = fix$p0) {
  drop(p0 %*% ctmc_reference_phi(fix$t_init, fix$t_end, fix))
}

ctmc_reference_state_value <- function(s, fix, values, p0 = fix$p0) {
  vapply(s, function(si) {
    m <- ctmc_reference_month(si, fix)
    drop(p0 %*% ctmc_reference_phi(fix$t_init, si, fix) %*% values[[m]])
  }, numeric(1))
}

ctmc_reference_state_value_cond <- function(s, fix, values, ei, p0 = fix$p0) {
  denom <- drop(p0 %*% ctmc_reference_phi(fix$t_init, fix$t_end, fix) %*% ei)
  vapply(s, function(si) {
    m <- ctmc_reference_month(si, fix)
    num <- p0 %*%
      ctmc_reference_phi(fix$t_init, si, fix) %*%
      diag(values[[m]], nrow = 2) %*%
      ctmc_reference_phi(si, fix$t_end, fix) %*%
      ei
    drop(num) / denom
  }, numeric(1))
}

ctmc_reference_integral_state_value <- function(fix, values, theta,
                                                sigma_bar,
                                                power = 1,
                                                p0 = fix$p0) {
  stats::integrate(
    function(s) {
      exp(-power * theta * (fix$t_end - s)) *
        sigma_bar(s)^power *
        ctmc_reference_state_value(s, fix, values, p0)
    },
    lower = fix$t_init,
    upper = fix$t_end,
    rel.tol = 1e-10,
    stop.on.error = FALSE
  )$value
}

ctmc_reference_integral_state_value_cond <- function(fix, values, theta,
                                                     sigma_bar, ei,
                                                     power = 1,
                                                     p0 = fix$p0) {
  stats::integrate(
    function(s) {
      exp(-power * theta * (fix$t_end - s)) *
        sigma_bar(s)^power *
        ctmc_reference_state_value_cond(s, fix, values, ei, p0)
    },
    lower = fix$t_init,
    upper = fix$t_end,
    rel.tol = 1e-10,
    stop.on.error = FALSE
  )$value
}

ctmc_reference_cross_mu <- function(s, u, fix, p0 = fix$p0) {
  if (s > u) {
    tmp <- s
    s <- u
    u <- tmp
  }
  m_s <- ctmc_reference_month(s, fix)
  m_u <- ctmc_reference_month(u, fix)
  drop(
    p0 %*%
      ctmc_reference_phi(fix$t_init, s, fix) %*%
      diag(fix$mu[[m_s]], nrow = 2) %*%
      ctmc_reference_phi(s, u, fix) %*%
      fix$mu[[m_u]]
  )
}

ctmc_reference_cross_mu_cond <- function(s, u, fix, ei, p0 = fix$p0) {
  if (s > u) {
    tmp <- s
    s <- u
    u <- tmp
  }
  denom <- drop(p0 %*% ctmc_reference_phi(fix$t_init, fix$t_end, fix) %*% ei)
  m_s <- ctmc_reference_month(s, fix)
  m_u <- ctmc_reference_month(u, fix)
  drop(
    p0 %*%
      ctmc_reference_phi(fix$t_init, s, fix) %*%
      diag(fix$mu[[m_s]], nrow = 2) %*%
      ctmc_reference_phi(s, u, fix) %*%
      diag(fix$mu[[m_u]], nrow = 2) %*%
      ctmc_reference_phi(u, fix$t_end, fix) %*%
      ei
  ) / denom
}

ctmc_reference_integral_cross_mu <- function(fix, theta, sigma_bar,
                                             conditional_ei = NULL,
                                             p0 = fix$p0) {
  stats::integrate(
    function(s_vec) {
      vapply(s_vec, function(s) {
        stats::integrate(
          function(u) {
            cross <- if (is.null(conditional_ei)) {
              vapply(u, ctmc_reference_cross_mu, numeric(1), s = s, fix = fix, p0 = p0)
            } else {
              vapply(u, ctmc_reference_cross_mu_cond, numeric(1),
                     s = s, fix = fix, ei = conditional_ei, p0 = p0)
            }
            exp(-theta * (2 * fix$t_end - s - u)) *
              sigma_bar(s) * sigma_bar(u) * cross
          },
          lower = fix$t_init,
          upper = fix$t_end,
          rel.tol = 1e-7,
          stop.on.error = FALSE
        )$value
      }, numeric(1))
    },
    lower = fix$t_init,
    upper = fix$t_end,
    rel.tol = 1e-7,
    stop.on.error = FALSE
  )$value
}

ctmc_reference_integral_Yt_Xt <- function(fix, theta, kappa, sigma_bar,
                                          sigma_X, rho, ei, p0 = fix$p0) {
  values <- purrr::map2(fix$sd, rho, ~.x * .y)
  stats::integrate(
    function(s) {
      exp(-theta * (fix$t_end - s)) *
        exp(-kappa * (fix$t_end - s)) *
        sigma_bar(s) *
        sigma_X *
        ctmc_reference_state_value_cond(s, fix, values, ei, p0)
    },
    lower = fix$t_init,
    upper = fix$t_end,
    rel.tol = 1e-10,
    stop.on.error = FALSE
  )$value
}

ctmc_expected_lengths_for_half_open_dates <- function(t_now, t_hor) {
  dates <- seq.Date(as.Date(t_now), as.Date(t_hor) - 1, by = "day")
  as.numeric(table(factor(lubridate::month(dates), levels = 1:12)))
}

ctmc_mock_ctmc_model <- function(t_now = "2022-01-15",
                                 p0 = c(0.4, 0.6),
                                 lambda = 0) {
  date_seq <- seq.Date(as.Date(t_now), as.Date(t_now) + 30, by = "day")
  alpha <- 0.05
  beta <- 0.90
  RY <- function(R, C) {
    z <- (R / C - (1 - alpha - beta)) / beta
    qlogis(pmin(pmax(z, 1e-8), 1 - 1e-8))
  }
  iRY <- function(Y, C) {
    C * (1 - alpha - beta + beta * plogis(Y))
  }
  CTMC <- list(
    params = list(
      mu = ctmc_mu_list(c(mu1 = 0.20, mu0 = -0.05)),
      sig = ctmc_sd_list(c(sd1 = 0.45, sd0 = 0.75)),
      Pm = ctmc_month_list(diag(1, 2, 2)),
      Qm = ctmc_month_list(matrix(0, nrow = 2, ncol = 2))
    ),
    alpha = cbind(
      alpha1 = rep(p0[1], length(date_seq)),
      alpha2 = rep(p0[2], length(date_seq))
    ),
    data = tibble::tibble(date = date_seq)
  )
  Ct_fun <- function(d) rep(100, length(d))
  Ybar_fun <- function(d) rep(0, length(d))
  data <- tibble::tibble(
    date = date_seq,
    Year = lubridate::year(date_seq),
    Month = lubridate::month(date_seq),
    Day = lubridate::day(date_seq),
    GHI = rep(50, length(date_seq)),
    Yt = RY(rep(50, length(date_seq)), Ct_fun(date_seq))
  )

  structure(list(
    model = list(
      data = data,
      spec = list(
        transform = list(alpha = alpha, beta = beta, RY = RY, iRY = iRY),
        mean.model = list(phi = exp(-0.1)),
        seasonal.variance = list(predict = function(d) rep(1, length(d)))
      )
    ),
    theta = 0.1,
    seasonal_variance = list(
      extra_params = list(seasonal_function = ctmc_sigma_bar_one)
    ),
    CTMC = CTMC,
    HMM = CTMC,
    Ct = Ct_fun,
    dCt = function(d) rep(0, length(d)),
    Yt_bar = Ybar_fun,
    dYt_bar = function(d) rep(0, length(d)),
    sigma_bar = function(d) rep(1, length(d)),
    lambda = lambda
  ), class = "radiationModel_CTMC")
}

ctmc_mock_joint_scenario_inputs <- function(t_now = "2022-01-15",
                                            t_hor = "2022-01-18",
                                            p0 = c(0.4, 0.6)) {
  date_seq <- seq.Date(as.Date(t_now), as.Date(t_hor), by = "day")
  model_Rt <- ctmc_mock_ctmc_model(t_now = t_now, p0 = p0)
  model_Rt$model$data <- dplyr::filter(model_Rt$model$data, date %in% date_seq)
  model_Rt$CTMC$data <- tibble::tibble(date = date_seq)
  model_Rt$CTMC$alpha <- cbind(
    alpha1 = rep(p0[1], length(date_seq)),
    alpha2 = rep(p0[2], length(date_seq))
  )

  model_Et <- list(
    data = tibble::tibble(
      date = date_seq,
      PUN = rep(100, length(date_seq)),
      Xt = log(rep(100, length(date_seq)))
    ),
    model = list(lambda = 0, kappa = 0.4, mu = log(100), sigma = 0.2),
    F_E = function(t, T, E0) rep(E0, length(T))
  )

  list(
    model_Et = model_Et,
    model_Rt = model_Rt,
    rho = ctmc_month_list(c(rho1 = 0, rho0 = 0))
  )
}
