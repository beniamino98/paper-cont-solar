# ---
#' @description
#' Test `radiationModel_CTMC_moments()`.
#'
#' @tests
#' - returns finite two-Gaussian moment rows for a minimal CTMC model;
#' - mixture probabilities and conditional standard deviations are valid;
#' - generated CDF closures are monotone with sensible tails;
#' - terminal probabilities agree with `Phi_C()`.
# ---

test_that("radiationModel_CTMC_moments returns stable mixture moments", {
  t_now <- "2022-01-15"
  t_hor <- as.Date(c("2022-01-17", "2022-01-19"))
  model <- ctmc_mock_ctmc_model(t_now = t_now, p0 = c(0.4, 0.6), lambda = 0)

  mom <- radiationModel_CTMC_moments(
    t_now = t_now,
    t_hor = t_hor,
    model_Rt = model,
    R0 = 50,
    maxEval = 20000,
    tol = 1e-8
  )

  expect_equal(nrow(mom), length(t_hor))

  core_numeric <- c(
    "e_Yt", "sd_Yt", "M_Y1", "M_Y0", "S_Y1", "S_Y0",
    "p1", "Ct", "GHI_bar", "alpha", "beta", "RT_min", "RT_max",
    "M_mu_1", "M_mu_0", "S2_sigma_1", "S2_sigma_0",
    "M_gamma_1", "M_gamma_0", "S2_mu_1", "S2_mu_0"
  )
  expect_true(all(is.finite(unlist(mom[core_numeric]))))
  expect_true(all(mom$p1 >= -1e-12 & mom$p1 <= 1 + 1e-12))
  expect_true(all(mom$S_Y1 > 0))
  expect_true(all(mom$S_Y0 > 0))
  expect_true(all(vapply(mom$pdf_Y, is.function, logical(1))))
  expect_true(all(vapply(mom$cdf_Y, is.function, logical(1))))
})

test_that("radiationModel_CTMC_moments CDF functions are monotone with sensible tails", {
  t_now <- "2022-01-15"
  t_hor <- as.Date(c("2022-01-17", "2022-01-19"))
  model <- ctmc_mock_ctmc_model(t_now = t_now, p0 = c(0.4, 0.6), lambda = 0)

  mom <- radiationModel_CTMC_moments(
    t_now = t_now,
    t_hor = t_hor,
    model_Rt = model,
    R0 = 50,
    maxEval = 20000,
    tol = 1e-8
  )

  for (idx in seq_len(nrow(mom))) {
    center <- mean(c(mom$M_Y1[idx], mom$M_Y0[idx]))
    scale <- max(c(mom$S_Y1[idx], mom$S_Y0[idx]))
    grid <- seq(center - 8 * scale, center + 8 * scale, length.out = 101)
    cdf_values <- mom$cdf_Y[[idx]](grid)

    expect_true(all(diff(cdf_values) >= -1e-10))
    expect_lt(cdf_values[1], 1e-6)
    expect_gt(cdf_values[length(cdf_values)], 1 - 1e-6)
  }
})

test_that("radiationModel_CTMC_moments p_T is consistent with Phi_C", {
  t_now <- "2022-01-15"
  t_hor <- as.Date(c("2022-01-17", "2022-01-19"))
  p0 <- c(0.4, 0.6)
  model <- ctmc_mock_ctmc_model(t_now = t_now, p0 = p0, lambda = 0)

  mom <- radiationModel_CTMC_moments(
    t_now = t_now,
    t_hor = t_hor,
    model_Rt = model,
    R0 = 50,
    maxEval = 20000,
    tol = 1e-8
  )

  Q <- if (!is.null(model$CTMC$params$Qm)) {
    model$CTMC$params$Qm
  } else {
    transition_list_to_generator_2state(model$CTMC$params$Pm)
  }
  for (idx in seq_along(t_hor)) {
    bounds <- create_bounds(t_now, t_hor[idx], Q)
    p_T <- drop(p0 %*% Phi_C(bounds$n[1], bounds$tau[1], bounds)[[1]])
    expect_equal(mom$p1[idx], p_T[1], tolerance = 1e-12)
  }
})

