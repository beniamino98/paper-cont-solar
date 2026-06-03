# ---
#' @description
#' Test `integral_E_Yt_Xt()`.
#'
#' @tests
#' - conditional radiation/electricity covariance integrals match independent
#'   Markov-bridge references;
#' - state-specific monthly correlations are used through the same month index
#'   convention as the CTMC transition products.
# ---

test_that("integral_E_Yt_Xt matches independent conditional covariance references", {
  fix <- ctmc_fixture()
  theta <- 0.1
  kappa <- 0.25
  sigma_X <- 0.8
  rho <- ctmc_month_list(c(rho1 = 0.2, rho0 = -0.35))

  got <- integral_E_Yt_Xt(
    fix$bounds, fix$p0, fix$sd, theta, kappa,
    ctmc_sigma_bar_seasonal, sigma_X, rho,
    maxEval = 20000, tol = 1e-8
  )

  expected <- c(
    ctmc_reference_integral_Yt_Xt(fix, theta, kappa, ctmc_sigma_bar_seasonal,
                                  sigma_X, rho, c(1, 0)),
    ctmc_reference_integral_Yt_Xt(fix, theta, kappa, ctmc_sigma_bar_seasonal,
                                  sigma_X, rho, c(0, 1))
  )

  expect_equal(as.numeric(got), expected, tolerance = 1e-7)
})

