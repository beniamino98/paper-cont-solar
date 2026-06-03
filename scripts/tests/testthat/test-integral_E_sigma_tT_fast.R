# ---
#' @description
#' Test `integral_E_sigma_tT_fast()`.
#'
#' @tests
#' - unconditional one-dimensional diffusion-variance integrals match
#'   independent `stats::integrate()` references;
#' - checks both constant and seasonal volatility scalings.
# ---

test_that("integral_E_sigma_tT_fast matches independent base-integrate references", {
  fix <- ctmc_fixture()

  for (theta in c(0, 0.1, 0.5)) {
    for (sigma_bar in list(ctmc_sigma_bar_one, ctmc_sigma_bar_seasonal)) {
      expect_equal(
        as.numeric(integral_E_sigma_tT_fast(fix$bounds, fix$p0, fix$sd, theta, sigma_bar)),
        ctmc_reference_integral_state_value(
          fix, purrr::map(fix$sd, ~.x^2), theta, sigma_bar, power = 2
        ),
        tolerance = 1e-7
      )
    }
  }
})

