# ---
#' @description
#' Test `integral_E_mu_tT_cond_fast()`.
#'
#' @tests
#' - conditional drift integrals match independent Markov-bridge references;
#' - probability-weighted conditional integrals recover the unconditional
#'   drift integral.
# ---

test_that("integral_E_mu_tT_cond_fast matches independent conditional references", {
  fix <- ctmc_fixture()
  theta <- 0.1
  sigma_bar <- ctmc_sigma_bar_seasonal

  got_mu <- integral_E_mu_tT_cond_fast(
    fix$bounds, fix$p0, fix$mu, theta, sigma_bar,
    maxEval = 20000, tol = 1e-8
  )

  expected_mu <- c(
    ctmc_reference_integral_state_value_cond(fix, fix$mu, theta, sigma_bar, c(1, 0), power = 1),
    ctmc_reference_integral_state_value_cond(fix, fix$mu, theta, sigma_bar, c(0, 1), power = 1)
  )

  expect_equal(as.numeric(got_mu), expected_mu, tolerance = 1e-7)
})

test_that("integral_E_mu_tT_cond_fast recovers the unconditional integral", {
  fix <- ctmc_fixture()
  theta <- 0.1
  sigma_bar <- ctmc_sigma_bar_one
  p_T <- drop(fix$p0 %*% Phi_C(fix$t_init, fix$t_end, fix$bounds)[[1]])

  mu_cond <- integral_E_mu_tT_cond_fast(
    fix$bounds, fix$p0, fix$mu, theta, sigma_bar,
    maxEval = 20000, tol = 1e-8
  )

  expect_equal(
    as.numeric(integral_E_mu_tT_fast(fix$bounds, fix$p0, fix$mu, theta, sigma_bar)),
    sum(p_T * mu_cond),
    tolerance = 1e-8
  )
})

