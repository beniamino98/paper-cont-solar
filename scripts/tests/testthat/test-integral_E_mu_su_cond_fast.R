# ---
#' @description
#' Test `integral_E_mu_su_cond_fast()`.
#'
#' @tests
#' - conditional two-dimensional drift cross integrals match independent
#'   Markov-bridge references;
#' - probability-weighted conditional values recover the unconditional
#'   cross integral;
#' - implied conditional variance contributions are non-negative.
# ---

test_that("integral_E_mu_su_cond_fast matches independent conditional references", {
  fix <- ctmc_fixture()
  theta <- 0.1
  sigma_bar <- ctmc_sigma_bar_one
  p_T <- drop(fix$p0 %*% Phi_C(fix$t_init, fix$t_end, fix$bounds)[[1]])

  got_uncond <- integral_E_mu_su_fast(
    fix$bounds, fix$p0, fix$mu, theta, sigma_bar,
    maxEval = 30000, tol = 1e-6
  )
  got_cond <- integral_E_mu_su_cond_fast(
    fix$bounds, fix$p0, fix$mu, theta, sigma_bar, p_T,
    maxEval = 30000, tol = 1e-6
  )

  expected_cond <- c(
    ctmc_reference_integral_cross_mu(fix, theta, sigma_bar, conditional_ei = c(1, 0)),
    ctmc_reference_integral_cross_mu(fix, theta, sigma_bar, conditional_ei = c(0, 1))
  )

  expect_equal(as.numeric(got_cond), expected_cond, tolerance = 1e-5)
  expect_equal(as.numeric(got_uncond), sum(p_T * got_cond), tolerance = 1e-7)
})

test_that("integral_E_mu_su_cond_fast variance contribution is non-negative", {
  fix <- ctmc_fixture()
  theta <- 0.1
  sigma_bar <- ctmc_sigma_bar_one
  p_T <- drop(fix$p0 %*% Phi_C(fix$t_init, fix$t_end, fix$bounds)[[1]])

  M_mu <- integral_E_mu_tT_cond_fast(
    fix$bounds, fix$p0, fix$mu, theta, sigma_bar,
    maxEval = 20000, tol = 1e-8
  )
  E_mu2 <- integral_E_mu_su_cond_fast(
    fix$bounds, fix$p0, fix$mu, theta, sigma_bar, p_T,
    maxEval = 30000, tol = 1e-6
  )

  expect_true(all(as.numeric(E_mu2) - as.numeric(M_mu)^2 >= -1e-8))
})

