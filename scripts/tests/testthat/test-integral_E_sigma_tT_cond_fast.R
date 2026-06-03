# ---
#' @description
#' Test `integral_E_sigma_tT_cond_fast()`.
#'
#' @tests
#' - conditional diffusion-variance integrals match independent Markov-bridge
#'   references;
#' - probability-weighted conditional integrals recover the unconditional
#'   diffusion-variance integral.
# ---

test_that("integral_E_sigma_tT_cond_fast matches independent conditional references", {
  fix <- ctmc_fixture()
  theta <- 0.1
  sigma_bar <- ctmc_sigma_bar_seasonal
  sd2 <- purrr::map(fix$sd, ~.x^2)

  got_sigma <- integral_E_sigma_tT_cond_fast(
    fix$bounds, fix$p0, fix$sd, theta, sigma_bar,
    maxEval = 20000, tol = 1e-8
  )

  expected_sigma <- c(
    ctmc_reference_integral_state_value_cond(fix, sd2, theta, sigma_bar, c(1, 0), power = 2),
    ctmc_reference_integral_state_value_cond(fix, sd2, theta, sigma_bar, c(0, 1), power = 2)
  )

  expect_equal(as.numeric(got_sigma), expected_sigma, tolerance = 1e-7)
})

test_that("integral_E_sigma_tT_cond_fast recovers the unconditional integral", {
  fix <- ctmc_fixture()
  theta <- 0.1
  sigma_bar <- ctmc_sigma_bar_one
  p_T <- drop(fix$p0 %*% Phi_C(fix$t_init, fix$t_end, fix$bounds)[[1]])

  sig_cond <- integral_E_sigma_tT_cond_fast(
    fix$bounds, fix$p0, fix$sd, theta, sigma_bar,
    maxEval = 20000, tol = 1e-8
  )

  expect_equal(
    as.numeric(integral_E_sigma_tT_fast(fix$bounds, fix$p0, fix$sd, theta, sigma_bar)),
    sum(p_T * sig_cond),
    tolerance = 1e-8
  )
})

