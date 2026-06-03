# ---
#' @description
#' Test `integral_E_mu_su_fast()`.
#'
#' @tests
#' - two-dimensional drift cross integrals match independent nested-integral
#'   references;
#' - swapped integration rectangles give the same value.
# ---

test_that("integral_E_mu_su_fast matches independent cross-integral references", {
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

  expected_uncond <- ctmc_reference_integral_cross_mu(fix, theta, sigma_bar)

  expect_equal(as.numeric(got_uncond), expected_uncond, tolerance = 1e-5)
  expect_equal(as.numeric(got_uncond), sum(p_T * got_cond), tolerance = 1e-7)
})

test_that("integral_E_mu_su_fast is symmetric through swapped rectangles", {
  fix <- ctmc_fixture()
  theta <- 0.1
  sigma_bar <- ctmc_sigma_bar_one
  eps <- 0.05
  s <- fix$t_init + 1.0
  u <- fix$t_init + 3.0

  got_su <- integral_E_mu_su_fast(
    fix$bounds, fix$p0, fix$mu, theta, sigma_bar,
    t0 = c(s - fix$t_init, u - fix$t_init),
    T0 = c(s + eps - fix$t_end, u + eps - fix$t_end),
    maxEval = 10000, tol = 1e-8
  )
  got_us <- integral_E_mu_su_fast(
    fix$bounds, fix$p0, fix$mu, theta, sigma_bar,
    t0 = c(u - fix$t_init, s - fix$t_init),
    T0 = c(u + eps - fix$t_end, s + eps - fix$t_end),
    maxEval = 10000, tol = 1e-8
  )

  expect_equal(as.numeric(got_su), as.numeric(got_us), tolerance = 1e-9)
})

