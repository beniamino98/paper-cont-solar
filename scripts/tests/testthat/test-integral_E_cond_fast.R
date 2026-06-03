# ---
#' @description
#' Test `integral_E_cond_fast()`.
#'
#' @tests
#' - output names match the production mathematical quantities;
#' - rows 1-2 are conditional drift integrals;
#' - rows 3-4 are conditional diffusion-variance contributions;
#' - rows 5-6 are conditional first-moment diffusion-scale integrals;
#' - optional C-backed transition kernels leave the aggregate result unchanged.
# ---

test_that("integral_E_cond_fast rows carry intended mathematical quantities", {
  fix <- ctmc_fixture()
  theta <- 0.1
  sigma_bar <- ctmc_sigma_bar_one
  sd2 <- purrr::map(fix$sd, ~.x^2)

  got <- integral_E_cond_fast(
    fix$bounds, fix$p0, fix$mu, fix$sd, theta, sigma_bar,
    maxEval = 20000, tol = 1e-8
  )

  expect_named(got, c(
    "E_mu_tT_1", "E_mu_tT_0",
    "V_sigma_tT_1", "V_sigma_tT_0",
    "E_sigma_tT_1", "E_sigma_tT_0"
  ))

  expected_by_position <- c(
    ctmc_reference_integral_state_value_cond(fix, fix$mu, theta, sigma_bar, c(1, 0), power = 1),
    ctmc_reference_integral_state_value_cond(fix, fix$mu, theta, sigma_bar, c(0, 1), power = 1),
    ctmc_reference_integral_state_value_cond(fix, sd2, theta, sigma_bar, c(1, 0), power = 2),
    ctmc_reference_integral_state_value_cond(fix, sd2, theta, sigma_bar, c(0, 1), power = 2),
    ctmc_reference_integral_state_value_cond(fix, fix$sd, theta, sigma_bar, c(1, 0), power = 1),
    ctmc_reference_integral_state_value_cond(fix, fix$sd, theta, sigma_bar, c(0, 1), power = 1)
  )
  expect_equal(as.numeric(got), unname(expected_by_position), tolerance = 1e-7)
})

test_that("integral_E_cond_fast is unchanged by optional C transition wrappers", {
  source(file.path(ctmc_radiation_dir, "ctmc-integrals-C.R"))
  fix <- ctmc_month_boundary_fixture("2022-12-31", "2023-01-02")
  theta <- 0.1
  sigma_bar <- ctmc_sigma_bar_seasonal
  p_T <- drop(fix$p0 %*% Phi_C(fix$t_init, fix$t_end, fix$bounds)[[1]])

  r_out <- integral_E_cond_fast(
    fix$bounds, fix$p0, fix$mu, fix$sd, theta, sigma_bar,
    p_T = p_T, maxEval = 20000, tol = 1e-8
  )

  source(file.path(ctmc_radiation_dir, "ctmc-integrals-C-wrappers.R"))
  c_out <- integral_E_cond_fast(
    fix$bounds, fix$p0, fix$mu, fix$sd, theta, sigma_bar,
    p_T = p_T, maxEval = 20000, tol = 1e-8
  )

  expect_equal(as.numeric(c_out), as.numeric(r_out), tolerance = 1e-7)
  expect_named(c_out, names(r_out))

  source(file.path(ctmc_radiation_dir, "ctmc-integrals-C.R"))
})

