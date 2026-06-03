# ---
#' @description
#' Test `e_mu_B_fast()`.
#'
#' @tests
#' - scalar and vector inputs match explicit CTMC transition-product references;
#' - state-dependent monthly drift parameters are selected correctly after a
#'   month boundary.
# ---

test_that("e_mu_B_fast matches explicit transition-product references", {
  fix <- ctmc_fixture()
  s_scalar <- fix$t_init + 2.5
  s_vector <- fix$t_init + c(0.1, 1.25, 3.75, 4.9)

  expect_equal(
    e_mu_B_fast(s_scalar, fix$bounds, fix$p0, fix$mu),
    ctmc_reference_state_value(s_scalar, fix, fix$mu),
    tolerance = 1e-9
  )
  expect_equal(
    e_mu_B_fast(s_vector, fix$bounds, fix$p0, fix$mu),
    ctmc_reference_state_value(s_vector, fix, fix$mu),
    tolerance = 1e-9
  )
})

test_that("e_mu_B_fast uses correct monthly parameters after crossing a month", {
  fix <- ctmc_month_boundary_fixture("2022-01-30", "2022-02-02")
  s <- fix$t_init + c(0.25, 1.25, 2.25)

  expect_equal(
    e_mu_B_fast(s, fix$bounds, fix$p0, fix$mu),
    ctmc_reference_state_value(s, fix, fix$mu),
    tolerance = 1e-9
  )
})

