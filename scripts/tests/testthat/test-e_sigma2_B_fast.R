# ---
#' @description
#' Test `e_sigma2_B_fast()`.
#'
#' @tests
#' - vector inputs match explicit CTMC transition-product references;
#' - state-dependent monthly variance parameters are selected correctly after a
#'   month boundary.
# ---

test_that("e_sigma2_B_fast matches explicit transition-product references", {
  fix <- ctmc_fixture()
  sd2 <- purrr::map(fix$sd, ~.x^2)
  s <- fix$t_init + c(0.1, 1.25, 3.75, 4.9)

  expect_equal(
    e_sigma2_B_fast(s, fix$bounds, fix$p0, sd2),
    ctmc_reference_state_value(s, fix, sd2),
    tolerance = 1e-9
  )
})

test_that("e_sigma2_B_fast uses correct monthly parameters after crossing a month", {
  fix <- ctmc_month_boundary_fixture("2022-01-30", "2022-02-02")
  sd2 <- purrr::map(fix$sd, ~.x^2)
  s <- fix$t_init + c(0.25, 1.25, 2.25)

  expect_equal(
    e_sigma2_B_fast(s, fix$bounds, fix$p0, sd2),
    ctmc_reference_state_value(s, fix, sd2),
    tolerance = 1e-9
  )
})

