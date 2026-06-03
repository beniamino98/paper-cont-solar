# ---
#' @description
#' Test `ctmc_phi_one()`.
#'
#' @tests
#' - returns identity on zero-length intervals;
#' - matches explicit products of monthly transition matrices;
#' - respects half-open month/year boundary products.
# ---

test_that("ctmc_phi_one matches explicit half-open transition products", {
  cases <- list(
    list(t_now = "2022-12-31", t_hor = "2023-01-01"),
    list(t_now = "2022-12-31", t_hor = "2023-01-02"),
    list(t_now = "2022-01-15", t_hor = "2022-01-20")
  )

  for (case in cases) {
    fix <- ctmc_month_boundary_fixture(case$t_now, case$t_hor)
    got <- ctmc_phi_one(fix$t_init, fix$t_end, fix$bounds)
    expected <- ctmc_reference_phi(fix$t_init, fix$t_end, fix)

    expect_equal(unname(got), unname(expected), tolerance = 1e-9)
  }
})

test_that("ctmc_phi_one returns identity on zero-length intervals", {
  fix <- ctmc_fixture()
  expect_equal(ctmc_phi_one(fix$t_init, fix$t_init, fix$bounds), diag(1, 2, 2))
})

