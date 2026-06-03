# ---
#' @description
#' Test `get_month_index_C()`.
#'
#' @tests
#' - returns the month active at the left endpoint of each half-open interval;
#' - agrees with an explicit calendar reference across month/year boundaries;
#' - optional C wrapper matches the production R version.
# ---

test_that("get_month_index_C follows the half-open month convention", {
  fix <- ctmc_month_boundary_fixture("2022-12-31", "2023-01-02")
  tau <- fix$t_init + c(0, 0.5, 1, 1.5)
  expected <- vapply(tau, ctmc_reference_month, numeric(1), fix = fix)

  expect_equal(get_month_index_C(tau, fix$bounds), expected)
})

test_that("get_month_index_C optional C wrapper matches production R", {
  source(file.path(ctmc_radiation_dir, "ctmc-integrals-C.R"))
  get_month_index_C_R <- get_month_index_C

  cases <- list(
    list(t_now = "2022-12-31", t_hor = "2023-01-02"),
    list(t_now = "2022-01-31", t_hor = "2022-02-02")
  )

  for (case in cases) {
    fix <- ctmc_month_boundary_fixture(case$t_now, case$t_hor)
    tau <- sort(unique(c(fix$bounds$n, fix$bounds$N,
                         seq(fix$t_init, fix$t_end, length.out = 7))))
    tau <- tau[tau >= fix$t_init & tau <= fix$t_end]
    expected <- get_month_index_C_R(tau, fix$bounds)

    source(file.path(ctmc_radiation_dir, "ctmc-integrals-C-wrappers.R"))
    expect_equal(get_month_index_C(tau, fix$bounds), expected)
    source(file.path(ctmc_radiation_dir, "ctmc-integrals-C.R"))
  }
})

