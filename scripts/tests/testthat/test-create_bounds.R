# ---
#' @description
#' Test `create_bounds()`.
#'
#' @tests
#' - uses the half-open interval convention `[t, T)`;
#' - handles month and year boundaries without adding an endpoint month;
#' - keeps monthly segment lengths consistent with explicit calendar counts.
# ---

test_that("create_bounds uses half-open calendar lengths at month boundaries", {
  cases <- list(
    list(t_now = "2022-12-31", t_hor = "2023-01-01",
         expected_months = c(12), expected_lengths = c(1)),
    list(t_now = "2022-12-31", t_hor = "2023-01-02",
         expected_months = c(12, 1), expected_lengths = c(1, 1)),
    list(t_now = "2022-01-31", t_hor = "2022-02-02",
         expected_months = c(1, 2), expected_lengths = c(1, 1))
  )

  for (case in cases) {
    fix <- ctmc_month_boundary_fixture(case$t_now, case$t_hor)
    positive_rows <- which((fix$bounds$N - fix$bounds$n) > 0)

    expect_equal(fix$bounds$Month[positive_rows], case$expected_months,
                 info = paste(case$t_now, "to", case$t_hor))
    expect_equal(fix$bounds$N[positive_rows] - fix$bounds$n[positive_rows],
                 case$expected_lengths,
                 info = paste(case$t_now, "to", case$t_hor))
  }
})

