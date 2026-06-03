# ---
#' @description
#' Test `transition_list_to_generator_2state()`.
#'
#' @tests
#' - converts every monthly transition matrix in a list;
#' - preserves list length and order;
#' - rejects non-list inputs.
# ---

test_that("transition_list_to_generator_2state converts monthly transitions", {
  Pm <- ctmc_month_list(
    matrix(c(0.90, 0.10,
             0.05, 0.95),
           2, byrow = TRUE)
  )

  Qm <- transition_list_to_generator_2state(Pm)

  expect_length(Qm, length(Pm))
  for (m in seq_along(Pm)) {
    expect_equal(unname(expm::expm(Qm[[m]])), unname(Pm[[m]]), tolerance = 1e-12)
    expect_equal(rowSums(Qm[[m]]), c(0, 0), tolerance = 1e-12)
  }
})

test_that("transition_list_to_generator_2state rejects non-list inputs", {
  expect_error(transition_list_to_generator_2state(diag(1, 2, 2)), "list")
})

