# ---
#' @description
#' Test `ctmc_as_generator_list()`.
#'
#' @tests
#' - returns generator lists unchanged;
#' - converts transition-matrix lists to generator lists;
#' - rejects mixed or invalid monthly lists.
# ---

test_that("ctmc_as_generator_list accepts generator lists and transition lists", {
  Qm <- ctmc_month_list(ctmc_generator(0.08, 0.03))
  expect_identical(ctmc_as_generator_list(Qm), Qm)

  Pm <- lapply(Qm, function(Q) expm::expm(Q))
  got <- ctmc_as_generator_list(Pm)

  for (m in seq_along(Pm)) {
    expect_equal(unname(expm::expm(got[[m]])), unname(Pm[[m]]), tolerance = 1e-12)
  }
})

test_that("ctmc_as_generator_list rejects mixed monthly matrix types", {
  mixed <- ctmc_month_list(ctmc_generator())
  mixed[[2]] <- expm::expm(ctmc_generator())

  expect_error(ctmc_as_generator_list(mixed), "either all")
  expect_error(ctmc_as_generator_list(diag(1, 2, 2)), "list")
})

