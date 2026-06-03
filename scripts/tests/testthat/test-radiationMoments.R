# ---
#' @description
#' Test `radiationMoments()`.
#'
#' @tests
#' - dispatches `radiationModel_CTMC` objects to
#'   `radiationModel_CTMC_moments()`;
#' - preserves the row-level moment output for the CTMC branch.
# ---

test_that("radiationMoments dispatches CTMC models to radiationModel_CTMC_moments", {
  t_now <- "2022-01-15"
  t_hor <- as.Date(c("2022-01-17", "2022-01-19"))
  model <- ctmc_mock_ctmc_model(t_now = t_now, p0 = c(0.4, 0.6), lambda = 0)

  got <- radiationMoments(
    t_now = t_now,
    t_hor = t_hor,
    model_Rt = model,
    R0 = 50,
    maxEval = 20000,
    tol = 1e-8
  )
  expected <- radiationModel_CTMC_moments(
    t_now = t_now,
    t_hor = t_hor,
    model_Rt = model,
    R0 = 50,
    maxEval = 20000,
    tol = 1e-8
  )

  expect_equal(got$p1, expected$p1, tolerance = 1e-12)
  expect_equal(got$M_Y1, expected$M_Y1, tolerance = 1e-12)
  expect_equal(got$M_Y0, expected$M_Y0, tolerance = 1e-12)
  expect_equal(got$S_Y1, expected$S_Y1, tolerance = 1e-12)
  expect_equal(got$S_Y0, expected$S_Y0, tolerance = 1e-12)
})

