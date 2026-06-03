# ---
#' @description
#' Test `scenarios_radiationModel_CT()`.
#'
#' @tests
#' - returns one row per simulation and time step;
#' - includes simulated transformed and physical radiation paths;
#' - honors the terminal half-open step count implied by `t_hor` and `dt`.
# ---

test_that("scenarios_radiationModel_CT returns coherent CTMC radiation paths", {
  model <- ctmc_mock_ctmc_model(t_now = "2022-01-15", p0 = c(0.4, 0.6))

  sim <- suppressMessages(scenarios_radiationModel_CT(
    model,
    t_now = "2022-01-15",
    t_hor = 2,
    nsim = 3,
    dt = 1,
    seed = 42
  ))

  expect_equal(nrow(sim), 3 * 3)
  expect_equal(sort(unique(sim$sim)), 1:3)
  expect_true(all(c("Yt", "Rt", "Rt_Yt", "gamma", "date_dt", "time") %in% names(sim)))
  expect_true(all(is.finite(sim$Yt)))
  expect_true(all(is.finite(sim$Rt)))
  expect_true(all(sim$Rt >= 0))
})

