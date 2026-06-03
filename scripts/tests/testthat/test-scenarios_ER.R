# ---
#' @description
#' Test `scenarios_ER()`.
#'
#' @tests
#' - preprocessing builds the joined daily electricity/radiation grid;
#' - residual simulation returns matrices with expected dimensions;
#' - the end-to-end wrapper returns radiation, electricity, and futures-net
#'   scenario matrices.
# ---

test_that("scenarios_ER preprocessing and end-to-end wrapper are coherent", {
  inputs <- ctmc_mock_joint_scenario_inputs(
    t_now = "2022-01-15",
    t_hor = "2022-01-18",
    p0 = c(0.4, 0.6)
  )

  preproc <- scenarios_ER_proproc(
    inputs$model_Et,
    inputs$model_Rt,
    inputs$rho,
    t_now = "2022-01-15",
    t_hor = "2022-01-18"
  )

  expect_equal(nrow(preproc$data_sim), 4)
  expect_true(all(c("GHI", "Yt", "PUN", "Xt", "Ct", "Yt_bar", "sigma_bar") %in% names(preproc$data_sim)))
  expect_equal(as.numeric(preproc$p0), c(0.4, 0.6), tolerance = 1e-12)

  residuals <- suppressMessages(scenarios_ER_residuals(nsim = 4, preproc = preproc, seed = 7))
  expect_equal(dim(residuals$dMt), c(4, 4))
  expect_equal(dim(residuals$dWt), c(4, 4))
  expect_equal(dim(residuals$B_t), c(4, 4))

  out <- suppressMessages(scenarios_ER(
    inputs$model_Et,
    inputs$model_Rt,
    inputs$rho,
    t_now = "2022-01-15",
    t_hor = "2022-01-18",
    nsim = 4,
    seed = 7
  ))

  expect_equal(dim(out$scenarios$RT), c(4, 4))
  expect_equal(dim(out$scenarios$ET_P), c(4, 4))
  expect_equal(dim(out$scenarios$ET_Q), c(4, 4))
  expect_equal(dim(out$scenarios$F_net_P), c(4, 4))
  expect_equal(dim(out$scenarios$F_net_Q), c(4, 4))
  expect_true(all(is.finite(out$scenarios$RT)))
  expect_true(all(is.finite(out$scenarios$ET_Q)))
})

