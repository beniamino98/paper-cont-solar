# ---
#' @description
#' Test `ctmc_monthly_scenarios()`.
#'
#' @tests
#' - simulated transition frequencies approximate `expm(dt * Q)`;
#' - transition frequencies depend on the previous state, preserving CTMC
#'   persistence.
# ---

test_that("ctmc_monthly_scenarios transition frequencies approximate expm(dt * Q)", {
  Q <- ctmc_generator(q12 = 0.25, q21 = 0.10)
  P <- expm::expm(Q)
  model <- ctmc_mock_ctmc_model(t_now = "2022-01-15", p0 = c(1, 0))
  model$CTMC$params$Qm <- ctmc_month_list(Q)

  set.seed(11)
  sim <- suppressMessages(ctmc_monthly_scenarios(
    model$CTMC,
    p0 = c(1, 0),
    month_idx = 1,
    nsim = 20000,
    dt = 1,
    nsteps = 1
  ))

  observed <- c(mean(sim$B_t[1, ] == 1), mean(sim$B_t[1, ] == 0))
  expect_equal(observed, unname(P[1, ]), tolerance = 0.015)
})

test_that("ctmc_monthly_scenarios persistence depends on previous state", {
  Q <- ctmc_generator(q12 = 0.08, q21 = 0.03)
  P <- expm::expm(Q)
  model <- ctmc_mock_ctmc_model(t_now = "2022-01-15", p0 = c(0.5, 0.5))
  model$CTMC$params$Qm <- ctmc_month_list(Q)

  set.seed(12)
  sim <- suppressMessages(ctmc_monthly_scenarios(
    model$CTMC,
    p0 = c(0.5, 0.5),
    month_idx = c(1, 1),
    nsim = 30000,
    dt = 1,
    nsteps = 2
  ))

  prev <- ifelse(sim$B_t[1, ] == 1, 1, 2)
  next_state <- ifelse(sim$B_t[2, ] == 1, 1, 2)
  freq_1 <- as.numeric(table(factor(next_state[prev == 1], levels = 1:2))) / sum(prev == 1)
  freq_2 <- as.numeric(table(factor(next_state[prev == 2], levels = 1:2))) / sum(prev == 2)

  expect_equal(freq_1, unname(P[1, ]), tolerance = 0.02)
  expect_equal(freq_2, unname(P[2, ]), tolerance = 0.02)
  expect_gt(freq_1[1], freq_2[1])
})

