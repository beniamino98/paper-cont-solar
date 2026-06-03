# ---
#' @description
#' Test `ctmc_is_transition_2state()`.
#'
#' @tests
#' - accepts valid two-state row-stochastic transition matrices;
#' - rejects non-matrix, wrong-dimension, non-finite, non-row-stochastic, and
#'   negative-probability candidates.
# ---

test_that("ctmc_is_transition_2state validates two-state transition matrices", {
  P <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, byrow = TRUE)

  expect_true(ctmc_is_transition_2state(P, quiet = TRUE))
  expect_false(ctmc_is_transition_2state(c(1, 2), quiet = TRUE))
  expect_false(ctmc_is_transition_2state(matrix(0, 3, 3), quiet = TRUE))
  expect_false(ctmc_is_transition_2state(matrix(c(0.9, 0.1, Inf, -Inf), 2), quiet = TRUE))
  expect_false(ctmc_is_transition_2state(matrix(c(0.5, 0.6, 0.2, 0.8), 2, byrow = TRUE), quiet = TRUE))
  expect_false(ctmc_is_transition_2state(matrix(c(1.1, -0.1, 0.2, 0.8), 2, byrow = TRUE), quiet = TRUE))
})

