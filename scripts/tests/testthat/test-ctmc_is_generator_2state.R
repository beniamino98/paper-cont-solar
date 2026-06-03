# ---
#' @description
#' Test `ctmc_is_generator_2state()`.
#'
#' @tests
#' - accepts valid two-state CTMC generators;
#' - rejects non-matrix, wrong-dimension, non-finite, non-zero-row-sum, and
#'   negative-off-diagonal candidates.
# ---

test_that("ctmc_is_generator_2state validates two-state generators", {
  expect_true(ctmc_is_generator_2state(ctmc_generator(), quiet = TRUE))
  expect_false(ctmc_is_generator_2state(c(1, 2), quiet = TRUE))
  expect_false(ctmc_is_generator_2state(matrix(0, 3, 3), quiet = TRUE))
  expect_false(ctmc_is_generator_2state(matrix(c(-1, 1, Inf, -Inf), 2), quiet = TRUE))
  expect_false(ctmc_is_generator_2state(matrix(c(-0.1, 0.2, 0.1, -0.1), 2, byrow = TRUE), quiet = TRUE))
  expect_false(ctmc_is_generator_2state(matrix(c(0.1, -0.1, 0.2, -0.2), 2, byrow = TRUE), quiet = TRUE))
})

