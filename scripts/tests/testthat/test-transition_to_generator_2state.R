# ---
#' @description
#' Test `transition_to_generator_2state()`.
#'
#' @tests
#' - embeds a regular two-state transition matrix into a CTMC generator;
#' - maps identity transitions to the zero generator;
#' - rejects non-stochastic or non-embeddable transition matrices.
# ---

test_that("transition_to_generator_2state embeds a one-step transition matrix", {
  P <- matrix(c(0.88, 0.12,
                0.04, 0.96),
              nrow = 2, byrow = TRUE)
  Q <- transition_to_generator_2state(P)

  expect_equal(unname(expm::expm(Q)), unname(P), tolerance = 1e-12)
  expect_equal(rowSums(Q), c(0, 0), tolerance = 1e-12)
  expect_true(Q[1, 2] >= 0)
  expect_true(Q[2, 1] >= 0)
})

test_that("transition_to_generator_2state maps identity to zero generator", {
  expect_equal(
    transition_to_generator_2state(diag(1, 2, 2)),
    matrix(0, 2, 2),
    tolerance = 1e-14
  )
})

test_that("transition_to_generator_2state rejects invalid inputs", {
  expect_error(
    transition_to_generator_2state(matrix(c(0.5, 0.6, 0.2, 0.8), 2, byrow = TRUE)),
    "valid transition"
  )
  expect_error(
    transition_to_generator_2state(matrix(c(0.1, 0.9, 0.3, 0.7), 2, byrow = TRUE)),
    "not embeddable"
  )
})

