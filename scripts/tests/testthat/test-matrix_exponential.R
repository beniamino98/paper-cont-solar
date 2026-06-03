# ---
#' @description
#' Test `matrix_exponential()`.
#'
#' @tests
#' - agrees with `expm::expm(k * Q)` for valid two-state CTMC generators;
#' - returns the identity at zero horizon;
#' - preserves transition-matrix constraints;
#' - optional C wrapper gives the same values as the production R version.
# ---

test_that("matrix_exponential agrees with expm and preserves probabilities", {
  Q <- ctmc_generator(q12 = 0.12, q21 = 0.04)

  for (k in c(0, 0.1, 1, 7, 30)) {
    P <- matrix_exponential(k, Q)
    expected <- expm::expm(k * Q)

    expect_equal(unname(P), unname(expected), tolerance = 1e-10)
    expect_equal(drop(rowSums(P)), c(1, 1), tolerance = 1e-12)
    expect_true(all(P >= -1e-12))
  }
})

test_that("matrix_exponential returns identity at zero horizon", {
  Q <- ctmc_generator(q12 = 0.12, q21 = 0.04)
  expect_equal(matrix_exponential(0, Q), diag(1, 2, 2), tolerance = 1e-14)
})

test_that("matrix_exponential optional C wrapper matches production R", {
  source(file.path(ctmc_radiation_dir, "ctmc-integrals-C.R"))
  matrix_exponential_R <- matrix_exponential

  source(file.path(ctmc_radiation_dir, "ctmc-integrals-C-wrappers.R"))
  Q <- ctmc_generator(q12 = 0.12, q21 = 0.04)

  for (k in c(0, 0.1, 1, 7, 30)) {
    expect_equal(
      unname(matrix_exponential(k, Q)),
      unname(matrix_exponential_R(k, Q)),
      tolerance = 1e-12
    )
  }

  source(file.path(ctmc_radiation_dir, "ctmc-integrals-C.R"))
})

