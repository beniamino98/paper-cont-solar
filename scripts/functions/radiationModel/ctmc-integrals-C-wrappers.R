#' Optional C Kernels for CTMC Integral Utilities
#'
#' Source this file after `scripts/s0-load.R` to replace the production R
#' implementations of `matrix_exponential()`, `get_month_index_C()`, and
#' `Phi_C()` with C-backed versions. If this file is not sourced, the production
#' R implementations remain active.
#'
#' @keywords internal

.ctmc_kernel_root <- normalizePath(getwd(), mustWork = TRUE)
repeat {
  .ctmc_kernel_c <- file.path(
    .ctmc_kernel_root,
    "scripts", "functions", "C", "ctmc_integrals_kernels.c"
  )
  if (file.exists(.ctmc_kernel_c)) {
    break
  }
  .ctmc_kernel_parent <- dirname(.ctmc_kernel_root)
  if (identical(.ctmc_kernel_parent, .ctmc_kernel_root)) {
    stop("Cannot find scripts/functions/C/ctmc_integrals_kernels.c from current working directory.")
  }
  .ctmc_kernel_root <- .ctmc_kernel_parent
}

.ctmc_kernel_dll <- "ctmc_integrals_kernels"
.ctmc_kernel_so <- file.path(tempdir(), paste0(.ctmc_kernel_dll, .Platform$dynlib.ext))

if (!.ctmc_kernel_dll %in% names(getLoadedDLLs())) {
  if (!file.exists(.ctmc_kernel_so) ||
      file.info(.ctmc_kernel_so)$mtime < file.info(.ctmc_kernel_c)$mtime) {
    .ctmc_kernel_old_wd <- getwd()
    .ctmc_kernel_tmp_c <- file.path(tempdir(), "ctmc_integrals_kernels.c")
    file.copy(.ctmc_kernel_c, .ctmc_kernel_tmp_c, overwrite = TRUE)
    setwd(tempdir())
    .ctmc_kernel_status <- system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "SHLIB", "-o", .ctmc_kernel_so, .ctmc_kernel_tmp_c),
      stdout = TRUE,
      stderr = TRUE
    )
    setwd(.ctmc_kernel_old_wd)
    if (!file.exists(.ctmc_kernel_so)) {
      stop(
        "Unable to compile CTMC C kernels.\n",
        paste(.ctmc_kernel_status, collapse = "\n")
      )
    }
  }
  dyn.load(.ctmc_kernel_so)
}

#' Matrix Exponential for a Two-State CTMC Generator
#'
#' C-backed replacement for the production R version. Computes `expm(k * Q)`
#' for a valid two-state CTMC generator.
#'
#' @param k Numeric scalar. Non-negative time horizon.
#' @param Q Numeric `2 x 2` CTMC generator.
#'
#' @return Numeric `2 x 2` transition matrix.
#' @export
matrix_exponential <- function(k, Q) {
  .Call(
    "ctmc_matrix_exponential_call",
    as.numeric(k),
    as.matrix(Q),
    PACKAGE = "ctmc_integrals_kernels"
  )
}

#' Identify the Active Monthly CTMC Interval
#'
#' C-backed replacement for the production R version. Uses the half-open
#' convention `[t, T)`.
#'
#' @param tau Numeric vector of time points in day-of-year coordinates.
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#'
#' @return Numeric vector of month indices.
#' @keywords internal
get_month_index_C <- function(tau, bounds) {
  .Call(
    "ctmc_get_month_index_call",
    as.numeric(tau),
    bounds,
    PACKAGE = "ctmc_integrals_kernels"
  )
}

#' CTMC Transition Product Over Possibly Vectorized Intervals
#'
#' C-backed replacement for the production R version. Computes transition
#' products over one or more half-open intervals `[a, b)`.
#'
#' @param a Numeric scalar or vector of start times.
#' @param b Numeric scalar or vector of end times.
#' @param bounds CTMC bounds tibble from `create_bounds()`.
#'
#' @return List of `2 x 2` transition matrices.
#' @export
Phi_C <- function(a, b, bounds) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  n_out <- max(length(a), length(b))
  if (length(a) == 1) {
    a <- rep(a, n_out)
  }
  if (length(b) == 1) {
    b <- rep(b, n_out)
  }
  if (length(a) != n_out || length(b) != n_out) {
    stop("a and b must have the same length, unless one is scalar.")
  }
  .Call("ctmc_phi_call", a, b, bounds, PACKAGE = "ctmc_integrals_kernels")
}

rm(list = intersect(
  c(
    ".ctmc_kernel_c",
    ".ctmc_kernel_dll",
    ".ctmc_kernel_parent",
    ".ctmc_kernel_root",
    ".ctmc_kernel_so",
    ".ctmc_kernel_tmp_c",
    ".ctmc_kernel_old_wd",
    ".ctmc_kernel_status"
  ),
  ls()
))
