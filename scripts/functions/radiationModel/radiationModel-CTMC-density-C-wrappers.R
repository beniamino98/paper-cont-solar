#' Optional C Kernel for CTMC Radiation Density Propagation
#'
#' Source this file after `scripts/s0-load.R` to replace the production R
#' implementation of `radiationModel_CTMC_density()` with a C-backed propagation
#' loop. If this file is not sourced, the production R implementation remains
#' active.
#'
#' @keywords internal

.ctmc_density_root <- normalizePath(getwd(), mustWork = TRUE)
repeat {
  .ctmc_density_c <- file.path(
    .ctmc_density_root,
    "scripts", "functions", "C", "ctmc_density_kernel.c"
  )
  if (file.exists(.ctmc_density_c)) {
    break
  }
  .ctmc_density_parent <- dirname(.ctmc_density_root)
  if (identical(.ctmc_density_parent, .ctmc_density_root)) {
    stop("Cannot find scripts/functions/C/ctmc_density_kernel.c from current working directory.")
  }
  .ctmc_density_root <- .ctmc_density_parent
}

.ctmc_density_dll <- "ctmc_density_kernel"
.ctmc_density_so <- file.path(tempdir(), paste0(.ctmc_density_dll, .Platform$dynlib.ext))

if (!.ctmc_density_dll %in% names(getLoadedDLLs())) {
  if (!file.exists(.ctmc_density_so) ||
      file.info(.ctmc_density_so)$mtime < file.info(.ctmc_density_c)$mtime) {
    .ctmc_density_old_wd <- getwd()
    .ctmc_density_tmp_c <- file.path(tempdir(), "ctmc_density_kernel.c")
    file.copy(.ctmc_density_c, .ctmc_density_tmp_c, overwrite = TRUE)
    setwd(tempdir())
    .ctmc_density_status <- system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "SHLIB", "-o", .ctmc_density_so, .ctmc_density_tmp_c),
      stdout = TRUE,
      stderr = TRUE
    )
    setwd(.ctmc_density_old_wd)
    if (!file.exists(.ctmc_density_so)) {
      stop(
        "Unable to compile CTMC density C kernel.\n",
        paste(.ctmc_density_status, collapse = "\n")
      )
    }
  }
  dyn.load(.ctmc_density_so)
}

#' Approximate CTMC Density for the Radiation DTMC
#'
#' C-backed replacement for the production R implementation. R handles model
#' extraction and grid construction; C handles the time-stepping Gaussian density
#' propagation.
#'
#' @param t_now Date or character scalar. Conditioning date.
#' @param t_hor Date or character scalar. Horizon date.
#' @param model_Rt A `radiationModel_CTMC` object.
#' @param R0 Optional numeric initial GHI at `t_now`.
#' @param y_grid Optional numeric grid for transformed radiation.
#' @param n_grid Integer. Grid size used when `y_grid` is not supplied.
#' @param dt Numeric. CTMC/kernel time step in days.
#' @param y_sd_mult Numeric. Width multiplier for automatic grid construction.
#' @param normalize Logical. If `TRUE`, renormalize density mass after each
#'   time step.
#'
#' @return List containing state densities, total density/CDF closures, grid,
#'   initial condition, and diagnostics.
#' @export
radiationModel_CTMC_density <- function(
    t_now,
    t_hor,
    model_Rt,
    R0 = NULL,
    y_grid = NULL,
    n_grid = 1201,
    dt = 1/30,
    y_sd_mult = 8,
    normalize = TRUE
) {
  stopifnot(any(class(model_Rt) %in% "radiationModel_CTMC"))
  t_now <- as.Date(t_now)
  t_hor <- as.Date(t_hor)
  tau <- as.numeric(difftime(t_hor, t_now, units = "days"))
  if (tau <= 0) {
    stop("t_hor must be after t_now.")
  }

  theta <- model_Rt$theta
  CTMC <- model_Rt$CTMC
  mu <- CTMC$params$mu
  sig <- CTMC$params$sig
  if (!is.null(CTMC$params$Qm)) {
    Q <- CTMC$params$Qm
  } else {
    Q <- transition_list_to_generator_2state(CTMC$params$Pm)
  }

  if (is.null(R0)) {
    R0 <- dplyr::filter(model_Rt$model$data, date == t_now)$GHI
  }

  C0 <- model_Rt$Ct(t_now)
  Y0 <- model_Rt$model$spec$transform$RY(R0, C0)
  p0 <- as.numeric(CTMC$alpha[CTMC$data$date == t_now, ])
  p0 <- p0 / sum(p0)

  s_grid <- seq(0, tau, by = dt)
  if (tail(s_grid, 1) < tau) {
    s_grid <- c(s_grid, tau)
  }

  n0 <- number_of_day(t_now)
  Ybar_fun <- function(s) model_Rt$Yt_bar(n0 + s)
  sigbar_fun <- function(s) model_Rt$sigma_bar(n0 + s)
  month_fun <- function(s) lubridate::month(t_now + floor(s + 1e-12))

  if (is.null(y_grid)) {
    s_tmp <- seq(0, tau, length.out = 200)
    Ybar_vals <- Ybar_fun(s_tmp)
    sigbar_vals <- sigbar_fun(s_tmp)
    max_sig <- max(unlist(sig), na.rm = TRUE)
    approx_sd <- sqrt(sum((sigbar_vals * max_sig)^2) * tau / length(s_tmp))
    y_grid <- seq(
      min(Y0, Ybar_vals, na.rm = TRUE) - y_sd_mult * approx_sd,
      max(Y0, Ybar_vals, na.rm = TRUE) + y_sd_mult * approx_sd,
      length.out = n_grid
    )
  }

  if (length(y_grid) < 2) {
    stop("y_grid must have at least two points.")
  }
  dy <- y_grid[2] - y_grid[1]
  if (!is.finite(dy) || dy <= 0 || max(abs(diff(y_grid) - dy)) > 1e-8) {
    stop("y_grid must be strictly increasing and evenly spaced.")
  }

  s0 <- head(s_grid, -1)
  s1 <- tail(s_grid, -1)
  ds <- s1 - s0
  month_idx <- month_fun(s0)
  P_list <- lapply(seq_along(ds), function(k) {
    matrix_exponential(ds[k], Q[[month_idx[k]]])
  })
  P_step <- t(vapply(P_list, function(P) {
    c(P[1, 1], P[1, 2], P[2, 1], P[2, 2])
  }, numeric(4)))
  mu_step <- t(vapply(seq_along(ds), function(k) {
    as.numeric(mu[[month_idx[k]]])
  }, numeric(2)))
  sig_step <- t(vapply(seq_along(ds), function(k) {
    as.numeric(sig[[month_idx[k]]])
  }, numeric(2)))

  f <- .Call(
    "ctmc_density_propagate_call",
    as.numeric(y_grid),
    as.numeric(Y0[1]),
    as.numeric(p0),
    as.numeric(ds),
    P_step,
    as.numeric(Ybar_fun(s0)),
    as.numeric(Ybar_fun(s1)),
    as.numeric(sigbar_fun(s1)),
    mu_step,
    sig_step,
    as.numeric(theta),
    as.logical(normalize),
    PACKAGE = "ctmc_density_kernel"
  )

  fY <- rowSums(f)
  cdfY <- cumsum(fY) * dy
  cdfY <- pmin(pmax(cdfY, 0), 1)

  list(
    y_grid = y_grid,
    f_state = f,
    f_Y = fY,
    pdf_Y = approxfun(y_grid, fY, yleft = 0, yright = 0, rule = 2),
    cdf_Y = approxfun(y_grid, cdfY, yleft = 0, yright = 1, rule = 2),
    mass = sum(fY) * dy,
    Y0 = Y0,
    p0 = p0,
    t_now = t_now,
    t_hor = t_hor,
    dt = dt
  )
}

rm(list = intersect(
  c(
    ".ctmc_density_c",
    ".ctmc_density_dll",
    ".ctmc_density_parent",
    ".ctmc_density_root",
    ".ctmc_density_so",
    ".ctmc_density_tmp_c",
    ".ctmc_density_old_wd",
    ".ctmc_density_status"
  ),
  ls()
))
