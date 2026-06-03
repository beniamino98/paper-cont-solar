# ---
#' @description
#' Test `radiationModel_CTMC_density()`.
#'
#' @tests
#' - production R density propagation returns valid density/CDF objects;
#' - density mass is normalized when requested;
#' - optional C density kernel matches production R on month boundaries;
#' - optional C density kernel matches production R with automatic grid
#'   construction.
# ---

test_that("radiationModel_CTMC_density returns valid production R density objects", {
  model <- ctmc_mock_ctmc_model(t_now = "2022-01-15", p0 = c(0.4, 0.6))
  y_grid <- seq(-5, 5, length.out = 121)

  dens <- radiationModel_CTMC_density(
    "2022-01-15", "2022-01-17", model,
    R0 = 50, y_grid = y_grid, dt = 0.5, normalize = TRUE
  )

  expect_equal(dens$y_grid, y_grid)
  expect_equal(dim(dens$f_state), c(length(y_grid), 2))
  expect_true(all(is.finite(dens$f_Y)))
  expect_true(all(dens$f_Y >= -1e-12))
  expect_equal(dens$mass, 1, tolerance = 1e-8)
  expect_true(all(diff(dens$cdf_Y(y_grid)) >= -1e-10))
  expect_true(all(dens$cdf_Y(y_grid) >= -1e-12 & dens$cdf_Y(y_grid) <= 1 + 1e-12))
})

test_that("radiationModel_CTMC_density optional C kernel matches production R on month boundaries", {
  source(file.path(ctmc_radiation_dir, "radiationModel_CTMC-R6.R"))
  density_R <- radiationModel_CTMC_density

  model <- ctmc_mock_ctmc_model(t_now = "2022-12-31", p0 = c(0.35, 0.65))
  model$sigma_bar <- function(d) 1 + 0.05 * sin(2 * pi * as.numeric(d) / 365)
  model$CTMC$params$Qm <- ctmc_month_list(
    ctmc_generator(0.04, 0.02),
    replacements = list(
      `1` = ctmc_generator(0.18, 0.05),
      `12` = ctmc_generator(0.07, 0.13)
    )
  )
  model$CTMC$params$mu <- ctmc_mu_list(
    c(mu1 = 0.20, mu0 = -0.05),
    replacements = list(
      `1` = c(mu1 = 0.40, mu0 = -0.15),
      `12` = c(mu1 = -0.25, mu0 = 0.30)
    )
  )
  model$CTMC$params$sig <- ctmc_sd_list(
    c(sd1 = 0.45, sd0 = 0.75),
    replacements = list(
      `1` = c(sd1 = 0.65, sd0 = 0.85),
      `12` = c(sd1 = 0.55, sd0 = 0.95)
    )
  )

  y_grid <- seq(-5, 5, length.out = 121)
  dens_R <- density_R(
    "2022-12-31", "2023-01-02", model,
    R0 = 50, y_grid = y_grid, dt = 0.5, normalize = TRUE
  )

  source(file.path(ctmc_radiation_dir, "radiationModel-CTMC-density-C-wrappers.R"))
  expect_false(identical(body(radiationModel_CTMC_density), body(density_R)))
  dens_C <- radiationModel_CTMC_density(
    "2022-12-31", "2023-01-02", model,
    R0 = 50, y_grid = y_grid, dt = 0.5, normalize = TRUE
  )

  expect_equal(dens_C$y_grid, dens_R$y_grid, tolerance = 0)
  expect_equal(dens_C$f_state, dens_R$f_state, tolerance = 1e-12)
  expect_equal(dens_C$f_Y, dens_R$f_Y, tolerance = 1e-12)
  expect_equal(dens_C$mass, dens_R$mass, tolerance = 1e-12)
  expect_equal(dens_C$pdf_Y(y_grid), dens_R$pdf_Y(y_grid), tolerance = 1e-12)
  expect_equal(dens_C$cdf_Y(y_grid), dens_R$cdf_Y(y_grid), tolerance = 1e-12)

  source(file.path(ctmc_radiation_dir, "radiationModel_CTMC-R6.R"))
})

test_that("radiationModel_CTMC_density optional C kernel matches production R with automatic grid", {
  source(file.path(ctmc_radiation_dir, "radiationModel_CTMC-R6.R"))
  density_R <- radiationModel_CTMC_density

  model <- ctmc_mock_ctmc_model(t_now = "2022-01-15", p0 = c(0.4, 0.6))
  model$sigma_bar <- function(d) rep(1, length(d))
  model$CTMC$params$Qm <- ctmc_month_list(ctmc_generator(0.12, 0.04))

  dens_R <- density_R(
    "2022-01-15", "2022-01-17", model,
    R0 = 50, n_grid = 101, dt = 0.5, normalize = FALSE
  )

  source(file.path(ctmc_radiation_dir, "radiationModel-CTMC-density-C-wrappers.R"))
  dens_C <- radiationModel_CTMC_density(
    "2022-01-15", "2022-01-17", model,
    R0 = 50, n_grid = 101, dt = 0.5, normalize = FALSE
  )

  probe <- seq(min(dens_R$y_grid), max(dens_R$y_grid), length.out = 25)
  expect_equal(dens_C$y_grid, dens_R$y_grid, tolerance = 0)
  expect_equal(dens_C$f_state, dens_R$f_state, tolerance = 1e-12)
  expect_equal(dens_C$f_Y, dens_R$f_Y, tolerance = 1e-12)
  expect_equal(dens_C$mass, dens_R$mass, tolerance = 1e-12)
  expect_equal(dens_C$pdf_Y(probe), dens_R$pdf_Y(probe), tolerance = 1e-12)
  expect_equal(dens_C$cdf_Y(probe), dens_R$cdf_Y(probe), tolerance = 1e-12)

  source(file.path(ctmc_radiation_dir, "radiationModel_CTMC-R6.R"))
})

