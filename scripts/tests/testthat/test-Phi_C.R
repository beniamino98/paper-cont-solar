# ---
#' @description
#' Test `Phi_C()`.
#'
#' @tests
#' - vectorized transition products match explicit half-open products;
#' - row sums and non-negativity remain valid;
#' - propagated probabilities agree with explicit recursion;
#' - optional C wrapper matches the production R version.
# ---

test_that("Phi_C matches explicit half-open monthly transition products", {
  cases <- list(
    list(t_now = "2022-12-31", t_hor = "2023-01-01"),
    list(t_now = "2022-12-31", t_hor = "2023-01-02"),
    list(t_now = "2022-01-15", t_hor = "2022-01-20"),
    list(t_now = "2022-01-31", t_hor = "2022-02-02")
  )

  for (case in cases) {
    fix <- ctmc_month_boundary_fixture(case$t_now, case$t_hor)
    P <- Phi_C(fix$t_init, fix$t_end, fix$bounds)[[1]]
    expected <- ctmc_reference_phi(fix$t_init, fix$t_end, fix)

    expect_equal(unname(P), unname(expected), tolerance = 1e-9,
                 info = paste(case$t_now, "to", case$t_hor))
    expect_equal(drop(rowSums(P)), c(1, 1), tolerance = 1e-10)
    expect_true(all(P >= -1e-12), info = paste(case$t_now, "to", case$t_hor))
  }
})

test_that("p0 %*% Phi_C matches explicit recursive products", {
  fix <- ctmc_month_boundary_fixture("2022-12-31", "2023-01-02")
  p0 <- c(1, 0)

  got <- drop(p0 %*% Phi_C(fix$t_init, fix$t_end, fix$bounds)[[1]])
  expected <- drop(p0 %*% ctmc_reference_phi(fix$t_init, fix$t_end, fix))

  expect_equal(got, expected, tolerance = 1e-9)
  expect_true(all(got >= -1e-12))
  expect_equal(sum(got), 1, tolerance = 1e-10)
})

test_that("Phi_C optional C wrapper matches production R", {
  source(file.path(ctmc_radiation_dir, "ctmc-integrals-C.R"))
  Phi_C_R <- Phi_C

  fix <- ctmc_month_boundary_fixture("2022-12-31", "2023-01-02")
  tau_grid <- sort(unique(c(fix$bounds$n, fix$bounds$N,
                            seq(fix$t_init, fix$t_end, length.out = 7))))
  tau_grid <- tau_grid[tau_grid >= fix$t_init & tau_grid <= fix$t_end]
  phi_R <- Phi_C_R(fix$t_init, tau_grid, fix$bounds)

  source(file.path(ctmc_radiation_dir, "ctmc-integrals-C-wrappers.R"))
  phi_C <- Phi_C(fix$t_init, tau_grid, fix$bounds)

  for (i in seq_along(phi_R)) {
    expect_equal(unname(phi_C[[i]]), unname(phi_R[[i]]), tolerance = 1e-10)
  }

  source(file.path(ctmc_radiation_dir, "ctmc-integrals-C.R"))
})

