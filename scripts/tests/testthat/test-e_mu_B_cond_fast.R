# ---
#' @description
#' Test `e_mu_B_cond_fast()`.
#'
#' @tests
#' - conditional Markov-bridge drift expectations match explicit bridge
#'   formulas;
#' - probability-weighted conditional expectations recover the unconditional
#'   expectation.
# ---

test_that("e_mu_B_cond_fast matches explicit bridge formulas", {
  fix <- ctmc_fixture()
  s <- fix$t_init + c(0.2, 2.0, 4.8)
  p_T <- ctmc_reference_pT(fix)

  for (state in seq_len(2)) {
    ei <- diag(1, 2, 2)[, state]

    expect_equal(
      e_mu_B_cond_fast(s, fix$bounds, fix$p0, fix$mu, ei = ei, p_T = p_T[state]),
      ctmc_reference_state_value_cond(s, fix, fix$mu, ei),
      tolerance = 1e-8
    )
  }
})

test_that("e_mu_B_cond_fast satisfies total expectation", {
  fix <- ctmc_fixture()
  s <- fix$t_init + c(0.2, 2.0, 4.8)
  p_T <- drop(fix$p0 %*% Phi_C(fix$t_init, fix$t_end, fix$bounds)[[1]])

  mu_cond_1 <- e_mu_B_cond_fast(s, fix$bounds, fix$p0, fix$mu,
                                ei = c(1, 0), p_T = p_T[1])
  mu_cond_0 <- e_mu_B_cond_fast(s, fix$bounds, fix$p0, fix$mu,
                                ei = c(0, 1), p_T = p_T[2])

  expect_equal(
    p_T[1] * mu_cond_1 + p_T[2] * mu_cond_0,
    e_mu_B_fast(s, fix$bounds, fix$p0, fix$mu),
    tolerance = 1e-9
  )
})

