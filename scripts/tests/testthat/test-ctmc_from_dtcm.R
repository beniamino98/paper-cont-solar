# ---
#' @description
#' Test `ctmc_from_dtcm()`.
#'
#' @tests
#' - preserves the fitted DTMC object structure;
#' - recomputes `params$Qm` from `params$Pm`;
#' - leaves the original transition matrices available for auditability.
# ---

test_that("ctmc_from_dtcm recomputes Qm from Pm", {
  Pm <- ctmc_month_list(
    matrix(c(0.89, 0.11,
             0.04, 0.96),
           2, byrow = TRUE)
  )
  dtmc <- list(params = list(Pm = Pm, Qm = ctmc_month_list(matrix(99, 2, 2))))

  ctmc <- ctmc_from_dtcm(dtmc)

  expect_length(ctmc$params$Qm, length(Pm))
  expect_equal(ctmc$params$Pm, Pm)
  expect_false(identical(ctmc$params$Qm[[1]], dtmc$params$Qm[[1]]))
  expect_equal(unname(expm::expm(ctmc$params$Qm[[1]])), unname(Pm[[1]]), tolerance = 1e-12)
})

