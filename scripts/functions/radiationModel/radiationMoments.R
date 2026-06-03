#' Dispatch Radiation Moment Computation
#'
#' Compute radiation moments using the HMM implementation when `model_Rt` is a
#' `radiationModelHMM`, otherwise use the IID-mixture `radiationModel`
#' implementation.
#'
#' @param t_now Date or character scalar. Conditioning date.
#' @param t_hor Date or character vector. Horizon date(s).
#' @param model_Rt A `radiationModel` or `radiationModelHMM` object.
#' @param R0 Optional numeric initial GHI at `t_now`.
#' @param maxEval Integer. Maximum number of integrand evaluations for HMM
#'   numerical integration.
#' @param tol Numeric. Numerical integration tolerance for HMM moments.
#'
#' @return A tibble returned by `radiationModel_moments()` or
#'   `radiationModelHMM_moments()`.
#' @export
radiationMoments <- function(t_now, t_hor, model_Rt, R0 = NULL, maxEval = 50000, tol = 0.0001){
  if (class(model_Rt)[1] == "radiationModel_CTMC") {
    mom <- radiationModel_CTMC_moments(t_now = t_now, t_hor = t_hor, model_Rt = model_Rt, R0 = R0, maxEval = maxEval, tol = tol)
  } else if (class(model_Rt)[1] == "radiationModel") {
    mom <- radiationModel_moments(t_now = t_now, t_hor = t_hor, model_Rt = model_Rt, R0 = R0)
  }
  return(mom)
}
