#' Exponential-Decay Flexible Probabilities
#'
#' Compute normalized exponential-decay weights around a target observation.
#'
#' @param t_bar Integer scalar. Number of observations.
#' @param tau_hl Numeric scalar. Half-life parameter.
#' @param t_star Numeric scalar. Target observation index.
#'
#' @return Numeric vector of probabilities summing to one.
#' @export
exp_decay_fp <- function(t_bar, tau_hl = 1, t_star = t_bar){
  t <- 1:t_bar
  p <- exp(-log(2)/tau_hl*abs(t_star - t))
  p <- p/sum(p)
  return(p)
}
