#' Mean-Variance Supply and Demand Functions
#'
#' Build buyer demand, seller supply, and equilibrium price functions for a
#' contract payoff under quadratic mean-variance preferences.
#'
#' @param M_Gamma Numeric. Expected contract payoff.
#' @param v_Gamma Numeric. Buyer payoff variance.
#' @param S_R_Gamma Numeric. Covariance between buyer exposure and contract
#'   payoff.
#' @param v_Gamma_seller Numeric. Seller payoff variance.
#' @param r Numeric. Daily risk-free rate.
#' @param tau Numeric. Maturity in days.
#' @param w_Gamma Numeric. Contract scaling factor.
#'
#' @return List of functions `demand(P0, nu_b)`, `supply(P0, nu_s)`, and
#'   `price(nu_s, nu_b)`.
#' @export
supplyDemand_mv <- function(M_Gamma, v_Gamma, S_R_Gamma, v_Gamma_seller, r = 0, tau = 365, w_Gamma = 1){
  force(M_Gamma); force(v_Gamma); force(S_R_Gamma); force(w_Gamma)
  # Discount factor
  DtT <- exp(-r * tau)
  # Demand of the buyers
  demand <- function(P0, nu_b){
    # Scaled risk-adversion
    w_Gamma * ((M_Gamma - P0 * DtT / w_Gamma) / (nu_b * w_Gamma * v_Gamma) - S_R_Gamma / v_Gamma)
  }
  # Supply of the seller
  supply <- function(P0, nu_s){
    # Recaled factor 
    w_Gamma * (P0  * DtT / w_Gamma  -  M_Gamma) / (nu_s * w_Gamma * v_Gamma_seller) 
  }
  # Equilibrium price 
  price <- function(nu_s, nu_b){
    # Recaled factor 
    w_Gamma * ((M_Gamma / DtT) - (nu_s * nu_b * v_Gamma_seller) / (nu_b * v_Gamma + nu_s * v_Gamma_seller) * (S_R_Gamma / DtT))
  }
  
  list(
    demand = demand,
    supply = supply, 
    price = price
  )
}
