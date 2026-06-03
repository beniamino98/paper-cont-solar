#' Summarise SoREdIDX Hedging Outcomes
#'
#' Compute theoretical and realized cashflow summaries for indexed SoREd
#' contracts under no, strip, or total futures hedging.
#'
#' @param sored SoREdIDX moment row/table.
#' @param scenario Scenario object returned by `scenarios_ER()`.
#' @param nu_b Numeric. Buyer risk aversion.
#' @param nu_s Numeric. Seller risk aversion.
#' @param w_Gamma Numeric. Contract scaling factor.
#' @param K_fun Strike function of day-of-year.
#' @param hedging Character. One of `"none"`, `"strip"`, or `"total"`.
#'
#' @return List with `pay`, `sim`, `day`, and `emp` summary tables.
#' @export
solarOption_hedging_soredidx <- function(sored, scenario, nu_b, nu_s, w_Gamma, K_fun, hedging = c("none", "strip", "total")){
  v_Gamma_hedged <- sored$v_E_Gamma
  S_ER_EGamma_hedged <- sored$S_ER_EGamma
  hedging <- match.arg(hedging, choices = c("none", "strip", "total"))
  if (hedging == "strip"){
    v_Gamma_hedged <- sored$v_E_Gamma_mid_E_strip
    S_ER_EGamma_hedged <- sored$S_ER_E_mid_E_strip
  } else if (hedging == "total"){
    v_Gamma_hedged <- sored$v_E_Gamma_mid_E_total
    S_ER_EGamma_hedged <- sored$S_ER_E_mid_E_total
  }
  
  # Supply and demands for SoREdIDX
  supply_demand_sored <- supplyDemand_mv(sored$M_E_Gamma, v_Gamma_hedged, S_ER_EGamma_hedged, v_Gamma_hedged, r = 0, tau = 365, w_Gamma = 1)
  # Equilibrium Price
  P_eq <- supply_demand_sored$price(nu_s, nu_b)
  # Equilibrium quantity
  Q_eq <- supply_demand_sored$supply(P_eq, nu_s)
  
  P_eq <- ifelse(Q_eq < 0, 0, P_eq)
  Q_eq <- ifelse(Q_eq < 0, 0, Q_eq)
  
  phi_seller <- phi_buyer_h <- phi_buyer_uh <- rep(0, length(Q_eq))
  if (hedging == "strip"){
    beta <- sored$beta_strip[[1]]
    alpha <- sored$alpha_strip[[1]]
    phi_seller <- Q_eq * beta
    phi_buyer_h <- alpha + phi_seller
    phi_buyer_uh <- alpha
  } else if (hedging == "total"){
    beta <- sored$beta_total[[1]]
    alpha <- sored$alpha_total[[1]]
    phi_seller <- Q_eq * beta
    phi_buyer_h <- alpha + phi_seller
    phi_buyer_uh <- alpha
  }
  
  # Realized radiation 
  RT <- scenario$data_sim$GHI[-1]
  # Realized electricity 
  ET <- scenario$data_sim$PUN[-1]
  # Strike price 
  K <- K_fun(scenario$data_sim$n[-1])
  # Realized payoff
  payoff <- ET * solarOption_payoff(RT, K)
  # Reference date
  dates <- seq.Date(sored$t_now+1, sored$t_hor, 1)
  # Daily premium for cash flows 
  M_E_Gamma <- apply(scenario$scenarios$ET_Q[-1,]*solarOption_payoff(scenario$scenarios$RT[-1,], K*matrix(1, nrow = nrow(scenario$scenarios$RT)-1, ncol = ncol(scenario$scenarios$RT))), 1, mean)
  Pt_eq <- ( M_E_Gamma/sum( M_E_Gamma)) * P_eq
  # Future prices
  FtT <- apply(scenario$scenarios$ET_Q[-1,], 1, mean)
  PnL <- ET - FtT 
  
  # 1) Derivative data 
  cf_pay <- tibble(
    Year = lubridate::year(dates[3]),
    payoff = w_Gamma * sum(payoff),
    Rt = w_Gamma * sum(RT),
    Et = w_Gamma * sum(ET),
    Rt_Et = w_Gamma * sum(ET*RT),
    P_eq = w_Gamma * P_eq,
    M_Gamma = w_Gamma * sored$M_E_Gamma,
    Q_eq = w_Gamma * Q_eq,
    premium = P_eq / M_Gamma
  )
  
  # 2) Moments of cash flows for buyers and sellers 
  cf_sim <- tibble(
    Year = lubridate::year(dates[1]),
    # Seller (unhedged)
    e_seller_uh = w_Gamma * Q_eq * (P_eq - sored$M_E_Gamma),
    v_seller_uh = (w_Gamma * Q_eq)^2 * sored$v_E_Gamma,
    mv_seller_uh = e_seller_uh - (nu_s/w_Gamma)/2 * v_seller_uh,
    # Seller (unhedged)
    e_seller_h = e_seller_uh,
    v_seller_h = (w_Gamma * Q_eq)^2 * sored$v_E_Gamma_mid_E_total,
    v_seller_strip = (w_Gamma * Q_eq)^2 * sored$v_E_Gamma_mid_E_strip,
    mv_seller_h = e_seller_h - (nu_s/w_Gamma)/2 * v_seller_h,
    # SPP (unhedged)
    e_buyer_uh = w_Gamma * sored$M_ER,
    v_buyer_uh = w_Gamma^2 * sored$v_ER,
    mv_buyer_uh = e_buyer_uh - (nu_b/w_Gamma)/2 * v_buyer_uh,
    # SPP (hedged)
    e_buyer_h = e_buyer_uh +  w_Gamma * Q_eq * (sored$M_E_Gamma - P_eq),
    v_buyer_h = v_buyer_uh + (w_Gamma * Q_eq)^2 * sored$v_E_Gamma  + 2 * w_Gamma^2 * Q_eq * sored$S_ER_EGamma,
    mv_buyer_h = e_buyer_h - (nu_b/w_Gamma)/2 * v_buyer_h,
    # Stats 
    mv_increase = (mv_buyer_h - mv_buyer_uh)/mv_buyer_uh*100,
    var_reduction = (v_buyer_h - v_buyer_uh)/v_buyer_uh*100
  )
  # Seasonal benchmark
  pi_seasonal <- w_Gamma * K * apply(scenario$scenarios$ET_Q[-1,], 1, mean)
  
  # 3) Realized daily cash flows 
  cf_day = tibble(
    date = dates,
    Year = lubridate::year(date),
    n = number_of_day(date),
    Gamma = w_Gamma * Q_eq * payoff,
    PnL = PnL, 
    Qt = w_Gamma * Q_eq,
    # Daily premium 
    Pt = w_Gamma * Q_eq * Pt_eq,
    # Premium 
    premium = P_eq/sored$M_E_Gamma-1,
    # Seller's cash flows 
    pi_seller_uh = Pt - w_Gamma * Q_eq * payoff,
    # Seller's cash flows 
    pi_seller_h = Pt - w_Gamma * Q_eq * payoff + w_Gamma * phi_seller * PnL,
    # Cumulated cash flows 
    cum_pi_seller_uh = cumsum(w_Gamma * Q_eq * (Pt_eq - payoff)),
    cum_ret_seller_uh = cum_pi_seller_uh / c(w_Gamma * Q_eq * P_eq),
    # Cumulated cash flows  
    cum_pi_seller_h = cumsum(w_Gamma * Q_eq * (Pt_eq-  payoff)) + cumsum(phi_seller * PnL),
    cum_ret_seller_h = cum_pi_seller_h/c(w_Gamma * Q_eq * P_eq),
    # Buyer's cash flows 
    pi_buyer_uh = w_Gamma * RT * ET - w_Gamma * phi_buyer_uh * PnL, 
    pi_buyer_h = w_Gamma * RT * ET + w_Gamma * Q_eq * (payoff - Pt_eq) - phi_buyer_h * PnL,
    # Sesonal production to benchmark ES 
    pi_seasonal = pi_seasonal
  )
  
  cf_emp <- cf_day %>%
    group_by(Year) %>%
    summarise(
      Pt = sum(Pt), 
      Qt = Q_eq,
      Gamma = sum(Gamma),
      premium = mean(premium),
      phi_buyer_h = sum(phi_buyer_h),
      phi_buyer_uh = sum(phi_buyer_uh),
      phi_seller = sum(phi_seller),
      Pi_seller_h = sum(pi_seller_h),
      Pi_seller_uh = sum(pi_seller_uh),
      Pi_buyer_uh = sum(pi_buyer_uh),
      Pi_buyer_h = sum(pi_buyer_h),
      # Variance buyer
      v_pi_buyer_uh = var(pi_buyer_uh),
      v_pi_buyer_h = var(pi_buyer_h),
      # Variance seller
      v_pi_seller_uh = var(pi_seller_uh),
      v_pi_seller_h = var(pi_seller_h),
      # Variance reduction 
      var_reduction_buyer = v_pi_buyer_h/v_pi_buyer_uh-1,
      var_reduction_seller = v_pi_seller_h/v_pi_seller_uh-1,
      # Expected shortfall when Pi is below seasonal level 
      ES_buyer_uh = mean(pmax(pi_seasonal - pi_buyer_uh, 0)),
      ES_buyer_h = mean(pmax(pi_seasonal - pi_buyer_h, 0)),
      # Expected shortfall decrease
      ES_reduction = ES_buyer_h/ES_buyer_uh - 1
    )
  
  structure(
    list(
      pay = cf_pay,
      sim = cf_sim, 
      day = cf_day,
      emp = cf_emp
    )
  )
}

#' Summarise SoRadIDX Hedging Outcomes
#'
#' Compute theoretical and realized cashflow summaries for indexed SoRad
#' contracts.
#'
#' @param sorad SoRadIDX moment row/table.
#' @param scenario Scenario object returned by `scenarios_ER()`.
#' @param nu_b Numeric. Buyer risk aversion.
#' @param nu_s Numeric. Seller risk aversion.
#' @param w_Gamma Numeric. Contract scaling factor.
#' @param K_fun Strike function of day-of-year.
#'
#' @return List with `pay`, `sim`, `day`, and `emp` summary tables.
#' @export
solarOption_hedging_soradidx <- function(sorad, scenario, nu_b, nu_s, w_Gamma = 1, K_fun){
  # Supply and demands for SoRadIDX
  supply_demand_sorad <- supplyDemand_mv(sorad$M_Gamma, sorad$v_Gamma, sorad$S_R_Gamma, sorad$v_Gamma, r = 0, tau = 365, w_Gamma = 1)
  # Equilibrium Price
  P_eq <- supply_demand_sorad$price(nu_s, nu_b)
  # Equilibrium quantity
  Q_eq <- supply_demand_sorad$supply(P_eq, nu_s)
  # Realized radiation 
  RT <- scenario$data_sim$GHI[-1]
  # Strike price 
  K <- K_fun(scenario$data_sim$n[-1])
  # Realized payoff
  payoff <- solarOption_payoff(RT, K)
  # Reference date
  dates <- seq.Date(sorad$t_now+1, sorad$t_hor, 1)
  # Daily premium for cash flows 
  M_Gamma <- apply(solarOption_payoff(scenario$scenarios$RT[-1,], K*matrix(1, nrow = nrow(scenario$scenarios$RT)-1, ncol = ncol(scenario$scenarios$RT))), 1, mean)
  # Compute a proxy for daily premiums rescaling the total premium according to M_Gamma
  Pt_eq <- (M_Gamma/sum(M_Gamma))* P_eq
  
  # 1) Derivative data 
  cf_pay <- tibble(
    Year = lubridate::year(dates[3]),
    payoff = w_Gamma * sum(payoff),
    Rt = w_Gamma * sum(RT),
    P_eq = w_Gamma * P_eq,
    M_Gamma = w_Gamma * sorad$M_Gamma,
    Q_eq = w_Gamma * Q_eq,
    premium = P_eq / M_Gamma
  )
  
  # 2) Moments of cash flows for buyers and sellers 
  cf_sim <- tibble(
    Year = lubridate::year(dates[4]),
    Pt = w_Gamma * P_eq,
    # Seller 
    e_seller = w_Gamma * Q_eq * (P_eq - sorad$M_Gamma),
    v_seller = (w_Gamma * Q_eq)^2 * sorad$v_Gamma,
    mv_seller = e_seller - (nu_s/w_Gamma)/2 * v_seller,
    # SPP (unhedged)
    e_buyer_uh = w_Gamma * sorad$M_R,
    v_buyer_uh = w_Gamma^2 * sorad$v_R,
    mv_buyer_uh = e_buyer_uh - (nu_b/w_Gamma)/2 * v_buyer_uh,
    # SPP (hedged)
    e_buyer_h = e_buyer_uh +  w_Gamma * Q_eq * (sorad$M_Gamma - P_eq),
    v_buyer_h = v_buyer_uh + (w_Gamma * Q_eq)^2 * sorad$v_Gamma  + 2 * w_Gamma^2 * Q_eq * sorad$S_R_Gamma,
    mv_buyer_h = e_buyer_h - (nu_b/w_Gamma)/2 * v_buyer_h,
    # Stats 
    mv_increase = (mv_buyer_h - mv_buyer_uh)/mv_buyer_uh*100,
    var_reduction = (v_buyer_h - v_buyer_uh)/v_buyer_uh*100
  )
  
  # 3) Realized daily cash flows 
  cf_day = tibble(
    date = dates,
    Year = lubridate::year(date),
    n = number_of_day(date),
    Gamma = w_Gamma * Q_eq * payoff,
    Qt = w_Gamma * Q_eq,
    # Daily premium 
    Pt = w_Gamma * Q_eq * Pt_eq,
    # Premium 
    premium = P_eq/sorad$M_Gamma-1,
    # Seller's cash flows 
    pi_seller = w_Gamma * Q_eq * (Pt_eq - payoff),
    # Cumulated cash flows 
    cum_pi_seller = w_Gamma * cumsum(Q_eq * (Pt_eq - payoff)),
    cum_ret_seller = cum_pi_seller/c(w_Gamma * Q_eq * P_eq),
    # Buyer's cash flows 
    pi_buyer_uh = w_Gamma * RT,
    pi_buyer_h = pi_buyer_uh + w_Gamma * Q_eq * (payoff - Pt_eq),
    pi_seasonal = w_Gamma * K
  )
  
  cf_emp <- cf_day %>%
    group_by(Year) %>%
    summarise(
      Pt = sum(Pt), 
      Qt = mean(Qt),
      Gamma = sum(w_Gamma * Q_eq * payoff),
      premium = mean(premium),
      Pi_seller = sum(pi_seller),
      Pi_buyer_uh = sum(pi_buyer_uh),
      Pi_buyer_h = sum(pi_buyer_h),
      # Variance reduction 
      v_pi_buyer_uh = var(pi_buyer_uh),
      v_pi_buyer_h = var(pi_buyer_h),
      var_reduction = v_pi_buyer_h/v_pi_buyer_uh-1,
      # Expected shortfalls 
      ES_buyer_uh = mean(pmax(pi_seasonal - pi_buyer_uh, 0)),
      ES_buyer_h = mean(pmax(pi_seasonal - pi_buyer_h, 0)),
      # Expected shortfall when Pi is below seasonal level 
      ES_buyer_uh = mean(pmax(pi_seasonal - pi_buyer_uh, 0)),
      ES_buyer_h = mean(pmax(pi_seasonal - pi_buyer_h, 0)),
      # Expected shortfall decrease
      ES_reduction = ES_buyer_h/ES_buyer_uh - 1
    )
  
  
  structure(
    list(
      pay = cf_pay,
      sim = cf_sim, 
      day = cf_day,
      emp = cf_emp
    )
  )
}

#' Summarise Daily SoRad Hedging Outcomes
#'
#' Compute theoretical and realized cashflow summaries for daily SoRad
#' contracts.
#'
#' @param sorad SoRad daily moment table.
#' @param scenario Scenario object returned by `scenarios_ER()`.
#' @param nu_b Numeric. Buyer risk aversion.
#' @param nu_s Numeric. Seller risk aversion.
#' @param w_Gamma Numeric. Contract scaling factor.
#' @param K_fun Strike function of day-of-year.
#'
#' @return List with `pay`, `sim`, `day`, and `emp` summary tables.
#' @export
solarOption_hedging_sorad <- function(sorad, scenario, nu_b, nu_s, w_Gamma = 1, K_fun){
  # Supply and demands for SoRadIDX
  supply_demand_sorad <- supplyDemand_mv(sorad$M_Gamma, sorad$v_Gamma, sorad$S_R_Gamma, sorad$v_Gamma, r = 0, tau = 365, w_Gamma = 1)
  # Equilibrium Price
  P_eq <- supply_demand_sorad$price(nu_s, nu_b)
  # Equilibrium quantity
  Q_eq <- supply_demand_sorad$supply(P_eq, nu_s)
  # Realized radiation 
  RT <- scenario$data_sim$GHI[-1]
  # Strike price 
  K <- K_fun(scenario$data_sim$n[-1])
  # Realized payoff
  Gamma <- solarOption_payoff(RT, K)
  # Reference date
  dates <- sorad$date
  # 1) Derivative data 
  cf_pay <- tibble(
    date = dates,
    Year = lubridate::year(date),
    n = number_of_day(date),
    payoff = Gamma,
    strike = K,
    Rt = RT,
    P_eq = P_eq,
    M_Gamma = sorad$M_Gamma,
    Q_eq = Q_eq,
    premium = P_eq / M_Gamma
  )
  # 2) Theoric moments for seller and buyer (hedged and unhedged)
  cf_sim <- tibble(
    date = dates,
    Year = lubridate::year(date),
    n = number_of_day(date),
    # Moments pi seller 
    e_seller = w_Gamma * Q_eq * (P_eq - sorad$M_Gamma),
    v_seller = (w_Gamma * Q_eq)^2 * sorad$v_Gamma,
    # Moments pi buyer (unhedged)
    e_buyer_uh = w_Gamma * sorad$M_R,
    v_buyer_uh = w_Gamma^2 * sorad$v_R,
    # Moments pi buyer (hedged)
    e_buyer_h = e_buyer_uh +  w_Gamma * Q_eq * (sorad$M_Gamma - P_eq),
    v_buyer_h = v_buyer_uh + (w_Gamma * Q_eq)^2 * sorad$v_Gamma  + 2 * w_Gamma^2 * Q_eq * sorad$S_R_Gamma,
    # Mean-variance utilities 
    mv_seller = e_seller - (nu_s/w_Gamma)/2 * v_seller,
    mv_buyer_uh = e_buyer_uh - (nu_b/w_Gamma)/2 * v_buyer_uh,
    mv_buyer_h = e_buyer_h - (nu_b/w_Gamma)/2 * v_buyer_h,
    # Increase in mean-variance 
    mv_increase = (mv_buyer_h - mv_buyer_uh)/mv_buyer_uh*100,
    # Decrease in variance (theoric)
    var_reduction = (v_buyer_h - v_buyer_uh)/v_buyer_uh*100
  )
  # 3) Realized daily cash flows for seller and buyer (hedged and unhedged)
  cf_day = tibble(
    date = dates,
    Year = lubridate::year(dates),
    n = number_of_day(dates),
    # Pre-scale payoff 
    payoff = w_Gamma * Q_eq * Gamma,
    Qt = w_Gamma * Q_eq,
    # Premium 
    premium = P_eq/sorad$M_Gamma-1,
    # Daily premium 
    Pt = w_Gamma * Q_eq * P_eq,
    # Seller's cash flows 
    pi_seller = Pt - w_Gamma * Q_eq * Gamma,
    # Cumulated cash flows 
    cum_pi_seller = cumsum(w_Gamma * Q_eq * (P_eq - Gamma)),
    # Cumulated net return 
    cum_ret_seller = cum_pi_seller / sum(w_Gamma * Q_eq * P_eq),
    # Buyer's cash flows (unhedged)
    pi_buyer_uh = w_Gamma * RT,
    # Buyer's cash flows (hedged)
    pi_buyer_h = pi_buyer_uh + w_Gamma * Q_eq * (Gamma - P_eq),
    # Seasonal cash flows for benchmark 
    pi_seasonal = w_Gamma * K,
  )
  # Aggregate daily data for each yer
  cf_emp <- cf_day %>%
    group_by(Year) %>%
    summarise(
      Pt = sum(Pt), 
      Qt = sum(Qt),
      Gamma = sum(payoff),
      premium = mean(premium),
      Pi_seller = sum(pi_seller),
      Pi_buyer_uh = sum(pi_buyer_uh),
      Pi_buyer_h = sum(pi_buyer_h),
      # Variance reduction 
      v_pi_buyer_uh = var(pi_buyer_uh),
      v_pi_buyer_h = var(pi_buyer_h),
      var_reduction = v_pi_buyer_h/v_pi_buyer_uh-1,
      # Expected shortfall when Pi is below seasonal level 
      ES_buyer_uh = mean(pmax(pi_seasonal - pi_buyer_uh, 0)),
      ES_buyer_h = mean(pmax(pi_seasonal - pi_buyer_h, 0)),
      # Expected shortfall decrease
      ES_reduction = ES_buyer_h/ES_buyer_uh - 1
    )
  
  structure(
    list(
      pay = cf_pay,
      sim = cf_sim, 
      day = cf_day,
      emp = cf_emp
    )
  )
}

#' Summarise Daily SoREd Hedging Outcomes
#'
#' Compute theoretical and realized cashflow summaries for daily SoREd
#' contracts with optional strip hedging.
#'
#' @param sored SoREd daily moment table.
#' @param scenario Scenario object returned by `scenarios_ER()`.
#' @param nu_b Numeric. Buyer risk aversion.
#' @param nu_s Numeric. Seller risk aversion.
#' @param w_Gamma Numeric. Contract scaling factor.
#' @param K_fun Strike function of day-of-year.
#' @param hedging Character. One of `"none"` or `"strip"`.
#'
#' @return List with `sim`, `day`, and `emp` summary tables.
#' @export
solarOption_hedging_sored <- function(sored, scenario, nu_b, nu_s, w_Gamma, K_fun, hedging = c("none", "strip")){
  print(sored$Year[1])
  # Type of hedging
  hedging <- match.arg(hedging, choices = c("none", "strip"))
  v_Gamma_hedged <- sored$v_E_Gamma
  S_ER_EGamma_hedged <- sored$S_ER_EGamma
  if (hedging == "strip"){
    v_Gamma_hedged <- sored$v_E_Gamma_mid_E
    S_ER_EGamma_hedged <- (sored$v_E * sored$S_ER_EGamma  - sored$S_EGamma_E * sored$S_ER_E)/(sored$v_E)
  }
  
  # Supply and demands for SoRadIDX
  supply_demand_sored <- supplyDemand_mv(sored$M_E_Gamma, v_Gamma_hedged, S_ER_EGamma_hedged, v_Gamma_hedged, r = 0, tau = 365, w_Gamma = 1)
  # Equilibrium Price
  P_eq <- supply_demand_sored$price(nu_s, nu_b)
  # Equilibrium quantity
  Q_eq <- supply_demand_sored$supply(P_eq, nu_s) 
  P_eq <- P_eq * ifelse(Q_eq < 0, 0, 1)
  Q_eq <- Q_eq * ifelse(Q_eq < 0, 0, 1)
  beta_buyer_uh <- beta_buyer_h <- rep(0, length(Q_eq))
  if (hedging == "strip"){
    # Future weights (buyer)
    beta_buyer_uh <- sored$S_ER_E / sored$v_E
    beta_buyer_h <- beta_buyer_uh + Q_eq * sored$S_EGamma_E / sored$v_E
  }
  
  # Realized radiation 
  RT <- scenario$data_sim$GHI[-1]
  # Realized electricity 
  ET <- scenario$data_sim$PUN[-1]
  # Strike price 
  K <- K_fun(scenario$data_sim$n[-1])
  # Realized payoff
  Gamma <- ET * solarOption_payoff(RT, K)
  # Reference date
  dates <- sored$date
  
  # 2) Moments of cash flows for buyers and sellers 
  cf_sim <- tibble(
    date = dates,
    Year = lubridate::year(dates),
    n = number_of_day(dates),
    # Seller (unhedged)
    e_seller_uh = w_Gamma * Q_eq * (P_eq - sored$M_E_Gamma),
    v_seller_uh = (w_Gamma * Q_eq)^2 * sored$v_E_Gamma,
    mv_seller_uh = e_seller_uh - (nu_s/w_Gamma)/2 * v_seller_uh,
    # Seller (unhedged)
    e_seller_h = e_seller_uh,
    v_seller_h = v_seller_uh + sored$beta^2 * sored$v_E,
    mv_seller_h = e_seller_h - (nu_s/w_Gamma)/2 * v_seller_h,
    # SPP (unhedged)
    e_buyer_uh = w_Gamma * sored$M_ER,
    v_buyer_uh = w_Gamma^2 * sored$v_ER + (w_Gamma * beta_buyer_uh)^2 * sored$v_E - 2 * beta_buyer_uh * w_Gamma *  sored$S_ER_E,
    mv_buyer_uh = e_buyer_uh - (nu_b/w_Gamma)/2 * v_buyer_uh,
    # SPP (hedged)
    e_buyer_h = w_Gamma * sored$M_ER +  w_Gamma * Q_eq * (sored$M_E_Gamma - P_eq),
    v_buyer_h = w_Gamma^2 * sored$v_ER + (w_Gamma * Q_eq)^2 * sored$v_E_Gamma_mid_E  + (w_Gamma * beta_buyer_h)^2 * sored$v_E - 2 * beta_buyer_h * w_Gamma * sored$S_ER_E  - 2 * beta_buyer_h * w_Gamma * Q_eq * sored$S_EGamma_E +  2 * w_Gamma^2 * Q_eq * sored$S_ER_EGamma ,
    mv_buyer_h = e_buyer_h - (nu_b/w_Gamma)/2 * v_buyer_h,
    # Stats 
    mv_increase = (mv_buyer_h - mv_buyer_uh)/mv_buyer_uh*100,
    var_reduction = (v_buyer_h - v_buyer_uh)/v_buyer_uh*100
  )
  
  # 3) Realized daily cash flows 
  cf_day = tibble(
    date = dates,
    Year = lubridate::year(date),
    n = number_of_day(date),
    payoff = w_Gamma * Q_eq * Gamma,
    Qt = w_Gamma * Q_eq,
    # Premium 
    premium = P_eq/sored$M_E_Gamma-1,
    # Daily premium 
    Pt = w_Gamma * Q_eq * P_eq,
    # Seller's cash flows 
    pi_seller_uh = Pt - w_Gamma * Q_eq * Gamma,
    # Cumulated cash flows 
    cum_pi_seller_uh = cumsum(w_Gamma * Q_eq * (P_eq-Gamma)),
    cum_ret_seller_uh = cum_pi_seller_uh/sum(Pt),
    # Seller's cash flows 
    pi_seller_h = Pt - w_Gamma * Q_eq * Gamma + w_Gamma * Q_eq * sored$PnL,
    # Cumulated cash flows 
    cum_pi_seller_h = cumsum(w_Gamma * Q_eq * (P_eq-Gamma)) + cumsum(w_Gamma * Q_eq * sored$PnL),
    cum_ret_seller_h = cum_pi_seller_h/sum(Pt),
    # Buyer's cash flows 
    pi_buyer_uh = w_Gamma * RT * ET - w_Gamma * beta_buyer_uh * sored$net_PnL,
    pi_buyer_h = w_Gamma * RT * ET + w_Gamma * Q_eq * Gamma - Pt - w_Gamma * beta_buyer_h * sored$net_PnL,
    # Sesonal production to benchmark ES 
    pi_seasonal = w_Gamma * K * sored$M_E 
  )
  
  cf_emp <- cf_day %>%
    group_by(Year) %>%
    summarise(
      Pt = sum(Pt), 
      Qt = sum(Qt),
      Gamma = sum(payoff),
      premium = mean(premium),
      phi_buyer_h = sum(beta_buyer_h),
      phi_buyer_uh = sum(beta_buyer_uh),
      phi_seller = sum(sored$beta),
      Pi_seller_h = sum(pi_seller_h),
      Pi_seller_uh = sum(pi_seller_uh),
      Pi_buyer_uh = sum(pi_buyer_uh),
      Pi_buyer_h = sum(pi_buyer_h),
      # Variance buyer
      v_pi_buyer_uh = var(pi_buyer_uh),
      v_pi_buyer_h = var(pi_buyer_h),
      # Variance seller
      v_pi_seller_uh = var(pi_seller_uh),
      v_pi_seller_h = var(pi_seller_h),
      # Variance reduction 
      var_reduction_buyer = v_pi_buyer_h/v_pi_buyer_uh-1,
      var_reduction_seller = v_pi_seller_h/v_pi_seller_uh-1,
      # Expected shortfall when Pi is below seasonal level 
      ES_buyer_uh = mean(pmax(pi_seasonal - pi_buyer_uh, 0)),
      ES_buyer_h = mean(pmax(pi_seasonal - pi_buyer_h, 0)),
      # Expected shortfall decrease
      ES_reduction = ES_buyer_h/ES_buyer_uh - 1
    )
  cf_emp
  
  structure(
    list(
      sim = cf_sim, 
      day = cf_day,
      emp = cf_emp
    )
  )
}



