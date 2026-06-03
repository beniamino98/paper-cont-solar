#' Preprocess Inputs for Joint Electricity-Radiation Scenarios
#'
#' Prepare the daily simulation grid, observed data, deterministic radiation
#' components, initial DTMC probabilities, and model references used by
#' `scenarios_ER()`.
#' 
#' @param model_Et Electricity model object.
#' @param model_Rt Radiation CTMC model object.
#' @param rho List with 12 monthly state-correlation vectors.
#' @param t_now Date or character scalar. Conditioning date.
#' @param t_hor Date or character scalar. Horizon date.
#'
#' @return A preprocessing list used internally by the scenario pipeline.
#' @keywords internal
scenarios_ER_proproc <- function(model_Et, model_Rt, rho, t_now, t_hor){
  t_now <- as.Date(t_now)
  t_hor <- as.Date(t_hor)
  # Grid in days relative to start (0, dt, 2dt, ..., t_hor)
  time <- seq.Date(t_now, t_hor, by = 1)
  nsteps <- length(time) - 1   # number of *steps* (increments)
  # Join true data only on daily dates
  data_sim <- dplyr::inner_join(
    dplyr::tibble(date = as.Date(time), n =  number_of_day(time)),
    dplyr::select(model_Rt$model$data, date, Year, Month, Day, GHI, Yt),
    by = "date") %>% 
    dplyr::inner_join(dplyr::select(model_Et$data, date, PUN, Xt), by = "date") 
  # Pre-compute deterministic variables 
  # Seasonal clear-sky 
  data_sim$Ct <- model_Rt$Ct(data_sim$n)
  # Seasonal mean Yt 
  data_sim$Yt_bar <- model_Rt$Yt_bar(data_sim$n)
  # Seasonal mean Rt
  data_sim$GHI_bar <- model_Rt$model$spec$transform$iRY(data_sim$Yt_bar, data_sim$Ct)
  # Seasonal std. deviation Yt 
  data_sim$sigma_bar <- sqrt(model_Rt$model$spec$seasonal.variance$predict(data_sim$n))
  # Initial probabilities 
  p0 <- model_Rt$CTMC$alpha[model_Rt$CTMC$data$date == t_now, ]
  
  list(
    model_Rt = model_Rt,
    model_Et = model_Et,
    rho = rho,
    p0 = p0, 
    data_sim = data_sim
  )
}

#' Simulate Joint Electricity and Radiation Residuals
#'
#' Simulate DTMC regimes, radiation shocks, and electricity shocks with
#' monthly/state-dependent correlations.
#' 
#' @param nsim Integer. Number of Monte Carlo paths.
#' @param preproc List returned by `scenarios_ER_proproc()`.
#' @param seed Integer random seed.
#' 
#' @return Updated preprocessing list containing simulated residual matrices.
#' @keywords internal
scenarios_ER_residuals <- function(nsim = 1, preproc, seed = 1){
  set.seed(seed)
  data_sim <- preproc$data_sim
  model_Rt <- preproc$model_Rt
  # CTMC parameters 
  mu <- model_Rt$CTMC$params$mu
  sd <- model_Rt$CTMC$params$sig
  Pm <- model_Rt$CTMC$params$Pm
  rho <- preproc$rho
  # Initialize a matrix of residuals  
  nsteps <- nrow(data_sim)
  B_t <- matrix(0, nrow = nsteps, ncol = nsim)
  # Residuals for solar model 
  dMt <- matrix(0, nrow = nsteps, ncol = nsim)
  # Residuals for electricity model 
  dWt <- matrix(0, nrow = nsteps, ncol = nsim)
  # Drift for solar model 
  mu_B <- matrix(0, nrow = nsteps, ncol = nsim)
  sigma_B <- matrix(0, nrow = nsteps, ncol = nsim)
  # Initial probabilities 
  p0 <- preproc$p0
  p1 <- matrix(0, nrow = nsteps, ncol = nsim)
  # Monthly index
  tm_t <- data_sim$Month
  sim <- 1
  # Simulation Bernoulli 
  for (sim in 1:nsim) {
    message(sim, "/", nsim, paste0(" (", round((sim/nsim)*100, 2), " %)"), "\r", appendLF = FALSE)
    # Initialize the probabilities 
    p_t <- p0
    p1[1, sim] <- p_t[1]
    # Initialize state 
    state <- sample(1:2, size = 1, prob = p0)
    t <- 2
    for(t in 2:nsteps){
      t_m <- tm_t[t]
      # Evolve probability depending on the state
      p_t <- Pm[[t_m]][state, ]
      p1[t, sim] <- p_t[1]
      # Next step state 
      state <- sample(1:2, size = 1, prob = p_t)
      # Evolve Markov Chain
      B_t[t, sim] <- ifelse(state == 1, 1, 0) 
      # Joint correlation matrix 
      Sigma_XY <- diag(1, 3, 3)
      # Cor(E, R1)
      Sigma_XY[1,2] <- Sigma_XY[2,1] <- rho[[t_m]][1]
      # Cor(E, R0)
      Sigma_XY[1,3] <- Sigma_XY[3,1] <- rho[[t_m]][2]
      # Simulate joint residuals 
      sim_joint <- mvtnorm::rmvnorm(1, sigma = Sigma_XY)
      dWt[t,sim] <- sim_joint[,1]
      dMt[t,sim] <- sim_joint[,2] * B_t[t,sim] + sim_joint[,3] * (1-B_t[t,sim])
    }
    # Simulated drift
    mu_B[,sim] <- purrr::map2_dbl(tm_t, B_t[,sim], ~mu[[.x]][ifelse(.y == 0, 2, .y)])
    sigma_B[,sim] <- purrr::map2_dbl(tm_t, B_t[,sim], ~sd[[.x]][ifelse(.y == 0, 2, .y)])
  }
  # Number of simulations 
  preproc$nsim <- nsim 
  # Store simulated bernoulli's 
  preproc$dMt <- dMt
  preproc$dWt <- dWt
  preproc$mu_B <- mu_B
  preproc$sigma_B <- sigma_B
  preproc$B_t <- B_t
  return(preproc)
}

#' Propagate Joint Electricity and Radiation Scenarios
#'
#' Convert simulated residuals into daily paths for radiation, electricity, and
#' forward-price deviations.
#' 
#' @param preproc List returned by `scenarios_ER_residuals()`.
#'
#' @return Updated preprocessing list with `scenarios$RT`, `scenarios$ET`, and
#'   `scenarios$FT` matrices.
#' @keywords internal
scenarios_ER_filter <- function(preproc){
  # Extract info 
  data_sim <- preproc$data_sim
  nsteps <- nrow(data_sim)
  dMt <- preproc$dMt
  dWt <- preproc$dWt
  mu_B <- preproc$mu_B
  sigma_B <- preproc$sigma_B
  nsim <- ncol(dMt)
  
  # Radiation model 
  model_Rt <- preproc$model_Rt
  theta <- exp(-model_Rt$theta)
  # Electricity model 
  model_Et <- preproc$model_Et
  lam <- model_Et$model$lambda
  kappa <- model_Et$model$kappa
  mu_X <- model_Et$model$mu
  sigma_X <- model_Et$model$sigma
  
  # Compute future prices for each maturity 
  F0T <- c(data_sim$PUN[1], model_Et$F_E(data_sim$date[1], data_sim$date[-1], data_sim$PUN[1]))
  
  # Initialize matrixes for scenarios 
  RT <- matrix(0, nsteps, ncol = nsim)
  ET_P <- matrix(0, nsteps, ncol = nsim)
  ET_Q <- matrix(0, nsteps, ncol = nsim)
  F_net_P   <- matrix(0, nsteps, ncol = nsim)
  F_net_Q   <- matrix(0, nsteps, ncol = nsim)
  # Simulations
  sim <- 1
  for (sim in 1:nsim) {
    Xt_P <- c(data_sim$Xt[1], numeric(nsteps-1))
    Xt_Q <- c(data_sim$Xt[1], numeric(nsteps-1))
    Yt <- c(data_sim$Yt[1], numeric(nsteps-1))
    i <- 1
    for (i in 1:(nsteps-1)) {
      # **************************************************
      # Drift of Yt 
      mu_Y <- data_sim$Yt_bar[i+1] + theta * (Yt[i] - data_sim$Yt_bar[i]) + data_sim$sigma_bar[i+1] * mu_B[i+1,sim]
      # Diffusion 
      sigma_Y <- data_sim$sigma_bar[i+1] * sigma_B[i+1, sim] * dMt[i+1, sim]
      # Increment 
      Yt[i+1] <- mu_Y + sigma_Y
      # Moments of X
      M_X_P <- Xt_P[i] * exp(-kappa) + 1/kappa * (kappa * mu_X - sigma_X^2/2) * (1 - exp(-kappa))
      M_X_Q <- Xt_Q[i] * exp(-kappa) + 1/kappa * (kappa * mu_X + lam * sigma_X - sigma_X^2/2) * (1 - exp(-kappa))
      S2_X  <- sigma_X^2 * (1 - exp(-2*kappa))/(2*kappa)
      # Increment 
      Xt_Q[i+1] <- M_X_Q + sqrt(S2_X) * dWt[i+1, sim]
      Xt_P[i+1] <- M_X_P + sqrt(S2_X) * dWt[i+1, sim]
      # **************************************************
      #M_X <-  Xt[i] * exp(-kappa * (nsteps - i):i  ) + 1/kappa * (kappa * mu_X + lam * sigma_X - sigma_X^2/2) * (1 - exp(-kappa * (nsteps - i):i))
      #S2_X <-  sigma_X^2 * (1 - exp(-2*kappa * (nsteps - i):i))/(2*kappa)
      #Ft[i+1] <- exp(M_X + 0.5*S2_X) * (1 + sigma_X * exp(-kappa) * dWt[i+1, sim])
    }
    # Next step value of Rt from exact solution  
    RT[,sim] <- model_Rt$model$spec$transform$iRY(Yt, data_sim$Ct) 
    ET_P[,sim] <- exp(Xt_P)
    ET_Q[,sim] <- exp(Xt_Q)
    F_net_Q[,sim] <- ET_Q[,sim] - F0T
    F_net_P[,sim] <- ET_P[,sim] - F0T
  }
  preproc$scenarios <- list(ET_Q = ET_Q, ET_P = ET_P, RT = RT,
                            F_net_Q = F_net_Q, F_net_P = F_net_P)
  preproc
}


#' Simulate Joint Electricity-Radiation Scenarios
#'
#' Generate daily Monte Carlo paths for electricity and solar radiation using
#' the fitted electricity model, radiation CTMC model, and monthly
#' state-dependent residual correlations.
#' 
#' @param model_Et Electricity model object.
#' @param model_Rt Radiation CTMC model object.
#' @param rho List with 12 monthly state-correlation vectors.
#' @param t_now Date or character scalar. Conditioning date.
#' @param t_hor Date or character scalar. Horizon date.
#' @param nsim Integer. Number of Monte Carlo paths.
#' @param seed Integer random seed.
#' 
#' @return A list containing the simulation grid, residuals, regimes, and
#'   scenario matrices.
#' @export
scenarios_ER <- function(model_Et, model_Rt, rho, t_now, t_hor, nsim = 1, seed = 1){
  # Preprocess 
  preproc <- scenarios_ER_proproc(model_Et, model_Rt, rho, t_now, t_hor)
  # Simulate residuals 
  preproc <- scenarios_ER_residuals(nsim = nsim, preproc, seed = seed)
  # Scenarios 
  preproc <- scenarios_ER_filter(preproc)
  return(preproc)
}

