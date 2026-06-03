#' Simulate IID Mixture Regime Paths for a Radiation Model
#'
#' Simulate Bernoulli mixture states and Brownian increments from the IID
#' mixture model.
#'
#' @param GM_model Fitted Gaussian mixture object.
#' @param month_idx Integer vector of monthly indices for each simulation step.
#' @param dt Numeric scalar. Time step in days.
#' @param nsteps Integer. Number of time steps.
#' @param nsim Integer. Number of Monte Carlo paths.
#'
#' @return List with matrices `mu_B`, `sigma_B`, `B_t`, and `dMt`.
#' @keywords internal
gm_monthly_scenarios <- function(GM_model, month_idx, dt = 1, nsteps = 10, nsim = 1){
  # Initialization 
  B_t <- matrix(0, nrow = nsteps, ncol = nsim)
  dMt <- matrix(0, nrow = nsteps, ncol = nsim)
  mu_B    <- matrix(0, nrow = nsteps, ncol = nsim)
  sigma_B <- matrix(0, nrow = nsteps, ncol = nsim)
  # Predicted probabilities 
  probs <- GM_model$prob$predict(month_idx)
  mu10   <- cbind(GM_model$mu1$predict(month_idx), GM_model$mu2$predict(month_idx))
  sd10   <- cbind(GM_model$sd1$predict(month_idx), GM_model$sd2$predict(month_idx))
  for (sim in 1:nsim) {
    message("Simulations IID Bernoulli regimes: ", sim, "/", nsim, "\r", appendLF = FALSE)
    # Simulated Bernoulli
    B_t[, sim] <- purrr::map_dbl(probs, ~rbinom(1, 1, .x))
    # Realized moments 
    mu_B[, sim] <- purrr::map_dbl(1:nsteps, ~ifelse(B_t[.x, sim] == 1, mu10[.x, 1], mu10[.x, 2]))
    sigma_B[, sim] <- purrr::map_dbl(1:nsteps, ~ifelse(B_t[.x, sim] == 1, sd10[.x, 1], sd10[.x, 2]))
    # Brownian simulations
    dW_1 <- rnorm(nsteps, 0, sqrt(dt))
    dW_0 <- rnorm(nsteps, 0, sqrt(dt))
    # Final simulation 
    dMt[, sim] <- dW_1 * B_t[, sim] + dW_0 * (1 - B_t[, sim])
  }
  list(
    mu_B = mu_B, 
    sigma_B = sigma_B,
    B_t = B_t, 
    dMt = dMt
  )
}

#' Simulate Continuous-Time Radiation Scenarios
#'
#' Simulate transformed and physical solar radiation paths under either the HMM
#' CTMC regime model or the IID mixture regime model.
#'
#' @param model_Rt A `radiationModel` or `radiationModelHMM` object.
#' @param t_now Date or character scalar. Conditioning date.
#' @param t_hor Numeric scalar. Horizon in days.
#' @param nsim Integer. Number of Monte Carlo paths.
#' @param dt Numeric scalar. Time step in days.
#' @param seed Integer random seed.
#'
#' @return Tibble containing simulated paths and deterministic diagnostics.
#' @export
scenarios_radiationModel_CT  <- function(model_Rt, t_now, t_hor = 365, nsim = 3, dt = 0.5, seed = 53){
  set.seed(seed)
  # Create the sequence of dates 
  date <- as.Date(t_now)
  # Daily sequence
  date_seq <- seq.Date(date, date + t_hor, 1)
  # Repeated sequence 
  date_dt <- rep(date_seq, 1/dt)
  # Time sequence 
  time <- seq(0, t_hor-dt, length.out = length(date_dt)) + number_of_day(date[1])
  # Dataset for simulated data 
  data_sim <- dplyr::tibble(date_dt = date_dt[order(date_dt)], time = time)
  data_sim$step <- rep(1:(1/dt), t_hor+1)
  # Number of the day of the year 
  data_sim$n <- number_of_day(data_sim$date_dt)
  # Daily dates 
  data_sim$date <- data_sim$date_dt
  #data_sim$date[data_sim$step != 1/dt] <- NA
  data_sim$date[data_sim$step != 1] <- NA
  # Add realized data 
  data_sim <- dplyr::left_join(data_sim, dplyr::select(model_Rt$model$data, date, GHI, Yt), by = "date")
  data_sim <- data_sim[1:max(which(!is.na(data_sim$date))),]
  
  # *************************************************************************************
  # Pre-compute deterministic variables 
  # Seasonal clear-sky 
  data_sim$Ct <- model_Rt$Ct(data_sim$time)
  # Differential Seasonal clear-sky 
  data_sim$dCt_dt <- model_Rt$dCt(data_sim$time)
  # Seasonal mean Yt 
  data_sim$Yt_bar <- model_Rt$Yt_bar(data_sim$time)
  # Differential Seasonal mean Yt 
  data_sim$dYt_dt <- model_Rt$dYt_bar(data_sim$time)
  # Seasonal mean Rt
  data_sim$GHI_bar <- model_Rt$model$spec$transform$iRY(data_sim$Yt_bar, data_sim$Ct)
  # Seasonal std. deviation Yt 
  data_sim$sigma_bar <- model_Rt$sigma_bar(data_sim$time)
  # Extract model's parameters
  beta <- model_Rt$model$spec$transform$beta
  # Mean reversion parameter for a time step dt
  theta <- model_Rt$model$spec$mean.model$phi^dt
  # Number of steps
  nsteps <- length(data_sim$time)
  # ************************************************************ 
  #                         Mixture Brownian
  # ************************************************************ 
  if (any(class(model_Rt) %in% "radiationModel_CTMC")) {
     # HMM model 
     CTMC <- model_Rt$CTMC
     # Initial probability
     p0 <- CTMC$alpha[CTMC$data$date == t_now, ]
     # Monthly indexes 
     month_idx <- lubridate::month(data_sim$date_dt)
     # Scenarios 
     Mt <- ctmc_monthly_scenarios(CTMC, p0, month_idx, dt, nsteps, nsim)
  } else {
    month_idx <- lubridate::month(data_sim$date_dt)
    # Generate scenarios for IID Bernoullis 
    GM_model <- model_Rt$model$spec$mixture.model
    # Scenarios 
    Mt <- gm_monthly_scenarios(GM_model, month_idx, dt, nsteps, nsim)
  }
  B_t <- Mt$B_t
  dMt <- Mt$dMt
  mu_B <- Mt$mu_B
  sigma_B <- Mt$sigma_B
  # ************************************************************ 
  #                          Simulations
  # ************************************************************ 
  # Simulations 
  scenarios <- list()
  # Simulation under P
  sim <- 1
  # Initial point 
  R0 <- filter(model_Rt$model$data, date == date_seq[1])$GHI
  for (sim in 1:nsim) {
    message("Filter sim:", sim, "/", nsim, "\r", appendLF = FALSE)
    # Starting point 
    Rt <- Rt_Yt <- R0
    Yt <-  model_Rt$model$spec$transform$RY(Rt[1], data_sim$Ct[1]) 
    i <- 1
    for (i in 1:nsteps) {
      # Simulate transformed variable 
      # Drift 
      mu_Y <- data_sim$dYt_dt[i] + theta * (data_sim$Yt_bar[i] - Yt[i]) + data_sim$sigma_bar[i] * mu_B[i,sim]
      # Diffusion 
      sigma_Y <- data_sim$sigma_bar[i] * sigma_B[i, sim]
      # Increment 
      dYt <- mu_Y * dt + sigma_Y * dMt[i, sim] 
      # ************************************
      # Simulate solar radiation 
      Yt_i <- model_Rt$model$spec$transform$RY(Rt_Yt[i], data_sim$Ct[i]) 
      # Drift 
      mu_Y_i <- data_sim$dYt_dt[i] + theta * (data_sim$Yt_bar[i] - Yt_i) + data_sim$sigma_bar[i] * mu_B[i,sim]
      # Diffusion 
      sigma_Y_i <- data_sim$sigma_bar[i] * sigma_B[i, sim]
      # Drift 
      mu_R <- Rt_Yt[i]/data_sim$Ct[i] * data_sim$dCt_dt[i]  + data_sim$Ct[i] * beta * exp(Yt_i - exp(Yt_i)) * (mu_Y_i + 0.5 * (1 - exp(Yt_i)) * sigma_Y_i^2)
      # Diffusion 
      sigma_R <- data_sim$Ct[i] * beta * exp(Yt_i - exp(Yt_i)) * sigma_Y_i
      # Increment 
      dRt <- mu_R * dt + sigma_R * dMt[i, sim] 
      # ************************************
      # Next step value of Yt 
      Yt[i+1] <- Yt[i] + dYt
      # Next step value of Rt from exact solution  
      Rt[i+1] <- model_Rt$model$spec$transform$iRY(Yt[i+1], data_sim$Ct[i+1]) 
      # Next step value of Rt from implied dynamic 
      Rt_Yt[i+1] <- Rt_Yt[i] + dRt
    }
    # Store the simulations 
    scenarios[[sim]] <- data_sim
    scenarios[[sim]]$i <- 1:nrow(data_sim)
    scenarios[[sim]]$sim <- sim
    scenarios[[sim]]$Yt <- Yt[-nsteps]
    scenarios[[sim]]$Rt <- Rt[-nsteps]
    scenarios[[sim]]$Rt_Yt <- Rt_Yt[-nsteps]
    # SoRad Payoff 
    scenarios[[sim]]$gamma <- ((data_sim$GHI_bar -  Rt[-nsteps]) * ifelse(data_sim$GHI_bar >  Rt[-nsteps], 1, 0))
  }
  dplyr::bind_rows(scenarios) 
}
