#' Electricity Price Model
#'
#' R6 wrapper around a fitted Log-OU electricity price model. The class stores
#' the original dataset, training split, flexible probability weights, fitted
#' log-price residuals, and methods for conditional moments under the selected
#' measure.
#'
#' @field date_train Dates used for training.
#'
#' @export
electricityModel <- R6::R6Class("electricityModel",
                                public = list(
                                  date_train = NA,
#' @description
#' Initialize an electricity model from a `PUN` price time series.
#'
#' @param data Data frame with at least `date` and `PUN`.
#' @param max_date Last training date.
#' @param tau_hl Optional half-life used for exponential-decay flexible
#'   probabilities.
                                  initialize = function(data, max_date = max(data$date), tau_hl){
                                    # Indicator for train data
                                    data <- dplyr::mutate(data, isTrain = ifelse(date <= max_date, TRUE, FALSE))
                                    # Add time references 
                                    data$Year <- lubridate::year(data$date)
                                    data$Month <- lubridate::month(data$date)
                                    data$Day <- lubridate::day(data$date)
                                    # Reorder the variables 
                                    data <- dplyr::select(data, date, Year, Month, Day, isTrain, dplyr::everything())
                                    # Compute log-value
                                    data$Xt <- log(data$PUN)
                                    # Train data
                                    data_train <- dplyr::filter(data, date <= max_date)
                                    # Store train date
                                    self$date_train <- data_train$date
                                    # Initialize and fit the model 
                                    ou_model <- logOUModel$new(data_train$PUN, tau_hl)
                                    # Flexible probability weights 
                                    weights <- ou_model$data$weights
                                    # Set weights of test set equal to zero 
                                    data$weights <- c(weights, rep(0, nrow(data) - length(weights)))
                                    # Compute fitted log-values 
                                    data$Xt_hat <- predict(ou_model$model, newdata = data.frame(x = data$Xt))
                                    # Compute fitted log-values 
                                    data$St_hat <- exp(data$Xt_hat)
                                    # Compute residuals on complete data
                                    data$z_E <- data$Xt - data$Xt_hat
                                    # Compute standardized residuals
                                    data$z_tilde_E <- data$z_E / ou_model$model$params.AR[["sigma_u"]]
                                    # *****************************************************
                                    # Store fitted model
                                    private$..model <- ou_model
                                    # Store complete data
                                    private$..data <- data
                                  },
                                  #' @description
                                  #' Change the probability measure from P to Q and viceversa. 
                                  #' Set the market risk premium for Q-measure.
                                  #' @param measure Character, probability measure "P" or "Q".
                                  #' @param lambda Numeric, market risk premium. 
                                  change_measure = function(measure, lambda){
                                    # If lambda is not missing is updated 
                                    if (!missing(lambda)) {
                                      private$..model$lambda <- lambda
                                    }
                                    private$..model$change_measure(measure)
                                  },
                                  #' @description
                                  #' Conditional expectation Normal variable 
                                  #' @param t_now Date, conditioning date.
                                  #' @param t_hor Date, horizon date.
                                  #' @param E0 Numeric, value of electricity prices at time `t_now`.
                                  M_X = function(t_now, t_hor, E0){
                                    t_now = as.Date(t_now)
                                    t_hor = as.Date(t_hor)
                                    # Compute time to maturity in days
                                    tau = as.numeric(difftime(t_hor, t_now, units = "days"))
                                    self$model$M_X(tau, E0)
                                  },
                                  #' @description
                                  #' Conditional variance Normal variable 
                                  #' @param t_now Date, conditioning date.
                                  #' @param t_hor Date, horizon date.
                                  V_X = function(t_now, t_hor){
                                    t_now = as.Date(t_now)
                                    t_hor = as.Date(t_hor)
                                    # Compute time to maturity in days
                                    tau = as.numeric(difftime(t_hor, t_now, units = "days"))
                                    self$model$V_X(tau)
                                  },
                                  #' @description
                                  #' Expected value of the log-normal variable
                                  #' @param t_now Date, conditioning date.
                                  #' @param t_hor Date, horizon date.
                                  #' @param E0 Numeric, value of electricity prices at time `t_now`.
                                  F_E = function(t_now, t_hor, E0){
                                    exp(self$M_X(t_now, t_hor, E0) + 0.5 * self$V_X(t_now, t_hor))
                                  },
                                  #' @description
                                  #' Variance of the log-normal variable
                                  #' @param t_now Date, conditioning date.
                                  #' @param t_hor Date, horizon date.
                                  #' @param E0 Numeric, value of electricity prices at time `t_now`.
                                  S_E = function(t_now, t_hor, E0){
                                    self$F_E(t_now, t_hor, E0)^2 * (exp(self$V_X(t_now, t_hor)) - 1)
                                  },
                                  #' @description
                                  #' `print` method for `electricityModel` class. 
                                  print = function(){
                                    cat("----------- Electricity Model ----------- \n")
                                    cat("- From: ", as.character(min(self$date_train)), " \n")
                                    cat("- To: ", as.character(max(self$date_train)), " \n")
                                    cat("- Nobs: ", length(self$date_train), " \n")
                                    cat("- Measure: ", self$model$measure, " \n")
                                    cat("- MRP: ", self$model$lambda, " \n")
                                  }
                                ),
                                private = list(
                                  ..model = NA,
                                  ..data = NA
                                ),
                                active = list(
                                  #' @field params Named numeric vector of Log-OU parameters.
                                  params = function(){
                                    self$model$params
                                  },
                                  #' @field model Underlying `logOUModel` object.
                                  model = function(){
                                    private$..model
                                  },
                                  #' @field data Full electricity dataset with fitted values and residuals.
                                  data = function(){
                                    private$..data
                                  }
                                ))

#' Evaluate an Electricity Change-of-Measure Calibration Moment
#'
#' Compute the weighted average realized-to-model future price ratio for a
#' candidate market risk premium.
#'
#' @param lambda Numeric scalar. Market risk premium.
#' @param tau Integer scalar. Time to maturity in days.
#' @param r Numeric scalar. Daily risk-free rate used for discounting.
#' @param model_Et An `electricityModel` object.
#' @param tau_hl Numeric scalar. Half-life for exponential-decay weights.
#'
#' @return Numeric scalar used as a calibration moment.
#' @export
dQdP_electricity <- function(lambda = 0, tau = 10, r = 0, model_Et, tau_hl = 65){
  model <- model_Et$clone(deep = TRUE)
  model$change_measure("Q", lambda)
  data_lag <- dplyr::filter(model$data, date %in% model$date_train)
  data_lag$probs <- exp_decay_fp(nrow(data_lag), tau_hl, nrow(data_lag)) # model$model$data$probs
  data_lag$L_PUN <- dplyr::lag(data_lag$PUN, tau)
  data_lag$t_now <- dplyr::lag(data_lag$date, tau)
  data_lag <- na.omit(data_lag)
  # Future price as Q-expectation
  FE_tT <- model$F_E(data_lag$t_now, data_lag$date, data_lag$L_PUN) * exp(-r * tau)
  # Loss function weighted
  sum((data_lag$PUN / FE_tT) * data_lag$probs)
}

#' Calibrate Electricity Market Risk Premium
#'
#' Choose the market risk premium so that the average realized-to-model future
#' price ratio is close to one across maturities.
#'
#' @param tau_max Integer scalar. Maximum maturity in days.
#' @param r Numeric scalar. Daily risk-free rate used for discounting.
#' @param model_Et An `electricityModel` object.
#' @param tau_hl Numeric scalar. Half-life for exponential-decay weights.
#' @param quiet Logical. If `TRUE`, suppress optimizer progress messages.
#'
#' @return A cloned `electricityModel` with the calibrated premium stored.
#' @export
calibrate_dQdP_electricity <- function(tau_max = 10, r = 0, model_Et, tau_hl = 65, quiet = FALSE){
  # Objective used to calibrate the market risk premium.
  loss_function <- function(lambda){
    avg <- c()
    for(tau in 2:tau_max){
      avg[tau-1] <- dQdP_electricity(lambda, tau, r, model_Et, tau_hl)
    }
    loss <- sum((avg - 1)^2, na.rm = TRUE)
    if(!quiet) message("Loss: ", loss, " Lambda: ", lambda, "\r", appendLF = FALSE)
    return(loss)
  }
  opt <- optim(par = c(0), loss_function, method = "Brent", lower = -1, upper = 1)
  message("Loss: ", opt$value, " Lambda: ", opt$par)
  model <- model_Et$clone(deep = TRUE)
  model$change_measure("P", opt$par)
  return(model)
}

