#' Log OU Model
#' @export
logOUModel <- R6::R6Class("logOUModel",
                          public = list(
                            #' @field kappa Numeric, mean reversion parameter.
                            kappa = NA,
                            #' @field mu Numeric, long-term mean log-prices.
                            mu = NA,
                            #' @field sigma Numeric, volatility of log-prices.
                            sigma = NA,
                            #' @field lambda Numeric, market risk premium.
                            lambda = 0,
                            #' @field lambda Numeric, half-life parameter.
                            tau_hl = 0,
                            #' @description
                            #' Initialize a Log-OU model
                            #' @param St Numeric, time series of prices.
                            #' @param tau_hl Integer, half-life parameter for exponential decay probability. If missing will be used constant probabilities.
                            initialize = function(St, tau_hl){
                              # Number of observations
                              N <- length(St)
                              # Log-prices
                              Xt <- log(St)
                              # Flexible probabilities
                              if (!missing(tau_hl)) {
                                # Half-life parameter
                                self$tau_hl <- tau_hl
                                #tau_hl <- N/tau_hl
                                # Exponential decay probabilities
                                probs <- exp_decay_fp(N, tau_hl)
                              } else {
                                # No weights (or uniform weights)
                                probs <- rep(1/N, N)
                              }
                              # Dataset for fit
                              data <- dplyr::tibble(St = St, Xt = Xt, probs = probs)
                              # Euler discretization for the SDE
                              AR_model <- fit_OU_model(Xt, probs)
                              # Add fitted log-prices
                              data$Xt_hat <- AR_model$fitted.values
                              # Add residuals on complete data
                              data$z_E <- data$Xt - data$Xt_hat
                              # Add standardized residuals
                              data$z_tilde_E <- data$z_E / AR_model$params.AR[["sigma_u"]]
                              # Store flexible weights
                              data$weights <- probs
                              # Add fitted exp-prices
                              data$St_hat <- exp(data$Xt_hat)
                              # OU parameters (dt = 1)
                              self$mu <- AR_model$params.OU[["mu"]]
                              self$kappa <- AR_model$params.OU[["kappa"]]
                              self$sigma <- AR_model$params.OU[["sigma"]]
                              # Add fitted model
                              private$..model <- AR_model
                              # Add complete data
                              private$..data <- data
                            },
                            #' @description
                            #' Change the probability measure
                            #' @param measure Character, reference measure, can be `P` or `Q`.
                            change_measure = function(measure){
                              measure <- match.arg(measure, choices = c("P", "Q"))
                              private$..measure <- measure
                            },
                            #' @description
                            #' Convert from log-prices into prices.
                            #' @param x Numeric, log-prices.
                            X_to_S = function(x){
                              exp(x)
                            },
                            #' @description
                            #' Convert from prices into log-prices
                            #' @param s Numeric, prices.
                            S_to_X = function(s){
                              log(s)
                            },
                            #' @description
                            #' Compute the expected value of the log variable.
                            #' @param tau Numeric. time to maturity in days.
                            #' @param S0 Numeric, price at time 0.
                            M_X = function(tau, S0){
                              # Parameters
                              kappa <- self$kappa
                              lambda <- ifelse(self$measure == "Q", self$lambda, 0)
                              sigma <- self$sigma
                              mu <- self$mu
                              # Exponential factor
                              exp_tau <- exp(-kappa * tau)
                              # Conditional mean of X
                              self$S_to_X(S0) * exp_tau + (1 / kappa) * (lambda * sigma + kappa * mu - sigma^2 / 2) * (1 - exp_tau)
                            },
                            #' @description
                            #' Compute the variance the log variable.
                            #' @param tau Numeric. time to maturity in days.
                            V_X = function(tau){
                              # Parameters
                              kappa <- self$kappa
                              sigma <- self$sigma
                              # Exponential factor
                              exp_2tau <- exp(-2 * kappa * tau)
                              # Conditional std. deviation of X
                              sigma^2 * (1 - exp_2tau) / (2 * kappa)
                            },
                            #' @description
                            #' Compute the expected value of the exponential variable.
                            #' @param tau Numeric. time to maturity in days.
                            #' @param S0 Numeric, price at time 0.
                            M_S = function(tau, S0){
                              exp(self$M_X(tau, S0) + 0.5 * self$V_X(tau))
                            },
                            #' @description
                            #' Compute the std deviation of the exponential variable.
                            #' @param tau Numeric. time to maturity in days.
                            #' @param S0 Numeric, price at time 0.
                            V_S = function(tau, S0){
                              self$M_S(tau, S0)^2 * (exp(self$V_X(tau)) - 1)
                            },
                            #' @description
                            #' Density function log-variable
                            #' @param x Numeric, value for computing the density.
                            #' @param tau Numeric. time to maturity in days.
                            #' @param S0 Numeric, price at time 0.
                            dX = function(x, tau, S0){
                              dnorm(x, mean = self$M_X(tau, S0), sd = sqrt(self$S_X(tau)))
                            },
                            #' @description
                            #' Distribution function log-variable
                            #' @param q Numeric, quantile for computing the distribution.
                            #' @param tau Numeric. time to maturity in days.
                            #' @param S0 Numeric, price at time 0.
                            pX = function(q, tau, S0){
                              pnorm(q, mean = self$M_X(tau, S0), sd = sqrt(self$S_X(tau)))
                            },
                            #' @description
                            #' Quantile function log-variable
                            #' @param p Numeric, probability in [0,1] for the quantile
                            #' @param tau Numeric. time to maturity in days.
                            #' @param S0 Numeric, price at time 0.
                            qX = function(p, tau, S0) {
                              qnorm(p, mean = self$M_X(tau, S0), sd = sqrt(self$S_X(tau)))
                            },
                            #' @description
                            #' Density function target variable
                            #' @param x Numeric, value for computing the density.
                            #' @param tau Numeric. time to maturity in days.
                            #' @param S0 Numeric, price at time 0.
                            dS = function(x, tau, S0){
                              (1 / x) * self$dX(self$S_to_X(x), tau, S0)
                            },
                            #' @description
                            #' Distribution function target variable
                            #' @param q Numeric, quantile for computing the distribution.
                            #' @param tau Numeric. time to maturity in days.
                            #' @param S0 Numeric, price at time 0.
                            pS = function(q, tau, S0){
                              self$pX(self$S_to_X(q), tau, S0)
                            },
                            #' @description
                            #' Quantile function target variable
                            #' @param p Numeric, probability in [0,1] for the quantile
                            #' @param tau Numeric. time to maturity in days.
                            #' @param S0 Numeric, price at time 0.
                            qS = function(p, tau, S0){
                              self$X_to_S(self$qX(p, tau, S0))
                            },
                            #' @description
                            #' Compute the grades given S0
                            #' @param x Numeric, value for computing the grades.
                            #' @param tau Numeric. time to maturity in days.
                            #' @param S0 Numeric, price at time 0.
                            grade = function(x, tau, S0){
                              z <- (self$S_to_X(x) - self$M_X(tau, S0)) / sqrt(self$S_X(tau))
                              pnorm(z)
                            },
                            #' @description
                            #' Compute the drift of the exponential variable.
                            #' @param S0 Price at time 0.
                            drift_S = function(S0){
                              drift_Q <- ifelse(self$measure == "Q", self$lambda * self$sigma, 0)
                              S0 * self$kappa * (self$mu - log(S0) + drift_Q)
                            },
                            #' @description
                            #' Compute the diffusion of the exponential variable.
                            #' @param S0 Price at time 0.
                            diffusion_S = function(S0){
                              S0 * self$sigma
                            },
                            #' @description
                            #' Compute the drift of the log variable.
                            #' @param X0 Log-price at time 0.
                            drift_X = function(X0){
                              drift_Q <- ifelse(self$measure == "Q", self$lambda * self$sigma, 0)
                              self$kappa * (self$mu - self$sigma^2 / (2*self$kappa) - X0) - drift_Q
                            },
                            #' @description
                            #' Compute the diffusion of the log variable.
                            diffusion_X = function(){
                              self$sigma
                            },
                            #' @description
                            #' Simulate X and S
                            #' @param S0 Price at time 0.
                            #' @param tau Time to maturity in days.
                            #' @param dt time step.
                            #' @param seed random seed
                            #' @param dW Brownian motions.
                            simulate = function(S0, tau, dt = 1, seed = 1, dW){
                              # Initial value
                              X <- c(log(S0))
                              # Time index
                              t_index = seq(0, tau, by = dt)
                              # Simulated Brownian
                              set.seed(seed)
                              if (missing(dW)) {
                                dW <- rnorm(length(t_index), 0, sqrt(dt))
                              }
                              for(i in 2:length(t_index)){
                                dX <- self$drift_X(X[i - 1]) * dt + self$diffusion_X() * dW[i]
                                X[i] <- X[i - 1] + dX
                              }
                              S <- self$X_to_S(X)
                              dplyr::tibble(seed = seed, t = t_index, X = X, S = S, dW = dW)
                            }
                          ),
                          private = list(
                            ..model = NA,
                            ..data = NA,
                            ..measure = "P"
                          ),
                          active = list(
                            #' @field params Numeric vector. OU parameters.
                            params = function(){
                              c(kappa = self$kappa, mu = self$mu, sigma = self$sigma)
                            },
                            #' @field model Fitted `lm` object containing the AR model.
                            model = function(){
                              private$..model
                            },
                            #' @field data Data used to fit the AR model
                            data = function(){
                              private$..data
                            },
                            #' @field measure Probability measure
                            measure = function(){
                              private$..measure
                            }
                          ))

#' Convert AR Parameters to Log-OU Parameters
#' 
#' @param phi Numeric vector. AR intercept and autoregressive slope.
#' @param sigma_u Numeric, standard deviation of the AR residuals.
#'
#' @return Named numeric vector with `mu`, `kappa`, and `sigma`.
#' @keywords internal
from_AR_to_OU <- function(phi, sigma_u){
  # Extract AR parameters
  phi_0 <- phi[1]
  phi_1 <- phi[2]
  # OU parameters (dt = 1)
  kappa <- -log(phi_1)
  sigma <- sigma_u * sqrt(- 2 * log(phi_1) / (1 - phi_1^2))
  mu <- phi_0 / (1 - exp(-kappa)) + (sigma^2) / (2 * kappa)
  c(mu = mu[[1]], kappa = kappa[[1]], sigma = sigma[[1]])
}

#' Jacobian of the AR-to-OU Reparameterization
#'
#' Compute the analytical Jacobian used to propagate AR parameter uncertainty
#' to continuous-time Log-OU parameters.
#'
#' @param phi Numeric vector. AR intercept and autoregressive slope.
#' @param sigma_u Numeric scalar. Standard deviation of AR residuals.
#'
#' @return A `3 x 3` matrix with rows `mu`, `kappa`, `sigma` and columns
#'   `phi_0`, `phi_1`, `sigma_u`.
#' @keywords internal
Jacobian_AR_to_OU <- function(phi, sigma_u){
  
  # To check
  # numDeriv::jacobian(function(x) from_AR_to_OU(x[1:2], x[3]), x = c(phi, sigma_u) )
  
  # Extract AR parameters
  phi_0 <- phi[1]
  phi_1 <- phi[2]
  # OU parameters 
  params.OU <- from_AR_to_OU(phi, sigma_u)
  kappa <- params.OU["kappa"]
  sigma <- params.OU["sigma"]
  
  # Initialization 
  J <- matrix(0, nrow = 3, ncol = 3)
  colnames(J) <- c("phi_0", "phi_1", "sigma_u")
  rownames(J) <- c("mu", "kappa", "sigma")
  
  # d mu_X / d phi_0 
  J[1,1] <- 1 / (1-phi_1)
  # d mu_X / d phi_1
  J[1,2] <- phi_0 / (1-phi_1)^2 + (2 * phi_1 * sigma_u^2) / (1 - phi_1^2)^2 
  # d mu_X / d sigma_u
  J[1,3] <- 2 * sigma_u / (1 - phi_1^2) 
  
  # d kappa / d phi_0 
  J[2,1] <- 0
  # d kappa / d phi_1
  J[2,2] <- - 1 / phi_1
  # d kappa / d sigma_u
  J[2,3] <- 0 
  
  # d sigma_X / d phi_0 
  J[3,1] <- 0
  # d sigma_X / d phi_1
  J[3,2] <- (sigma_u^2 * (2 * (1 - exp(-2*kappa)) - 4 * kappa * exp(-2*kappa))) / (2 * sigma * phi_1 * (1 - exp(-2*kappa))^2)
  # d sigma_X / d sigma_u
  J[3,3] <- sigma / sigma_u
  return(J)
}

#' Fit a Discrete AR Model and Map It to Log-OU Parameters
#'
#' Fit the AR representation of log-prices with optional weights and append
#' continuous-time Log-OU parameters and standard errors to the fitted `lm`
#' object.
#'
#' @param x Numeric vector. Log-price observations.
#' @param weights Optional numeric vector of observation weights.
#'
#' @return An `lm` object augmented with AR/OU parameters, standard errors,
#'   fitted values, and residuals.
#' @keywords internal
fit_OU_model <- function(x, weights){
  
  # Number of observations
  n <- length(x)
  # Custom weights 
  if (missing(weights)) {
    weights = rep(1, n)
  }
  
  AR_model <- lm(x ~ I(dplyr::lag(x, 1)), weights = weights)
  # Fitted values
  x_hat <- c(x[1], AR_model$fitted.values)
  # Residuals 
  eps <- x - x_hat
  # Estimated parameters
  phi <- AR_model$coefficients
  sigma_u <- sd(eps[-1])
  # Jacobian 
  J_AR_to_OU <- Jacobian_AR_to_OU(phi, sigma_u)
  # Matrix for the AR std. errors 
  vcov.AR <- matrix(0, 3, 3)
  colnames(vcov.AR) <- rownames(vcov.AR) <- colnames(J_AR_to_OU)
  diag(vcov.AR) <- c(broom::tidy(AR_model)$std.error^2, sigma_u^2/(2*length(eps)))
  # Matrix for the OU std. errors 
  vcov.OU <- J_AR_to_OU %*% vcov.AR %*% t(J_AR_to_OU)
  # Standard errors 
  std.errors.AR <- sqrt(diag(vcov.AR))
  std.errors.OU <- sqrt(diag(vcov.OU))
  
  AR_model$params.AR = c(phi_0 = phi[[1]], phi_1 = phi[[2]], sigma_u = sigma_u[[1]])
  AR_model$std.errors.AR = std.errors.AR
  AR_model$params.OU = from_AR_to_OU(phi, sigma_u)
  AR_model$std.errors.OU = std.errors.OU
  
  AR_model$fitted.values <- x_hat
  AR_model$residuals <- eps
  
  return(AR_model)
}





















