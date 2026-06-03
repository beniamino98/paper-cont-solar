#' Preprocess Radiation Data for Direct Monthly CTMC Estimation
#'
#' Wrapper around `dtmc_monthly_preprocess()` used by the direct CTMC estimator.
#' The preprocessing of residuals, seasonal components, and emissions is kept
#' identical to the existing DTMC estimator; only the transition M-step changes.
#'
#' @param Y Numeric vector. Transformed radiation observations.
#' @param dates Date vector aligned with `Y`.
#' @param weights Optional numeric/logical training weights.
#' @param params List containing initial DTMC/CTMC parameters.
#' @param eps Numeric tolerance used to stabilize probabilities.
#' @param update_emissions Logical. If `TRUE`, update starting emission
#'   parameters from residuals.
#'
#' @return A list with preprocessed `data` and completed `params`.
#' @keywords internal
dtmc_monthly_preprocess_direct <- function(Y, dates, weights, params, eps = 1e-8, update_emissions = FALSE) {
  dtmc_monthly_preprocess(
    Y = Y,
    dates = dates,
    weights = weights,
    params = params,
    eps = eps,
    update_emissions = update_emissions
  )
}

#' Run One Direct Monthly CTMC EM Update
#'
#' Perform one forward-backward E-step as in `dtmc_monthly_EM()`, but update
#' monthly transition dynamics by estimating the position-indexed CTMC rates
#' `q12_m` and `q21_m` directly. With state order `(1, 0)`, these correspond
#' to transitions `1 -> 0` and `0 -> 1`, respectively. The one-day transition
#' matrix used in the HMM recursion is then obtained from the corresponding
#' two-state CTMC transition formula.
#'
#' @param data Preprocessed data returned by `dtmc_monthly_preprocess_direct()`.
#' @param params List of CTMC/DTMC parameters. Must contain `pi`, `Pm`, `mu`,
#'   `sig`, `a`, `theta`, and `b`; may also contain `Qm`.
#' @param update_emissions Logical. If `TRUE`, update monthly emission means
#'   and standard deviations.
#' @param update_transitions Logical. If `TRUE`, update monthly CTMC rates.
#' @param eps Numeric tolerance used to stabilize probabilities.
#'
#' @return A list containing updated parameters, filtered/smoothed
#'   probabilities, transition posteriors, log-likelihood, and optimizer
#'   diagnostics for the direct CTMC transition M-step.
#' @keywords internal
dtmc_monthly_EM_direct <- function(data, params,
                                   update_emissions = TRUE,
                                   update_transitions = TRUE,
                                   eps = 1e-8) {
  Tn <- nrow(data)
  month_idx <- data$Month
  dates     <- data$date
  weights   <- data$weights
  r_t       <- data$r_t
  eps_small <- eps

  logSumExp <- function(x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
  }

  ctmc_P_2state <- function(q12, q21, delta = 1) {
    lambda <- q12 + q21
    if (!is.finite(lambda) || lambda <= 0) {
      return(diag(2))
    }
    move <- 1 - exp(-lambda * delta)
    p12 <- q12 / lambda * move
    p21 <- q21 / lambda * move
    matrix(c(1 - p12, p12,
             p21, 1 - p21), 2, 2, byrow = TRUE)
  }

  pi <- as.numeric(params$pi)
  if (length(pi) != 2) stop("params$pi must have length 2")
  pi <- pmax(pi, eps_small)
  pi <- pi / sum(pi)

  Pm <- params$Pm
  if (length(Pm) != 12) stop("params$Pm must be a list of length 12")
  for (m in 1:12) {
    if (!is.matrix(Pm[[m]]) || any(dim(Pm[[m]]) != c(2, 2))) {
      stop("Each Pm[[m]] must be a 2x2 matrix")
    }
    Pm[[m]] <- pmax(Pm[[m]], eps_small)
    Pm[[m]] <- Pm[[m]] / rowSums(Pm[[m]])
  }

  Qm <- params$Qm
  if (is.null(Qm)) {
    Qm <- transition_list_to_generator_2state(Pm, delta = 1)
  }

  mu <- params$mu
  mu_1 <- purrr::map_dbl(1:Tn, ~mu[[month_idx[.x]]][1])
  mu_2 <- purrr::map_dbl(1:Tn, ~mu[[month_idx[.x]]][2])
  sig <- params$sig
  sd_1 <- purrr::map_dbl(1:Tn, ~sig[[month_idx[.x]]][1])
  sd_2 <- purrr::map_dbl(1:Tn, ~sig[[month_idx[.x]]][2])

  logf <- cbind(
    stats::dnorm((r_t - mu_1) / sd_1, log = TRUE) - log(sd_1),
    stats::dnorm((r_t - mu_2) / sd_2, log = TRUE) - log(sd_2)
  )

  logalpha <- matrix(NA_real_, Tn, 2)
  logbeta  <- matrix(NA_real_, Tn, 2)

  logalpha[1, ] <- log(pi + 1e-12) + logf[1, ]
  for (t in 2:Tn) {
    m_prev <- month_idx[t - 1]
    for (j in 1:2) {
      trans_log <- logalpha[t - 1, ] + log(pmax(Pm[[m_prev]][, j], eps_small))
      logalpha[t, j] <- logf[t, j] + logSumExp(trans_log)
    }
  }

  logbeta[Tn, ] <- 0
  for (t in (Tn - 1):1) {
    m_curr <- month_idx[t]
    for (i in 1:2) {
      trans_log <- log(pmax(Pm[[m_curr]][i, ], eps_small)) + logf[t + 1, ] + logbeta[t + 1, ]
      logbeta[t, i] <- logSumExp(trans_log)
    }
  }

  loglik <- logSumExp(logalpha[Tn, ])

  gamma <- matrix(NA_real_, Tn, 2)
  for (t in 1:Tn) {
    z <- logalpha[t, ] + logbeta[t, ]
    gamma[t, ] <- exp(z - logSumExp(z))
  }

  # Joint probabilities in row-stochastic code order: xi_t(previous, next).
  xi <- array(NA_real_, dim = c(Tn - 1, 2, 2))
  for (t in 1:(Tn - 1)) {
    m_curr <- month_idx[t]
    Z <- matrix(NA_real_, 2, 2)
    for (i in 1:2) {
      for (j in 1:2) {
        Z[i, j] <- logalpha[t, i] +
          log(pmax(Pm[[m_curr]][i, j], eps_small)) +
          logf[t + 1, j] +
          logbeta[t + 1, j]
      }
    }
    z_norm <- logSumExp(as.vector(Z))
    xi[t, , ] <- exp(Z - z_norm)
  }

  pi_new <- gamma[1, ]
  pi_new <- pmax(pi_new, eps_small)
  pi_new <- pi_new / sum(pi_new)

  Pm_new <- Pm
  Qm_new <- Qm
  transition_optim <- lapply(seq_along(Pm), function(m) {
    list(
      updated = FALSE,
      convergence = NA_integer_,
      value = NA_real_,
      q12 = Qm[[m]][1, 2],
      q21 = Qm[[m]][2, 1],
      qsum = Qm[[m]][1, 2] + Qm[[m]][2, 1],
      p12 = Pm[[m]][1, 2],
      p21 = Pm[[m]][2, 1],
      det = 1 - Pm[[m]][1, 2] - Pm[[m]][2, 1]
    )
  })
  if (update_transitions) {
    for (m in 1:12) {
      idx_t <- which(month_idx[1:(Tn - 1)] == m)
      if (length(idx_t) == 0) next

      xi_m <- xi[idx_t, , , drop = FALSE]
      w_m <- weights[idx_t]
      if (sum(w_m) <= 0) next

      q_start <- c(Qm[[m]][1, 2], Qm[[m]][2, 1])
      q_start <- pmax(q_start, 1e-6)

      objective <- function(eta) {
        q12 <- exp(eta[1])
        q21 <- exp(eta[2])
        P <- pmax(ctmc_P_2state(q12, q21, delta = 1), eps_small)
        ll_t <- xi_m[, 1, 1] * log(P[1, 1]) +
          xi_m[, 1, 2] * log(P[1, 2]) +
          xi_m[, 2, 1] * log(P[2, 1]) +
          xi_m[, 2, 2] * log(P[2, 2])
        out <- -sum(w_m * ll_t)
        if (!is.finite(out)) Inf else out
      }

      fit <- try(
        stats::optim(
          par = log(q_start),
          fn = objective,
          method = "BFGS",
          control = list(maxit = 200)
        ),
        silent = TRUE
      )

      if (!inherits(fit, "try-error") && is.finite(fit$value)) {
        q_hat <- exp(fit$par)
        Qm_new[[m]] <- matrix(c(-q_hat[1], q_hat[1],
                                q_hat[2], -q_hat[2]), 2, 2, byrow = TRUE)
        Pm_new[[m]] <- ctmc_P_2state(q_hat[1], q_hat[2], delta = 1)
        transition_optim[[m]] <- list(
          updated = TRUE,
          convergence = fit$convergence,
          value = fit$value,
          q12 = q_hat[1],
          q21 = q_hat[2],
          qsum = sum(q_hat),
          p12 = Pm_new[[m]][1, 2],
          p21 = Pm_new[[m]][2, 1],
          det = 1 - Pm_new[[m]][1, 2] - Pm_new[[m]][2, 1]
        )
      } else {
        transition_optim[[m]] <- list(
          updated = FALSE,
          convergence = NA_integer_,
          value = NA_real_,
          q12 = q_start[1],
          q21 = q_start[2],
          qsum = sum(q_start),
          p12 = Pm[[m]][1, 2],
          p21 = Pm[[m]][2, 1],
          det = 1 - Pm[[m]][1, 2] - Pm[[m]][2, 1],
          warning = "optim failed; previous transition retained"
        )
      }
    }
  }

  mu_new  <- mu
  sig_new <- sig
  if (update_emissions) {
    for (m in 1:12) {
      idx_m <- which(month_idx == m)
      if (length(idx_m) == 0) next
      for (i in 1:2) {
        w <- gamma[idx_m, i]
        wsum <- sum(weights[idx_m] * w)
        if (wsum <= 0) next
        mu_i <- sum(w * weights[idx_m] * r_t[idx_m]) / wsum
        var_i <- sum(w * weights[idx_m] * (r_t[idx_m] - mu_i)^2) / wsum
        mu_new[[m]][i]  <- mu_i
        sig_new[[m]][i] <- sqrt(max(var_i, eps_small))
      }
    }
  }

  alpha <- matrix(NA_real_, nrow(gamma), 2)
  lSE <- function(z) {
    m <- max(z)
    m + log(sum(exp(z - m)))
  }
  for (t in 1:nrow(alpha)) {
    z <- logalpha[t, ]
    alpha[t, ] <- exp(z - lSE(z))
  }

  w_pred <- matrix(NA_real_, nrow(alpha), 2)
  w_pred[1, ] <- pi_new / sum(pi_new)
  for (t in 2:nrow(alpha)) {
    m_prev <- as.integer(format(as.Date(dates[t - 1]), "%m"))
    w_pred[t, ] <- as.numeric(alpha[t - 1, ] %*% Pm_new[[m_prev]])
    w_pred[t, ] <- pmax(w_pred[t, ], eps_small)
    w_pred[t, ] <- w_pred[t, ] / sum(w_pred[t, ])
  }

  colnames(alpha) <- paste0("alpha", 1:ncol(alpha))
  colnames(gamma) <- paste0("gamma", 1:ncol(gamma))
  colnames(w_pred) <- paste0("w_pred", 1:ncol(w_pred))

  mu_new <- purrr::map(mu_new, ~setNames(.x, c("mu1", "mu2")))
  sig_new <- purrr::map(sig_new, ~setNames(.x, c("sd1", "sd2")))

  list(
    data = data,
    params = list(
      pi = pi_new,
      Pm = Pm_new,
      Qm = Qm_new,
      mu = mu_new,
      sig = sig_new,
      a = params$a,
      theta = params$theta,
      b = params$b
    ),
    alpha = alpha,
    gamma = gamma,
    w_pred = w_pred,
    xi = xi,
    loglik = loglik,
    transition_optim = transition_optim
  )
}

#' Fit a Monthly Two-State CTMC by Direct Rate Estimation
#'
#' Fit the same monthly two-state regime model as `dtmc_monthly_fit()`, but in
#' each transition M-step optimize the position-indexed CTMC rates `q12_m` and
#' `q21_m` directly.
#'
#' @param Y Numeric vector. Transformed radiation observations.
#' @param dates Date vector aligned with `Y`.
#' @param weights Numeric vector of observation weights.
#' @param model Fitted discrete radiation model used for starting values.
#' @param p0 Numeric length-two initial state probability.
#' @param maxit Integer. Maximum EM iterations.
#' @param tol Numeric. Log-likelihood convergence tolerance.
#'
#' @return Fitted CTMC list returned by `dtmc_monthly_EM_direct()`.
#' @keywords internal
dtmc_monthly_fit_direct <- function(Y, dates, weights, model, p0 = c(0.5, 0.5), maxit = 1000, tol = 0.01) {
  NM_model <- model$spec$mixture.model
  params <- list()
  params[["mu"]] <- purrr::map(1:12, ~unlist(NM_model$means[.x, ]))
  params[["sig"]] <- purrr::map(1:12, ~unlist(NM_model$sd[.x, ]))

  build_P_from_probs <- function(pi_vec, lambda0 = 0.3) {
    p1 <- pi_vec[1]
    p2 <- pi_vec[2]
    a <- lambda0 * p2
    b <- lambda0 * p1
    matrix(c(1 - a, a,
             b, 1 - b), nrow = 2, byrow = TRUE)
  }

  probs <- purrr::map(1:12, ~unlist(NM_model$p[.x, ]))
  params[["Pm"]] <- purrr::map(1:12, ~build_P_from_probs(probs[[.x]], lambda0 = 0.3))
  params[["Qm"]] <- transition_list_to_generator_2state(params[["Pm"]], delta = 1)
  params[["pi"]] <- p0
  params$a <- model$spec$seasonal.mean$coefficients
  params$theta <- -log(model$spec$mean.model$phi)
  params$b <- seasonalModel_params_to_zeta(model$spec$seasonal.variance$coefficients)

  data <- dtmc_monthly_preprocess_direct(Y, dates, weights, params, eps = 1e-8)

  loglik <- -Inf
  em <- NULL
  loglik_trace <- numeric(0)
  for (iter in seq_len(maxit)) {
    em <- dtmc_monthly_EM_direct(
      data = data$data,
      params = params,
      update_emissions = TRUE,
      update_transitions = TRUE,
      eps = 1e-8
    )
    print(em$loglik)
    loglik_trace <- c(loglik_trace, em$loglik)

    if (!is.finite(em$loglik)) break
    if (iter > 1 && abs(loglik - em$loglik) <= tol) break

    params <- em$params
    loglik <- em$loglik
  }

  em$iter <- length(loglik_trace)
  em$loglik_trace <- loglik_trace
  em
}

#' Fit the Continuous-Time Radiation CTMC by Direct Rate Estimation
#'
#' Variant of `radiationModel_CTMC_fit()` that estimates monthly CTMC
#' generators directly during the transition M-step. The returned object has
#' the same broad structure as the existing CTMC fit and includes both `Qm`
#' and the implied one-day transition matrices `Pm`.
#'
#' @param model Fitted discrete radiation model.
#' @param p0 Numeric length-two initial state probability.
#' @param maxit Integer. Maximum EM iterations.
#' @param tol Numeric. Log-likelihood convergence tolerance.
#'
#' @return Fitted CTMC list with parameters, filtered probabilities, and data.
#' @keywords internal
radiationModel_CTMC_fit_direct <- function(model, p0 = c(0.5, 0.5), maxit = 1000, tol = 0.01) {
  data_train <- dplyr::filter(model$data, isTrain)

  em_train <- dtmc_monthly_fit_direct(
    Y = data_train$Yt,
    dates = data_train$date,
    weights = data_train$weights,
    model = model,
    p0 = p0,
    maxit = maxit,
    tol = tol
  )

  data <- model$data
  preproc <- dtmc_monthly_preprocess_direct(
    Y = data$Yt,
    dates = data$date,
    params = em_train$params,
    eps = 1e-8
  )

  em_full <- dtmc_monthly_EM_direct(
    data = preproc$data,
    params = preproc$params,
    update_emissions = FALSE,
    update_transitions = FALSE,
    eps = 1e-8
  )

  em_full$transition_optim <- em_train$transition_optim
  em_full$training_loglik <- em_train$loglik
  em_full$training_iter <- em_train$iter
  em_full$training_loglik_trace <- em_train$loglik_trace
  em_full
}
