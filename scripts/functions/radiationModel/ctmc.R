#' Convert a Two-State Transition Matrix to a CTMC Generator
#'
#' Recover the unique regular two-state CTMC generator embedding a one-step
#' transition matrix over interval `delta`.
#'
#' @param P Numeric `2 x 2` row-stochastic transition matrix.
#' @param delta Numeric scalar. Time interval represented by `P`.
#' @param tol Numeric tolerance for validation and embeddability checks.
#'
#' @return Numeric `2 x 2` CTMC generator matrix.
#' @export
transition_to_generator_2state <- function(P, delta = 1, tol = 1e-10) {
  P <- as.matrix(P)
  condition <- ctmc_is_transition_2state(P, tol)
  # Check is 2x2 matrix 
  if (!condition) {
    stop("P is not a valid transition matrix.")
  }
  P[P < 0 & P > -tol] <- 0
  p12 <- P[1, 2]
  p21 <- P[2, 1]
  s <- p12 + p21
  if (s < tol) {
    return(matrix(c(0, 0, 0, 0), 2, 2, byrow = TRUE))
  }
  detP <- 1 - s
  # Check CTMC condition 
  if (detP <= tol) {
    stop(
      "Transition matrix is not embeddable as a regular two-state CTMC over this delta: ",
      "need 1 - p12 - p21 > 0."
    )
  }
  # Lambda 
  lambda <- -log(detP) / delta
  # With state order (1, 0), q12 is B = 1 -> 0 and q21 is B = 0 -> 1.
  q12 <- lambda * p12 / s
  q21 <- lambda * p21 / s
  # Generator matrix Q
  matrix(c(-q12, q12,
           q21, -q21), 2, 2, byrow = TRUE)
}

#' Convert a List of Two-State Transition Matrices to Generators
#'
#' Apply `transition_to_generator_2state()` to a monthly list of transition
#' matrices.
#'
#' @param Pm List of `2 x 2` transition matrices.
#' @param delta Numeric scalar. Time interval represented by each matrix.
#' @param tol Numeric tolerance passed to `transition_to_generator_2state()`.
#'
#' @return List of `2 x 2` CTMC generators.
#' @export
transition_list_to_generator_2state <- function(Pm, delta = 1, tol = 1e-10) {
  if (!is.list(Pm)) {
    stop("Pm must be a list of transition matrices.")
  }
  lapply(Pm, transition_to_generator_2state, delta = delta, tol = tol)
}

#' Test Whether a Matrix Is a Two-State CTMC Generator
#'
#' @param Q Candidate numeric matrix.
#' @param tol Numeric tolerance.
#' @param quiet Logical. If `TRUE`, suppress validation warnings.
#'
#' @return Logical scalar.
#' @keywords internal
ctmc_is_generator_2state <- function(Q, tol = 1e-8, quiet = FALSE) {
  condition <- TRUE
  # Check matrix type.
  if (!is.matrix(Q)) {
    if (!quiet) warning("Q must be a matrix.")
    condition <- FALSE
    return(condition)
  }
  # Check matrix dimensions.
  if (!all(dim(Q) == c(2, 2))) {
    if (!quiet) warning("Q must be a 2x2 matrix.")
    condition <- FALSE
    return(condition)
  }
  if (any(!is.finite(Q))) {
    if (!quiet) warning("Q contains non-finite values.")
    condition <- FALSE
    return(condition)
  }
  if (max(abs(rowSums(Q))) > tol) {
    if (!quiet) warning("Rows of Q must sum to 0.")
    condition <- FALSE
    return(condition)
  }

  if (Q[1, 2] < -tol) {
    if (!quiet) warning("Q[1,2] is lower than -tol.")
    condition <- FALSE
    return(condition)
  }
  if (Q[2, 1] < -tol) {
    if (!quiet) warning("Q[2,1] is lower than -tol.")
    condition <- FALSE
    return(condition)
  }
  condition
}

#' Test Whether a Matrix Is a Two-State Transition Matrix
#'
#' @param P Candidate numeric matrix.
#' @param tol Numeric tolerance.
#' @param quiet Logical. If `TRUE`, suppress validation warnings.
#'
#' @return Logical scalar.
#' @keywords internal
ctmc_is_transition_2state <- function(P, tol = 1e-8, quiet = FALSE) {
  condition <- TRUE
  # Check matrix type.
  if (!is.matrix(P)) {
    if (!quiet) warning("P must be a matrix.")
    condition <- FALSE
    return(condition)
  }
  # Check matrix dimensions.
  if (!all(dim(P) == c(2, 2))) {
    if (!quiet) warning("P must be a 2x2 transition matrix.")
    condition <- FALSE
    return(condition)
  }
  if (any(!is.finite(P))) {
    if (!quiet) warning("P contains non-finite values.")
    condition <- FALSE
    return(condition)
  }
  if (max(abs(rowSums(P) - 1)) > 1e-8) {
    if (!quiet) warning("Rows of P must sum to 1.")
    condition <- FALSE
    return(condition)
  }
  if (any(P < -tol)) {
    if (!quiet) warning("P contains negative probabilities.")
    condition <- FALSE
    return(condition)
  }
  condition
}

#' Coerce Monthly CTMC Inputs to Generator Form
#'
#' Accept either a list of generators or a list of transition matrices and
#' return generators.
#'
#' @param Q_or_P List of monthly `2 x 2` matrices.
#' @param delta Numeric scalar used when matrices are transitions.
#' @param tol Numeric tolerance.
#'
#' @return List of monthly CTMC generators.
#' @keywords internal
ctmc_as_generator_list <- function(Q_or_P, delta = 1, tol = 1e-10) {
  if (!is.list(Q_or_P)) {
    stop("Monthly transitions must be supplied as a list.")
  }
  is_generator <- vapply(Q_or_P, ctmc_is_generator_2state, logical(1), quiet = TRUE)
  is_transition <- vapply(Q_or_P, ctmc_is_transition_2state, logical(1),quiet = TRUE)

  if (all(is_generator)) {
    return(Q_or_P)
  }
  if (all(is_transition)) {
    return(transition_list_to_generator_2state(Q_or_P, delta = delta, tol = tol))
  }

  stop("Monthly transitions must be either all 2x2 CTMC generators or all 2x2 transition matrices.")
}

#' Convert DTMC in CTMC 
#'
#' @param dtmc Fitted dtmc, see `dtmc_monthly_fit()`. 
#' @param delta Numeric scalar used when matrices are transitions.
#' @param tol Numeric tolerance.
#'
#' @return CTMC model. 
#' @keywords internal
ctmc_from_dtcm <- function(dtmc, delta = 1, tol = 1e-10){
  # Clone object
  ctmc <- dtmc
  # Convert P in Q
  ctmc$params$Qm <- purrr::map(ctmc$params$Pm, ~transition_to_generator_2state(.x, delta = delta, tol = tol))
  ctmc
}

#' Simulate CTMC Regime Paths for a Radiation HMM
#'
#' Simulate regime indicators and Brownian increments using monthly CTMC
#' transition matrices recovered from fitted HMM transition probabilities.
#'
#' @param CTMC Fitted radiation HMM object/list.
#' @param p0 Numeric length-two initial filtered probability.
#' @param month_idx Integer vector of monthly indices for each simulation step.
#' @param dt Numeric scalar. Time step in days.
#' @param nsteps Integer. Number of time steps.
#' @param nsim Integer. Number of Monte Carlo paths.
#'
#' @return List with matrices `mu_B`, `sigma_B`, `B_t`, and `dMt`.
#' @keywords internal
ctmc_monthly_scenarios <- function(CTMC, p0, month_idx, dt = 1, nsteps = 10, nsim = 1){
  # Initialization 
  B_t <- matrix(0, nrow = nsteps, ncol = nsim)
  dMt <- matrix(0, nrow = nsteps, ncol = nsim)
  mu_B    <- matrix(0, nrow = nsteps, ncol = nsim)
  sigma_B <- matrix(0, nrow = nsteps, ncol = nsim)
  # Initial probability
  p0 <- as.numeric(p0)
  p0 <- p0 / sum(p0)
  # Recover the CTMC generator from the fitted one-step transition matrices.
  CTMC_mu <- CTMC$params$mu
  CTMC_sd <- CTMC$params$sig
  if (!is.null(CTMC$params$Qm)){
    CTMC_Qm <- CTMC$params$Qm
  } else {
    CTMC_Qm <- transition_list_to_generator_2state(CTMC$params$Pm)
  }
  for (sim in 1:nsim) {
    message("Simulations CTMC regimes: ", sim, "/", nsim, "\r", appendLF = FALSE)
    # Initial state sampled from filtered probability at t_now
    state <- sample(1:2, size = 1, prob = p0)
    # Evolve the Markov chain 
    for (t in 1:nsteps) {
      # Next step probability 
      P_dt <- matrix_exponential(dt, CTMC_Qm[[month_idx[t]]])
      # Next step state 
      state <- sample(1:2, size = 1, prob = P_dt[state, ])
      # Store state as 1/0 to preserve convention (B_t = 1 corresponds to state 1)
      B_t[t, sim] <- ifelse(state == 1, 1, 0)
      # Realized moments 
      mu_B[t, sim] <- CTMC_mu[[month_idx[t]]][state]
      sigma_B[t, sim] <- CTMC_sd[[month_idx[t]]][state]
    }
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
