#' Create Lagged Date Sequences
#'
#' Build paired current/horizon date sequences separated by a fixed lag.
#'
#' @param t_now Date or character scalar. Start date.
#' @param t_hor Date or character scalar. End date.
#' @param tau Integer. Lag in days between `now` and `hor`.
#'
#' @return List with date vectors `now` and `hor`.
#' @keywords internal
create_t_seq <- function(t_now, t_hor, tau = 1){
  t_now <- as.Date(t_now)
  t_hor <- as.Date(t_hor)
  t_seq <- seq.Date(t_now, t_hor, 1)
  structure(
    list(
      now = t_seq - tau,
      hor = t_seq
    )
  )
}
