#' Print Command-Line Script Arguments
#'
#' Print named script arguments in the standard project style.
#'
#' @param ... Named values to print.
#'
#' @return Invisibly returns `NULL`; called for console output.
#' @keywords internal
print_script_args <- function(...){
  .l <- list(...)
  .names <- names(.l)
  msg <- "\n"
  for(i in 1:length(.l)){
    param_i <- .l[[i]]
    if(is.logical(param_i)){
      msg[i] <- paste0(.names[i], ": ", ifelse(param_i, crayon::green$bold(param_i), crayon::red$bold(param_i)), "\n")
    } else {
      msg[i] <- paste0(.names[i], ": ", crayon::blue$bold(param_i), "\n")
    }
  }
  msg <- paste0(msg, collapse = "")
  cat(crayon::red$bold(paste0(c(rep("-", 50), "\n"), collapse = "")))
  cat(msg)
  cat(crayon::red$bold(paste0(c(rep("-", 50), "\n"), collapse = "")))
}



