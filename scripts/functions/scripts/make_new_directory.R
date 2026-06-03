#' Create a Directory If Needed
#'
#' Create an output directory when it does not already exist and optionally
#' print a CLI confirmation message.
#'
#' @param dir_output Character scalar. Directory path to create.
#' @param quiet Logical. If `TRUE`, suppress the confirmation message.
#'
#' @return Invisibly returns `NULL`; called for its filesystem side effect.
#' @keywords internal
make_new_directory <- function(dir_output, quiet = FALSE){
  # Initialize a directory 
  if (!file.exists(dir_output)) {
    system(paste0("mkdir ", dir_output))
    if (!quiet) cli::cli_alert_success(paste0("New folder: ", "\033[1;32m", dir_output, "\033[0m", " created!\n"))
  }
}
