#' Save an R Object With a Project Message
#'
#' Save one or more R objects to disk using `save()` and print a standardized
#' success message.
#'
#' @param dir_output Character scalar or `NULL`. Output directory. If `NULL`,
#'   the current working directory is used.
#' @param file.name Character scalar. File name without extension.
#' @param file.format Character scalar. File extension, usually `"RData"` or
#'   `"rds"` where appropriate.
#' @param quiet Logical. If `TRUE`, suppress the confirmation message.
#' @param ... Objects passed to `save()`.
#'
#' @return Invisibly returns `NULL`; called for its filesystem side effect.
#' @keywords internal
save_new_file <- function(dir_output = NULL, file.name = "", file.format = "RData", quiet = FALSE, ...){
  file.format <- ifelse(missing(file.format), "RData", file.format)
  # Create a custom filename
  filename <- paste0(file.name, ".", file.format)
  # Output 
  dir_file <- filename
  if (!is.null(dir_output)) {
    dir_file <- file.path(dir_output, dir_file)
  } else {
    dir_output <- tail(strsplit(getwd(), "/")[[1]],1)
  }
  # Save the model
  save(..., file = dir_file)
  msg <- paste0("File: ", "\033[1;35m", filename, "\033[0m", " saved in ", "\033[1;34m", dir_output, "\033[0m", "\n")
  #msg <- paste0("File: ", crayon::magenta(filename), " saved in ", crayon::blue(dir_output), "\n")
  if (!quiet) cli::cli_alert_success(msg, "\n")
}

