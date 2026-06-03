#' Describe Time Since Last Modification
#'
#' Format the elapsed time since a file was last modified.
#'
#' @param last.modified POSIXct timestamp, typically from `file.info()$mtime`.
#'
#' @return Character scalar describing elapsed minutes, hours, or days.
#' @keywords internal
count_time_passed <- function(last.modified) {
  time.now <- Sys.time()
  # Count time passed from last modified 
  time.passed <- as.numeric(difftime(time.now, last.modified, units = "hours"))
  if (time.passed > 24) {
    time.passed <- as.numeric(difftime(time.now, last.modified, units = "days"))
    time.passed <- paste0(round(time.passed), " days ago.")
  } else if (time.passed < 1) {
    time.passed <- as.numeric(difftime(time.now, last.modified, units = "mins"))
    time.passed <- paste0(round(time.passed), " minutes ago.")
  } else {
    time.passed <- paste0(round(time.passed), " hours ago.")
  }
  time.passed
}

#' Load a Project Data File
#'
#' Load an `.RData` file from a project directory and optionally print file
#' metadata.
#'
#' @param dir_ Character scalar or `NULL`. Directory containing the file.
#' @param filename Character scalar. File name without extension.
#' @param file.format Character scalar. File extension.
#' @param quiet Logical. If `TRUE`, suppress console output.
#' @param envir Environment where loaded objects are assigned.
#'
#' @return Character vector returned by `load()`.
#' @keywords internal
load_data <- function(dir_ = getwd(), filename = "", file.format = "RData", quiet = FALSE, envir = .GlobalEnv){
  if (!is.null(dir_)) {
    dir_input <- dir_
    filepath <- file.path(dir_, paste0(filename, ".", file.format))
  } else {
    dir_input <- tail(strsplit(getwd(), "/")[[1]],1)
    filepath <- paste0(filename, ".", file.format)
  }
  # Get info 
  last.modified <- file.info(filepath)$mtime 
  h <- lubridate::hour(last.modified)
  m <- lubridate::minute(last.modified)
  s <- lubridate::second(last.modified)
  # Count time passed from last modified 
  time.passed <- count_time_passed(last.modified)
  last.modified <- as.character(as.Date(last.modified) )
  last.modified <- paste0(last.modified, " ", h, ":", m, ":", round(s), " (", time.passed, ")")
  
  if (!quiet) cat(crayon::red$bold(paste0(c(rep("-", 30), "\n"), collapse = "-")))
  if (!quiet) cat(paste0(" Loading Data: ", crayon::yellow$bold(filename) , "\n",
                         " Directory: ", crayon::blue$bold(dir_input), " \n", 
                         " Last modified: ", crayon::red$bold(last.modified), "\n"))
  load(filepath, envir = envir)
  if (!quiet) cat(crayon::red$bold(paste0(c(rep("-", 30), "\n"), collapse = "-")))
} 
