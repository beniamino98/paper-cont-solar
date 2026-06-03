#' Print Script Header Information
#'
#' Read the roxygen-style header delimited by `# ---` in a script and print it
#' with the current last-modified time.
#'
#' @param dir_script Character scalar. Path to the script.
#'
#' @return Invisibly returns `NULL`; called for console output.
#' @keywords internal
print_script_info <- function(dir_script){
  txt <- suppressWarnings(readLines(dir_script))
  info <- file.info(dir_script)
  # Detect description
  idx <- which(stringr::str_detect(txt, "---"))
  if (purrr::is_empty(idx) | is.na(idx[2])){
    return(invisible(NULL))
  }
  txt <- txt[(idx[1]+1):(idx[2]-1)]
  txt <- purrr::map_chr(txt, ~paste0(.x, "\n"))
  time.passed <- count_time_passed(info$mtime)
  date_modified <- as.Date(info$mtime)
  txt[stringr::str_detect(txt, "#' @modified ")] <- paste0("#' @modified ", date_modified, " (", time.passed, ")", "\n")
  tags <- unlist(stringr::str_extract_all(txt, "@[a-zA-Z]+"))
  for(tag in tags){
    txt <- stringr::str_replace_all(txt, tag, crayon::yellow(tag))
  }
  txt[1] <- paste0(" ", txt[1])
  cat(txt)
}
