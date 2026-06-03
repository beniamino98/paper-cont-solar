#' Save a ggplot Figure With a Project Message
#'
#' Save a figure with `ggplot2::ggsave()` and optionally print a standardized
#' success message.
#'
#' @param dir_output Character scalar. Output directory.
#' @param fig A ggplot object.
#' @param fig.name Character scalar. File name without extension.
#' @param file.format Character scalar. Output extension, for example `"jpg"`
#'   or `"pdf"`.
#' @param quiet Logical. If `TRUE`, suppress the confirmation message.
#' @param dpi Numeric. Resolution passed to `ggsave()`.
#' @param width Numeric. Figure width passed to `ggsave()`.
#' @param height Numeric. Figure height passed to `ggsave()`.
#'
#' @return Invisibly returns `NULL`; called for its filesystem side effect.
#' @keywords internal
save_new_fig <- function(dir_output, fig, fig.name = "", file.format = "jpg", quiet = FALSE, dpi = 500, width = 15, height = 10){
  # File.format
  file.format <- ifelse(missing(file.format), "jpg", file.format)
  # Create a custom filename
  filename <- paste0(fig.name, ".", file.format)
  # Output 
  filepath <- file.path(dir_output, filename)
  # Save the figure
  ggsave(plot = fig, filename = filepath, dpi = dpi, width = width, height = height)
  if (!quiet) cli::cli_alert_success(paste0("Figure: ", "\033[1;35m", filename, "\033[0m", " saved in ", "\033[1;34m", dir_output, "\033[0m", "\n"))
}
