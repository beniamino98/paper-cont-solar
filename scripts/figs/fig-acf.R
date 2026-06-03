# ---
#' @description
#' Generate the autocorrelation diagnostic figure for CTMC radiation-model residuals.
#'
#' @section `main`
#' @label `fig-acf`
#' @name `fig-acf`
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse, mixtools
#'
#' @arguments
#'   - param[1] (place): reference location, e.g. "Bologna".
#'   - param[2] (save_plot): write the generated figure ("TRUE" or "FALSE").
#'   - param[3] (nyear): training year used for the diagnostic figure.
#'
#' @example
#' Rscript scripts/figs/fig-acf.R "Bologna" "TRUE" "2022"
#' for place in Bologna Palermo Roma; do for year in {2013..2022}; do Rscript scripts/figs/fig-acf.R "$place" "TRUE" "$year"; done; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData` (CTMC radiation models)
#'
#' @outputs
#'   - `figs/fig-acf-{nyear}.pdf` for the paper figure
#'   - `figs/models/radiation/fig-acf/{place}/fig-acf-{nyear}.pdf`
#'
#' @depends
#'   - `scripts/data/s2b-models-radiation-P-CTMC-place.R`
#'
#' @tags
#'   - figures
#'   - main
# ---
# Load the required functions
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "figs", "fig-acf.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Bologna"
# Save the plot 
save_plot <- TRUE
# Reference year
nyear <- 2022
# ***************************** Fixed Arguments ********************************
# Figures base name
common_name <- "fig-acf"
# ******************************************************************************
#                             Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly = TRUE)
if (!purrr::is_empty(args)) {
  # Reference location
  place <- ifelse(is.na(args[1]), place, args[1])
  # Save the plot 
  save_plot <- ifelse(is.na(args[2]), save_plot, args[2])
  # Reference year
  nyear <- as.character(ifelse(is.na(args[3]), nyear, args[3]))
  print_script_args(place = place, save_plot = save_plot, nyear = nyear)
}
# ******************************************************************************
#                               Load data 
# ******************************************************************************
# Load CTMC models
load_data(file.path(outputs$dir$data$models$radiation$P, place), "radiation_models_CTMC")
model_Rt <- radiation_models_CTMC[[as.character(nyear)]]
# ******************************************************************************
#                           Generating functions  
# ******************************************************************************
fig_acf <- function(x, lag_max = 365, label = "x", caption = NULL, ci = 0.05, ci_color = "blue", limits){
  # residuals
  lag.index <- seq.default(0, lag_max, 1)
  x_breaks <- seq.default(0, lag_max, length.out = 5)
  # Autocorrelation ut
  acf_x <- acf(x, lag.max = lag_max, plot = FALSE)
  acf_x_bounds <- qnorm(1 - ci/2)/sqrt(acf_x$n.used)
  plot <- ggplot() +
    geom_segment(aes(x = lag.index, xend = lag.index,
                     y = acf_x$acf[,,1],
                     yend = 0)) +
    geom_point(aes(x = lag.index, y = acf_x$acf[,,1]), size = 0.2) +
    geom_line(aes(lag.index, 0)) +
    geom_line(aes(lag.index, acf_x_bounds), color = ci_color,
              linetype = "dashed") +
    geom_line(aes(lag.index, -acf_x_bounds), color = ci_color,
              linetype = "dashed") +
    scale_x_continuous(breaks = x_breaks,
                       labels = round(x_breaks)) +
    labs(x = NULL, y = latex2exp::TeX(paste0("$", label, "$"))) +
    theme_bw()+
    #figure_theme+
    theme(# Title
      plot.title  = element_text(face = "bold", size = 30),
      # Subtitle
      plot.subtitle = element_text(size = 24),
      # Caption
      plot.caption = element_text(face = "italic"),
      # Axis-x
      axis.title.x = element_text(face = "bold", size = 20),
      axis.text.x = element_text(face = "bold", size = 17),
      axis.ticks.x = element_line(linewidth = 0.2),
      axis.line.x = element_line(),
      # Grid x-axis
      panel.grid.minor.x = element_line(),
      panel.grid.major.x = element_line(),
      # Axis-y
      axis.title.y = element_text(size = 25),
      axis.text.y = element_text(size = 20),
      axis.ticks.y = element_line(linewidth = 0.2),
      axis.line.y = element_line(),
      # Grid x-axis
      panel.grid.minor.y = element_line(),
      panel.grid.major.y = element_line(),
      # Legend
      legend.title = element_text(face = "bold", size = 25),
      legend.text = element_text(face = "italic", size = 25),
      legend.box.background = element_rect())
  
  if(!missing(limits)){
    plot <- plot +
      scale_y_continuous(limits = limits)
  }
  plot
}
# Add z tilde
add_z_tilde_to_data <- function(model_Rt, nyear = NULL){
  # Extract data
  data <- model_Rt$model$data
  if (!is.null(nyear)) {
    data <- filter(data, Year <= (as.numeric(nyear)+1))
  }
  # Index month 
  t_month <- data$Month
  # Residuals 
  u_tilde <- data$u_tilde
  # CTMC object
  CTMC <- model_Rt$CTMC
  # Reference dates 
  ref_dates <- data$date
  # Predictive probabilities
  alpha_ctmc <- CTMC$alpha[CTMC$data$date %in% ref_dates,]
  
  z_tilde <- c()
  for(i in 1:length(t_month)){
    # Monthly parameters
    mu_m <- CTMC$params$mu[[t_month[i]]]
    sd_m <- CTMC$params$sig[[t_month[i]]]
    pT_m <- drop(alpha_ctmc[i-1,] %*% CTMC$params$Pm[[t_month[i]]])
    # Moments 
    mom <- GM_moments(mu_m, sd_m, pT_m)
    # Standardized residuals 
    z_tilde[i] <- (u_tilde[i] - mom$mean) / sqrt(mom$variance) 
  }
  data$z_tilde <- z_tilde
  data
}
# ******************************************************************************
#                            Generate Figure  
# ******************************************************************************
# Data 
data <- add_z_tilde_to_data(model_Rt)
limits <- range(acf(data$Yt, lag.max = 365, plot = FALSE)$acf[,,1])*1.01
# Acf in mean
fig_1 <- fig_acf(data$Yt, label = "Y_t", limits = limits)
fig_2 <- fig_acf(data$Yt_tilde, label = "\\widetilde{Y}_t",  limits = limits)
fig_3 <- fig_acf(data$eps, label = "\\epsilon_t",  limits = limits)
fig_4 <- fig_acf(data$z_tilde[-1], label = "\\tilde{z}_t",  limits = limits)
# Acf in variance
fig_12 <- fig_acf(data$Yt^2, label = "Y_t^2")
fig_22 <- fig_acf(data$Yt_tilde^2, label = "\\widetilde{Y}_t^2")
fig_32 <- fig_acf(data$eps^2, label = "\\epsilon_t^2")
fig_42 <- fig_acf(data$z_tilde[-1]^2, label = "\\tilde{z}_t^2")
# Labels 
yleft = gridtext::richtext_grob("ACF", rot = 90, gp = grid::gpar(fontsize = 25))
bottom = gridtext::richtext_grob(text = "Lag", gp = grid::gpar(fontsize = 25))
# Composite figure 
fig <- gridExtra::grid.arrange(fig_1, fig_2, fig_3, fig_4, 
                               fig_12, fig_22, fig_32, fig_42, nrow = 2, 
                               left = yleft, bottom = bottom)
# ******************************************************************************
#                             Save Figure  
# ******************************************************************************
if (save_plot){
  # Initialize a folder to store the figures
  dir_output <- file.path(outputs$dir$figs$main, common_name)
  make_new_directory(dir_output)
  # Initialize a folder to store the figures
  dir_output <- file.path(dir_output, place)
  make_new_directory(dir_output)
  control <- outputs$dir$figs$control
 # Save figure 
  save_new_fig(dir_output, fig, fig.name = paste0(common_name, "-", nyear), 
               file.format = stringr::str_remove_all(control$format, "\\."), 
               quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
  # Save the paper figure in main 
  if (place == "Bologna" & nyear == "2022"){
    save_new_fig(dir_output = outputs$dir$figs$main, fig, fig.name = common_name, 
                file.format = stringr::str_remove_all(control$format, "\\."), 
                quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
  }
}


