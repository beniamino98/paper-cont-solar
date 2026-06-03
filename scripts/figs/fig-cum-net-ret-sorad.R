# ---
#' @description
#' Generate cumulative net-return figures for SoRad contracts.
#'
#' @section `main`
#' @label `fig-cum-net-ret-sorad`
#' @name `fig-cum-net-ret-sorad`
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse, mixtools
#'
#' @arguments
#'   - param[1] (save_plot): write the generated figures ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/figs/fig-cum-net-ret-sorad.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/solarOptions/mv/Bologna/data_hedging.RData`
#'   - `data/solarOptions/mv/Roma/data_hedging.RData`
#'   - `data/solarOptions/mv/Palermo/data_hedging.RData`
#'
#' @outputs
#'   - `figs/fig-cum-net-ret-sorad.pdf`
#'
#' @depends
#'   - `scripts/data/s8-solarOptions-mv-place.R`
#'
#' @tags
#'   - figures
#'   - main
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "figs", "fig-cum-net-ret-sorad.R"))
# ******************************************************************************
#                                Inputs 
# ******************************************************************************
# Reference places 
places <- c(Bologna = "Bologna", Rome = "Roma", Palermo = "Palermo")
# Figure names
common_name <- "fig-cum-net-ret-sorad"
# Reference year
nyear <- "2018"
# Save plot
save_plot <- FALSE
# ******************************************************************************
#                             Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly=TRUE)
if (!purrr::is_empty(args)) {
  # Save figure
  save_plot <- ifelse(is.na(args[1]), save_plot, args[1])
  print_script_args(save_plot = save_plot)
}
# ******************************************************************************
#                                Load data
# ******************************************************************************
data_sorad <- data_soradidx <- list()
for(place in places){
  # Load scenarios 
  load_data(file.path("data/solarOptions/mv", place), "data_hedging")
  data_sorad[[place]] <- data_hedging$sorad
  data_soradidx[[place]] <- data_hedging$soradidx
}
# ******************************************************************************
#                             Generating functions
# ******************************************************************************
fig_cum_cash_flow_sorad <- function(data, place, y_lab = NULL, legend = FALSE){
  # Breaks and labels y-axis
  y_breaks <- seq(-2, 2, 0.1)
  y_labels <- paste0(format(y_breaks*100, digits = 1, scientific = FALSE), "%")
  # Breaks and labels x-axis
  x_breaks <- seq(0, 365, length.out = 6)
  x_labels <- round(x_breaks, 0)
  plt <- data %>%
    ggplot()+
    geom_line(aes(n, cum_ret_seller, group = Year, color = "ret"), alpha = 0.5, linetype = "dotted")+
    geom_line(aes(n, e_cum_ret_seller, group = Year, color = "avg"), alpha = 0.5)+
    labs(x = NULL, y = NULL, subtitle = place, color = NULL)+
    scale_y_continuous(breaks = y_breaks, labels = y_labels, 
                       sec.axis = dup_axis(
                         name = y_lab,
                         breaks = NULL,
                         labels = NULL),)+
    scale_x_continuous(breaks = x_breaks, labels = x_labels)+
    scale_color_manual(values = c(avg = "black", ret = "black"), 
                       labels = c(avg = "Avg.", ret = "Yearly")) +
    theme_bw()+
    figure_theme+
    theme(legend.position = "none")
  if (legend){
    plt +
      theme(legend.position = "top")
  } else {
    plt
  }
}
# ******************************************************************************
#                             Generate all figures
# ******************************************************************************
# Place[1]
fig_11 <- fig_cum_cash_flow_sorad(data_sorad[[1]], names(places[1]), legend = TRUE, y_lab = NULL)
fig_12 <- fig_cum_cash_flow_sorad(data_soradidx[[1]], NULL, y_lab = NULL)
# Place[2]
fig_21 <- fig_cum_cash_flow_sorad(data_sorad[[2]], names(places[2]), y_lab = NULL)
fig_22 <- fig_cum_cash_flow_sorad(data_soradidx[[2]],NULL, y_lab = NULL)
# Place[3]
fig_31 <- fig_cum_cash_flow_sorad(data_sorad[[3]], names(places[3]), y_lab = "SoRad")
fig_32 <- fig_cum_cash_flow_sorad(data_soradidx[[3]],NULL, y_lab = "SoRadIDX")
# ******************************************************************************
#                         Generate Figure: SoRad
# ******************************************************************************
# Build the full figure 
glist <- gridExtra::arrangeGrob(
  fig_11+theme(legend.position = "none"), fig_21, fig_31, 
  fig_12, fig_22, fig_32, nrow = 2)
# Common labels y and x axis 
yleft = gridtext::richtext_grob("Cum. net return (%)", rot = 90, gp = grid::gpar(fontsize = 25))
bottom = gridtext::richtext_grob(text = "Day of the year", gp = grid::gpar(fontsize = 25))
# Figure 
fig <- gridExtra::grid.arrange(g_legend(fig_11+theme(legend.position = "top")), glist, 
                               nrow=2,heights=c(1, 10), left = yleft, bottom = bottom)
# ******************************************************************************
#                              Save Figure  
# ******************************************************************************
if (save_plot){
  control <- outputs$dir$figs$control
  save_new_fig(dir_output = outputs$dir$figs$main, fig, fig.name = common_name[1], 
               file.format = stringr::str_remove_all(control$format, "\\."), 
               quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
}
