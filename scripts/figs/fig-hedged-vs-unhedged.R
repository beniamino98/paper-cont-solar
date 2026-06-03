# ---
#' @description
#' Generate hedged-versus-unhedged cash-flow figures for SoRad and SoREd contracts.
#'
#' @section `main`
#' @label `fig-hedged-vs-unhedged-sorad`, `fig-hedged-vs-unhedged-sored`
#' @name `fig-hedged-vs-unhedged-sorad`, `fig-hedged-vs-unhedged-sored`
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
#' Rscript scripts/figs/fig-hedged-vs-unhedged.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/solarOptions/mv/Bologna/data_hedging.RData`
#'   - `data/solarOptions/mv/Roma/data_hedging.RData`
#'   - `data/solarOptions/mv/Palermo/data_hedging.RData`
#'
#' @outputs
#'   - `figs/fig-hedged-vs-unhedged-sorad.pdf`
#'   - `figs/fig-hedged-vs-unhedged-sored.pdf`
#'
#' @depends
#'   - `scripts/data/s8-solarOptions-mv-place.R`
#'
#' @tags
#'   - figures
#'   - main
# ---
# Load the required functions
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "figs", "fig-hedged-vs-unhedged.R"))
# ******************************************************************************
#                                Inputs 
# ******************************************************************************
# Reference places 
places <- c(Bologna = "Bologna", Rome = "Roma", Palermo = "Palermo")
# Figure names
common_name <- c("fig-hedged-vs-unhedged-sorad", "fig-hedged-vs-unhedged-sored")
# Reference year
nyear <- "2018"
# Save plot
save_plot <- FALSE
# ******************************************************************************
#                                Load data
# ******************************************************************************
data_sorad <- data_soradidx <- list()
data_sored <- data_sored_strip <- data_soredidx <- list()
for(place in places){
  # Load scenarios 
  load_data(file.path("data/solarOptions/mv", place), "data_hedging")
  data_sorad[[place]] <- data_hedging$sorad
  data_soradidx[[place]] <- data_hedging$soradidx
  data_sored[[place]] <- data_hedging$sored
  data_sored_strip[[place]] <- data_hedging$sored_strip
  data_soredidx[[place]] <- data_hedging$soredidx_strip
}
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
#                             Generating functions
# ******************************************************************************
fig_buyer_cash_flows <- function(data, nyear, y_limits, y_lab = NULL, subtitle = NULL, seasonal = FALSE, legend = FALSE){
  if (is.null(data$Gamma[1])){
    data$Gamma <- data$payoff
  }
  data <- dplyr::filter(data, Year == nyear)
  # Limits for y-axis 
  if (missing(y_limits)){
    # Range of the cash flows 
    y_limits <- range(data$pi_buyer_uh)
  } 
  # Breaks and labels y-axis
  y_breaks <- seq(y_limits[1], y_limits[2], length.out = 4)
  y_labels <- round(y_breaks, 2)
  
  # Breaks and labels x-axis
  x_breaks <- seq(0, 365, length.out = 6)
  x_labels <- round(x_breaks, 0)
  
  plt <- ggplot(data)+
    # Unhedged cash flows             
    geom_point(aes(n, pi_buyer_uh, color = "uh"), alpha = 0.5, shape = 4, size = 1.3)+
    # Hedged cash flows (black)
    geom_point(aes(n, pi_buyer_h), size = 1.5)+
    # Hedged cash flows (color)
    geom_point(aes(n, pi_buyer_h, color = Gamma > 0), size = 1)
  
  # Seasonal cash flows 
  if (seasonal) {
    plt <- plt+
      geom_line(aes(n, pi_seasonal, color = "k"))
  }
   plt <- plt + 
    # Scale 
    scale_color_manual(values = c(`TRUE` = "green", `FALSE` = "red", uh = "black", k = "black"), 
                       labels = c(`TRUE` = latex2exp::TeX("Hedged($\\Gamma > 0$)"), 
                                  `FALSE` = latex2exp::TeX("Hedged($\\Gamma = 0$)"),
                                  k = "Strike",
                                  uh = "Unhedged"))+
    scale_x_continuous(breaks = x_breaks, labels = x_labels)+
    scale_y_continuous(breaks = y_breaks, labels = y_labels, limits = y_limits,
                        sec.axis = dup_axis(
                          name = y_lab,
                          breaks = NULL,
                          labels = NULL))+
    # Labels 
    labs(x = NULL, y = NULL, color = NULL, subtitle = subtitle)+
    # Theme
    theme_bw()+
    theme(legend.position = "top", panel.border = element_blank(),
          panel.spacing = element_blank(),
          strip.background = element_rect(colour = "black", fill = "white"),
          panel.background = element_blank(),
          strip.text = element_text(angle = 0, face = "bold", size = 15))+
    figure_theme
  
  if(!legend){
    plt + theme(legend.position = "none")
  } else {
    plt
  }
}
# ******************************************************************************
#                         Generate all figure: cash flows 
# ******************************************************************************
# Place[1]
fig_11 <- fig_buyer_cash_flows(data_sorad[[1]], nyear = nyear, seasonal = TRUE, subtitle = names(places[1]), legend = TRUE)
fig_12 <- fig_buyer_cash_flows(data_soradidx[[1]],  seasonal = TRUE, nyear = nyear)
fig_13 <- fig_buyer_cash_flows(data_sored[[1]], nyear = nyear, legend = TRUE)
fig_14 <- fig_buyer_cash_flows(data_sored_strip[[1]], nyear = nyear)
fig_15 <- fig_buyer_cash_flows(data_soredidx[[1]], nyear = nyear)
# Place[2]
fig_21 <- fig_buyer_cash_flows(data_sorad[[2]], seasonal = TRUE, nyear = nyear, subtitle = names(places[2]))
fig_22 <- fig_buyer_cash_flows(data_soradidx[[2]], seasonal = TRUE, nyear = nyear)
fig_23 <- fig_buyer_cash_flows(data_sored[[2]], nyear = nyear)
fig_24 <- fig_buyer_cash_flows(data_sored_strip[[2]], nyear = nyear)
fig_25 <- fig_buyer_cash_flows(data_soredidx[[2]], nyear = nyear)
# Place[3]
fig_31 <- fig_buyer_cash_flows(data_sorad[[3]], y_lab = "SoRad", seasonal = TRUE, nyear = nyear, subtitle = names(places[3]))
fig_32 <- fig_buyer_cash_flows(data_soradidx[[3]], y_lab = "SoRadIDX", seasonal = TRUE, nyear = nyear)
fig_33 <- fig_buyer_cash_flows(data_sored[[3]], y_lab = "SoREd", nyear = nyear)
fig_34 <- fig_buyer_cash_flows(data_sored_strip[[3]], y_lab = "SoREd+futures", nyear = nyear)
fig_35 <- fig_buyer_cash_flows(data_soredidx[[3]], y_lab = "SoREdIDX+futures", nyear = nyear)
# ******************************************************************************
#                        Generate figure: cash flows SoRad 
# ******************************************************************************
# Full figure
glist <- gridExtra::arrangeGrob(
  fig_11+ theme(legend.position = "none"), fig_21, fig_31,
  fig_12, fig_22, fig_32,
  nrow = 2, ncol = 3)
# Labels 
yleft = gridtext::richtext_grob("Daily P&L (Eur)", rot = 90, gp = grid::gpar(fontsize = 25))
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
# ******************************************************************************
#                         Generate figure: SoREd
# ******************************************************************************
glist <- gridExtra::arrangeGrob(
  fig_13+labs(subtitle = names(places[1]))+ theme(legend.position = "none"), 
  fig_23+labs(subtitle = names(places[2])), 
  fig_33+labs(subtitle = names(places[3])),
  fig_14, fig_24, fig_34,
  fig_15, fig_25, fig_35,
  nrow = 3, ncol = 3)
# labels
yleft = gridtext::richtext_grob("Daily P&L (Eur)", rot = 90, gp = grid::gpar(fontsize = 25))
bottom = gridtext::richtext_grob(text = "Day of the year", gp = grid::gpar(fontsize = 25))
fig <- gridExtra::grid.arrange(g_legend(fig_13+theme(legend.position = "top")), glist, 
                               nrow=2,heights=c(1, 10), left = yleft, bottom = bottom)
# ******************************************************************************
#                              Save Figure  
# ******************************************************************************
if (save_plot){
  control <- outputs$dir$figs$control
  save_new_fig(dir_output = outputs$dir$figs$main, fig, fig.name = common_name[2], 
               file.format = stringr::str_remove_all(control$format, "\\."), 
               quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
}
