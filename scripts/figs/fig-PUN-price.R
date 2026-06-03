# ---
#' @description
#' Generate the PUN electricity-price time-series figure.
#'
#' @section `main`
#' @label `fig-PUN-price`
#' @name `fig-PUN-price`
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse, mixtools
#'
#' @arguments
#'   - param[1] (save_plot): write the generated figure ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/figs/fig-PUN-price.R "TRUE"
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/GME_daily.csv` (PUN price data)
#'
#' @outputs
#'   - `figs/fig-PUN-price.pdf`
#'
#' @depends
#'   - none
#'
#' @tags
#'   - figures
#'   - main
# ---
load("outputs.RData")
source(file.path("scripts", "s0-load.R"))
print_script_info(file.path("scripts", "figs", "fig-PUN-price.R"))
# ******************************************************************************
#                                  Inputs
# ******************************************************************************
# Logical value to denote if the plot should be saved
save_plot <- FALSE
# ***************************** Fixed Arguments ********************************
# Figures base name
common_name <- "fig-PUN-price"
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
#                                  Load data 
# ******************************************************************************
# PUN data
data <- readr::read_csv(outputs$dir_GME, show_col_types = FALSE)  %>%
  dplyr::select(date, PUN) %>%
  dplyr::mutate(Year = lubridate::year(date), PUN = as.numeric(PUN)/1000)
# ******************************************************************************
#                                Generate data 
# ******************************************************************************
# Labels for x-axis
x_breaks <- seq.Date(as.Date("2005-01-01"), as.Date("2024-01-01"), by = "4 year")
x_breaks <- unique(c(x_breaks, as.Date("2024-01-01")))
# Compute cumulated annual averages 
nyears <- 2005:2023
PUN_avg <- purrr::map_dbl(nyears, ~mean(dplyr::filter(data, Year == .x)$PUN))
data_avg <- dplyr::tibble(Year = nyears, t = paste0(Year, "-06-01"), PUN_avg = PUN_avg)
data <- dplyr::left_join(data, data_avg, by = "Year")
# Highlight structural break in 2021-2023
data_break <- dplyr::filter(data, date >= "2021-01-01" & date <= "2023-01-01")
# Labels for y-axis
y_breaks <- c(seq(min(data$PUN), max(data$PUN), length.out = 4), mean(data$PUN))
y_labels <- paste0(format(y_breaks, digits = 1), "")
# ******************************************************************************
#                                Generate figure 
# ******************************************************************************
fig <- ggplot(data)+
  geom_area(data = data_break, aes(date, y = max(data$PUN)), fill = "darkgray", alpha = 0.2)+
  geom_line(aes(date, PUN), linewidth = 0.2, alpha = 0.8)+
  geom_line(data = filter(data, date <= "2023-06-01"), aes(date, PUN_avg), color = "red", linetype = "solid", linewidth = 1)+
  geom_point(data = data_avg, aes(as.Date(t), PUN_avg), color = "red",  size = 3)+
  geom_point(data = data_avg, aes(as.Date(t), PUN_avg), color = "black", size = 3.7, shape = 1)+
  theme_bw()+
  labs(x = NULL, y = "PUN (Eur/kWh)")+
  scale_y_continuous(breaks = y_breaks, labels = y_labels)+
  scale_x_date(breaks = x_breaks, date_labels = "%Y")+
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank())+
  figure_theme
print(fig)
# ******************************************************************************
#                                Save figure 
# ******************************************************************************
if (save_plot){
  dir_output <- outputs$dir$figs$main
  control <- outputs$dir$figs$control
  save_new_fig(outputs$dir$figs$main, fig, fig.name = common_name, file.format = stringr::str_remove_all(control$format, "\\."), 
               quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
}
