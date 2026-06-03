# ---
#' @description
#' Generate the simulated solar-radiation trajectory figure.
#'
#' @section `main`
#' @label `fig-GHI-sim`
#' @name `fig-GHI-sim`
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr, tidyverse, mixtools
#'
#' @arguments
#'   - param[1] (place): reference location, e.g. "Bologna".
#'   - param[2] (save_plot): write the generated figure ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/figs/fig-GHI-sim.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/figs/fig-GHI-sim.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData` (CTMC radiation models)
#'
#' @outputs
#'   - `figs/fig-GHI-sim.pdf` for the paper figure
#'   - `figs/models/radiation/fig-GHI-sim/{place}/fig-GHI-sim.pdf`
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
print_script_info(file.path("scripts", "figs", "fig-GHI-sim.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Bologna"
# Save the plot 
save_plot <- TRUE
# ***************************** Fixed Arguments ********************************
# Reference year
nyear <- 2014
# Figures base name
common_name <- "fig-GHI-sim"
# Number of simulations
nsim <- 1
# Horizon time 
t_hor <- 365*10
# Time step 
dt <- 0.05
# Random seed
seed <- 1
# ******************************************************************************
#                             Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly=TRUE)
if (!purrr::is_empty(args)) {
  # Reference location
  place <- ifelse(is.na(args[1]), place, args[1])
  # Save figure
  save_plot <- ifelse(is.na(args[2]), save_plot, args[2])
  print_script_args(place = place, save_plot = save_plot)
}
# ******************************************************************************
#                                  Load data 
# ******************************************************************************
# Load radiation models
load_data(file.path(outputs$dir$data$models$radiation$P, place), "radiation_models_CTMC")
# Reference model
model_Rt <- radiation_models_CTMC[[as.character(nyear)]]$clone(TRUE)
# ******************************************************************************
#                                Generate data 
# ******************************************************************************
# Today date 
t_now <- as.Date(paste0(nyear-1, "-12-31"))
# Simulation in continuous time 
df_sim_CT_intraday <- scenarios_radiationModel_CT(model_Rt, t_now, t_hor = t_hor, nsim = nsim, dt = dt, seed = seed)
df_sim_CT_day <- scenarios_radiationModel_CT(model_Rt, t_now, t_hor = t_hor, nsim = nsim, dt = 1, seed = seed)
# ******************************************************************************
#                                Generate figure
# ******************************************************************************
plot_sim_intraday <- filter(df_sim_CT_intraday, !is.na(date)) %>%
  ggplot()+
  geom_point(aes(date, Rt, group = sim), alpha = 0.7, size = 0.8)+
  geom_line(aes(date, Ct, group = sim), alpha = 1, color = "blue")+
  facet_wrap(~paste0("Simulated (dt = ", dt, ")"))+
  theme_bw()+
  figure_theme+
  theme(strip.text = element_text(angle = 0, face = "bold", size = 15), 
        axis.text.x = element_blank())+
  labs(x = NULL, y = NULL)

plot_sim_day <- ggplot()+
  geom_point(data = df_sim_CT_day, aes(date, Rt, group = sim), alpha = 1, size = 0.8)+
  geom_line(data = df_sim_CT_day, aes(date, Ct, group = sim), alpha = 1, color = "blue")+
  facet_wrap(~paste0("Simulated (dt = ", 1, ")"))+
  theme_bw()+
  figure_theme+
  theme(strip.text = element_text(angle = 0, face = "bold", size = 15),
        axis.text.x = element_blank())+
  labs(x = NULL, y = NULL)

plot_emp_day <- df_sim_CT_day %>%
  ggplot()+
  geom_point(aes(date, GHI), color = "black", size = 0.8)+
  geom_line(aes(date, Ct), alpha = 1, color = "blue")+
  facet_wrap(~"Realized (daily)")+
  theme_bw()+
  figure_theme+
  theme(strip.text = element_text(angle = 0, face = "bold", size = 15))+
  labs(x = NULL, y = NULL)

yleft = gridtext::richtext_grob("GHI (kWh/m²)", rot = 90, gp = grid::gpar(fontsize = 25))
bottom = gridtext::richtext_grob(text = '', gp = grid::gpar(fontsize = 15))
fig <- gridExtra::grid.arrange(plot_sim_intraday, plot_sim_day, plot_emp_day, 
                               nrow = 3, ncol = 1, left = yleft, bottom = bottom)
# ******************************************************************************
#                                Save figure 
# ******************************************************************************
if (save_plot){
  # Initialize a folder to store the figures
  dir_output <- file.path(outputs$dir$figs$models$radiation, common_name)
  make_new_directory(dir_output)
  # Initialize a folder to store the figures
  dir_output <- file.path(dir_output, place)
  make_new_directory(dir_output)
  control <- outputs$dir$figs$control
  # Save figure 
  save_new_fig(dir_output, fig, fig.name = common_name, file.format = stringr::str_remove_all(control$format, "\\."), 
               quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
  # Save the paper figure in main 
  if (place == "Bologna"){
    save_new_fig(outputs$dir$figs$main, fig, fig.name = common_name, file.format = stringr::str_remove_all(control$format, "\\."), 
                 quiet = FALSE, dpi = control$dpi, width = control$width, height = control$height)
  }
}
