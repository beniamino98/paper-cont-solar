# ---
#' @description
#' Generate mean-variance summaries and hedging data for one location.
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr
#'
#' @arguments
#'   - param[1] (place): reference location, e.g. "Bologna".
#'   - param[2] (save_data): write generated data ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/data/s8-solarOptions-mv-place.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/data/s8-solarOptions-mv-place.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/models/radiation/P/{place}/radiation_models_CTMC.RData`
#'   - `data/solarOptions/moments/{place}/moments_index.RData`
#'   - `data/solarOptions/moments/{place}/moments_strip.RData`
#'   - `data/scenarios/{place}/scenarios.RData`
#'
#' @outputs
#'   - `data/solarOptions/mv/{place}/data_mv.RData` (mean-variance summaries)
#'   - `data/solarOptions/mv/{place}/data_hedging.RData` (hedging summaries)
#'
#' @depends
#'   - `scripts/data/s2b-models-radiation-P-CTMC-place.R`
#'   - `scripts/data/s7-solarOptions-moments-place.R`
#'   - `scripts/data/s6-simulate-Rt-Et-place.R`
#'
#' @tags
#'   - data
#'   - main
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
print_script_info(file.path("scripts", "data", "s8-solarOptions-mv-place.R"))
#summarise_soredidx_hedging <- solarOption_hedging_soredidx
#summarise_sored_hedging <- solarOption_hedging_sored
#summarise_soradidx_hedging <- solarOption_hedging_soradidx
#summarise_sorad_hedging <- solarOption_hedging_sorad
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Roma"
# Save the output
save_data <- TRUE
# ***************************** Fixed Arguments ********************************
# Reference names
common_name <- c("data_mv", "data_hedging")
# Reference years
nyears <- as.character(outputs$nyears)
# Risk aversion for SoRadIDX 
nu_sorad <- c(nu_b = 0.01, nu_s = 0.01)
# Risk aversion for SoREdIDX 
nu_sored <- c(nu_b = 0.1, nu_s = 0.1) 
# ******************************************************************************
#                             Inputs (Command line)
# ******************************************************************************
# Supply arguments from command line
args <- commandArgs(trailingOnly=TRUE)
if (!purrr::is_empty(args)) {
  # Reference location
  place <- ifelse(is.na(args[1]), place, args[1])
  # Save output 
  save_data <- ifelse(is.na(args[2]), save_data, args[2])
  print_script_args(place = place, save_data = save_data)
}
# ******************************************************************************
#                                 Load data 
# ******************************************************************************
# Load radiation models
load_data(file.path(outputs$dir$data$models$radiation$P, place), "radiation_models_CTMC")
# Same strike
K_fun <- function(n) radiation_models_CTMC[[1]]$Rt_bar(n)
# Load moments strip
load_data(file.path(outputs$dir$data$main, "solarOptions", "moments", place), "moments_strip")
# Load moments index
load_data(file.path(outputs$dir$data$main, "solarOptions", "moments", place), "moments_index")
# Load scenarios 
load_data(file.path(outputs$dir$data$main, "scenarios", place), "scenarios")
# ******************************************************************************
#                         SoRad and SoRadIDX
# ******************************************************************************
nu_b <- nu_sorad[1]
nu_s <- nu_sorad[2]
sorad_moments <- purrr::map_df(moments_strip, ~.x$sorad_day) %>%
  mutate(Year = as.character(Year))
sorad_hedging  <-purrr::map_df(nyears, ~solarOption_hedging_sorad(moments_strip[[.x]]$sorad, scenarios[[.x]], nu_b, nu_s, 1, K_fun)$emp)%>%
  mutate(Year = as.character(Year))
soradidx_moments <- purrr::map_df(moments_index, ~.x$sorad) %>%
  mutate(Year = as.character(lubridate::year(t_hor))) %>%
  select(-t_hor, -t_now) %>%
  select(Year, everything()) 
soradidx_hedging <- purrr::map_df(nyears, ~solarOption_hedging_soradidx(moments_index[[.x]]$sorad, scenarios[[.x]], nu_b, nu_s, 1, K_fun)$emp)%>%
  mutate(Year = as.character(Year))
# ******************************************************************************
#                         SoREd and SoREdIDX
# ******************************************************************************
nu_b <- nu_sored[1]
nu_s <- nu_sored[2]
sored_moments <- purrr::map_df(moments_strip, ~.x$sored_day) %>%
  mutate(Year = as.character(Year)) 
sored_hedging <- purrr::map_df(nyears, ~solarOption_hedging_sored(moments_strip[[.x]]$sored, scenarios[[.x]], nu_b, nu_s, 1, K_fun, hedging = "none")$emp) %>%
  mutate(Year = as.character(Year))
sored_hedging_strip <- purrr::map_df(nyears, ~solarOption_hedging_sored(moments_strip[[.x]]$sored, scenarios[[.x]], nu_b, nu_s, 1, K_fun, hedging = "strip")$emp) %>%
  mutate(Year = as.character(Year))
soredidx_moments <- purrr::map_df(moments_index, ~.x$sored) %>%
  mutate(Year = as.character(lubridate::year(t_hor))) %>%
  select(-t_hor, -t_now) %>%
  select(Year, everything()) 
soredidx_hedging <- purrr::map_df(nyears, ~solarOption_hedging_soredidx(moments_index[[.x]]$sored, scenarios[[.x]], nu_b, nu_s, 1, K_fun, hedging = "none")$emp) %>%
  mutate(Year = as.character(Year))
soredidx_hedging_strip <- purrr::map_df(nyears, ~solarOption_hedging_soredidx(moments_index[[.x]]$sored, scenarios[[.x]], nu_b, nu_s, 1, K_fun, hedging = "strip")$emp) %>%
  mutate(Year = as.character(Year))
# Structure output
data_mv <- list(
  sorad = list(moments = sorad_moments, hedging = sorad_hedging),
  soradidx = list(moments = soradidx_moments, hedging = soradidx_hedging),
  sored = list(moments = sored_moments, hedging = sored_hedging, hedging_strip = sored_hedging_strip),
  soredidx = list(moments = soredidx_moments, hedging = soredidx_hedging, hedging_strip = soredidx_hedging_strip)
)
# ******************************************************************************
#                              Initialize output
# ******************************************************************************
# Output data
data_hedging <- list(
  sorad = NA,
  soradidx = NA,
  sored = NA,
  sored_strip = NA,
  soredidx = NA,
  soredidx_strip = NA,
  soredidx_tot = NA
)
# ******************************************************************************
#                            Generate Data: SoRad
# ******************************************************************************
nu_b <- nu_sorad[1]
nu_s <- nu_sorad[2]
df <- list()
for(nyear in nyears){
  nyear <- as.character(nyear)
  # Inputs
  sorad <- moments_strip[[nyear]]$sorad
  scenario <- scenarios[[nyear]]
  # Same strike
  K_fun <- function(n) radiation_models_CTMC[[1]]$Rt_bar(n)
  hedge <- solarOption_hedging_sorad(sorad, scenario, nu_b, nu_s, 1, K_fun)
  df[[nyear]] <- hedge$day
}
# Save data 
data_hedging$sorad <- bind_rows(df)%>%
  group_by(n) %>%
  filter(n != 366)%>%
  mutate(e_cum_ret_seller = mean(cum_ret_seller))
# ******************************************************************************
#                            Generate Data: SoRadIDX
# ******************************************************************************
df <- list()
for(nyear in nyears){
  nyear <- as.character(nyear)
  # Inputs
  sorad <- moments_index[[nyear]]$sorad
  scenario <- scenarios[[nyear]]
  # Same strike
  K_fun <- function(n) radiation_models_CTMC[[1]]$Rt_bar(n)
  hedge <- solarOption_hedging_soradidx(sorad, scenario, nu_b, nu_s, 1, K_fun)
  df[[nyear]] <- hedge$day
}
# Save data 
data_hedging$soradidx <- bind_rows(df)%>%
  group_by(n) %>%
  filter(n != 366)%>%
  mutate(
    e_cum_ret_seller = mean(cum_ret_seller))
# ******************************************************************************
#                            Generate Data: SoREd
# ******************************************************************************
nu_b <- nu_sored[1]
nu_s <- nu_sored[2]
df <- list()
for(nyear in nyears){
  nyear <- as.character(nyear)
  # Inputs
  sored <- moments_strip[[nyear]]$sored
  scenario <- scenarios[[nyear]]
  # Same strike
  K_fun <- function(n) radiation_models_CTMC[[1]]$Rt_bar(n)
  hedge <- solarOption_hedging_sored(sored, scenario, nu_b, nu_s, 1, K_fun)
  df[[nyear]] <- hedge$day
}
# Save data 
data_hedging$sored <- bind_rows(df)%>%
  group_by(n) %>%
  filter(n != 366)%>%
  mutate(
    e_cum_ret_seller_uh = mean(cum_ret_seller_uh),
    e_cum_ret_seller_h = mean(cum_ret_seller_h))
# ******************************************************************************
#                        Generate Data: SoREd (hedged)
# ******************************************************************************
df <- list()
for(nyear in nyears){
  nyear <- as.character(nyear)
  # Inputs
  sored <- moments_strip[[nyear]]$sored
  scenario <- scenarios[[nyear]]
  # Same strike
  K_fun <- function(n) radiation_models_CTMC[[1]]$Rt_bar(n)
  hedge <- solarOption_hedging_sored(sored, scenario, nu_b, nu_s, 1, K_fun, hedging = "strip")
  df[[nyear]] <- hedge$day
}
# Save data 
data_hedging$sored_strip <- bind_rows(df)%>%
  group_by(n) %>%
  filter(n != 366)%>%
  mutate(
    e_cum_ret_seller_uh = mean(cum_ret_seller_uh),
    e_cum_ret_seller_h = mean(cum_ret_seller_h))
# ******************************************************************************
#                            Generate Data: SoREdIDX
# ******************************************************************************
df <- list()
for(nyear in nyears){
  nyear <- as.character(nyear)
  # Inputs
  sored <- moments_index[[nyear]]$sored
  scenario <- scenarios[[nyear]]
  # Same strike
  K_fun <- function(n) radiation_models_CTMC[[1]]$Rt_bar(n)
  hedge <- solarOption_hedging_soredidx(sored, scenario, nu_b, nu_s, 1, K_fun, hedging = "none")
  df[[nyear]] <- hedge$day
}
# Save data 
data_hedging$soredidx <- bind_rows(df)%>%
  group_by(n) %>%
  filter(n != 366)%>%
  mutate(
    e_cum_ret_seller_uh = mean(cum_ret_seller_uh),
    e_cum_ret_seller_h = mean(cum_ret_seller_h))
# ******************************************************************************
#                   Generate Data: SoREdIDX (hedged, strip)
# ******************************************************************************
df <- list()
for(nyear in nyears){
  nyear <- as.character(nyear)
  # Inputs
  sored <- moments_index[[nyear]]$sored
  scenario <- scenarios[[nyear]]
  # Same strike
  K_fun <- function(n) radiation_models_CTMC[[1]]$Rt_bar(n)
  hedge <- solarOption_hedging_soredidx(sored, scenario, nu_b, nu_s, 1, K_fun, hedging = "strip")
  df[[nyear]] <- hedge$day
}
# Save data 
data_hedging$soredidx_strip <- bind_rows(df)%>%
  group_by(n) %>%
  filter(n != 366)%>%
  mutate(
    e_cum_ret_seller_uh = mean(cum_ret_seller_uh),
    e_cum_ret_seller_h = mean(cum_ret_seller_h))
# ******************************************************************************
#                             Save data 
# ******************************************************************************
if (save_data) {
  # Initialize output directory 
  dir_output <- file.path(outputs$dir$data$main, "solarOptions")
  make_new_directory(dir_output)
  # Initialize output directory 
  dir_output <- file.path(dir_output, "mv")
  make_new_directory(dir_output)
  # Initialize output directory 
  dir_output <- file.path(dir_output, place)
  make_new_directory(dir_output)
  # **********************************************************************
  # Save outputs
  save_new_file(dir_output, file.name = common_name[1],
                file.format = "RData", quiet = FALSE, data_mv)
  # Save outputs
  save_new_file(dir_output, file.name = common_name[2],
                file.format = "RData", quiet = FALSE, data_hedging)
}
