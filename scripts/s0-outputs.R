# ---
#' @description
#' Initialize the project output registry and directory structure.
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies base R
#'
#' @example
#' Rscript scripts/s0-outputs.R
#'
#' @inputs
#'   - none
#'
#' @outputs
#'   - `outputs.RData`
#'   - project data, table, figure, animation, and TeX output directories
#'
#' @depends
#'   - none
#'
#' @tags
#'   - pipeline
#'   - setup
# ---
library(solarr)
# Initialize an outputs.RData file to store the outputs
outputs <- list()
outputs$dir <- list()
# Benchmark years
outputs$nyears <- c(2013:2022)
# Directory for GME data
outputs$dir_GME <- "data/GME_daily.csv"
outputs$places <- setNames( c("Bologna", "Roma", "Palermo"),  c("Bologna", "Rome", "Palermo"))
# ======================================================================================
#                                 Directory: data
# ======================================================================================
# Main: directory data
outputs$dir$data <- list()
outputs$dir$data$main <- "data"
system(paste0("mkdir ", outputs$dir$data$main))
# ======================================================================================
#                                 Directory: data/models
# ======================================================================================
# data
#     models
## Subdir: data/models
outputs$dir$data$models <- list()
outputs$dir$data$models$main <- "data/models"
system(paste0("mkdir ", outputs$dir$data$models$main))
### Sub-sub-dir: data/models/rho
outputs$dir$data$models$rho <- file.path(outputs$dir$data$models$main, "rho")
system(paste0("mkdir ", outputs$dir$data$models$rho))
# ======================================================================================
#                                 Directory: data/models/electricity
# ======================================================================================
# data
#     models
#          electricity
#                     P
#                     Q
### Sub-sub-dir: data/models/electricity
outputs$dir$data$models$electricity <- list()
outputs$dir$data$models$electricity$main <- file.path(outputs$dir$data$models$main, "electricity")
system(paste0("mkdir ", outputs$dir$data$models$electricity$main))
#### Sub-sub-sub-dir: data/models/electricity/P
outputs$dir$data$models$electricity$P <- file.path(outputs$dir$data$models$electricity$main, "P")
system(paste0("mkdir ", outputs$dir$data$models$electricity$P))
#### Sub-sub-sub-dir: data/models/electricity/Q
outputs$dir$data$models$electricity$Q <- file.path(outputs$dir$data$models$electricity$main, "Q")
system(paste0("mkdir ", outputs$dir$data$models$electricity$Q))
# ======================================================================================
#                                 Directory: data/models/radiation
# ======================================================================================
# data
#     models
#          radiation
#                   P
#                   Q
#### Sub-sub-dir: data/models/radiation
outputs$dir$data$models$radiation <- list()
outputs$dir$data$models$radiation$main <- file.path(outputs$dir$data$models$main, "radiation")
system(paste0("mkdir ", outputs$dir$data$models$radiation$main))
#### Sub-sub-sub-dir: data/models/radiation/P
outputs$dir$data$models$radiation$P <- file.path(outputs$dir$data$models$radiation$main, "P")
#### Sub-sub-sub-dir: data/models/radiation/Q
outputs$dir$data$models$radiation$Q <- file.path(outputs$dir$data$models$radiation$main, "Q")
# ======================================================================================
#                                 Directory: data/SoRad
# ======================================================================================
# data
#     SoRad
#          yearly
#          monthly
## Subdir: data/SoRad
outputs$dir$data$sorad <- list()
outputs$dir$data$sorad$main <- file.path(outputs$dir$data$main, "SoRad")
system(paste0("mkdir ", outputs$dir$data$sorad$main))
### Sub-sub-dir: data/SoRad/yearly
outputs$dir$data$sorad$yearly <- file.path(outputs$dir$data$sorad$main, "yearly")
system(paste0("mkdir ", outputs$dir$data$sorad$yearly))
### Sub-sub-dir: data/SoRad/monthly
outputs$dir$data$sorad$monthly <- file.path(outputs$dir$data$sorad$main, "monthly")
system(paste0("mkdir ", outputs$dir$data$sorad$monthly))
# ======================================================================================
#                                 Directory: data/SoREd
# ======================================================================================
# data
#     SoREd
#          yearly
#          monthly
## Subdir: data/SoREd
outputs$dir$data$sored <- list()
outputs$dir$data$sored$main <- file.path(outputs$dir$data$main, "SoREd")
system(paste0("mkdir ", outputs$dir$data$sored$main))
### Sub-sub-dir: data/SoREd/yearly
outputs$dir$data$sored$yearly <- file.path(outputs$dir$data$sored$main, "yearly")
system(paste0("mkdir ", outputs$dir$data$sored$yearly))
### Sub-sub-dir: data/SoREd/monthly
outputs$dir$data$sored$monthly <- file.path(outputs$dir$data$sored$main, "monthly")
system(paste0("mkdir ", outputs$dir$data$sored$monthly))
# ======================================================================================
#                                 Directory: figs
# ======================================================================================
# Main: directory figs
outputs$dir$figs <- list()
outputs$dir$figs$main <- "figs"
system(paste0("mkdir ", outputs$dir$figs$main))
# Parameters for figures outputs
outputs$dir$figs$control <- list(dpi = 600, width = 15, height = 10, format = ".pdf")
outputs$dir$figs$models <- list()
outputs$dir$figs$models$main <- file.path(outputs$dir$figs$main, "models")
system(paste0("mkdir ", outputs$dir$figs$models$main))
outputs$dir$figs$models$electricity <- file.path(outputs$dir$figs$models$main, "electricity")
system(paste0("mkdir ", outputs$dir$figs$models$electricity))
outputs$dir$figs$models$radiation <- file.path(outputs$dir$figs$models$main, "radiation")
system(paste0("mkdir ", outputs$dir$figs$models$radiation))
outputs$dir$figs$models$rho <- file.path(outputs$dir$figs$models$main, "rho")
system(paste0("mkdir ", outputs$dir$figs$models$rho))
# ======================================================================================
#                                 Directory: figs/SoREd
# ======================================================================================
# figs
#     SoREd
#          yearly
#          monthly
## Sub-dir: figs/SoREd
outputs$dir$figs$sored <- list()
outputs$dir$figs$sored$main <- file.path(outputs$dir$figs$main, "SoREd")
system(paste0("mkdir ", outputs$dir$figs$sored$main))
### Sub-sub-dir: figs/SoREd/yearly
outputs$dir$figs$sored$yearly <- file.path(outputs$dir$figs$sored$main, "yearly")
system(paste0("mkdir ", outputs$dir$figs$sored$yearly))
#### Sub-sub-sub-dir: figs/SoREd/monthly
outputs$dir$figs$sored$monthly <- file.path(outputs$dir$figs$sored$main, "monthly")
system(paste0("mkdir ", outputs$dir$figs$sored$monthly))
# ======================================================================================
#                                 Directory: figs/SoRad
# ======================================================================================
# figs
#     SoRad
#          yearly
#          monthly
## Sub-dir: figs/SoRad
outputs$dir$figs$sorad <- list()
outputs$dir$figs$sorad$main <- file.path(outputs$dir$figs$main, "SoRad")
system(paste0("mkdir ", outputs$dir$figs$sorad$main))
### Sub-sub-dir: figs/SoRad/yearly
outputs$dir$figs$sorad$yearly <- file.path(outputs$dir$figs$sorad$main, "yearly")
system(paste0("mkdir ", outputs$dir$figs$sorad$yearly))
#### Sub-sub-sub-dir: figs/SoRad/monthly
outputs$dir$figs$sorad$monthly <- file.path(outputs$dir$figs$sorad$main, "monthly")
system(paste0("mkdir ", outputs$dir$figs$sorad$monthly))
# ======================================================================================
#                                 Directory: figs/animations
# ======================================================================================
#### Sub-dir: figs/animations
outputs$dir$figs$animations <- list()
outputs$dir$figs$animations$main <- file.path(outputs$dir$figs$main, "animations")
system(paste0("mkdir ", outputs$dir$figs$animations$main))
# ======================================================================================
#                                 Directory: figs/animations/SoREd
# ======================================================================================
#### Sub-sub-dir: figs/animations/SoREd
outputs$dir$figs$animations$sored <- list()
outputs$dir$figs$animations$sored$main <- file.path(outputs$dir$figs$animations$main, "SoREd")
system(paste0("mkdir ", outputs$dir$figs$animations$sored$main))
#### Sub-sub-sub-dir: figs/animations/SoREd/yearly
outputs$dir$figs$animations$sored$yearly <- file.path(outputs$dir$figs$animations$sored$main, "yearly")
system(paste0("mkdir ", outputs$dir$figs$animations$sored$yearly))
#### Sub-sub-sub-dir: figs/animations/SoREd/monthly
outputs$dir$figs$animations$sored$monthly <- file.path(outputs$dir$figs$animations$sored$main, "monthly")
system(paste0("mkdir ", outputs$dir$figs$animations$sored$monthly))
# ======================================================================================
#                                 Directory: figs/animations/SoRad
# ======================================================================================
#### Sub-sub-dir: figs/animations/SoRad
outputs$dir$figs$animations$sorad <- list()
outputs$dir$figs$animations$sorad$main <- file.path(outputs$dir$figs$animations$main, "SoRad")
system(paste0("mkdir ", outputs$dir$figs$animations$sorad$main))
#### Sub-sub-sub-dir: figs/animations/SoREd/yearly
outputs$dir$figs$animations$sorad$yearly <- file.path(outputs$dir$figs$animations$sorad$main, "yearly")
system(paste0("mkdir ", outputs$dir$figs$animations$sorad$yearly))
#### Sub-sub-sub-dir: figs/animations/SoREd/monthly
outputs$dir$figs$animations$sorad$monthly <- file.path(outputs$dir$figs$animations$sorad$main, "monthly")
system(paste0("mkdir ", outputs$dir$figs$animations$sorad$monthly))
# ======================================================================================
#                                 Save the output file
# ======================================================================================
outputs$table <- list()
outputs$tex <- list()
# Hedging parameters 
outputs$control_hedging <- control_solarHedging(n_panels = 10000, efficiency = 0.15, PUN = 0.08, tick = 0.1, n_contracts = 1)
# Save the output file
save(outputs, file = "outputs.RData")
