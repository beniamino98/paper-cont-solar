# ---
#' @description
#' Generate location-level empirical application tables.
#'
#' @section `main`
#' @label `tbl-avg-moments-sorad`, `tbl-hedging-sorad`, `tbl-hedging-soradidx`,
#'        `tbl-avg-moments-sored`, `tbl-hedging-sored-uh`, `tbl-hedging-sored-h`, `tbl-hedging-soredidx-strip`
#' @name `tbl_avg_moments_sorad`, `tbl_hedging_sorad`, `tbl_hedging_soradidx`,
#'        `tbl_avg_moments_sored`, `tbl_hedging_sored_uh`, `tbl_hedging_sored_h`, `tbl_hedging_soredidx_strip`
#'
#' @author Beniamino Sartini
#' @created 2025-12-25
#' @modified 2026-06-03
#' @dependencies solarr
#'
#' @arguments
#'   - param[1] (place): reference location, e.g. "Bologna".
#'   - param[2] (save_data): write generated tables ("TRUE" or "FALSE").
#'
#' @example
#' Rscript scripts/tables/tbl-mv-place.R "Bologna" "TRUE"
#' for place in Bologna Palermo Roma; do Rscript scripts/tables/tbl-mv-place.R "$place" "TRUE"; done
#'
#' @inputs
#'   - `outputs.RData`
#'   - `data/solarOptions/mv/{place}/data_mv.RData`
#'
#' @outputs
#'   - `outputs$table[["tbl_avg_moments_sorad"]][[place]]`
#'   - `outputs$table[["tbl_hedging_sorad"]][[place]]`
#'   - `outputs$table[["tbl_hedging_soradidx"]][[place]]`
#'   - `outputs$table[["tbl_avg_moments_sored"]][[place]]`
#'   - `outputs$table[["tbl_hedging_sored_uh"]][[place]]`
#'   - `outputs$table[["tbl_hedging_sored_h"]][[place]]`
#'   - `outputs$table[["tbl_hedging_soredidx_strip"]][[place]]`
#'
#' @depends
#'   - `scripts/data/s8-solarOptions-mv-place.R`
#'   - `scripts/tables/zzz.R`
#'
#' @tags
#'   - tables
#'   - main
# ---
suppressMessages(source(file.path("scripts", "s0-load.R")))
load("outputs.RData")
source(file.path("scripts", "tables", "zzz.R"))
print_script_info(file.path("scripts", "tables", "tbl-mv-place.R"))
# ******************************************************************************
#                                  Inputs 
# ******************************************************************************
# Reference location
place <- "Bologna"
# Save the output
save_data <- TRUE
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
# Load moments strip
load_data(file.path(outputs$dir$data$main, "solarOptions", "mv", place), "data_mv")
# ******************************************************************************
#                      Table: tbl-avg-moments-sorad
# ******************************************************************************
# Table name
tbl_name <- "tbl_avg_moments_sorad"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# Moments 
df_mom <- data_mv$sorad$moments %>%
  dplyr::select(Year, M_R, M_R_emp, v_R, v_R_emp, M_Gamma, M_Gamma_emp, v_Gamma, v_Gamma_emp, cr_R_Gamma, cr_R_Gamma_emp)
# Format
data_kable <- bind_rows(df_mom, bind_cols(Year = "2014-2023", summarise_all(df_mom[,-1], mean))) %>%
  dplyr::mutate(cr_R_Gamma = format_perc(cr_R_Gamma),
                cr_R_Gamma_emp = format_perc(cr_R_Gamma_emp)) %>%
  dplyr::mutate_if(is.numeric, round, digits = 3)%>%
  dplyr::mutate_if(is.numeric, as.character)
data_kable 
colnames(data_kable) <- c("Year", 
                          "$\\bar{M}_{\\text{R}}$", "$\\bar{M}_{\\text{R}}^{\\text{emp}}$", 
                          "$\\bar{S}^2_{\\text{R}}$", "$\\bar{S}_{\\text{R}}^{2\\text{emp}}$", 
                          "$\\bar{M}_{\\Gamma}$", "$\\bar{M}_{\\Gamma}^{\\text{emp}}$",
                          "$\\bar{S}^2_{\\Gamma}$", "$\\bar{S}^{2\\text{emp}}_{\\Gamma}$", 
                          "$\\bar{\\mathbb{C}r}_{\\text{R}, \\Gamma}$", "$\\bar{\\mathbb{C}r}^{\\text{emp}}_{\\text{R}, \\Gamma}$")
# Store R object 
outputs$table[[tbl_name]][[place]]  <- data_kable
# ******************************************************************************
# Generate LaTex table 
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(9, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 7)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
# ******************************************************************************
#                          Table: tbl-hedging-sorad
# ******************************************************************************
# Table name
tbl_name <- "tbl_hedging_sorad"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# Data 
df_mom <- data_mv$sorad$hedging
# Kable 
data_kable <- bind_rows(df_mom, bind_cols(Year = "2014-2023", summarise_all(df_mom[,-1], mean)))%>%
  mutate(var_reduction = format_perc(var_reduction),
         premium = format_perc(premium, 2),
         ES_reduction = format_perc(ES_reduction))  %>%
  mutate(
    v_pi_buyer_uh = round(v_pi_buyer_uh, digits = 2)%>% as.character(),
    v_pi_buyer_h = round(v_pi_buyer_h, digits = 2)%>% as.character(),
    ES_buyer_uh = round(ES_buyer_uh, digits = 2)%>% as.character(),
    ES_buyer_h = round(ES_buyer_h, digits = 2)%>% as.character(),
    Pi_seller = round(Pi_seller, digits = 2) %>% as.character(),
    v_pi_buyer_h = paste0(v_pi_buyer_h, "(", var_reduction, ")"),
    ES_buyer_h = paste0(ES_buyer_h, " (", ES_reduction, ")")
  ) %>%
  mutate_if(is.numeric, round, digits = 0)%>%
  mutate_if(is.numeric, as.character) %>%
  select(-ES_reduction, -var_reduction) %>%
  select(Year, Pt, Qt, Gamma, premium, 
         Pi_seller, Pi_buyer_uh, Pi_buyer_h, 
         v_pi_buyer_uh, v_pi_buyer_h, ES_buyer_uh, ES_buyer_h)
data_kable
colnames(data_kable) <- c("Year", 
                          "$\\sum V_t$", "$\\sum q_t$", "$\\sum \\Gamma$", "Premium",
                          "$\\Pi^s_{t \\to T}$", "$\\Pi^{b,u}_{t \\to T}$", "$\\Pi^{b,h}_{t \\to T}$",
                          "$\\mathbb{V}\\{\\pi^{u}\\}$", "$\\mathbb{V}\\{\\pi^{h}\\}$", 
                          "$\\mathbb{EL}\\{\\pi^{u}\\}$", "$\\mathbb{EL}\\{\\pi^{h}\\}$")
# Store R object 
outputs$table[[tbl_name]][[place]]  <- data_kable
# ******************************************************************************
# Generate LaTex table 
outputs$tex[[tbl_name]][[place]] <- data_kable  %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(8, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(10, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 7)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
# ******************************************************************************
#                        Table: tbl-hedging-soradidx
# ******************************************************************************
# Table name
tbl_name <- "tbl_hedging_soradidx"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# Data 
df_mom <- data_mv$soradidx$hedging
# Kable
data_kable <- bind_rows(mutate(df_mom, Year = as.character(Year)), bind_cols(Year = "2014-2023", summarise_all(df_mom[,-1], mean)))%>%
  mutate(var_reduction = format_perc(var_reduction),
         premium = format_perc(premium, 2),
         ES_reduction = format_perc(ES_reduction))  %>%
  mutate(
    Qt = round(Qt, digits = 2) %>% as.character(),
    v_pi_buyer_uh = round(v_pi_buyer_uh, digits = 2)%>% as.character(),
    v_pi_buyer_h = round(v_pi_buyer_h, digits = 2)%>% as.character(),
    ES_buyer_uh = round(ES_buyer_uh, digits = 2)%>% as.character(),
    ES_buyer_h = round(ES_buyer_h, digits = 2)%>% as.character(),
    Pi_seller = round(Pi_seller, digits = 2) %>% as.character(),
    v_pi_buyer_h = paste0(v_pi_buyer_h, "(", var_reduction, ")"),
    ES_buyer_h = paste0(ES_buyer_h, " (", ES_reduction, ")")
  ) %>%
  mutate_if(is.numeric, round, digits = 0) %>%
  select(-ES_reduction, -var_reduction) %>%
  mutate_if(is.numeric, as.character)%>%
  select(Year, Pt, Qt, Gamma, premium, 
         Pi_seller, Pi_buyer_uh, Pi_buyer_h, 
         v_pi_buyer_uh, v_pi_buyer_h, ES_buyer_uh, ES_buyer_h)
data_kable
colnames(data_kable) <- c("Year", 
                          "$V_t$", "$q_t$", "$\\Gamma$", "Premium",
                          "$\\Pi^s_{t \\to T}$", "$\\Pi^{b,u}_{t \\to T}$", "$\\Pi^{b,h}_{t \\to T}$",
                          "$\\mathbb{V}\\{\\pi^{u}\\}$", "$\\mathbb{V}\\{\\pi^{h}\\}$",
                          "$\\mathbb{EL}\\{\\pi^{u}\\}$", "$\\mathbb{EL}\\{\\pi^{h}\\}$")
# Store R object 
outputs$table[[tbl_name]][[place]]  <- data_kable
# ******************************************************************************
# Generate LaTex table 
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(8, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(10, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 6)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
# ******************************************************************************
#                        Table: tbl-avg-moments-sored
# ******************************************************************************
# Table name
tbl_name <- "tbl_avg_moments_sored"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# Moments 
df_mom <- data_mv$sored$moments  %>%
  select(Year, M_ER, M_ER_emp, v_ER, v_ER_emp, M_E_Gamma, M_E_Gamma_emp, v_E_Gamma, v_E_Gamma_emp, cr_ER_EGamma, cr_ER_EGamma_emp)

data_kable <- bind_rows(df_mom, bind_cols(Year = "2014-2023", summarise_all(df_mom[,-1], mean)))%>%
  mutate(cr_ER_EGamma = paste0(round(cr_ER_EGamma*100, digits = 1), "\\%"),
         cr_ER_EGamma_emp = paste0(round(cr_ER_EGamma_emp*100, digits = 1), "\\%")) %>%
  mutate_if(is.numeric, format, digits = 2, scientific = FALSE)
data_kable 
colnames(data_kable) <- c("Year", 
                          "$\\bar{M}_{\\text{ER}}$", "$\\bar{M}_{\\text{ER}}^{\\text{emp}}$", 
                          "$\\bar{S}^2_{\\text{ER}}$", "$\\bar{S}_{\\text{ER}}^{2\\text{emp}}$", 
                          "$\\bar{M}_{\\Gamma^{\\text{E}}}$", "$\\bar{M}_{\\Gamma^{\\text{E}}}^{\\text{emp}}$",
                          "$\\bar{S}^2_{\\Gamma^{\\text{E}}}$", "$\\bar{S}^{2\\text{emp}}_{\\Gamma^{\\text{E}}}$", 
                          "$\\bar{\\mathbb{C}r}_{\\text{ER}, \\Gamma^{\\text{E}}}$", "$\\bar{\\mathbb{C}r}^{\\text{emp}}_{\\text{ER}, \\Gamma^{\\text{E}}}$")
# Store R object 
outputs$table[[tbl_name]][[place]]  <- data_kable
# ******************************************************************************
# Generate LaTex table 
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(9, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 6)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
# ******************************************************************************
#                      Table: tbl-hedging-sored-uh
# ******************************************************************************
# Table name
tbl_name <- "tbl_hedging_sored_uh"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# No hedging
df_mom <- data_mv$sored$hedging %>%
  select(-Pi_seller_h, -v_pi_seller_h, -var_reduction_seller, -phi_buyer_h, -phi_buyer_uh, -phi_seller)
data_kable <- bind_rows(df_mom, bind_cols(Year = "2014-2023", summarise_all(df_mom[,-1], mean)))%>%
  mutate(var_reduction_buyer = format_perc(var_reduction_buyer),
         ES_reduction = format_perc(ES_reduction),
         premium = format_perc(premium, 2),
         Qt = round(Qt, digits = 2) %>% as.character(),
         v_pi_buyer_uh = format(v_pi_buyer_uh*100, digits = 2)%>% as.character(),
         v_pi_buyer_h = format(v_pi_buyer_h*100, digits = 2)%>% as.character(),
         ES_buyer_uh = format(ES_buyer_uh*100, digits = 2)%>% as.character(),
         ES_buyer_h = format(ES_buyer_h*100, digits = 2)%>% as.character(),
         Pi_seller_uh = format(Pi_seller_uh, digits = 2) %>% as.character(),
         v_pi_buyer_h = paste0(v_pi_buyer_h, "(", var_reduction_buyer, ")"),
         ES_buyer_h = paste0(ES_buyer_h, " (", ES_reduction, ")"))  %>%
  mutate_if(is.numeric, format, digits = 2, scientific = FALSE) %>%
  select(-ES_reduction, - var_reduction_buyer, -v_pi_seller_uh)%>%
  select(Year, Pt, Qt, Gamma, premium, 
         Pi_seller_uh, Pi_buyer_uh, Pi_buyer_h, 
         v_pi_buyer_uh, v_pi_buyer_h, ES_buyer_uh, ES_buyer_h)
colnames(data_kable) <- c("Year", 
                          "$\\sum V_t$", "$\\sum q_t$", "$\\sum \\Gamma$",  "Premium",
                          "$\\Pi_{t\\to T}^s$", "$\\Pi_{t\\to T}^{b,u}$", "$\\Pi_{t\\to T}^{b,h}$",
                          "$\\mathbb{V}\\{\\pi^{u}\\}$", "$\\mathbb{V}\\{\\pi^{h}\\}$", 
                          "$\\mathbb{EL}\\{\\pi^{u}\\}$", "$\\mathbb{EL}\\{\\pi^{h}\\}$")
# Store R object 
outputs$table[[tbl_name]][[place]]  <- data_kable
# ******************************************************************************
# Generate LaTex table 
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(8, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(11, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 6)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
# ******************************************************************************
#                       Table: tbl-hedging-sored-h
# ******************************************************************************
# Table name
tbl_name <- "tbl_hedging_sored_h"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
# Hedging
df_mom <- data_mv$sored$hedging_strip %>%
  select(-Pi_seller_uh, -v_pi_seller_uh, -v_pi_seller_h, -var_reduction_seller)
data_kable <- bind_rows(mutate(df_mom, Year = as.character(Year)), bind_cols(Year = "2014-2023", summarise_all(df_mom[,-1], mean)))%>%
  mutate(premium = format_perc(premium, 2),
         ES_reduction = format_perc(ES_reduction))  %>%
  mutate(
    Qt = round(Qt, digits = 0) %>% as.character(),
    phi_buyer_h = round(phi_buyer_h, 0),
    phi_buyer_uh = round(phi_buyer_uh, 0),
    phi_seller = round(phi_seller, 0),
    v_pi_buyer_uh = round(v_pi_buyer_uh*100, digits = 2)%>% as.character(),
    v_pi_buyer_h = round(v_pi_buyer_h*100, digits = 2)%>% as.character(),
    ES_buyer_uh = round(ES_buyer_uh*100, digits = 2)%>% as.character(),
    ES_buyer_h = round(ES_buyer_h*100, digits = 2)%>% as.character(),
    Pi_seller_h = round(Pi_seller_h, digits = 2) %>% as.character(),
    Pi_buyer_uh = round(Pi_buyer_uh, digits = 2) %>% as.character(),
    Pi_buyer_h = round(Pi_buyer_h, 0),
    v_pi_buyer_h = paste0(v_pi_buyer_h, "(", format_perc(var_reduction_buyer), ")"),
    ES_buyer_h = paste0(ES_buyer_h, " (", ES_reduction, ")")
  ) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  select(-ES_reduction, -var_reduction_buyer) %>%
  mutate_if(is.numeric, as.character) %>%
  select(Year, Pt, Qt, Gamma, premium, 
         phi_seller, phi_buyer_uh, phi_buyer_h, 
         Pi_seller_h, Pi_buyer_uh, Pi_buyer_h, 
         v_pi_buyer_uh, v_pi_buyer_h, ES_buyer_uh, ES_buyer_h)

colnames(data_kable) <- c("Year", 
                          "$\\sum \\widetilde{V}_{t}$", "$\\sum \\tilde{q}_{t}$", "$\\sum \\Gamma$", "Premium",
                          "$\\sum \\phi_{t}^{s}$", "$\\sum \\phi_{t}^{b,u}$", "$\\sum \\phi_t^{b,h}$",
                          "$\\widetilde{\\Pi}_{t\\to T}^{s,h}$", "$\\widetilde{\\Pi}_{t \\to T}^{b,u}$", "$\\widetilde{\\Pi}_{t \\to T}^{b,h}$",
                          "$\\mathbb{V}\\{\\tilde{\\pi}^{u}\\}$",  "$\\mathbb{V}\\{\\tilde{\\pi}^{h}\\}$", 
                          "$\\mathbb{EL}\\{\\tilde{\\pi}^{u}\\}$", "$\\mathbb{EL}\\{\\tilde{\\pi}^{h}\\}$")
# Store R object 
outputs$table[[tbl_name]][[place]]  <- data_kable
# ******************************************************************************
# Generate LaTex table 
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(8, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(10, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(13, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 4.5)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")
# ******************************************************************************
#                       Table: tbl-hedging-soredidx-strip
# ******************************************************************************
# Table name
tbl_name <- "tbl_hedging_soredidx_strip"
cli::cli_alert_info(paste0("Generating table: ", "\033[1;35m", tbl_name, "\033[0m", " (", place, ") \n"))
# ******************************************************************************
df_mom <- data_mv$soredidx$hedging_strip %>%
  select(-Pi_seller_uh, -v_pi_seller_h, -v_pi_seller_uh, -var_reduction_seller)
data_kable <- bind_rows(mutate(df_mom, Year = as.character(Year)), bind_cols(Year = "2014-2023", summarise_all(df_mom[,-1], mean))) %>%
  mutate(var_reduction_buyer = format_perc(var_reduction_buyer),
         premium = format_perc(premium, 2),
         ES_reduction = format_perc(ES_reduction),
         Qt = round(Qt, digits = 2) %>% as.character(),
         v_pi_buyer_uh = round(v_pi_buyer_uh*100, digits = 2)%>% as.character(),
         v_pi_buyer_h = round(v_pi_buyer_h*100, digits = 2)%>% as.character(),
         ES_buyer_uh = round(ES_buyer_uh*100, digits = 2)%>% as.character(),
         ES_buyer_h = round(ES_buyer_h*100, digits = 2)%>% as.character(),
         Pi_seller_h = round(Pi_seller_h, digits = 2) %>% as.character(),
         v_pi_buyer_h = paste0(v_pi_buyer_h, "(", var_reduction_buyer, ")"),
         ES_buyer_h = paste0(ES_buyer_h, " (", ES_reduction, ")")) %>%
  mutate_if(is.numeric, format, digits = 2, scientific = FALSE) %>%
  select(-var_reduction_buyer, -ES_reduction) %>%
  select(Year, Pt, Qt, Gamma, premium, 
         phi_seller, phi_buyer_uh, phi_buyer_h, 
         Pi_seller_h, Pi_buyer_uh, Pi_buyer_h, 
         v_pi_buyer_uh, v_pi_buyer_h, ES_buyer_uh, ES_buyer_h)
colnames(data_kable) <- c("Year", 
                          "$\\widetilde{V}_{t}$", "$\\tilde{q}_{t}$", "$\\Gamma$", "Premium",
                          "$\\boldsymbol{\\phi}_{t\\to T}^{s}$", "$\\boldsymbol{\\phi}_{t\\to T}^{b,u}$", "$\\boldsymbol{\\phi}_{t \\to T}^{b,h}$",
                          "$\\widetilde{\\Pi}_{t\\to T}^{s,h}$", "$\\widetilde{\\Pi}_{t \\to T}^{b,u}$", "$\\widetilde{\\Pi}_{t \\to T}^{b,h}$",
                          "$\\mathbb{V}\\{\\tilde{\\pi}^{u}\\}$",  "$\\mathbb{V}\\{\\tilde{\\pi}^{h}\\}$", 
                          "$\\mathbb{EL}\\{\\tilde{\\pi}^{u}\\}$", "$\\mathbb{EL}\\{\\tilde{\\pi}^{h}\\}$")

# Store R object 
outputs$table[[tbl_name]][[place]]  <- data_kable
# ******************************************************************************
# Generate LaTex table 
outputs$tex[[tbl_name]][[place]] <- data_kable %>%
  dplyr::mutate_if(is.numeric, format, scientific = FALSE, digits = 3) %>%
  knitr::kable(booktabs = TRUE,
               longtable = FALSE,
               label = NA,
               format = "latex",
               linesep = "",
               caption = NA,
               escape = FALSE) %>% 
  kableExtra::column_spec(1, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(5, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(8, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(10, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::column_spec(13, border_left = FALSE, border_right = TRUE) %>%
  kableExtra::row_spec(10, hline_after = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE, font_size = 4.5)
outputs$tex[[tbl_name]][[place]][1] <- stringr::str_replace(outputs$tex[[tbl_name]][[place]][1], "\\[!h\\]", "")
# ******************************************************************************
save(outputs, file = "outputs.RData")

