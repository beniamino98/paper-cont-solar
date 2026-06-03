fig_mv_demand_supply_sorad <- function(supply_demand_sorad, nu_b, nu_s, r_min_max = c(0.97, 1.03), digits = 3){
  
  # 1) Equilibrium price 
  P_eq <- supply_demand_sorad$price(nu_s, nu_b)
  # 2) Equilibrium quantities
  Q_eq <- supply_demand_sorad$supply(P_eq, nu_s)
  # ************************************
  P0 <- P_eq * seq(r_min_max[1], r_min_max[2], 0.001)
  # ************************************
  # Demand around the eq. price
  w_b <- supply_demand_sorad$demand(P0, nu_b)
  # Supply around the eq. price
  w_s <- supply_demand_sorad$supply(P0, nu_s)
  # Scales 
  x_breaks <- c(min(P0), P_eq, max(P0))
  x_labels <- round(x_breaks, digits = digits)
  y_breaks <- c(min(c(w_b, w_s)), Q_eq,  max(c(w_b, w_s)))
  y_labels <- round(y_breaks, digits = digits)
  
  ggplot()+
    geom_line(aes(P0, w_b, color = "demand"))+
    geom_line(aes(P0, w_s, color = "supply"))+
    geom_segment(aes(x = P_eq, xend = P_eq, y = min(c(w_s, w_b)), yend = Q_eq), linetype = "dashed")+
    geom_segment(aes(x = min(P0), xend = P_eq, y = Q_eq, yend = Q_eq), linetype = "dashed")+
    geom_point(aes(x = P_eq, y = Q_eq, color = "eq"), shape = 10, size = 4)+
    scale_x_continuous(breaks = x_breaks, labels = x_labels)+
    scale_y_continuous(breaks = y_breaks, labels = y_labels)+
    scale_color_manual(values = c(demand = "green", supply = "red", eq = "black"), 
                       labels = c(demand = latex2exp::TeX("$q_t^{b}$"), 
                                  supply = latex2exp::TeX("$q_t^{s}$"), 
                                  eq = latex2exp::TeX("$V_t$")))+
    theme_bw()+
    theme(legend.position = "top")+
    labs(color = NULL, x = "Price", y = "Quantity")+
    figure_theme
}

fig_demand_supply_sorad <- function(supply_demand_sorad, nu_b, nu_s, r_min_max = c(0.97, 1.03)){
  
  nu_b_seq <- nu_b * c(1, 2, 1, 3)
  nu_s_seq <- nu_s * c(1, 1, 2, 0.75)
  # Subtitles 
  sub_1 <- latex2exp::TeX(paste0("$\\nu_b = ", round(nu_b_seq[1], 3), ",\\; \\nu_s = ", round(nu_s_seq[1], 3), "$"))
  sub_2 <- latex2exp::TeX(paste0("$\\nu_b = ", round(nu_b_seq[2], 3), ",\\;  \\nu_s = ", round(nu_s_seq[2], 3), "$"))
  sub_3 <- latex2exp::TeX(paste0("$\\nu_b = ", round(nu_b_seq[3], 3), ",\\;  \\nu_s = ", round(nu_s_seq[3], 3), "$"))
  sub_4 <- latex2exp::TeX(paste0("$\\nu_b = ", round(nu_b_seq[4], 3), ",\\;  \\nu_s = ", round(nu_s_seq[4], 3), "$"))
  # Figures
  fig_1 <- fig_mv_demand_supply_sorad(supply_demand_sorad, nu_b_seq[1], nu_s_seq[1], r_min_max) + labs(x = NULL, y = NULL, subtitle = sub_1)
  fig_2 <- fig_mv_demand_supply_sorad(supply_demand_sorad, nu_b_seq[2], nu_s_seq[2], r_min_max) + labs(x = NULL, y = NULL, subtitle = sub_2)
  fig_3 <- fig_mv_demand_supply_sorad(supply_demand_sorad, nu_b_seq[3], nu_s_seq[3], r_min_max) + labs(x = NULL, y = NULL, subtitle = sub_3)
  fig_4 <- fig_mv_demand_supply_sorad(supply_demand_sorad, nu_b_seq[4], nu_s_seq[4], r_min_max) + labs(x = NULL, y = NULL, subtitle = sub_4)
  
  yleft = gridtext::richtext_grob("Quantity", rot = 90, gp = grid::gpar(fontsize = 25))
  bottom = gridtext::richtext_grob(text = 'Price', gp = grid::gpar(fontsize = 25))
  glist <- gridExtra::arrangeGrob(
    fig_1+theme(legend.position = "none"),
    fig_2+theme(legend.position = "none"),
    fig_3+theme(legend.position = "none"),
    fig_4+theme(legend.position = "none"), nrow = 2, ncol = 2)
  fig <- gridExtra::grid.arrange(g_legend(fig_1+theme(legend.position = "top")), glist, 
                                 nrow=2,heights=c(1, 10),
                                 left = yleft, bottom = bottom)
  fig
}

fig_mv_demand_supply_sored <- function(supply_demand_sored_h, supply_demand_sored_uh, nu_b, nu_s, r_min_max = c(0.97, 1.03), digits = 4){
  
  # 1) Equilibrium prices
  P_eq_uh <- supply_demand_sored_uh$price(nu_s, nu_b)
  P_eq_h <- supply_demand_sored_h$price(nu_s, nu_b)
  # 2) Equilibrium quantities
  Q_eq_uh <- supply_demand_sored_uh$supply(P_eq_uh, nu_s)
  Q_eq_h <- supply_demand_sored_h$supply(P_eq_h, nu_s)
  # ************************************
  P0 <- mean(c(P_eq_uh, P_eq_h)) * seq(r_min_max[1], r_min_max[2], 0.001)
  # ************************************
  # Demand around the eq. price
  w_b_h <- supply_demand_sored_h$demand(P0, nu_b)
  w_b_uh <- supply_demand_sored_uh$demand(P0, nu_b)
  # Supply around the eq. price
  w_s_h <- supply_demand_sored_h$supply(P0, nu_s)
  w_s_uh <- supply_demand_sored_uh$supply(P0, nu_s)
  
  # Scales 
  x_breaks <- c(min(P0), P_eq_uh, P_eq_h, max(P0))
  x_labels <- round(x_breaks, digits = digits)
  y_minmax <- range(c(w_b_h, w_b_uh, w_s_h, w_s_uh))
  y_breaks <- c(y_minmax[1], Q_eq_uh, Q_eq_h, y_minmax[2])
  y_labels <- round(y_breaks, digits = digits)
  
  ggplot()+
    geom_line(aes(P0, w_b_uh, color = "demand_uh"))+
    geom_line(aes(P0, w_b_h, color = "demand_h"))+
    geom_line(aes(P0, w_s_h, color = "supply_h"))+
    geom_line(aes(P0, w_s_uh, color = "supply_uh"))+
    # Hedged equilibrium 
    geom_segment(aes(x = P_eq_h, xend = P_eq_h, y = min(c(w_s_h, w_s_uh)), yend = Q_eq_h), linetype = "dashed")+
    geom_segment(aes(x = min(P0), xend = P_eq_h, y = Q_eq_h, yend = Q_eq_h), linetype = "dashed")+
    geom_point(aes(x = P_eq_h, y = Q_eq_h, color = "eq_h"), shape = 9, size = 4)+
    # UnHedged equilibrium 
    geom_segment(aes(x = P_eq_uh, xend = P_eq_uh, y = min(c(w_s_h, w_s_uh)), yend = Q_eq_uh), linetype = "dashed")+
    geom_segment(aes(x = min(P0), xend = P_eq_uh, y = Q_eq_uh, yend = Q_eq_uh), linetype = "dashed")+
    geom_point(aes(x = P_eq_uh, y = Q_eq_uh, color = "eq_uh"), shape = 10, size = 4)+
    scale_x_continuous(breaks = x_breaks, labels = x_labels)+
    scale_y_continuous(breaks = y_breaks, labels = y_labels)+
    scale_color_manual(breaks = c("demand_uh", "supply_uh", "eq_uh",
                                  "demand_h", "supply_h", "eq_h"),
                       values = c(demand_uh = "darkgreen", supply_uh = "purple", eq_uh = "black",
                                  demand_h = "green", supply_h = "red", eq_h = "black"), 
                       labels = c(demand_uh = latex2exp::TeX("$q_t^{b}$"), 
                                  supply_uh = latex2exp::TeX("$q_t^{s}$"), 
                                  eq_uh = latex2exp::TeX("$V_t$"), 
                                  demand_h = latex2exp::TeX("$\\tilde{q}_{t}^{b}$"), 
                                  supply_h = latex2exp::TeX("$\\tilde{q}_{t}^{s}$"), 
                                  eq_h = latex2exp::TeX("$\\widetilde{V}_t$")))+
    theme_bw()+
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(legend.position = "top")+
    labs(color = NULL, x = "Price", y = "Quantity")+
    figure_theme
}

fig_demand_supply_sored <- function(supply_demand_sored_h, supply_demand_sored_uh, nu_b, nu_s, r_min_max = c(0.97, 1.03), digits = 3){
  
  nu_b_seq <- nu_b * c(1, 2, 1, 3)
  nu_s_seq <- nu_s * c(1, 1, 2, 0.75)
  # Subtitles 
  sub_1 <- latex2exp::TeX(paste0("$\\nu_b = ", round(nu_b_seq[1], 3), ",\\; \\nu_s = ", round(nu_s_seq[1], 3), "$"))
  sub_2 <- latex2exp::TeX(paste0("$\\nu_b = ", round(nu_b_seq[2], 3), ",\\; \\nu_s = ", round(nu_s_seq[2], 3), "$"))
  sub_3 <- latex2exp::TeX(paste0("$\\nu_b = ", round(nu_b_seq[3], 3), ",\\; \\nu_s = ", round(nu_s_seq[3], 3), "$"))
  sub_4 <- latex2exp::TeX(paste0("$\\nu_b = ", round(nu_b_seq[4], 3), ",\\; \\nu_s = ", round(nu_s_seq[4], 3), "$"))
  # Figures
  fig_1 <- fig_mv_demand_supply_sored(supply_demand_sored_h, supply_demand_sored_uh, nu_b_seq[1], nu_s_seq[1], r_min_max, digits) + labs(x = NULL, y = NULL, subtitle = sub_1)
  fig_2 <- fig_mv_demand_supply_sored(supply_demand_sored_h, supply_demand_sored_uh, nu_b_seq[2], nu_s_seq[2], r_min_max, digits) + labs(x = NULL, y = NULL, subtitle = sub_2)
  fig_3 <- fig_mv_demand_supply_sored(supply_demand_sored_h, supply_demand_sored_uh, nu_b_seq[3], nu_s_seq[3], r_min_max, digits) + labs(x = NULL, y = NULL, subtitle = sub_3)
  fig_4 <- fig_mv_demand_supply_sored(supply_demand_sored_h, supply_demand_sored_uh, nu_b_seq[4], nu_s_seq[4], r_min_max, digits) + labs(x = NULL, y = NULL, subtitle = sub_4)
  
  yleft = gridtext::richtext_grob("Quantity", rot = 90, gp = grid::gpar(fontsize = 20))
  bottom = gridtext::richtext_grob(text = 'Price', gp = grid::gpar(fontsize = 20))
  glist <- gridExtra::arrangeGrob(
    fig_1+theme(legend.position = "none"),
    fig_2+theme(legend.position = "none"),
    fig_3+theme(legend.position = "none"),
    fig_4+theme(legend.position = "none"), nrow = 2, ncol = 2)
  fig <- gridExtra::grid.arrange(g_legend(fig_1+theme(legend.position = "top")), glist, 
                                 nrow=2,heights=c(1, 10),
                                 left = yleft, bottom = bottom)
  fig
}
