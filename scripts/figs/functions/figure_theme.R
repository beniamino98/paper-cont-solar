figure_theme <- theme(# Title
  plot.title  = element_text(face = "bold", size = 30),
  # Subtitle
  plot.subtitle = element_text(size = 24),
  # Caption
  plot.caption = element_text(face = "italic"),
  # Axis-x
  axis.title.x = element_text(face = "bold", size = 19),
  axis.text.x = element_text(face = "bold", size = 14),
  axis.ticks.x = element_line(linewidth = 0.2),
  axis.line.x = element_line(),
  # Grid x-axis
  panel.grid.minor.x = element_line(),
  panel.grid.major.x = element_line(),
  # Axis-y
  axis.title.y = element_text(size = 20),
  axis.text.y = element_text(size = 15),
  axis.ticks.y = element_line(linewidth = 0.2),
  axis.line.y = element_line(),
  # Grid x-axis
  panel.grid.minor.y = element_line(),
  panel.grid.major.y = element_line(),
  # Legend
  legend.title = element_text(face = "bold", size = 25),
  legend.text = element_text(face = "italic", size = 20),
  legend.box.background = element_rect())



# extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Function to create y-labels 
create_y_breaks_return <- function(x, n = 7, mean_x = NA, digits = 0) {
  if (!is.na(mean_x)) {
    if (mean_x < 0){
      # Upper breaks
      y_breaks_hi <- seq(0, max(x), length.out = trunc(n/3*2+1))[-1]
      # Lower breaks
      y_breaks_lo1 <- seq(min(x), mean_x, length.out = n/3)
      y_breaks_lo2 <- seq(mean_x, 0, length.out = n/3)[-1]
      y_breaks_lo <- c(y_breaks_lo1, y_breaks_lo2)
    } else {
      # Lower breaks
      y_breaks_lo <- seq(min(x), 0, length.out = trunc(n/3*2))
      # Upper breaks
      y_breaks_hi1 <- seq(0, mean_x, length.out = n/3)
      y_breaks_hi2 <- seq(mean_x, max(x), length.out = n/3)[-1]
      y_breaks_hi <- c(y_breaks_hi1, y_breaks_hi2)[-1]
    }
  } else {
    # Lower breaks
    y_breaks_lo <- seq(min(x), 0, length.out = n/2)
    # Upper breaks
    y_breaks_hi <- seq(0, max(x), length.out = n/2)[-1]
  }
  
  # y-breaks 
  y_breaks <- c(y_breaks_lo, y_breaks_hi)
  # y-labels 
  y_labels <- paste0(round(y_breaks*100, digits), " %")
  
  list(
    breaks = y_breaks,
    labels = y_labels
  )
}
