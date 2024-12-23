# load required libraries
library(dplyr)    # v1.1.2
library(readxl)   # v1.4.2
library(ggplot2)  # v3.4.2

plot_toposcore_density <- function(scores, clin, var = 'TOPO', color = 'OS12bis', 
                                 xlim = NULL, grey = NULL) {
  
  # Join data and create OS category
  data <- scores %>% 
    left_join(clin, by = 'Sample_id') %>%
    mutate(OS_category = ifelse(OS < 12, "OS<12", "OS≥12"))
  
  # Create the plot
  plt <- ggplot(data, aes(x = .data[[var]], fill = OS_category)) +
    # Add density plots for each OS category
    geom_density(alpha = 0.8) +
    # Add vertical lines for cutoffs
    geom_vline(xintercept = grey, 
              color = "black", 
              linewidth = 0.5,
              linetype = "solid") +
    # Set colors for OS categories
    scale_fill_manual(values = c(
      "OS<12" = "#000000",    # black for OS<12
      "OS≥12" = "#99FF99"       # Light green for OS≥12
    )) +
    # Customize theme
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.position = "top",
      legend.title = element_blank(),
      # Add more space at the bottom for labels
      plot.margin = margin(t = 5, r = 5, b = 40, l = 5, unit = "pt")
    ) +
    # Expand y axis limits to make room for labels
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    # Labels
    labs(x = "Score",
         y = "Density") +
    # Add region labels with percentages
    annotate("text", x = 0.3, y = 1.8, 
            label = "SIG1\n68%OS<12", 
            color = "black",
            size = 5) +
    annotate("text", x = 0.66, y = 2.5, 
            label = "Gray zone", 
            color = "black", 
            size = 5) +
    annotate("text", x = 0.9, y = 1.8, 
            label = "SIG2\n21%OS<12", 
            color = "darkgreen", 
            size = 5)
  
  # Add limits if specified
  if (!is.null(xlim)) plt <- plt + xlim(xlim)
  
  plt
}

plot_toposcoreb01_density <- function(scores, clin, resp = 'OS12', lims = c(0.535, 0.791)) {
  plot_toposcore_density(scores, clin, 
                        var = 'TOPOB01', 
                        xlim = c(0, 1), 
                        grey = lims)
}