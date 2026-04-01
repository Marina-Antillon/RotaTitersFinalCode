# Cross-Correlation Plot with Quadrant Interpretation 
# (not part of the analysis; just a graph for the methods)

library(ggplot2)
library(dplyr)

# Create the "out" dir if not available
source("./functions/fn_create_dir_if_nec.R")
create_dir_if_necessary("./out/")

# Sample data structure (replace with your actual final_results)
final_results <- data.frame(
  lag = -6:6,
  iga_median = c(-0.1, -0.2, -0.3, -0.4, -0.57, -0.2, 0.1, 0.25, 0.38, 0.3, 0.2, 0.1, 0.05),
  iga_q025 = c(-0.3, -0.4, -0.5, -0.6, -0.71, -0.4, -0.1, 0.05, 0.18, 0.1, 0.0, -0.1, -0.15),
  iga_q975 = c(0.1, 0.0, -0.1, -0.2, -0.42, 0.0, 0.3, 0.45, 0.58, 0.5, 0.4, 0.3, 0.25),
  igg_median = c(-0.05, -0.1, -0.2, -0.3, -0.55, 0.41, 0.2, 0.15, 0.1, 0.05, 0.0, -0.05, -0.1),
  igg_q025 = c(-0.25, -0.3, -0.4, -0.5, -0.7, 0.21, 0.0, -0.05, -0.1, -0.15, -0.2, -0.25, -0.3),
  igg_q975 = c(0.15, 0.1, 0.0, -0.1, -0.4, 0.61, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1),
  iga_significant = c(0.3, 0.4, 0.6, 0.7, 0.8, 0.5, 0.4, 0.6, 0.7, 0.6, 0.4, 0.3, 0.2),
  igg_significant = c(0.2, 0.3, 0.5, 0.6, 0.99, 0.95, 0.4, 0.3, 0.2, 0.2, 0.1, 0.1, 0.1)
)

# Prepare data for plotting
plot_data <- final_results %>%
  dplyr::select(lag, 
         iga_median, iga_q025, iga_q975, 
         igg_median, igg_q025, igg_q975, 
         iga_significant, igg_significant) %>%
  reshape2::melt(id.vars = c("lag", "iga_significant", "igg_significant"), 
                measure.vars = c("iga_median", "iga_q025", "iga_q975", 
                               "igg_median", "igg_q025", "igg_q975")) %>%
  mutate(
    antibody = ifelse(grepl("iga", variable), "IgA", "IgG"),
    statistic = case_when(
      grepl("median", variable) ~ "median",
      grepl("q025", variable) ~ "q025", 
      grepl("q975", variable) ~ "q975"
    ),
    significant_prop = ifelse(antibody == "IgA", iga_significant, igg_significant)
  ) %>%
  reshape2::dcast(lag + antibody + significant_prop ~ statistic, value.var = "value")

# Create quadrant annotations
quadrant_data <- data.frame(
  xmin = c(0, -6, -6, 0),
  xmax = c(6, 0, 0, 6),
  ymin = c(0, 0, -1, -1),
  ymax = c(1, 1, 0, 0),
  quadrant = c("I: Positive Lag\nPositive Correlation\n(Boosting)", 
               "II: Biologically implausible\n(High antibodies would predict increased disease)", 
               "III: Negative Lag\nNegative Correlation\n(Susceptible Buildup)", 
               "IV: Positive Lag\nNegative Correlation\n(Depletion)"),
  alpha_val = c(0.15, 0.15, 0.15, 0.15),
  fill_color = c("#ff9999", "#99ff99", "#9999ff", "#ffff99")
)


  # Add quadrant backgrounds
  # Enhanced plot with quadrants
  p_ccf_quadrants <- ggplot() +
  

  # Add quadrant backgrounds
    # geom_rect(data=quadrant_data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill_color),
    #           alpha=0.15, show.legend = F) + 
  # Add quadrant labels
    geom_text(data=quadrant_data, aes(x=(xmax+xmin)/2, y=(ymin+ymax)/1.25, label=quadrant), 
              size = 3, hjust = 0.5, vjust = 0.5, fontface = "bold") +  
    
  # Confidence intervals
  geom_ribbon(data=plot_data, aes(x = lag, ymin = q025, ymax = q975, fill = antibody),
              alpha = 0.3, color=NA) +
  
  # Median correlation lines
  geom_line(data=plot_data, aes(x = lag, color = antibody, y = median), size = 1.5) +
  geom_point(data=plot_data, aes(x = lag, color = antibody, y = median, size = significant_prop), alpha = 0.8) +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.7, size = 1) +
  
  # Highlight key points
  geom_point(data = data.frame(lag = c(-2.05, -1.95), median = c(-0.57, -0.55), 
                              antibody = c("IgA", "IgG")), 
             aes(x = lag, y = median, color = antibody), 
             size = 6, shape = 21, fill = "white", stroke = 3) +
  
  scale_color_manual(values = c("IgA" = "red", "IgG" = "blue")) +
  scale_fill_manual(values = c("IgA" = "red", "IgG" = "blue")) +
  scale_size_continuous(name = "Proportion\nSignificant", 
                       range = c(2, 5), 
                       labels = scales::percent,
                       guide = guide_legend(override.aes = list(color = "black"))) +
  scale_x_continuous(breaks = -6:6, limits = c(-6.5, 6.5)) +
  scale_y_continuous(limits = c(-1.25, 1.25)) +
  
  labs(
    title = "Cross-Correlation Function: Quadrant Interpretation",
    # subtitle = "Individual-Level Analysis with Bootstrap Uncertainty",
    x = "Lag (months)", 
    y = "Cross-Correlation",
    color = "Antibody Type",
    fill = "Antibody Type",
    caption = "Shaded areas = 95% bootstrap CI\nPoint size = proportion significant\nLarge outlined points = key findings: Depletion (lag +1), Protection (lag -1), Boosting (lag +4)"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90")
  ) + 
    annotate("text", x = -3, y = 1.2, label = "← Antibodies change before incidence", 
             size = 3, hjust = 0.5, vjust = 0.5, fontface = "italic", color = "gray20") +
    annotate("text", x = 3, y = 1.2, label = "Incidence changes before antibodies →", 
             size = 3, hjust = 0.5, vjust = 0.5, fontface = "italic", color = "gray20") +
    annotate("text", x = -3, y = 1.1, label = "↑ Antibodies & incidence move TOGETHER", 
             size = 3, hjust = 0.5, vjust = 0.5, fontface = "italic", color = "gray20") +
    annotate("text", x = -3, y = -1.1, label = "↓ Antibodies & incidence move in OPPOSITION", 
             size = 3, hjust = 0.5, vjust = 0.5, fontface = "italic", color = "gray20") 
    

print(p_ccf_quadrants)
ggsave("./out/ccf_quadrants.pdf", width=7, height=7, units="in", bg=NULL)

