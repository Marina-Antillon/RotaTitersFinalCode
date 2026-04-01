# CCF Validation with Simple Single-Peak Scenarios
# Testing cross-correlation functions with clear, interpretable patterns
# (not part of the analysis; just a graph for the methods)

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(cowplot)

# Create the "out" dir if not available
source("./functions/fn_create_dir_if_nec.R")
create_dir_if_necessary("./out/")

# =============================================================================
# 1. WRAP-AROUND FUNCTION
# =============================================================================

get_incidence_wraparound <- function(birth_month, lag, incidence_vec) {
  incidence_position <- birth_month  # No +2 offset needed for synthetic data (both start in Jan)
  target_position <- incidence_position - lag
  while (target_position <= 0) {
    target_position <- target_position + 12
  }
  while (target_position > 12) {
    target_position <- target_position - 12
  }
  return(incidence_vec[target_position])
}

# =============================================================================
# 2. INDIVIDUAL-LEVEL CCF FUNCTION
# =============================================================================

individual_ccf <- function(titer_data, incidence_vec, max_lag = 6) {
  lags <- -max_lag:max_lag
  results <- data.frame(
    lag = lags,
    correlation = NA,
    pvalue = NA
  )
  
  for (i in seq_along(lags)) {
    lag <- lags[i]
    
    # Get incidence values for each individual's birth month at this lag
    individual_incidence <- sapply(titer_data$month_num, function(bm) {
      get_incidence_wraparound(bm, lag, incidence_vec)
    })
    
    # Calculate correlation using all individuals
    test_result <- cor.test(titer_data$antibody, individual_incidence, 
                           method = "pearson", use = "complete.obs")
    
    results$correlation[i] <- test_result$estimate
    results$pvalue[i] <- test_result$p.value
  }
  
  return(results)
}

# =============================================================================
# 3. SIMPLE SCENARIO GENERATION
# =============================================================================

# Generate synthetic individuals
generate_individuals <- function(n_per_month = 30) {
  months <- rep(1:12, each = n_per_month)
  data.frame(
    month_num = months,
    month = month.abb[months],
    individual_id = 1:length(months)
  )
}

# Create very simple scenarios with single peaks
create_simple_scenario <- function(individuals, scenario_name, 
                                  incidence_peak_month, antibody_peak_month,
                                  base_incidence = 2, peak_incidence = 20,
                                  base_antibody = 100, peak_antibody = 500) {  # More realistic antibody range
  
  # Create simple incidence pattern: low everywhere except one peak/trough month
  incidence_monthly <- rep(base_incidence, 12)
  incidence_monthly[incidence_peak_month] <- peak_incidence
  
  # Create simple antibody pattern: one level everywhere except one peak/trough month
  individuals$antibody <- ifelse(individuals$month_num == antibody_peak_month, 
                                peak_antibody, base_antibody)
  
  # Add realistic log-normal noise (antibodies are on exponential scale)
  individuals$antibody <- 10^(log10(individuals$antibody) + rnorm(nrow(individuals), 0, 0.015))
  
  # Calculate expected lag and correlation
  expected_lag <- antibody_peak_month - incidence_peak_month
  if (expected_lag > 6) expected_lag <- expected_lag - 12
  if (expected_lag < -6) expected_lag <- expected_lag + 12
  
  # Expected correlation sign based on the actual relationship
  if (scenario_name == "Depletion (Lag +2)") {
    # High incidence (June) → Low antibodies (July): negative correlation
    expected_correlation_sign <- -1
  } else if (scenario_name == "Susceptible Buildup (Lag -2)") {
    # Low antibodies (May) → High incidence (June): negative correlation  
    expected_correlation_sign <- -1
  } else if (scenario_name == "Boosting (Lag +2)") {
    # High incidence (June) → High antibodies (August): positive correlation
    expected_correlation_sign <- 1
  } else if (scenario_name == "Protection (Lag -2)") {
    # High antibodies (June) → Low incidence (July): negative correlation
    expected_correlation_sign <- -1
  } else {
    # No relationship
    expected_correlation_sign <- 0
  }
  
  individuals$scenario <- scenario_name
  
  return(list(
    individuals = individuals,
    incidence_monthly = incidence_monthly,
    incidence_peak_month = incidence_peak_month,
    antibody_peak_month = antibody_peak_month,
    expected_lag = expected_lag,
    expected_correlation_sign = expected_correlation_sign
  ))
}

# =============================================================================
# 4. CREATE SIMPLE TEST SCENARIOS
# =============================================================================

set.seed(123)
base_individuals <- generate_individuals(n_per_month = 30)

# Let's use month 6 (June) as incidence peak for all scenarios
incidence_peak <- 6

scenarios <- list(
  # Depletion: Antibodies drop 2 months AFTER incidence peak
  # Incidence peaks in June (6), antibodies are low in Aug (8)
  # Expected lag: +2 (antibodies follow incidence)
  create_simple_scenario(base_individuals, "Depletion (Lag +2)", 
                        incidence_peak_month = 6, antibody_peak_month = 8,
                        base_antibody = 500, peak_antibody = 50),  # Low antibodies at peak month
  
  # Susceptible Buildup: Antibodies drop 2 months BEFORE incidence peak  
  # Incidence peaks in June (6), antibodies are low in April (4)
  # Expected lag: -2 (antibodies predict incidence)
  create_simple_scenario(base_individuals, "Susceptible Buildup (Lag -2)", 
                        incidence_peak_month = 6, antibody_peak_month = 4,
                        base_antibody = 500, peak_antibody = 50),  # Low antibodies at peak month
  
  # Boosting: Antibodies peak 2 months AFTER incidence peak
  # Incidence peaks in June (6), antibodies peak in August (8)  
  # Expected lag: +2 (antibodies follow incidence)
  create_simple_scenario(base_individuals, "Boosting (Lag +2)", 
                        incidence_peak_month = 6, antibody_peak_month = 8,
                        base_antibody = 50, peak_antibody = 500) #,  # High antibodies at peak month
  
)

# =============================================================================
# 5. RUN CCF ANALYSIS FOR ALL SCENARIOS
# =============================================================================

cat("=== RUNNING CCF ANALYSIS ON SIMPLE SCENARIOS ===\n")

ccf_results <- list()
time_series_data <- list()

for (i in seq_along(scenarios)) {
  scenario <- scenarios[[i]]
  scenario_name <- scenario$individuals$scenario[1]
  
  cat("\nAnalyzing:", scenario_name, "\n")
  cat("Incidence peak month:", scenario$incidence_peak_month, "\n")
  cat("Antibody peak month:", scenario$antibody_peak_month, "\n") 
  cat("Expected lag:", scenario$expected_lag, "\n")
  
  # Run CCF analysis
  ccf_result <- individual_ccf(scenario$individuals, scenario$incidence_monthly, max_lag = 6)
  ccf_result$scenario <- scenario_name
  ccf_result$expected_lag <- scenario$expected_lag
  ccf_result$expected_correlation_sign <- scenario$expected_correlation_sign
  
  ccf_results[[i]] <- ccf_result
  
  # Prepare time series data for plotting
  monthly_summary <- scenario$individuals %>%
    group_by(month_num, month) %>%
    summarise(
      mean_antibody = mean(antibody),
      se_antibody = sd(antibody)/sqrt(n()),
      .groups = 'drop'
    ) %>%
    mutate(
      incidence = scenario$incidence_monthly,
      scenario = scenario_name
    ) %>%
    arrange(month_num)
  
  time_series_data[[i]] <- monthly_summary
  
  # Show the actual patterns
  cat("Incidence pattern:", paste(scenario$incidence_monthly, collapse = ", "), "\n")
  cat("Antibody pattern:", paste(round(monthly_summary$mean_antibody, 1), collapse = ", "), "\n")
}

# Combine results
all_ccf <- do.call(rbind, ccf_results)
all_timeseries <- do.call(rbind, time_series_data)

# =============================================================================
# 6. CREATE CLEAN FACETED VISUALIZATION
# =============================================================================

# Prepare data for time series faceted plot
ts_facet_data <- all_timeseries %>%
  mutate(
    incidence_norm = incidence/max(incidence),
    antibody_norm = mean_antibody/max(mean_antibody),
    scenario_clean = case_when(
      scenario == "Depletion (Lag +2)" ~ "Depletion",
      scenario == "Susceptible Buildup (Lag -2)" ~ "Susceptible Buildup", 
      scenario == "Boosting (Lag +2)" ~ "Boosting",
      TRUE ~ scenario
    )
  ) %>%
  # Convert to long format for easier plotting
  pivot_longer(cols = c(incidence_norm, antibody_norm), 
               names_to = "variable", values_to = "value") %>%
  mutate(variable = case_when(
    variable == "incidence_norm" ~ "Incidence",
    variable == "antibody_norm" ~ "Antibodies",
    TRUE ~ variable
  )) %>%
  mutate(scenario_clean = factor(scenario_clean, 
                                 levels = c("Boosting", "Susceptible Buildup", "Depletion")))

# Create time series faceted plot
p_timeseries <- ggplot(ts_facet_data, aes(x = month_num, y = value, color = variable)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Incidence" = "black", "Antibodies" = "blue")) +
  scale_x_continuous(breaks = c(1, 3, 6, 9, 12), 
                     labels = c("Jan", "Mar", "Jun", "Sep", "Dec")) +
  facet_grid(scenario_clean ~ ., scales = "free_y") +
  labs(x = "Month", y = "Normalized Values", color = NULL) +
  theme_linedraw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text.y = element_blank(), # element_text(angle = 0, hjust = 0),
    legend.key.size = unit(0.5, "cm")
  ) +
  ylim(0, 1.1)

# Prepare data for CCF faceted plot
ccf_facet_data <- all_ccf %>%
  mutate(
    scenario_clean = case_when(
      scenario == "Depletion (Lag +2)" ~ "Q4: Depletion",
      scenario == "Susceptible Buildup (Lag -2)" ~ "Q3: Susceptible Buildup",
      scenario == "Boosting (Lag +2)" ~ "Q1: Boosting", 
      TRUE ~ scenario
    ),
    significant = pvalue < 0.05
  ) %>%
  group_by(scenario) %>%
  mutate(
    optimal_lag = lag[which.max(abs(correlation))],
    optimal_corr = correlation[which.max(abs(correlation))]
  ) %>%
  ungroup() %>%
  mutate(scenario_clean = factor(scenario_clean, 
                                 levels = c("Q1: Boosting", "Q3: Susceptible Buildup", "Q4: Depletion")))

# Create CCF faceted plot  
p_ccf <- ggplot(ccf_facet_data, aes(x = lag, y = correlation)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 1, linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 1, linewidth = 1) +
  geom_line(color = "gray50", linewidth = 1) +
  geom_point(aes(size = significant, color = significant), alpha = 0.8) +
  scale_size_manual(values = c("TRUE" = 4, "FALSE" = 2), guide = "none") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50"), guide = "none") +
  scale_x_continuous(breaks = seq(-6, 6, 2)) +
  facet_grid(scenario_clean ~ ., scales = "free_y") +
  labs(x = "Lag (months)", y = "Correlation") +
  theme_linedraw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(angle = 270, face="bold", size=12) # element_blank()  # Remove strip text since it's redundant with left plot
  ) +
  ylim(-1.1, 1.1)

# Combine the two plots with column labels
final_plot <- plot_grid(
  # Add column headers
  plot_grid(
    ggdraw() + draw_label("Time Series", fontface = "bold", size = 12),
    ggdraw() + draw_label("Cross-Correlation Functions", fontface = "bold", size = 12),
    ncol = 2, rel_widths = c(0.6, 0.4)
  ),
  # Add the main plots
  plot_grid(p_timeseries, p_ccf, ncol = 2, rel_widths = c(0.6, 0.4), align = "h", axis = "tb"),
  ncol = 1, rel_heights = c(0.05, 0.95)
)

print(final_plot)

# =============================================================================
# 7. VALIDATION SUMMARY
# =============================================================================

cat("\n=== SIMPLE SCENARIO VALIDATION SUMMARY ===\n")

validation_summary <- all_ccf %>%
  group_by(scenario) %>%
  summarise(
    max_correlation = correlation[which.max(abs(correlation))],
    optimal_lag = lag[which.max(abs(correlation))],
    expected_lag = first(expected_lag),
    expected_correlation_sign = first(expected_correlation_sign),
    detected_correlation_sign = sign(max_correlation),
    lag_match = (expected_lag == optimal_lag),
    sign_match = (expected_correlation_sign == detected_correlation_sign),
    validation_success = lag_match & sign_match,
    .groups = 'drop'
  )

print(validation_summary)

cat("\nDetailed Validation Results:\n")
for (i in 1:nrow(validation_summary)) {
  row <- validation_summary[i,]
  lag_result <- ifelse(row$lag_match, "✓", "✗")
  sign_result <- ifelse(row$sign_match, "✓", "✗")
  overall <- ifelse(row$validation_success, "✓ PASS", "✗ FAIL")
  
  cat(sprintf("%s:\n", row$scenario))
  cat(sprintf("  Expected lag: %d, Found lag: %d %s\n", 
              row$expected_lag, row$optimal_lag, lag_result))
  cat(sprintf("  Expected correlation sign: %d, Found sign: %d %s\n", 
              row$expected_correlation_sign, row$detected_correlation_sign, sign_result))
  cat(sprintf("  Max correlation: %.3f\n", row$max_correlation))
  cat(sprintf("  Overall: %s\n\n", overall))
}

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

ggsave("./out/ccf_simple_validation.pdf", final_plot, 
       width = 10, height = 8, units = "in", bg = "white")

save(scenarios, all_ccf, all_timeseries, validation_summary,
     file = "./out/ccf_simple_validation_results.RData")

cat("Simple validation complete!\n")
cat("Results saved to 'ccf_simple_validation.pdf' and 'ccf_simple_validation_results.RData'\n")
