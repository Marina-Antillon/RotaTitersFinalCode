# INDIVIDUAL-LEVEL CROSS-CORRELATION ANALYSIS WITH INCIDENCE UNCERTAINTY -----

# Uses all individual maternal titers and incorporates incidence CI uncertainty

# This file produces: 
# base_ccf is the CCF with mean incidence values (no uncertainty on incidence, but using all 329 Ig values). 
# It is done with the function 'individual_ccf'.
# uncertainty_impact summary of whether uncertainty in the incidence meaningfully changes the conclusions compared to base_ccf
# gamma_params: shape and rate parameters for incidence of each month, when modelled by gamma distribution
# validation_data: data frame comparing the incidence data to the modelled gamma distributions that describe uncertainty in the CCF analysis
# titer_data: titer data from mothers
# monthly_incidence: the incidence summaries from GEMS-1a, MSD cases only.
# bootstrap_results gives the correlation coefficient and the pvalue for each time 
  # lag from -6 to 6. Dimensions are iterations x lags x (which Ig and which result)
  # Use this one in the master .Rmd to combine with other results
# final_results is the ggplot with the CCF. Draws heavily from bootstrap_results.
  # Use this one in the master .Rmd to combine with other results
# final_results is saved as a pdf called "ccf_data_uncertainty.pdf"

library(dplyr)
library(ggplot2)
library(readr)
library(reshape2)
library(boot) # For bootstrap confidence intervals
library(MASS) # fitdistr

# Load the gamma parameter finding function
# source("./functions/fn_find_gamma.R")

# Create the "out" dir if not available
create_dir_if_necessary("./out/")

cat("=== INDIVIDUAL-LEVEL CCF WITH INCIDENCE UNCERTAINTY ===\n")

# *****************************************************************************
## 1. Data preparation --------------------------------------------------------
# *****************************************************************************

# Read individual-level data
# Read and process data
titer_data = readr::read_csv("./data/synthetic_rotavirus_titers-final.csv") %>%
  mutate(DOB=paste0("2012-",sprintf("%02d", month_num),"-",sprintf("%02d", day_num))) %>%
  mutate(DOB=base::as.Date(DOB)) %>% 
  filter(iga>7.5)

cat("Processing individual-level titer data...\n")
cat("Total individuals:", nrow(titer_data), "\n")
print(table(titer_data$month))

# Method of moments for gamma distribution
# Mean = shape/rate, Variance = shape/rate^2
# Therefore: rate = mean/variance, shape = mean * rate
# rate <- mean / approx_var
# shape <- mean * rate

monthly_incidence_lsd <- read.csv("./data/epi_desc_1a_lsd.csv")
monthly_incidence_msd <- read.csv("./data/epi_desc_1a_msd.csv")

monthly_incidence_lsd = monthly_incidence_lsd[4:15,] %>% 
  mutate(lower = ifelse(lower<1, 0.5, lower)) %>% 
  mutate(mean = ifelse(mean<1, lower+0.25, mean)) %>% 
  mutate(upper = ifelse(upper<1, mean+0.25, upper)) %>% 
  # rowwise() %>% 
  # mutate(shape = find_gamma_pars(lower, upper)[1], 
  #        rate = 1/find_gamma_pars(lower, upper)[2]) %>% 
  mutate(approx_sd = (upper - lower) / (2 * 1.96)) %>% 
  mutate(approx_var = approx_sd^2) %>% 
  mutate(shape = mean*mean / approx_var, 
         rate = mean / approx_var) %>% 
  ungroup %>%
  mutate(Month = c(month.abb[11:12], month.abb[1:10])) %>%
  mutate(Month = factor(Month, levels=month.abb)) %>% 
  relocate(Month)

# Apr, May, Jul need change
monthly_incidence_lsd$shape[monthly_incidence_lsd$Month=="Mar"] = monthly_incidence_lsd$shape[monthly_incidence_lsd$Month=="Feb"]
monthly_incidence_lsd$rate[monthly_incidence_lsd$Month=="Mar"] = monthly_incidence_lsd$rate[monthly_incidence_lsd$Month=="Feb"]

monthly_incidence_lsd$shape[monthly_incidence_lsd$Month=="Apr"] = monthly_incidence_lsd$shape[monthly_incidence_lsd$Month=="Feb"]
monthly_incidence_lsd$rate[monthly_incidence_lsd$Month=="Apr"] = monthly_incidence_lsd$rate[monthly_incidence_lsd$Month=="Feb"]

monthly_incidence_lsd$shape[monthly_incidence_lsd$Month=="May"] = monthly_incidence_lsd$shape[monthly_incidence_lsd$Month=="Feb"]
monthly_incidence_lsd$rate[monthly_incidence_lsd$Month=="May"] = monthly_incidence_lsd$rate[monthly_incidence_lsd$Month=="Feb"]

monthly_incidence_lsd$shape[monthly_incidence_lsd$Month=="Jun"] = monthly_incidence_lsd$shape[monthly_incidence_lsd$Month=="Aug"]
monthly_incidence_lsd$rate[monthly_incidence_lsd$Month=="Jun"] = monthly_incidence_lsd$rate[monthly_incidence_lsd$Month=="Aug"]

monthly_incidence_lsd$shape[monthly_incidence_lsd$Month=="Jul"] = monthly_incidence_lsd$shape[monthly_incidence_lsd$Month=="Aug"]
monthly_incidence_lsd$rate[monthly_incidence_lsd$Month=="Jul"] = monthly_incidence_lsd$rate[monthly_incidence_lsd$Month=="Aug"]

monthly_incidence_msd = monthly_incidence_msd[4:15,] %>% 
  mutate(lower = ifelse(lower<1, 0.5, lower)) %>% 
  mutate(mean = ifelse(mean<1, lower+0.25, mean)) %>% 
  mutate(upper = ifelse(upper<1, mean+0.25, upper)) %>% 
  # rowwise() %>% 
  # mutate(shape = find_gamma_pars(lower, upper)[1], 
  #        rate = 1/find_gamma_pars(lower, upper)[2]) %>% 
  mutate(approx_sd = (upper - lower) / (2 * 1.96)) %>% 
  mutate(approx_var = approx_sd^2) %>% 
  mutate(shape = mean*mean / approx_var, 
         rate = mean / approx_var) %>% 
  ungroup %>%
  mutate(Month = c(month.abb[11:12], month.abb[1:10])) %>%
  mutate(Month = factor(Month, levels=month.abb)) %>% 
  relocate(Month)

# Mar, June, Jul, Aug need change
monthly_incidence_msd$shape[monthly_incidence_msd$Month=="Mar"] = monthly_incidence_msd$shape[monthly_incidence_msd$Month=="Feb"]
monthly_incidence_msd$rate[monthly_incidence_msd$Month=="Mar"] = monthly_incidence_msd$rate[monthly_incidence_msd$Month=="Feb"]

monthly_incidence_msd$shape[monthly_incidence_msd$Month=="Jun"] = monthly_incidence_msd$shape[monthly_incidence_msd$Month=="May"]
monthly_incidence_msd$rate[monthly_incidence_msd$Month=="Jun"] = monthly_incidence_msd$rate[monthly_incidence_msd$Month=="May"]

monthly_incidence_msd$shape[monthly_incidence_msd$Month=="Jul"] = monthly_incidence_msd$shape[monthly_incidence_msd$Month=="May"]
monthly_incidence_msd$rate[monthly_incidence_msd$Month=="Jul"] = monthly_incidence_msd$rate[monthly_incidence_msd$Month=="May"]

monthly_incidence_msd$shape[monthly_incidence_msd$Month=="Jul"] = monthly_incidence_msd$shape[monthly_incidence_msd$Month=="May"]
monthly_incidence_msd$rate[monthly_incidence_msd$Month=="Jul"] = monthly_incidence_msd$rate[monthly_incidence_msd$Month=="May"]

monthly_lsd_sim <- mapply(function(shape, rate) {rgamma(10000, shape = shape, rate = rate)}, 
  shape = monthly_incidence_lsd$shape, 
  rate = monthly_incidence_lsd$rate)
monthly_msd_sim <- mapply(function(shape, rate) {rgamma(10000, shape = shape, rate = rate)}, 
  shape = monthly_incidence_msd$shape, 
  rate = monthly_incidence_msd$rate)
monthly_all_sim = monthly_lsd_sim + monthly_msd_sim

# Fit gamma to each column (month)
monthly_all_fits <- apply(monthly_all_sim, 2, function(x) {
  fitdistr(x, "gamma")
})

# Extract shape and rate parameters
monthly_all_mean <- apply(monthly_all_sim, 2, mean)
monthly_all_lower <- apply(monthly_all_sim, 2, quantile, 0.025)
monthly_all_upper <- apply(monthly_all_sim, 2, quantile, 0.975)
monthly_all_shape <- sapply(monthly_all_fits, function(fit) fit$estimate["shape"])
monthly_all_rate <- sapply(monthly_all_fits, function(fit) fit$estimate["rate"])

monthly_incidence_total = monthly_incidence_lsd %>% 
  dplyr::select("Month", "Group") %>%
  mutate(mean = monthly_all_mean,
         lower = monthly_all_lower, upper = monthly_all_upper,
         shape = monthly_all_shape, rate = monthly_all_rate)

cat("Incidence data prepared for gamma uncertainty:\n")

# *****************************************************************************
## 2. Functions ----------------------------------------------------------------
# *****************************************************************************

### A. Wrap-around function (same as catalytic model) ---------------------------

get_incidence_wraparound <- function(birth_month, lag, incidence_vec) {
  incidence_position <- birth_month + 2 # because it starts in November
  target_position <- incidence_position - lag
  while (target_position <= 0) {
    target_position <- target_position + 12
  }
  while (target_position > 12) {
    target_position <- target_position - 12
  }
  return(incidence_vec[target_position])
}

### B. Individual-level cross-correlation function ------------------------------

# Function to calculate individual-level correlation with proper uncertainty
individual_ccf <- function(titer_data, incidence_mean_vec, max_lag = 6) {
  lags <- -max_lag:max_lag
  results <- data.frame(
    lag = lags,
    iga_correlation = NA,
    igg_correlation = NA,
    iga_pvalue = NA,
    igg_pvalue = NA,
    n_individuals = NA
  )
  
  for (i in seq_along(lags)) {
    lag <- lags[i]
    
    # Get incidence values for each individual's birth month at this lag
    individual_incidence <- sapply(titer_data$month_num, function(bm) {
      get_incidence_wraparound(bm, lag, incidence_mean_vec)
    })
    
    # Calculate correlations using all individuals
    iga_test <- cor.test(titer_data$iga, individual_incidence, 
                         method="spearman", exact=FALSE, use = "complete.obs")
    igg_test <- cor.test(titer_data$igg, individual_incidence, 
                         method="spearman", exact=FALSE, use = "complete.obs")
    
    results$iga_correlation[i] <- iga_test$estimate
    results$igg_correlation[i] <- igg_test$estimate
    results$iga_pvalue[i] <- iga_test$p.value
    results$igg_pvalue[i] <- igg_test$p.value
    results$n_individuals[i] <- sum(complete.cases(titer_data$iga, individual_incidence))
  }
  
  return(results)
}

### C. Bootstrap ccf functions --------------------------------------------------

# Function to sample from incidence uncertainty using gamma distributions
bootstrap_ccf_with_gamma_uncertainty <- function(titer_data, monthly_incidence, 
                                                max_lag = 6, n_bootstrap = 1000) {
  
  cat("Setting up gamma distributions for incidence uncertainty...\n")

  # Calculate gamma parameters for each month
  gamma_params <- matrix(NA, nrow = 12, ncol = 2)
  colnames(gamma_params) <- c("shape", "rate")

  for (month in 1:12) {

    params <- monthly_incidence[monthly_incidence$Month==month.abb[month],c("shape", "rate")] %>% as.numeric()
    gamma_params[month, ] = params

    cat(sprintf("Month %d: mean=%.4f, CI=[%.4f, %.4f] -> Gamma(%.3f, %.3f)\n",
                month,
                monthly_incidence$mean[monthly_incidence$Month==month.abb[month]],
                monthly_incidence$lower[monthly_incidence$Month==month.abb[month]],
                monthly_incidence$upper[monthly_incidence$Month==month.abb[month]],
                params[1], params[2]))
  }
  
  cat("Running bootstrap analysis with", n_bootstrap, "samples...\n")
  
  # Storage for bootstrap results
  bootstrap_results <- array(NA, dim = c(n_bootstrap, 2*max_lag + 1, 4),
                           dimnames = list(NULL, paste0("lag_", -max_lag:max_lag), 
                                         c("iga_corr", "igg_corr", "iga_pval", "igg_pval")))
  
  # Progress tracking
  progress_interval <- max(1, floor(n_bootstrap / 10))
  
  for (boot_iter in 1:n_bootstrap) {
    if (boot_iter %% progress_interval == 0) {
      cat("Bootstrap iteration:", boot_iter, "/", n_bootstrap, "\n")
    }
    
    # Sample incidence values from gamma distributions
    sampled_incidence <- numeric(12)
    for (month in 1:12) {
      shape <- gamma_params[month, "shape"]
      rate <- gamma_params[month, "rate"]
      
      # Sample from gamma distribution
      sampled_incidence[month] <- rgamma(1, shape = shape, rate = rate)
    }
    
    # Calculate CCF for this bootstrap sample
    ccf_result <- individual_ccf(titer_data, sampled_incidence[c(11:12, 1:10)], max_lag)
    
    # Store results
    bootstrap_results[boot_iter, , "iga_corr"] <- ccf_result$iga_correlation
    bootstrap_results[boot_iter, , "igg_corr"] <- ccf_result$igg_correlation
    bootstrap_results[boot_iter, , "iga_pval"] <- ccf_result$iga_pvalue
    bootstrap_results[boot_iter, , "igg_pval"] <- ccf_result$igg_pvalue
  }
  
  # Return both results and gamma parameters for validation
  return(list(
    bootstrap_results = bootstrap_results,
    gamma_params = gamma_params
  ))
}

### D. Lag contribution function ------------------------------------------------

# Add this analysis after your main CCF calculation
analyze_lag_contributions <- function(titer_data, incidence_vec, lag) {
  # Get incidence for each individual at this lag
  individual_incidence <- sapply(titer_data$month_num, function(bm) {
    get_incidence_wraparound(bm, lag, incidence_vec)
  })
  
  # Calculate correlation contribution by birth month
  month_contributions <- titer_data %>%
    mutate(incidence_lag = individual_incidence) %>%
    group_by(month_num, month) %>%
    summarise(
      n = n(),
      mean_iga = mean(iga),
      mean_incidence = first(incidence_lag),
      correlation_iga = cor(iga, incidence_lag),
      .groups = 'drop'
    ) %>%
    arrange(month_num)
  
  return(month_contributions)
}

# *****************************************************************************
## 3. Analysis -----------------------------------------------------------------
# *****************************************************************************

### A. Main analysis ---------------------------------------------------------

maxlag = 3

for(i in 1:3){ # 3 situations: LSD, MSD, both
  
  if(i==1){
    inc_level = "lsd"
    monthly_incidence = monthly_incidence_lsd
  } else if(i==2){
    inc_level = "msd"
    monthly_incidence = monthly_incidence_msd
  } else {
    inc_level = "total"
    monthly_incidence = monthly_incidence_total
  }

  cat("\n=== RUNNING MAIN ANALYSIS ===\n")
  
  # First, calculate CCF with mean incidence values
  base_ccf <- individual_ccf(titer_data, monthly_incidence$mean, max_lag = maxlag)
  
  cat("Base CCF results (using mean incidence):\n")
  print(round(base_ccf, 4))
  
  # Analyze lag +1 specifically
  lag1_breakdown <- analyze_lag_contributions(titer_data, monthly_incidence$mean, lag = 1)
  print("Lag +1 breakdown by month:")
  print(lag1_breakdown)
  
  # Run bootstrap analysis with gamma uncertainty
  set.seed(123)  # For reproducibility
  n_bootstrap <- 2000  # Adjust based on computational constraints
  
  cat("\nStarting bootstrap analysis with gamma incidence uncertainty...\n")
  bootstrap_output <- bootstrap_ccf_with_gamma_uncertainty(
    titer_data, monthly_incidence, max_lag = maxlag, n_bootstrap = n_bootstrap
  )
  
  bootstrap_results <- bootstrap_output$bootstrap_results
  gamma_params <- bootstrap_output$gamma_params
  
  # *****************************************************************************
  ### B. Validate gamma distributions -------------------------------------------
  # *****************************************************************************
  
  cat("\n=== VALIDATING GAMMA DISTRIBUTIONS ===\n")
  
  # Create validation plots for gamma distributions
  validation_data <- data.frame()
  
  for (month in 1:12) {
    shape <- gamma_params[month, "shape"]
    rate <- gamma_params[month, "rate"]
    
    observed_mean <- monthly_incidence$mean[monthly_incidence$Month==month.abb[month]]
    observed_lower <- monthly_incidence$lower[monthly_incidence$Month==month.abb[month]]
    observed_upper <- monthly_incidence$upper[monthly_incidence$Month==month.abb[month]]
    
    # Calculate gamma distribution statistics
    gamma_mean <- shape / rate
    gamma_q025 <- qgamma(0.025, shape = shape, rate = rate)
    gamma_q975 <- qgamma(0.975, shape = shape, rate = rate)
    
    validation_data <- rbind(validation_data, data.frame(
      month = month,
      period_label = monthly_incidence$Month[monthly_incidence$Month==month.abb[month]],
      observed_mean = observed_mean,
      observed_lower = observed_lower,
      observed_upper = observed_upper,
      gamma_mean = gamma_mean,
      gamma_q025 = gamma_q025,
      gamma_q975 = gamma_q975,
      shape = shape,
      rate = rate
    ))
  }
  
  cat("Gamma distribution validation:\n")
  print(validation_data %>% 
        dplyr::select(period_label, observed_mean, gamma_mean, observed_lower, gamma_q025, 
               observed_upper, gamma_q975) %>%
        mutate(across(where(is.numeric), ~round(.x, 4))))
  
  # Check how well gamma distributions match the original CIs
  mean_difference <- mean(abs(validation_data$observed_mean - validation_data$gamma_mean))
  lower_difference <- mean(abs(validation_data$observed_lower - validation_data$gamma_q025))
  upper_difference <- mean(abs(validation_data$observed_upper - validation_data$gamma_q975))
  
  cat("\nGamma distribution fit quality:\n")
  cat("Mean absolute difference in means:", round(mean_difference, 4), "\n")
  cat("Mean absolute difference in lower bounds:", round(lower_difference, 4), "\n")
  cat("Mean absolute difference in upper bounds:", round(upper_difference, 4), "\n")
  
  # *****************************************************************************
  ### C. Summarize bootstrap results --------------------------------------------
  # *****************************************************************************
  
  cat("\n=== Summarize bootstrap results ===\n")
  
  # Calculate summary statistics from bootstrap
  lags <- (-maxlag):maxlag
  final_results <- data.frame(
    lag = lags,
    
    # IgA results
    iga_median = apply(bootstrap_results[, , "iga_corr"], 2, median, na.rm = TRUE),
    iga_mean = apply(bootstrap_results[, , "iga_corr"], 2, mean, na.rm = TRUE),
    iga_q025 = apply(bootstrap_results[, , "iga_corr"], 2, quantile, 0.025, na.rm = TRUE),
    iga_q975 = apply(bootstrap_results[, , "iga_corr"], 2, quantile, 0.975, na.rm = TRUE),
    iga_significant = apply(bootstrap_results[, , "iga_pval"], 2, function(x) mean(x < 0.05, na.rm = TRUE)),
    
    # IgG results  
    igg_median = apply(bootstrap_results[, , "igg_corr"], 2, median, na.rm = TRUE),
    igg_mean = apply(bootstrap_results[, , "igg_corr"], 2, mean, na.rm = TRUE),
    igg_q025 = apply(bootstrap_results[, , "igg_corr"], 2, quantile, 0.025, na.rm = TRUE),
    igg_q975 = apply(bootstrap_results[, , "igg_corr"], 2, quantile, 0.975, na.rm = TRUE),
    igg_significant = apply(bootstrap_results[, , "igg_pval"], 2, function(x) mean(x < 0.05, na.rm = TRUE))
  )
  
  cat("Final results with uncertainty:\n")
  print(round(final_results, 3))
  
  # Find optimal lags
  iga_optimal_idx <- which.max(abs(final_results$iga_median))
  igg_optimal_idx <- which.max(abs(final_results$igg_median))
  
  iga_optimal_lag <- final_results$lag[iga_optimal_idx]
  igg_optimal_lag <- final_results$lag[igg_optimal_idx]
  iga_optimal_corr <- final_results$iga_median[iga_optimal_idx]
  igg_optimal_corr <- final_results$igg_median[igg_optimal_idx]
  
  cat("\n=== OPTIMAL LAGS WITH UNCERTAINTY ===\n")
  cat("IgA optimal lag:", iga_optimal_lag, "months, correlation:", round(iga_optimal_corr, 3), "\n")
  cat("IgG optimal lag:", igg_optimal_lag, "months, correlation:", round(igg_optimal_corr, 3), "\n")
  
  # *****************************************************************************
  ### D. Visualization with gamma uncertainty -----------------------------------
  # *****************************************************************************
  
  cat("\n=== Visualizations ===\n")
  
  # Prepare data for plotting with confidence intervals
  plot_data <- final_results %>%
    dplyr::select(lag, iga_median, iga_q025, iga_q975, igg_median, igg_q025, igg_q975, 
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
  
  # Create a dataset for adding labels where significance > 0.3
  label_data <- plot_data %>%
    filter(significant_prop > 0.3) %>%
    mutate(
      # Position based on correlation sign and antibody
      label_y = case_when(
        antibody == "IgA" & lag==-3 ~ median - 0.025,
        antibody == "IgA" & lag==-2 ~ median + 0.05,
        antibody == "IgA" & lag==-1 ~ median - 0.04,
        antibody == "IgG" & lag==-2 ~ median - 0.03,
        antibody == "IgG" & lag==-1 ~ median + 0.04,
        antibody == "IgG" & lag==3 ~ median + 0.03,
      ),
      # Format the percentage labels
      label_text = paste0(round(significant_prop * 100), "%")
    )
  
  label_data$label_y[label_data$lag>0] = 0.15
  
  # Main plot with uncertainty bands
  p_ccf_uncertainty <- ggplot(plot_data, aes(x = lag, color = antibody, fill = antibody)) +
    # Confidence intervals
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.2, color = NA, show.legend = F) +
    
    # Median correlation lines
    geom_line(aes(y = median), size = 1.75) +
    geom_point(aes(y = median, size = significant_prop), alpha = 0.8) +
    
    # Add significance percentage labels for points where >30% significant
    geom_text(data = label_data, 
              aes(x = lag, y = label_y, label = label_text, color = antibody),
              size = 3, 
              fontface = "bold",
              show.legend = FALSE) +
    
    # Reference lines
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    
    scale_color_manual(values = c("IgA" = "red", "IgG" = "blue"), 
                       guide = guide_legend(override.aes = list(linetype = 1, shape = NA))) +
    scale_fill_manual(values = c("IgA" = "red", "IgG" = "blue"), guide = "none") +
    scale_size_continuous(name = "Proportion\nSignificant", 
                         range = c(2, 6), 
                         labels = scales::percent,
                         guide = guide_legend(override.aes = list(color = "black"))) +
    scale_x_continuous(breaks = (-maxlag):maxlag) +
    scale_y_continuous(limits=c(-0.51, 0.51), breaks=round(seq(-0.50, 0.50, by=0.1), 1)) + 
    labs(
      title = paste0("CCF between titers and rotavirus-attributable ", ifelse(i==1, "LSD", ifelse(i==2, "MSD", "diarrhea"))),
      # subtitle = paste0("N = ", nrow(titer_data), " maternal anti-rotavirus titers at birth | Bootstrap samples = ", n_bootstrap, " for incidence"),
      x = "Lag (months)", 
      y = "Spearman's rank correlation",
      color = "Antibody Type",
      # fill = "Antibody Type",
      caption = "Shaded areas = 95% bootstrap CI\nPoint size = proportion of bootstrap samples with p < 0.05\nPercentages shown where >30% of iterations significant"
    ) +
    theme_bw() + 
    theme(
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(hjust=0.5, size=12)
    ) + 
    annotate("text", x = -1.5, y = -0.05, label = "Antibodies fall before incidence rises\n(Susceptible Buildup)", 
             size = 3, hjust = 0.5, vjust = 0.5, fontface = "italic", color = "gray20") +
    annotate("text", x = 1.5, y = 0.35, label = "Incidence rises before antibodies rise\n(Boosting)", 
             size = 3, hjust = 0.5, vjust = 0.5, fontface = "italic", color = "gray20") +
    # annotate("text", x = -3, y = 1.1, label = "↑ Antibodies & incidence move TOGETHER", 
    #          size = 3, hjust = 0.5, vjust = 0.5, fontface = "italic", color = "gray20") +
    annotate("text", x = 1.5, y = -0.25, label = "Incidence rises before antibodies fall\n(Depletion)", 
             size = 3, hjust = 0.5, vjust = 0.5, fontface = "italic", color = "gray20") 
  
  print(p_ccf_uncertainty)
  final_results
  ggsave(paste0("./out/ccf_data_uncertainty_", inc_level, ".pdf"), width=7.5, height=5, units="in", bg=NULL, device = cairo_pdf)
  save(p_ccf_uncertainty, final_results, plot_data, label_data, file=paste0("./out/ccf_data_uncertainty_", inc_level, "_plot.RData"))
  
  # *****************************************************************************
  ### E. Significance testing and interpretation --------------------------------
  # *****************************************************************************
  
  cat("\n=== SIGNIFICANCE ANALYSIS ===\n")
  
  # Identify consistently significant lags
  significant_threshold <- 0.75  # At least 75% of bootstrap samples significant
  
  iga_significant_lags <- final_results$lag[final_results$iga_significant >= significant_threshold]
  igg_significant_lags <- final_results$lag[final_results$igg_significant >= significant_threshold]
  
  cat("IgA significantly correlated lags (", significant_threshold*100, "% threshold):", 
      paste(iga_significant_lags, collapse = ", "), "\n")
  cat("IgG significantly correlated lags (", significant_threshold*100, "% threshold):", 
      paste(igg_significant_lags, collapse = ", "), "\n")
  
  # Check if uncertainty meaningfully changes conclusions
  uncertainty_impact <- final_results %>%
    mutate(
      iga_ci_width = iga_q975 - iga_q025,
      igg_ci_width = igg_q975 - igg_q025,
      iga_base_in_ci = (base_ccf$iga_correlation >= iga_q025) & (base_ccf$iga_correlation <= iga_q975),
      igg_base_in_ci = (base_ccf$igg_correlation >= igg_q025) & (base_ccf$igg_correlation <= igg_q975)
    )
  
  cat("\nUncertainty impact assessment:\n")
  cat("Mean IgA CI width:", round(mean(uncertainty_impact$iga_ci_width), 3), "\n")
  cat("Mean IgG CI width:", round(mean(uncertainty_impact$igg_ci_width), 3), "\n")
  cat("Base IgA estimates within bootstrap CI:", mean(uncertainty_impact$iga_base_in_ci) * 100, "%\n")
  cat("Base IgG estimates within bootstrap CI:", mean(uncertainty_impact$igg_base_in_ci) * 100, "%\n")
  
  # *****************************************************************************
  ### F. Summary ----------------------------------------------------------------
  # *****************************************************************************
  
  cat("\n=== Summary ===\n")
  cat("Individual-level cross-correlation analysis with gamma incidence uncertainty:\n\n")
  
  cat("Methods:\n")
  cat("- Individual maternal samples: N =", nrow(titer_data), "\n")
  cat("- Gamma distributions for incidence uncertainty (respects non-negativity)\n")
  cat("- Bootstrap iterations: N =", n_bootstrap, "\n")
  cat("- Uses fn_find_gamma for parameter estimation from confidence intervals\n\n")
  
  cat("Primary findings:\n")
  interpret_lag <- function(lag, correlation, antibody_type) {
    if (lag > 0) {
      direction <- paste("peaks", lag, "months AFTER incidence")
    } else if (lag < 0) {
      direction <- paste("peaks", abs(lag), "months BEFORE incidence")
    } else {
      direction <- "peaks simultaneously with incidence"
    }
    cat(sprintf("- %s %s (median r = %.3f)\n", antibody_type, direction, correlation))
  }
  
  interpret_lag(iga_optimal_lag, iga_optimal_corr, "IgA")
  interpret_lag(igg_optimal_lag, igg_optimal_corr, "IgG")
  
  cat("\nUncertainty:\n")
  cat("- Gamma distributions respect non-negative nature of incidence data\n")
  cat("- Handles months with lower bounds at 0.01/1000 appropriately\n")
  cat("- Mean gamma fit error:", round(mean_difference, 4), "\n")
  cat("- Analysis accounts for asymmetric uncertainty distributions\n\n")
  
  cat("Robustness:\n")
  cat("- Analysis uses", nrow(titer_data), "individual measurements (not monthly means)\n")
  cat("- Bootstrap procedure uses appropriate gamma uncertainty distributions\n")
  cat("- Results show", round(mean(uncertainty_impact$iga_base_in_ci + uncertainty_impact$igg_base_in_ci)/2 * 100), 
      "% consistency between point estimates and uncertainty-adjusted estimates\n")
  
  cat("\nStatistical significance:\n")
  cat("- Significant IgA correlations at lags:", paste(iga_significant_lags, collapse = ", "), "\n")
  cat("- Significant IgG correlations at lags:", paste(igg_significant_lags, collapse = ", "), "\n")
  
  # Save comprehensive results
  save(final_results, bootstrap_results, base_ccf, uncertainty_impact, 
       gamma_params, validation_data, titer_data, monthly_incidence, 
       file = paste0("./out/ccf_gamma_uncertainty_", inc_level, ".RData"))
  
  cat(paste0("\nAnalysis complete. Results saved to './out/ccf_gamma_uncertainty_", inc_level, ".RData'"))

}
