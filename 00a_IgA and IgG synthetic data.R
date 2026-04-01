# Load necessary packages

# TODO: remove kurtosis? It doesn't seem that anything is made with that.
# TODO: add these data (mat-age and birthweight) to the fake dataset.

library(tidyverse)   # For data manipulation and visualization
library(ggplot2)
library(gridExtra)

# Function to generate synthetic data with specified moments
generate_synthetic_data <- function(n_samples, geom_mean, geom_sd, skewness) { # , kurtosis
  # For antibody titers, we often use log-normal distribution
  # Convert geometric mean and SD to parameters of log-normal distribution
  mu <- log(geom_mean)
  sigma <- log(geom_sd)
  
  # Generate data from normal distribution
  log_data <- rnorm(n_samples, mean = mu, sd = sigma)
  
  # Transform to match target skewness
  # This is a simplified approach - we adjust the generated normal distribution
  # to better match skewness, while maintaining the log-normal nature
  if (skewness != 0) {
    # Add skewness by using a skew-normal transformation
    # This is a simplification but works reasonably well for many biological data
    alpha <- skewness * 2  # Approximate transformation
    u <- runif(n_samples)
    v <- runif(n_samples)
    z <- sqrt(-2 * log(u)) * cos(2 * pi * v)  # Box-Muller transform
    skew_factor <- pmax(-3, pmin(3, alpha * z))  # Constrain to reasonable range
    
    # Apply skewness adjustment (small adjustment to avoid extreme values)
    log_data <- log_data + skew_factor * 0.2 * sigma
  }
  
  # Convert back to original scale
  data <- exp(log_data)
  
  return(data)
}

# Create synthetic monthly data for rotavirus antibody titers
create_monthly_synthetic_data <- function() {
  # Define months
  months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  
  # Example statistical properties for IgA by month
  # Based on data presented in data: overall GMT of 413, peak in Jan (714), lowest in Oct (274)
  # Format: list(n_samples, geom_mean, geom_sd, skewness, kurtosis)
  iga_stats <- list(
    "Jan" = list(30, 714, 1.9, 0.8),  # Peak month 
    "Feb" = list(25, 650, 1.9, 0.9), 
    "Mar" = list(28, 670, 1.8, 0.8), 
    "Apr" = list(26, 600, 1.8, 0.7), 
    "May" = list(27, 550, 1.8, 0.6),
    "Jun" = list(25, 500, 1.7, 0.6), 
    "Jul" = list(28, 450, 1.7, 0.5), 
    "Aug" = list(30, 400, 1.7, 0.5), 
    "Sep" = list(26, 350, 1.7, 0.4), 
    "Oct" = list(27, 274, 1.6, 0.4), # Trough month 
    "Nov" = list(28, 350, 1.7, 0.5),
    "Dec" = list(29, 500, 1.8, 0.7) 
  )
  
  # Similar stats for IgG 
  # Based on paper: IgG GMT of 731 overall, peak in Jan (1370), low in Oct (508)
  igg_stats <- list(
    "Jan" = list(30, 1370, 1.9, 0.9),  # Peak month 
    "Feb" = list(25, 1200, 1.9, 0.9), 
    "Mar" = list(28, 1100, 1.8, 0.8), 
    "Apr" = list(26, 1000, 1.8, 0.7), 
    "May" = list(27, 900, 1.8, 0.6), 
    "Jun" = list(25, 800, 1.7, 0.6), 
    "Jul" = list(28, 750, 1.7, 0.5), 
    "Aug" = list(30, 700, 1.7, 0.5),
    "Sep" = list(26, 600, 1.7, 0.4), 
    "Oct" = list(27, 508, 1.6, 0.4), # Trough month 
    "Nov" = list(28, 600, 1.7, 0.5), 
    "Dec" = list(29, 900, 1.8, 0.7) 
  )
  
  iga_lloq = c(1,2,4,1,0,1,1,1,1,0,1,0)
  names(iga_lloq) = months
  
  # Create empty data frame to store results
  synthetic_data <- data.frame()
  
  # Generate data for each month
  for (month in months) {
    # Get stats for current month
    iga_params <- iga_stats[[month]]
    igg_params <- igg_stats[[month]]
    
    n_samples <- iga_params[[1]]
    
    # Generate IgA and IgG values
    iga_values <- generate_synthetic_data(
      n_samples, 
      iga_params[[2]], # geom_mean
      iga_params[[3]], # geom_sd
      iga_params[[4]]  # skewness
    )
    
    igg_values <- generate_synthetic_data(
      n_samples, 
      igg_params[[2]], # geom_mean
      igg_params[[3]], # geom_sd
      igg_params[[4]]  # skewness
    )
    
    # Create data frame for this month
    month_data <- data.frame(
      month = rep(month, n_samples),
      month_num = rep(which(months == month), n_samples),
      day_num = 1:n_samples, 
      iga = pmax(1, iga_values),  # Ensure positive values
      igg = pmax(1, igg_values)   # Ensure positive values
    )
    
    if(iga_lloq[month]>0){
      n_samples_lloq = iga_lloq[month]
      lloq_data <- data.frame(
        month = rep(month, n_samples_lloq),
        month_num = rep(which(months == month), n_samples_lloq),
        day_num = 1:n_samples_lloq, 
        iga = 7.5,  # Ensure positive values
        igg = pmax(1, sample(igg_values, n_samples_lloq)) # Ensure positive values
        )   
      # Append to main data frame
      synthetic_data <- rbind(synthetic_data, month_data, lloq_data)
    } else {
      # Append to main data frame
      synthetic_data <- rbind(synthetic_data, month_data)
    }

  }
  
  return(synthetic_data)
}

# Generate the synthetic dataset
set.seed(123)  # For reproducibility
synthetic_df = create_monthly_synthetic_data()
synthetic_df = synthetic_df %>% 
  mutate(mat_age = rnorm(342, 25, 5.5), 
        bw = rnorm(342, 2977, 300), 
        baby_sex = rbinom(342, 1, 0.45)+1) # 1: Boys, 2: Girls

synthetic_df$mat_age[synthetic_df$mat_age<14] = 14
synthetic_df$mat_age[synthetic_df$mat_age>42] = 42

summary(synthetic_df$mat_age)
summary(synthetic_df$bw)
summary(synthetic_df$baby_sex)

# Calculate geometric means by month to verify the synthetic data behaves like the real one
monthly_stats <- synthetic_df %>%
  group_by(month, month_num) %>%
  summarize(
    iga_gmt = exp(mean(log(iga))),
    igg_gmt = exp(mean(log(igg))),
    n = n(),
    .groups = 'drop'
  ) %>%
  arrange(month_num)

print(monthly_stats)

# Plot the data to visualize
# Set up the plotting area
par(mfrow=c(2,1), mar=c(4,4,2,1))

# Plot IgA
plot(synthetic_df$month_num, synthetic_df$iga, 
     pch=19, col=rgb(0,0,1,0.3), 
     xlab="Month", ylab="IgA Titer",
     main="Synthetic IgA Titers by Month",
     xaxt="n")
axis(1, at=1:12, labels=unique(synthetic_df$month))
lines(monthly_stats$month_num, monthly_stats$iga_gmt, 
      col="red", lwd=2)

# Plot IgG
plot(synthetic_df$month_num, synthetic_df$igg, 
     pch=19, col=rgb(0,0,1,0.3), 
     xlab="Month", ylab="IgG Titer",
     main="Synthetic IgG Titers by Month",
     xaxt="n")
axis(1, at=1:12, labels=unique(synthetic_df$month))
lines(monthly_stats$month_num, monthly_stats$igg_gmt, 
      col="red", lwd=2)

# Save the synthetic dataset
write.csv(synthetic_df, "./data/synthetic_rotavirus_titers-final.csv", row.names = FALSE)

cat("Generated synthetic dataset with", nrow(synthetic_df), "samples\n")
cat("Statistical summary of synthetic data:\n")
print(summary(synthetic_df))

# Create a more comprehensive summary of the synthetic data
summary_stats <- synthetic_df %>%
  group_by(month) %>%
  summarize(
    n = n(),
    iga_gmt = exp(mean(log(iga))),
    iga_gsd = exp(sd(log(iga))),
    iga_95ci_lower = exp(mean(log(iga)) - 1.96 * sd(log(iga))/sqrt(n())),
    iga_95ci_upper = exp(mean(log(iga)) + 1.96 * sd(log(iga))/sqrt(n())),
    igg_gmt = exp(mean(log(igg))),
    igg_gsd = exp(sd(log(igg))),
    igg_95ci_lower = exp(mean(log(igg)) - 1.96 * sd(log(igg))/sqrt(n())),
    igg_95ci_upper = exp(mean(log(igg)) + 1.96 * sd(log(igg))/sqrt(n())),
    .groups = 'drop'
  ) %>%
  arrange(match(month, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))

print(summary_stats)

# Optional: create nicer plots with ggplot2

# IgA Plot
p1 <- ggplot(synthetic_df, aes(x = factor(month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")), 
                               y = iga)) +
  geom_jitter(alpha = 0.3, width = 0.2, color = "blue") +
  geom_point(data = monthly_stats, aes(x = factor(month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")), 
                                       y = iga_gmt), color = "red", size = 3) +
  geom_line(data = monthly_stats, aes(x = factor(month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")), 
                                      y = iga_gmt, group = 1), color = "red") +
  labs(title = "Synthetic IgA Titers by Month", 
       x = "Month", 
       y = "IgA Titer") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_log10()

# IgG Plot
p2 <- ggplot(synthetic_df, aes(x = factor(month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")), 
                               y = igg)) +
  geom_jitter(alpha = 0.3, width = 0.2, color = "blue") +
  geom_point(data = monthly_stats, aes(x = factor(month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")), 
                                       y = igg_gmt), color = "red", size = 3) +
  geom_line(data = monthly_stats, aes(x = factor(month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")), 
                                      y = igg_gmt, group = 1), color = "red") +
  labs(title = "Synthetic IgG Titers by Month", 
       x = "Month", 
       y = "IgG Titer") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_log10()

# Combine plots
combined_plot <- grid.arrange(p1, p2, nrow = 2)

# Save the combined plot
ggsave("./data/synthetic_data_visualization.png", combined_plot, width = 10, height = 8)