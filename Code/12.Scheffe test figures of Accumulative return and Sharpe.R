library(ggplot2)
library(reshape2)
library(scales)  
library(dplyr)
library(readxl)

hs_300 = read_xlsx('/Users/mac/Desktop/GitHub-English version/data/HS300 dataset.xlsx')

hs_300_returns <- as.numeric(hs_300[[2]])/100

hs_300_dates <- as.Date(hs_300[[1]])

r_daily <- numeric(length(hs_300_returns))
r_daily[1] <- hs_300_returns[1] 


for (t in 2:length(hs_300_returns)) {
  r_daily[t] <- (1 + hs_300_returns[t]) / (1 + hs_300_returns[t - 1]) - 1
}



Mean = as.matrix(read.csv('/Users/mac/Desktop/GitHub-English version/data/Daily return.csv'))


col_names <- c("SCAD(tau=0.01)", 
               "SCAD(tau=0.05)", 
               "SCAD(tau=0.1)", 
               "MVP", "EWP", "MVO")

colnames(Mean) <- col_names



calc_rolling_risk <- function(returns, window_size = 30, roll_step = 1, alpha = 0.05, rf_daily = 0) {
  n <- length(returns)
  roll_indices <- seq(1, n - window_size + 1, by = roll_step)
  
  results <- data.frame(
    window_start = roll_indices,
    window_end   = roll_indices + window_size - 1,
    cum_return   = NA,
    VaR          = NA,
    CVaR         = NA,
    max_drawdown = NA,
    RAROC        = NA,
    Sharpe       = NA,   
    RAROC_7d_MA  = NA    
  )
  
  for (i in seq_along(roll_indices)) {
    idx_start <- roll_indices[i]
    idx_end   <- idx_start + window_size - 1
    r_win     <- returns[idx_start:idx_end]
    

    cum_r <- prod(1 + r_win) - 1
    

    port_loss <- -r_win
    VaR_val   <- quantile(port_loss, alpha, type = 8)
    CVaR_val  <- mean(port_loss[port_loss >= VaR_val])
    

    wealth_index <- cumprod(1 + r_win)
    cum_max      <- cummax(wealth_index)
    drawdowns    <- (wealth_index - cum_max) / cum_max
    max_dd       <- min(drawdowns, na.rm = TRUE)
    

    raroc_val <- ifelse(CVaR_val != 0, cum_r / CVaR_val, NA)
    

    excess_r <- r_win - rf_daily
    sharpe_val <- ifelse(sd(excess_r, na.rm = TRUE) != 0,
                         mean(excess_r, na.rm = TRUE) / sd(excess_r, na.rm = TRUE),
                         NA)
    

    results$cum_return[i]   <- cum_r
    results$VaR[i]          <- VaR_val
    results$CVaR[i]         <- CVaR_val
    results$max_drawdown[i] <- max_dd
    results$RAROC[i]        <- raroc_val
    results$Sharpe[i]       <- sharpe_val
  }
  
  return(results)
}

window_size <- 60    # or window_size <-90
roll_step <- 1     
alpha <- 0.99       # this value makes no difference here 

risk_results <- lapply(1:ncol(Mean), function(j) {
  calc_rolling_risk(Mean[, j], window_size = window_size, roll_step = roll_step, alpha = alpha)
})
r_daily_risk_results <- calc_rolling_risk(r_daily, window_size = window_size, roll_step = roll_step, alpha = alpha)


library(dplyr)
library(tidyr)
library(agricolae)
library(ggplot2)

# ------------------------------------------------------------

# ------------------------------------------------------------
cum_returns <- lapply(risk_results, function(x) x$cum_return)
benchmark_cum_return <- r_daily_risk_results$cum_return

df_cum <- as.data.frame(do.call(cbind, cum_returns))
colnames(df_cum) <- c("SCAD(tau=0.01)", "SCAD(tau=0.05)", "SCAD(tau=0.1)", "MVP", "EWP", "MVO")
df_cum$benchmark <- benchmark_cum_return

df_cum_long <- pivot_longer(df_cum,
                            cols = everything(),
                            names_to = "strategy",
                            values_to = "cum_return"
) %>% drop_na()

# ------------------------------------------------------------

# ------------------------------------------------------------
anova_cum <- aov(cum_return ~ strategy, data = df_cum_long)
summary(anova_cum)

scheffe_cum <- scheffe.test(anova_cum, "strategy")
group_info_cum <- scheffe_cum$groups %>% tibble::rownames_to_column("strategy")

# ------------------------------------------------------------

# ------------------------------------------------------------
summary_cum <- df_cum_long %>%
  group_by(strategy) %>%
  summarise(
    mean = mean(cum_return, na.rm = TRUE),
    sd = sd(cum_return, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    lower = mean - qt(0.975, df = n - 1) * se,
    upper = mean + qt(0.975, df = n - 1) * se
  ) %>%
  arrange(mean) %>%
  mutate(strategy = factor(strategy, levels = strategy))

plot_data_cum <- summary_cum %>%
  left_join(group_info_cum, by = "strategy")

# ------------------------------------------------------------

# ------------------------------------------------------------
benchmark_info_cum <- summary_cum %>%
  dplyr::filter(.data$strategy == "benchmark")

benchmark_upper <- benchmark_info_cum$upper
benchmark_lower <- benchmark_info_cum$lower
benchmark_mean <- benchmark_info_cum$mean

# ------------------------------------------------------------

# ------------------------------------------------------------
plot_data_cum <- plot_data_cum %>%
  mutate(sig_higher = ifelse(lower > benchmark_upper, TRUE, FALSE))

# ------------------------------------------------------------

# ------------------------------------------------------------
ggplot(plot_data_cum, aes(x = mean, y = strategy)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = sig_higher), height = 0.25) +
  geom_point(aes(color = sig_higher), size = 3) +
  geom_text(aes(x = upper + 0.1, y = strategy, label = groups), size = 5, fontface = "bold") +
  geom_vline(xintercept = benchmark_mean, linetype = "dashed", color = "blue") +
  annotate(
    "rect",
    xmin = benchmark_lower, xmax = benchmark_upper,
    ymin = -Inf, ymax = Inf,
    alpha = 0.1, fill = "blue"
  ) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "gray40")) +
  labs(
    title = "Scheffé Comparison of Cumulative Return Means",
    subtitle = "Mean ± 95% CI",
    x = "Cumulative Return Mean ± 95% CI",
    y = "Strategy",
    color = "Significantly Higher\nthan Benchmark?"
  ) +
  theme_minimal(base_size = 14)

# ------------------------------------------------------------

# ------------------------------------------------------------
df_plot_cum <- df_cum_long %>%
  left_join(dplyr::select(plot_data_cum, strategy, sig_higher, groups), by = "strategy")

df_plot_cum$strategy <- factor(df_plot_cum$strategy,
                               levels = plot_data_cum$strategy)

ggplot(df_plot_cum, aes(x = strategy, y = cum_return)) +
  geom_boxplot(aes(fill = sig_higher),
               outlier.size = 0.8, outlier.alpha = 0.5) +
  geom_text(data = plot_data_cum,
            aes(x = strategy,
                y = max(df_plot_cum$cum_return, na.rm = TRUE) + 1,
                label = groups),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "gray70")) +
  labs(title = "Rolling Cumulative Return Boxplot with Scheffé Grouping",
       x = "Strategy",
       y = "Cumulative Return",
       fill = "Significantly Higher\nthan Benchmark?") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")



# ================================================

# ------------------------------------------------------------

# ------------------------------------------------------------
sharpe_values <- lapply(risk_results, function(x) x$Sharpe)
benchmark_sharpe <- r_daily_risk_results$Sharpe

df_sharpe <- do.call(cbind, sharpe_values)
colnames(df_sharpe) <- c("SCAD(tau=0.01)", "SCAD(tau=0.05)", "SCAD(tau=0.1)", "MVP", "EWP", "MVO")
df_sharpe <- cbind(df_sharpe, benchmark = benchmark_sharpe)

df_long_sharpe <- pivot_longer(as.data.frame(df_sharpe),
                               cols = everything(),
                               names_to = "strategy",
                               values_to = "Sharpe"
) %>% drop_na()

# ------------------------------------------------------------

# ------------------------------------------------------------
anova_sharpe <- aov(Sharpe ~ strategy, data = df_long_sharpe)
summary(anova_sharpe)

scheffe_sharpe <- scheffe.test(anova_sharpe, "strategy")
group_info_sharpe <- scheffe_sharpe$groups %>% tibble::rownames_to_column("strategy")

# ------------------------------------------------------------

# ------------------------------------------------------------
summary_sharpe <- df_long_sharpe %>%
  group_by(strategy) %>%
  summarise(
    mean = mean(Sharpe, na.rm = TRUE),
    sd = sd(Sharpe, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    lower = mean - qt(0.975, df = n - 1) * se,
    upper = mean + qt(0.975, df = n - 1) * se
  ) %>%
  arrange(mean) %>%
  mutate(strategy = factor(strategy, levels = strategy))

plot_data_sharpe <- summary_sharpe %>%
  left_join(group_info_sharpe, by = "strategy")

# ------------------------------------------------------------

# ------------------------------------------------------------
benchmark_info_sharpe <- summary_sharpe %>%
  dplyr::filter(.data$strategy == "benchmark")

benchmark_upper_sharpe <- benchmark_info_sharpe$upper
benchmark_lower_sharpe <- benchmark_info_sharpe$lower
benchmark_mean_sharpe <- benchmark_info_sharpe$mean

# ------------------------------------------------------------

# ------------------------------------------------------------
plot_data_sharpe <- plot_data_sharpe %>%
  mutate(sig_higher = ifelse(lower > benchmark_upper_sharpe, TRUE, FALSE))

# ------------------------------------------------------------

# ------------------------------------------------------------
ggplot(plot_data_sharpe, aes(x = mean, y = strategy)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = sig_higher), height = 0.25) +
  geom_point(aes(color = sig_higher), size = 3) +
  geom_text(aes(x = upper + 0.1, y = strategy, label = groups), size = 5, fontface = "bold") +
  geom_vline(xintercept = benchmark_mean_sharpe, linetype = "dashed", color = "blue") +
  annotate(
    "rect",
    xmin = benchmark_lower_sharpe, xmax = benchmark_upper_sharpe,
    ymin = -Inf, ymax = Inf,
    alpha = 0.1, fill = "blue"
  ) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "gray40")) +
  labs(
    title = "Scheffé-style Comparison of Sharpe Ratios",
    subtitle = "Mean ± 95% CI",
    x = "Sharpe Ratio Mean ± 95% CI",
    y = "Strategy",
    color = "Significantly Higher\nthan Benchmark?"
  ) +
  theme_minimal(base_size = 14)

# ------------------------------------------------------------

# ------------------------------------------------------------
df_plot_sharpe <- df_long_sharpe %>%
  left_join(dplyr::select(plot_data_sharpe, strategy, sig_higher, groups), by = "strategy")

df_plot_sharpe$strategy <- factor(df_plot_sharpe$strategy,
                                  levels = plot_data_sharpe$strategy)

ggplot(df_plot_sharpe, aes(x = strategy, y = Sharpe)) +
  geom_boxplot(aes(fill = sig_higher),
               outlier.size = 0.8, outlier.alpha = 0.5) +
  geom_text(data = plot_data_sharpe,
            aes(x = strategy,
                y = max(df_plot_sharpe$Sharpe, na.rm = TRUE) + 1,
                label = groups),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "gray70")) +
  labs(title = "Sharpe Ratio Boxplot with Benchmark and Scheffé Grouping",
       x = "Strategy",
       y = "Sharpe Ratio",
       fill = "Significantly Higher\nthan Benchmark?") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")


