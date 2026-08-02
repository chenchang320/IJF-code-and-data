library(ggplot2)
library(reshape2)
library(scales) 
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


CumReturn_compound <- apply(Mean, 2, function(col) cumprod(1 + col) - 1)


df_1 <- data.frame(
  Date = hs_300_dates[1:974],  
  EWP = CumReturn_compound[, 5],
  MVP = CumReturn_compound[, 4],
  MVO = CumReturn_compound[, 6],
  SCAD_tau_0.01 = CumReturn_compound[, 1],
  SCAD_tau_0.05 = CumReturn_compound[, 2],
  SCAD_tau_0.1 = CumReturn_compound[, 3],
  HS300 = hs_300_returns[1:974]
)


df_long_1 <- reshape2::melt(df_1, id.vars = "Date", variable.name = "Series", value.name = "Value")


df_1$Date <- as.Date(df_1$Date) 

'======================================='
# Accumulated Return
ggplot(df_long_1, aes(x = Date, y = Value, color = Series)) + 
  geom_line(linewidth = 1.2, alpha = 0.9) + 
  labs(
    title = "Cumulative Return Comparison",
    subtitle = "Comparison of EWP, MVP, MVO, SCAD Strategies and HS300 Benchmark",
    x = "Date",
    y = "Cumulative Return"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 14, color = "gray40"),
    axis.title = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(
    values = c(
      "EWP" = "#1f77b4",
      "MVP" = "#ff7f0e",
      "MVO" = "#2ca02c",
      "SCAD_tau_0.01" = "#d62728",
      "SCAD_tau_0.1" = "#9467bd",
      "SCAD_tau_0.05" = "#8c564b",
      "HS300" = "#e377c2"
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.01),
    breaks = pretty(df_long_1$Value, n = 8)
  ) +
  scale_x_date(
    date_labels = "%Y-%m",       
    date_breaks = "6 months"    
   
  )



final_return <- CumReturn_compound[nrow(CumReturn_compound), ]

annual_return <- (1 + final_return)^(252 / nrow(CumReturn_compound)) - 1


rf <- 0.02  

daily_sd <- apply(Mean, 2, sd, na.rm = TRUE)

annual_sd <- daily_sd * sqrt(252)

sharpe_ratio <- (annual_return - rf) / annual_sd

max_drawdown <- apply(CumReturn_compound, 2, function(x) {
  cummax_x <- cummax(1 + x)
  drawdown <- (1 + x) / cummax_x - 1
  min(drawdown)
})




final_hs300 <- hs_300_returns[length(hs_300_returns)]

annual_hs300 <- (1 + final_hs300)^(252 / length(hs_300_returns)) - 1

annual_sd_hs300 <- 0.191

sharpe_hs300 <- (annual_hs300 - rf) / annual_sd_hs300

max_drawdown_hs300 <- min((hs_300_returns - cummax(hs_300_returns)) / (1 + cummax(hs_300_returns)))





'======================================='
# Rolling Accumulated Return
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

window_size <- 90   # or window_size <-60
roll_step <- 1     
alpha <- 0.90   # this value makes no difference here，No adjustment is required   

risk_results <- lapply(1:ncol(Mean), function(j) {
  calc_rolling_risk(Mean[, j], window_size = window_size, roll_step = roll_step, alpha = alpha)
})


r_daily_risk_results <- calc_rolling_risk(r_daily, window_size = window_size, roll_step = roll_step, alpha = alpha)



strategy_names <- c("SCAD(tau=0.01)", "SCAD(tau=0.05)", "SCAD(tau=0.1)", "MVP", "EWP", "MVO")
colors <- c("blue","darkgreen",  "red", "purple", "orange", "brown")
ltypes <- c("solid", "solid", "solid", "dashed", "dashed", "dashed")
n_win <- length(risk_results[[1]]$CVaR)


data_all <- data.frame(
  Window_Index = rep(seq_len(n_win), times = length(strategy_names)),
  CVaR = unlist(lapply(risk_results, function(x) x$CVaR)),
  cum_return = unlist(lapply(risk_results, function(x) x$cum_return)),
  RAROC = unlist(lapply(risk_results, function(x) x$RAROC)),
  Strategy = factor(rep(strategy_names, each = n_win), levels = strategy_names)
)


plot_metric <- function(metric, ylab, title) {
  ggplot(data_all, aes(x = Window_Index, y = .data[[metric]], color = Strategy, linetype = Strategy)) +
    geom_line(size = 1) +
    labs(title = title, y = ylab, x = "Rolling Window Index") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = ltypes)
}


plot_metric("cum_return", "Rolling Accumulated Return", "Rolling Accumulated Return Comparison")





'======================================='
# t test
p_values_ttest <- apply(Mean, 2, function(x) {
  t.test(x, r_daily, paired = TRUE, alternative="greater")$p.value
})
print(p_values_ttest)

# Wilcoxon 
p_values_wilcox <- apply(Mean, 2, function(x) {
  wilcox.test(x, r_daily, paired = TRUE, alternative = "greater")$p.value
})
print(p_values_wilcox)

# delta Sharpe
sharpes <- apply(Mean - r_daily, 2, function(d) {
  mean(d) / sd(d)
})
print(sharpes)

