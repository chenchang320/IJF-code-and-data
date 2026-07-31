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

Mean_0.01 = Mean[,c(1,4,5,6)]
Mean_0.05 = Mean[,c(2,4,5,6)]
Mean_0.1 = Mean[,c(3,4,5,6)]


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

# adjust window_size、alpha、Mean_， calc_rolling_risk gets different result with different input
window_size <- 90   # or 60
roll_step <- 1      
alpha <- 0.90      # or 0.99、or 0.95

risk_results <- lapply(1:ncol(Mean_0.1), function(j) {
  calc_rolling_risk(Mean_0.1[, j], window_size = window_size, roll_step = roll_step, alpha = alpha)
})


strategy_names <- c("SCAD(tau=0.1)", "MVP", "EWP", "MVO")
colors <- c("red", "purple", "orange", "brown")
ltypes <- c("solid", "dashed", "dashed", "dashed")


n_win <- length(risk_results[[1]]$CVaR)

data_all <- data.frame(
  Window_Index = rep(seq_len(n_win), times = length(strategy_names)),
  CVaR         = unlist(lapply(risk_results, function(x) x$CVaR)),
  cum_return   = unlist(lapply(risk_results, function(x) x$cum_return)),
  RAROC        = unlist(lapply(risk_results, function(x) x$RAROC)),
  Strategy     = factor(rep(strategy_names, each = n_win), levels = strategy_names)
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


plot_metric("CVaR", "Rolling CVaR", "Rolling CVaR Comparison")
plot_metric("RAROC", "RAROC", "RAROC Comparison")


'========================================'
# RAROC statistics
raroc_summary <- data_all %>%
  group_by(Strategy) %>%
  summarise(
    Mean_RAROC   = mean(RAROC, na.rm = TRUE),
    Max_RAROC    = max(RAROC, na.rm = TRUE),
    Min_RAROC    = min(RAROC, na.rm = TRUE),
    SD_RAROC     = sd(RAROC, na.rm = TRUE),
    Median_RAROC = median(RAROC, na.rm = TRUE),
    Positive_Ratio = mean(RAROC > 0, na.rm = TRUE),
    N = n()
  ) %>%
  arrange(desc(Mean_RAROC))



'========================================'
n_win <- length(risk_results[[1]]$CVaR)
data_all <- data.frame(
  Window_Index = rep(seq_len(n_win), times = length(strategy_names)),
  CVaR         = unlist(lapply(risk_results, function(x) x$CVaR)),
  cum_return   = unlist(lapply(risk_results, function(x) x$cum_return)),
  RAROC        = unlist(lapply(risk_results, function(x) x$RAROC)),
  Sharpe       = unlist(lapply(risk_results, function(x) x$Sharpe)), 
  Strategy     = factor(rep(strategy_names, each = n_win), levels = strategy_names)
)


# Sharpe statistics
sharpe_summary <- data_all %>%
  group_by(Strategy) %>%
  summarise(
    Mean_Sharpe    = mean(Sharpe, na.rm = TRUE),
    Max_Sharpe     = max(Sharpe, na.rm = TRUE),
    Min_Sharpe     = min(Sharpe, na.rm = TRUE),
    SD_Sharpe      = sd(Sharpe, na.rm = TRUE),
    Median_Sharpe  = median(Sharpe, na.rm = TRUE),
    Positive_Ratio = mean(Sharpe > 0, na.rm = TRUE),
    N              = n()
  ) %>%
  arrange(desc(Mean_Sharpe))


