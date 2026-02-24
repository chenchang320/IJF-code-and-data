library(ggplot2)
library(reshape2)
library(scales)  # 用于格式化纵轴

hs_300 = read_xlsx('/Users/mac/Desktop/沪深300指数收益.xlsx')
# 假设第二列是沪深300日收益率
hs_300_returns <- as.numeric(hs_300[[2]])/100

hs_300_dates <- as.Date(hs_300[[1]])

r_daily <- numeric(length(hs_300_returns))
r_daily[1] <- hs_300_returns[1]  # 第一日收益假定等于第一日累计收益

# 从第二天开始逐步恢复日收益
for (t in 2:length(hs_300_returns)) {
  r_daily[t] <- (1 + hs_300_returns[t]) / (1 + hs_300_returns[t - 1]) - 1
}



VAR_1 = as.matrix(read.csv('/Users/mac/Desktop/期刊结果/表格tau=0.01/VAR.csv'))
CVAR_1 = as.matrix(read.csv('/Users/mac/Desktop/期刊结果/表格tau=0.01/CVAR.csv'))
Mean_1 = as.matrix(read.csv('/Users/mac/Desktop/期刊结果/表格tau=0.01/Mean.csv'))
VAR_2 = as.matrix(read.csv('/Users/mac/Desktop/期刊结果/表格tau=0.05/VAR.csv'))
CVAR_2 = as.matrix(read.csv('/Users/mac/Desktop/期刊结果/表格tau=0.05/CVAR.csv'))
Mean_2 = as.matrix(read.csv('/Users/mac/Desktop/期刊结果/表格tau=0.05/Mean.csv'))
VAR = cbind(VAR_1, VAR_2)
CVAR = cbind(CVAR_1, CVAR_2)
Mean = cbind(Mean_1, Mean_2)

col_names <- c("SCAD(tau=0.01)", 
               "SCAD(tau=0.1)", 
               "SCAD(tau=0.05)", 
               "MVP", "EWP", "MVO")

# 赋予列名
colnames(VAR)  <- col_names
colnames(CVAR) <- col_names
colnames(Mean) <- col_names


CumReturn_compound <- apply(Mean, 2, function(col) cumprod(1 + col) - 1)

# 构造数据框时使用日期作为X轴
df_1 <- data.frame(
  Date = hs_300_dates[1:974],  # 使用日期作为X轴，确保长度一致
  EWP = CumReturn_compound[, 5],
  MVP = CumReturn_compound[, 4],
  MVO = CumReturn_compound[, 6],
  SCAD_tau_0.01 = CumReturn_compound[, 1],
  SCAD_tau_0.1 = CumReturn_compound[, 2],
  SCAD_tau_0.05 = CumReturn_compound[, 3],
  HS300 = hs_300_returns[1:974]
)

# 宽转长格式
df_long_1 <- reshape2::melt(df_1, id.vars = "Date", variable.name = "Series", value.name = "Value")

# 直接使用数据中的日期（确保已经是 Date 类型）
df_1$Date <- as.Date(df_1$Date)  # 如果尚未转换，确保转为 Date 类型

# 绘图（自动适应数据中的日期范围）
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
    date_labels = "%Y-%m",       # 显示格式：YYYY-MM
    date_breaks = "6 months"     # 每 6 个月显示一个标签
    # 不再需要 limits = date_range，自动使用数据范围
  )


# 最终累计收益
final_return <- CumReturn_compound[nrow(CumReturn_compound), ]
# 计算年化收益
annual_return <- (1 + final_return)^(252 / nrow(CumReturn_compound)) - 1

# 无风险利率
rf <- 0.02  
# 每列的日波动率
daily_sd <- apply(Mean, 2, sd, na.rm = TRUE)
# 年化波动率
annual_sd <- daily_sd * sqrt(252)
# 夏普比率
sharpe_ratio <- (annual_return - rf) / annual_sd
# 最大回撤
max_drawdown <- apply(CumReturn_compound, 2, function(x) {
  cummax_x <- cummax(1 + x)
  drawdown <- (1 + x) / cummax_x - 1
  min(drawdown)
})



# 1. 最终累计收益（直接取最后一个值）
final_hs300 <- hs_300_returns[length(hs_300_returns)]
# 2. 计算年化收益（需从累计收益反推日收益率，再年化）
annual_hs300 <- (1 + final_hs300)^(252 / length(hs_300_returns)) - 1
# 3. 计算波动率（需先还原日收益率）
annual_sd_hs300 <- 0.191
# 4. 夏普比率（无风险利率rf = 0.02）
sharpe_hs300 <- (annual_hs300 - rf) / annual_sd_hs300
# 5. 最大回撤（直接基于累计收益计算）
max_drawdown_hs300 <- min((hs_300_returns - cummax(hs_300_returns)) / (1 + cummax(hs_300_returns)))







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
    Sharpe       = NA,   # ✅ 新增：Sharpe
    RAROC_7d_MA  = NA    # （可后期再计算）
  )
  
  for (i in seq_along(roll_indices)) {
    idx_start <- roll_indices[i]
    idx_end   <- idx_start + window_size - 1
    r_win     <- returns[idx_start:idx_end]
    
    # 累计收益
    cum_r <- prod(1 + r_win) - 1
    
    # VaR, CVaR
    port_loss <- -r_win
    VaR_val   <- quantile(port_loss, alpha, type = 8)
    CVaR_val  <- mean(port_loss[port_loss >= VaR_val])
    
    # 最大回撤
    wealth_index <- cumprod(1 + r_win)
    cum_max      <- cummax(wealth_index)
    drawdowns    <- (wealth_index - cum_max) / cum_max
    max_dd       <- min(drawdowns, na.rm = TRUE)
    
    # RAROC
    raroc_val <- ifelse(CVaR_val != 0, cum_r / CVaR_val, NA)
    
    # Sharpe
    excess_r <- r_win - rf_daily
    sharpe_val <- ifelse(sd(excess_r, na.rm = TRUE) != 0,
                         mean(excess_r, na.rm = TRUE) / sd(excess_r, na.rm = TRUE),
                         NA)
    
    # 存结果
    results$cum_return[i]   <- cum_r
    results$VaR[i]          <- VaR_val
    results$CVaR[i]         <- CVaR_val
    results$max_drawdown[i] <- max_dd
    results$RAROC[i]        <- raroc_val
    results$Sharpe[i]       <- sharpe_val
  }
  
  return(results)
}

window_size <- 90   # 例如 30 天持有期
roll_step <- 1      # 每天滚动一次窗口
alpha <- 0.90       # 5% 下行风险水平

risk_results <- lapply(1:ncol(Mean), function(j) {
  calc_rolling_risk(Mean[, j], window_size = window_size, roll_step = roll_step, alpha = alpha)
})

# 计算 r_daily 的滚动风险指标
r_daily_risk_results <- calc_rolling_risk(r_daily, window_size = window_size, roll_step = roll_step, alpha = alpha)


# ---------- 公共设置 ----------
strategy_names <- c("SCAD(tau=0.01)", "SCAD(tau=0.1)", "SCAD(tau=0.05)", "MVP", "EWP", "MVO")
colors <- c("blue", "red", "darkgreen", "purple", "orange", "brown")
ltypes <- c("solid", "solid", "solid", "dashed", "dashed", "dashed")
n_win <- length(risk_results[[1]]$CVaR)

# ---------- 构建统一数据框 ----------
data_all <- data.frame(
  Window_Index = rep(seq_len(n_win), times = length(strategy_names)),
  CVaR = unlist(lapply(risk_results, function(x) x$CVaR)),
  cum_return = unlist(lapply(risk_results, function(x) x$cum_return)),
  RAROC = unlist(lapply(risk_results, function(x) x$RAROC)),
  Strategy = factor(rep(strategy_names, each = n_win), levels = strategy_names)
)

# ---------- 画图函数 ----------
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

# ---------- 生成图 ----------
plot_metric("cum_return", "Rolling Accumulated Return", "Rolling Accumulated Return Comparison")






# 超额收益均值进行配对 t 检验
p_values_ttest <- apply(Mean, 2, function(x) {
  t.test(x, r_daily, paired = TRUE, alternative="greater")$p.value
})
print(p_values_ttest)

# Wilcoxon 符号秩检验：检验超额收益中位数是否显著大于0
p_values_wilcox <- apply(Mean, 2, function(x) {
  wilcox.test(x, r_daily, paired = TRUE, alternative = "greater")$p.value
})
print(p_values_wilcox)

# 超额收益的“夏普比率”
sharpes <- apply(Mean - r_daily, 2, function(d) {
  mean(d) / sd(d)
})
print(sharpes)

