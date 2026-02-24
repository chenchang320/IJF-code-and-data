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




Mean_1 = as.matrix(read.csv('/Users/mac/Desktop/期刊结果/表格tau=0.01/Mean.csv'))

Mean_2 = as.matrix(read.csv('/Users/mac/Desktop/期刊结果/表格tau=0.05/Mean.csv'))

Mean = cbind(Mean_1, Mean_2)

col_names <- c("SCAD(tau=0.01)", 
               "SCAD(tau=0.1)", 
               "SCAD(tau=0.05)", 
               "MVP", "EWP", "MVO")

# 赋予列名
colnames(Mean) <- col_names

Mean_0.01 = Mean[,c(1,4,5,6)]
Mean_0.05 = Mean[,c(3,4,5,6)]
Mean_0.1 = Mean[,c(2,4,5,6)]

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

#alpha与对应的策略分为点取一致
window_size <- 90   # 例如 30 天持有期
roll_step <- 1      # 每天滚动一次窗口
alpha <- 0.9       # 5% 下行风险水平

risk_results <- lapply(1:ncol(Mean_0.1), function(j) {
  calc_rolling_risk(Mean_0.1[, j], window_size = window_size, roll_step = roll_step, alpha = alpha)
})

# 计算 r_daily 的滚动风险指标
r_daily_risk_results <- calc_rolling_risk(r_daily, window_size = window_size, roll_step = roll_step, alpha = alpha)






library(dplyr)
library(tidyr)
library(agricolae)
library(ggplot2)

# ------------------------------------------------------------
# 1️⃣ 提取 RAROC（策略 + benchmark）
# ------------------------------------------------------------
raroc_values <- lapply(risk_results, function(x) x$RAROC)
benchmark_raroc <- r_daily_risk_results$RAROC

# 👉 只保留 Mean_0.01 对应的列
df_raroc <- do.call(cbind, raroc_values)
colnames(df_raroc) <- c("SCAD(tau=0.1)", "MVP", "EWP", "MVO")
df_raroc <- cbind(df_raroc, benchmark = benchmark_raroc)

# 转成长表
df_long <- df_raroc %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "strategy",
               values_to = "RAROC") %>%
  drop_na()

# ------------------------------------------------------------
# 2️⃣ ANOVA + Scheffé
# ------------------------------------------------------------
anova_result <- aov(RAROC ~ strategy, data = df_long)
summary(anova_result)

scheffe_result <- scheffe.test(anova_result, "strategy")
group_info <- scheffe_result$groups %>%
  tibble::rownames_to_column("strategy")

# ------------------------------------------------------------
# 3️⃣ 平均值 ± 95% CI
# ------------------------------------------------------------
summary_raroc <- df_long %>%
  group_by(strategy) %>%
  summarise(
    mean  = mean(RAROC, na.rm = TRUE),
    sd    = sd(RAROC, na.rm = TRUE),
    n     = n(),
    se    = sd / sqrt(n),
    lower = mean - qt(0.95, df = n - 1) * se,
    upper = mean + qt(0.95, df = n - 1) * se
  ) %>%
  arrange(mean) %>%
  mutate(strategy = factor(strategy, levels = strategy))

plot_data_raroc <- summary_raroc %>%
  left_join(group_info, by = "strategy")

# ------------------------------------------------------------
# 4️⃣ benchmark CI
# ------------------------------------------------------------
benchmark_info <- summary_raroc %>%
  dplyr::filter(.data$strategy == "benchmark")
benchmark_upper <- benchmark_info$upper
benchmark_lower <- benchmark_info$lower
benchmark_mean  <- benchmark_info$mean

# ------------------------------------------------------------
# 5️⃣ 标记显著大于 benchmark 的策略
# ------------------------------------------------------------
plot_data_raroc <- plot_data_raroc %>%
  mutate(sig_higher = ifelse(lower > benchmark_upper, TRUE, FALSE))

# ------------------------------------------------------------
# 6️⃣ 【均值 ± CI 图】
# ------------------------------------------------------------
ggplot(plot_data_raroc, aes(x = mean, y = strategy)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = sig_higher),
                 height = 0.25) +
  geom_point(aes(color = sig_higher), size = 3) +
  geom_text(aes(x = upper + 0.001, y = strategy, label = groups),
            size = 5, fontface = "bold") +
  geom_vline(xintercept = benchmark_mean, linetype = "dashed", color = "blue") +
  annotate("rect",
           xmin = benchmark_lower, xmax = benchmark_upper,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "blue") +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "gray40")) +
  labs(title = "Scheffé-style Comparison of RAROC Means",
       x = "RAROC Mean ± 95% CI",
       y = "Strategy",
       color = "Significantly Higher\nthan Benchmark?") +
  theme_minimal(base_size = 14)

# ------------------------------------------------------------
# 7️⃣ 【箱线图】Scheffé 分组
# ------------------------------------------------------------
df_plot_raroc <- df_long %>%
  left_join(dplyr::select(plot_data_raroc, strategy, sig_higher, groups), by = "strategy") %>%
  mutate(strategy = factor(strategy, levels = plot_data_raroc$strategy))

ggplot(df_plot_raroc, aes(x = strategy, y = RAROC)) +
  geom_boxplot(aes(fill = sig_higher),
               outlier.size = 0.8, outlier.alpha = 0.5) +
  geom_text(data = plot_data_raroc,
            aes(x = strategy,
                y = max(df_plot_raroc$RAROC, na.rm = TRUE) + 0.001,
                label = groups),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "gray70")) +
  labs(title = "RAROC Boxplot with Benchmark and Scheffé Grouping",
       x = "Strategy",
       y = "RAROC",
       fill = "Significantly Higher\nthan Benchmark?") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")



