library(readxl)
library(dplyr)
library(rugarch)
library(fGarch)
library(timeDate)
library(timeSeries)
library(fBasics)
library(mistr)
library(quantmod)
library(copula)
library(VineCopula)
library(kdecopula)
library(shiny)
library(network)
library(TSP)
library(rvinecopulib)
library(quantreg)
library(ncvreg)
library(MASS)
library(evir)
library(quadprog)  
library(fPortfolio)

greedy_reduce_condition_number <- function(r_matrix, cond_threshold = 100, min_assets = 5, verbose = TRUE) {
  
  current_indices <- 1:ncol(r_matrix)
  
  compute_condition <- function(mat) {
    X_centered <- scale(mat, center = TRUE, scale = FALSE)
    kappa <- kappa(X_centered)
    return(kappa)
  }
  
  current_cond <- compute_condition(r_matrix[, current_indices])
  
  if (verbose) {
    cat(sprintf("初始资产数: %d, 初始条件数: %.2f\n", length(current_indices), current_cond))
  }
  
  while (current_cond > cond_threshold && length(current_indices) > min_assets) {
    best_cond <- Inf
    worst_idx <- NA
    
    for (j in current_indices) {
      candidate_indices <- setdiff(current_indices, j)
      cond_j <- compute_condition(r_matrix[, candidate_indices])
      if (cond_j < best_cond) {
        best_cond <- cond_j
        worst_idx <- j
      }
    }
    
    if (is.na(worst_idx)) break
    
    current_indices <- setdiff(current_indices, worst_idx)
    current_cond <- best_cond
    
    if (verbose) {
      cat(sprintf("移除资产 %d，剩余 %d 个，当前条件数: %.2f\n", worst_idx, length(current_indices), current_cond))
    }
  }
  
  return(list(
    r_reduced = r_matrix[, current_indices, drop = FALSE],
    selected_indices = current_indices,
    final_condition = current_cond
  ))
}

optimize_portfolio_nonnegative <- function(r, tau = 0.05, gamma = 1, lambda = 1, a = 3.7) {
  library(nloptr)
  
  n <- nrow(r)
  p <- ncol(r)
  
  # ---------- SCAD 惩罚函数 ----------
  scad_penalty <- function(beta, lambda, a = 3.7) {
    sapply(abs(beta), function(x) {
      if (x <= lambda) lambda * x
      else if (x <= a * lambda) (-x^2 + 2 * a * lambda * x - lambda^2) / (2 * (a - 1))
      else (lambda^2 * (a + 1)) / 2
    }) |> sum()
  }
  
  # ---------- 目标函数 ----------
  objective <- function(beta) {
    # 若 beta 中有非法值，立即返回极大惩罚值
    if (any(is.na(beta)) || any(is.nan(beta)) || any(is.infinite(beta))) return(1e10)
    
    port_returns <- as.vector(r %*% beta)
    
    # 若乘积结果中也出现非法值，返回极大惩罚值
    if (any(is.na(port_returns)) || any(is.nan(port_returns)) || any(is.infinite(port_returns))) return(1e10)
    
    VaR_est <- quantile(port_returns, tau, type = 8)
    losses <- pmax(VaR_est - port_returns, 0)
    CVaR <- VaR_est + mean(losses) / tau
    
    mean_return <- mean(port_returns)
    penalty <- scad_penalty(beta, lambda, a)
    
    return(CVaR - gamma * mean_return + penalty)
  }
  
  
  # ---------- 等式约束 ∑β = 1 ----------
  equality_constraint <- function(beta) {
    sum(beta) - 1
  }
  
  # ---------- 初始点 ----------
  init_beta <- rep(1 / p, p)
  
  # ---------- 调用 nloptr ----------
  res <- nloptr(
    x0 = init_beta,
    eval_f = objective,
    eval_g_eq = equality_constraint,
    lb = rep(0, p),  # 非负约束
    ub = rep(1, p),
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 1e-6,
      maxeval = 5000
    )
    
  )
  
  return(list(
    weights = res$solution,
    objective = res$objective,
    status = res$status,
    message = res$message
  ))
}


# num随机挑选r中的一次运行
r = as.matrix(read.csv(paste0('/Users/mac/Desktop/期刊结果/r/r_', num, ".csv")))
log_diff_data_final = read.csv(paste0("/Users/mac/Desktop/期刊结果/第一次剔除后日收益/收益_", num, ".csv"))[,-1]
result <- greedy_reduce_condition_number(r, cond_threshold = 100, min_assets = 80)
r_reduced <- result$r_reduced
selected_indices <- result$selected_indices
final_cond <- result$final_condition
log_diff_data_final <- log_diff_data_final[selected_indices, drop = FALSE]
# scad
lambda_seq <- seq(0.3, 1.5, by = 0.1)  # lambda 从0.1到1.0，每步0.1
all_weights <- list()  # 存储每个lambda下的结果

# 更改tau值对应不同的分位点策略
for (lambda_val in lambda_seq) {
  gamma_val <- 5
  min_nonzero_weights <- 25
  threshold <- 0.001
  weight <- NULL
  
  repeat {
    g_2 <- optimize_portfolio_nonnegative(r_reduced, tau = 0.05, gamma = gamma_val, lambda = lambda_val)
    weight <- g_2$weights
    nonzero_count <- sum(abs(weight) > threshold)
    
    if (nonzero_count >= min_nonzero_weights || gamma_val <= 1) {
      break
    }
    gamma_val <- gamma_val - 1
  }
  
  # 记录每个lambda下的结果，包括lambda、gamma、非零个数和权重
  all_weights[[paste0("lambda_", lambda_val)]] <- list(
    lambda = lambda_val,
    gamma = gamma_val,
    nonzero_count = nonzero_count,
    weights = weight
  )
  print(lambda_val)
}

# 可选：转换为数据框用于查看或输出
library(tibble)
weight_df <- do.call(rbind, lapply(names(all_weights), function(name) {
  data.frame(lambda = all_weights[[name]]$lambda,
             gamma = all_weights[[name]]$gamma,
             nonzero_count = all_weights[[name]]$nonzero_count,
             t(all_weights[[name]]$weights))
}))


weight_matrix <- do.call(rbind, lapply(all_weights, function(x) {
  x$weights
}))

# 给行名加上 lambda 标识，方便识别
rownames(weight_matrix) <- paste0("lambda_", sapply(all_weights, function(x) x$lambda))

# 将绝对值小于 0.01 的权重置为 0
weight_matrix[abs(weight_matrix) < 0.01] <- 0

# 提取 lambda = 1.0 对应的行
target_weights <- weight_matrix["lambda_1", ]  # 或者 weight_matrix["lambda_1.0", ]，视具体行名

# 计算非零权重数量（已经处理过小于 0.01 为 0 了）
nonzero_count <- sum(target_weights != 0)

cat("在 lambda = 1.0 时，非零权重的资产数量为：", nonzero_count, "\n")


# 假设 weight_matrix 行是 lambda，不含 lambda 列
# 行名形如 "lambda_0.1", "lambda_0.2" 等

lambda_vals <- as.numeric(gsub("lambda_", "", rownames(weight_matrix)))
weight_df <- data.frame(lambda = lambda_vals, weight_matrix)

library(tidyr)
library(dplyr)

long_data <- weight_df %>%
  pivot_longer(cols = -lambda, names_to = "variable", values_to = "value")

library(ggplot2)



ggplot(long_data, aes(x = lambda, y = value, color = variable)) +
  geom_line() +
  geom_vline(xintercept = 1.0, lty = 2, lwd = 1, col = 'black') +
  annotate('text', x = 1.0, y = -0.2, label = 'lambda* = 1.0\n(Number of remaining shares = 36)') +
  labs(
    title = "Portfolio weights solution paths (tau = 0.05)",
    x = "Lambda",
    y = "Coefficients"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.title.position = "plot",
    legend.position = "none"
  )

