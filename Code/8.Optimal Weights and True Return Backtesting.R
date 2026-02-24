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




VAR = matrix(0,974,4)
CVAR = matrix(0,974,4)
Mean=matrix(0,974,4)
for (num in 1:974) {
  r = as.matrix(read.csv(paste0('/Users/mac/Desktop/期刊结果/r/r_', num, ".csv")))
  log_diff_data_final = read.csv(paste0("/Users/mac/Desktop/期刊结果/第一次剔除后日收益/收益_", num, ".csv"))[,-1]
  result <- greedy_reduce_condition_number(r, cond_threshold = 100, min_assets = 80)
  r_reduced <- result$r_reduced
  selected_indices <- result$selected_indices
  final_cond <- result$final_condition
  log_diff_data_final <- log_diff_data_final[selected_indices, drop = FALSE]
  # scad
  gamma_val <- 5
  min_nonzero_weights <- 25
  threshold <- 0.001
  weight <- NULL  # 初始化
  
  repeat {
    # 优化组合
    g_2 <- optimize_portfolio_nonnegative(r_reduced, tau = 0.05, gamma = gamma_val, lambda = 1)
    weight <- g_2$weights
    
    # 判断非零权重数量
    nonzero_count <- sum(abs(weight) > threshold)
    
    if (nonzero_count >= min_nonzero_weights || gamma_val <= 1) {
      break
    }
    
    # 减小 gamma 并重试
    gamma_val <- gamma_val - 1
  }
  sim_returns_scad=r_reduced%*%weight
  VaR_95_scad <- quantile(sim_returns_scad, probs = 0.05)
  CVaR_95_scad <- mean(sim_returns_scad[sim_returns_scad <= VaR_95_scad])
  
  
  #min-var
  mo=minvariancePortfolio(as.timeSeries(r_reduced), spec = portfolioSpec(), constraints = "LongOnly")
  qo=getWeights(mo)
  weights_minvar=matrix(qo,ncol(r_reduced),1)
  sim_returns_minvar <- r_reduced %*% weights_minvar
  VaR_95_minvar <- quantile(sim_returns_minvar, probs = 0.05)
  CVaR_95_minvar <- mean(sim_returns_minvar[sim_returns_minvar <= VaR_95_minvar])
  
  #mean-var
  Spec = portfolioSpec()
  setTargetReturn(Spec) = mean(colMeans(r_reduced))
  mo_1=efficientPortfolio(as.timeSeries(r_reduced), Spec , constraints = "LongOnly")
  qo_1=getWeights(mo_1)
  weights_meanvar=matrix(qo_1,ncol(r_reduced),1)
  sim_returns_meanvar <- r_reduced %*% weights_meanvar
  VaR_95_meanvar <- quantile(sim_returns_meanvar, probs = 0.05)
  CVaR_95_meanvar <- mean(sim_returns_meanvar[sim_returns_meanvar <= VaR_95_meanvar])
  
  #ewp
  n_assets <- ncol(r_reduced)
  weights_ewp <- rep(1 / n_assets, n_assets)  
  simulated_returns_ewp <- r_reduced %*% weights_ewp  # 得到10000个组合收益
  VaR_95_ewp <- quantile(simulated_returns_ewp, probs = 0.05)
  CVaR_95_ewp <- mean(simulated_returns_ewp[simulated_returns_ewp <= VaR_95_ewp])
  
  VAR[num,1]=VaR_95_scad
  VAR[num,2]=VaR_95_minvar
  VAR[num,3]=VaR_95_ewp
  VAR[num,4]=VaR_95_meanvar
  CVAR[num,1] =CVaR_95_scad
  CVAR[num,2]=CVaR_95_minvar
  CVAR[num,3]=CVaR_95_ewp
  CVAR[num,4]=CVaR_95_meanvar
  Mean[num,1]=log_diff_data_final%*%weight
  Mean[num,2]=log_diff_data_final%*%weights_minvar
  Mean[num,3]=log_diff_data_final%*%weights_ewp
  Mean[num,4]=log_diff_data_final%*%weights_meanvar
  print(num)
  write.csv(VAR, file = '/Users/mac/Desktop/期刊结果/表格tau=0.05/VAR.csv', row.names = FALSE)
  write.csv(CVAR, file = '/Users/mac/Desktop/期刊结果/表格tau=0.05/CVAR.csv', row.names = FALSE)
  write.csv(Mean, file = '/Users/mac/Desktop/期刊结果/表格tau=0.05/Mean.csv', row.names = FALSE)
}






#tau=0.01和tau=0.1放在一起计算
VAR = as.matrix(read.csv('/Users/mac/Desktop/期刊结果/表格tau=0.01/VAR.csv'))
CVAR = as.matrix(read.csv('/Users/mac/Desktop/期刊结果/表格tau=0.01/CVAR.csv'))
Mean = as.matrix(read.csv('/Users/mac/Desktop/期刊结果/表格tau=0.01/Mean.csv'))
for (num in 1:974) {
  r = as.matrix(read.csv(paste0('/Users/mac/Desktop/期刊结果/r/r_', num, ".csv")))
  log_diff_data_final = read.csv(paste0("/Users/mac/Desktop/期刊结果/第一次剔除后日收益/收益_", num, ".csv"))[,-1]
  result <- greedy_reduce_condition_number(r, cond_threshold = 100, min_assets = 80)
  r_reduced <- result$r_reduced
  selected_indices <- result$selected_indices
  final_cond <- result$final_condition
  log_diff_data_final <- log_diff_data_final[selected_indices, drop = FALSE]
  # scad
  gamma_val <- 5
  min_nonzero_weights <- 25
  threshold <- 0.001
  weight <- NULL  # 初始化
  
  repeat {
    # 优化组合
    g_2 <- optimize_portfolio_nonnegative(r_reduced, tau = 0.01, gamma = gamma_val, lambda = 1)
    weight <- g_2$weights
    
    # 判断非零权重数量
    nonzero_count <- sum(abs(weight) > threshold)
    
    if (nonzero_count >= min_nonzero_weights || gamma_val <= 1) {
      break
    }
    
    # 减小 gamma 并重试
    gamma_val <- gamma_val - 1
  }
  sim_returns_scad=r_reduced%*%weight
  VaR_95_scad <- quantile(sim_returns_scad, probs = 0.05)
  CVaR_95_scad <- mean(sim_returns_scad[sim_returns_scad <= VaR_95_scad])
  
  gamma_val_1 <- 5
  min_nonzero_weights_1 <- 25
  threshold_1 <- 0.001
  weight_1 <- NULL  # 初始化
  repeat {
    # 优化组合
    g_3 <- optimize_portfolio_nonnegative(r_reduced, tau = 0.1, gamma = gamma_val_1, lambda = 1)
    weight_1 <- g_3$weights
    
    # 判断非零权重数量
    nonzero_count_1 <- sum(abs(weight_1) > threshold_1)
    
    if (nonzero_count_1 >= min_nonzero_weights_1 || gamma_val_1 <= 1) {
      break
    }
    
    # 减小 gamma 并重试
    gamma_val_1 <- gamma_val_1 - 1
  }
  sim_returns_scad_1 = r_reduced%*%weight_1
  VaR_95_scad_1 <- quantile(sim_returns_scad_1, probs = 0.05)
  CVaR_95_scad_1 <- mean(sim_returns_scad_1[sim_returns_scad_1 <= VaR_95_scad_1])
  
  VAR[num,1]=VaR_95_scad
  VAR[num,2]=VaR_95_scad_1
  CVAR[num,1] =CVaR_95_scad
  CVAR[num,2]=CVaR_95_scad_1 
  Mean[num,1]=log_diff_data_final%*%weight
  Mean[num,2]=log_diff_data_final%*%weight_1
  print(num)
  write.csv(VAR, file = '/Users/mac/Desktop/期刊结果/表格tau=0.01/VAR.csv', row.names = FALSE)
  write.csv(CVAR, file = '/Users/mac/Desktop/期刊结果/表格tau=0.01/CVAR.csv', row.names = FALSE)
  write.csv(Mean, file = '/Users/mac/Desktop/期刊结果/表格tau=0.01/Mean.csv', row.names = FALSE)
}







#最后一次循环画的权重稀疏图
log_diff_name = read.csv(paste0("/Users/mac/Desktop/期刊结果/第一次剔除后日收益/收益_", num, ".csv"))[,1]
log_diff_name <- log_diff_name[selected_indices, drop = FALSE]
# 去掉前缀 X
group <- gsub("^X", "", log_diff_name)

# 假设 weight 是你的权重向量
values <- weight
values[abs(values) < 0.01] <- 0  # 小于0.001的权重归零


# 构建数据框
data_ggp <- data.frame(group = as.factor(group), values = values)

# 绘制柱状图
library(ggplot2)
ggplot(data_ggp, aes(x = group, y = values)) +
  geom_bar(stat = "identity") +
  labs(title = "Optimal portfolio weights (tau = 0.01)", 
       x = "Stock Code", 
       y = "Weight") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.title.position = "plot")

# ==================================

# 假设 weight 是你的权重向量
values <- weight_1
values[abs(values) < 0.01] <- 0  # 小于0.001的权重归零


# 构建数据框
data_ggp <- data.frame(group = as.factor(group), values = values)

# 绘制柱状图
library(ggplot2)
ggplot(data_ggp, aes(x = group, y = values)) +
  geom_bar(stat = "identity") +
  labs(title = "Optimal portfolio weights (tau = 0.1)", 
       x = "Stock Code", 
       y = "Weight") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.title.position = "plot")
