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

data = read.csv('/Users/mac/Desktop/沪深131支个股股票数据.csv')
data = data[c(-1,-2), -1]
# 将第一列作为行索引，并移除该列
rownames(data) <- data[[1]]  # 将第一列的值设置为行索引
data <- data[,-1]  # 删除第一列，保留其他列
# 删除第一行值为0的列
data_1 <- data[, !(as.numeric(data[1, ]) == 0)]
# 确保所有数据为数值型
data_1[] <- lapply(data_1, as.numeric)
# 对每一列进行 log-diff（对数差分）计算
log_diff_data <- apply(data_1, 2, function(x) diff(log(x)))
# 2. 转换为 data.frame（如果是 matrix）
log_diff_data <- as.data.frame(log_diff_data)
# 3. 检查每一列是否含有 NaN、Inf 或 -Inf
bad_cols <- sapply(log_diff_data, function(col) {
  any(!is.finite(col))  # 检查是否存在 NaN, Inf, -Inf
})
# 4. 删除含有异常值的列
log_diff_data_clean <- log_diff_data[, !bad_cols]
log_diff_data_clean = as.matrix(log_diff_data_clean)

VAR = matrix(0,974,4)
CVAR = matrix(0,974,4)
Mean=matrix(0,974,4)

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

rq.fit.scad.with.return <- function(x, y, tau = 0.5, alpha = 3.2, lambda = 1, gamma = 10,
                                    start = "rq", beta = 0.9995, eps = 1e-06) {
  n <- length(y)
  p <- ncol(x)
  
  if (n != nrow(x))
    stop("x and y don't match n")
  if (tau < eps || tau > 1 - eps)
    stop("tau must be in (0,1)")
  if (length(lambda) == 1)
    lambda <- c(0, rep(lambda, p - 1))  # 不惩罚截距项
  if (length(lambda) != p)
    stop(paste("lambda must be either of length ", p, " or length one"))
  if (any(lambda < 0))
    stop("negative lambdas disallowed")
  
  # 惩罚矩阵 R 和 r 向量
  R <- diag(lambda, nrow = length(lambda))
  R <- R[which(lambda != 0), , drop = FALSE]
  r <- rep(0, nrow(R))
  
  # 拼接 X, Y
  X <- rbind(x, R)
  Y <- c(y, r)
  N <- length(Y)
  
  # ==== ⬇️ 添加收益项（期望收益） ====
  mean_return <- colMeans(x)  # 各资产的平均收益
  rhs <- (1 - tau) * colSums(x) + colSums(R) - gamma * mean_return
  # ==== ⬆️ 添加收益项 ====
  
  # SCAD导数
  dscad <- function(x, a = 3.7, lambda = 2) {
    lambda * sign(x) * (abs(x) <= lambda) +
      sign(x) * (a * lambda - abs(x)) / (a - 1) * (abs(x) <= a * lambda & abs(x) > lambda)
  }
  
  # 初始化
  binit <- switch(start,
                  rq = rq.fit.fnb(x, y, tau = tau)$coef[-1],
                  lasso = rq.fit.lasso(x, y, tau = tau, lambda = lambda)$coef[-1]
  )
  coef <- rep(.Machine$double.xmax, p)
  vscad <- rhs - c(0, dscad(binit, a = alpha, lambda = lambda[-1]) * sign(binit))
  it <- 0
  
  # 主循环
  while (sum(abs(binit - coef[-1])) > eps) {
    it <- it + 1
    d <- rep(1, N)
    u <- rep(1, N)
    wn <- rep(0, 10 * N)
    wn[1:N] <- c(rep((1 - tau), n), rep(0.5, nrow(R)))
    
    vrhs <- rhs - vscad
    binit <- coef[-1]
    
    z <- .Fortran("rqfnb",
                  as.integer(N), as.integer(p),
                  a = as.double(t(as.matrix(X))), c = as.double(-Y),
                  vrhs = as.double(vrhs), d = as.double(d), as.double(u),
                  beta = as.double(beta), eps = as.double(eps), wn = as.double(wn),
                  wp = double((p + 3) * p), it.count = integer(3),
                  info = integer(1)
    )
    
    coef <- -z$wp[1:p]
    vscad <- c(0, dscad(coef[2:p], a = alpha, lambda = lambda[-1]) * sign(coef[2:p]))
  }
  
  if (z$info != 0)
    stop(paste("Error info =", z$info, "in stepy2: singular design"))
  
  coefficients <- -z$wp[1:p]
  names(coefficients) <- dimnames(x)[[2]]
  residuals <- y - x %*% coefficients
  it.count <- z$it.count
  
  return(list(
    coefficients = coefficients,
    residuals = residuals,
    tau = tau,
    lambda = lambda,
    gamma = gamma,
    it = it.count
  ))
}

# 开始做循环
for (i in 957:974) {
  f = log_diff_data_clean[i:(1698+i), ]
  
  # 初始化所有存储对象
  garch <- vector("list", 131)
  error_cols <- c()
  warn_cols <- c()
  non_converge_cols <- c()
  sigma <- array(0, dim = c(1699, 131))
  mu_forecast <- array(0, dim = c(1, 131))
  sigma_forecast <- array(0, dim = c(1, 131))
  z <- array(0, dim = c(1699, 131))
  CDF <- vector("list", 131)
  inverse_CDF <- vector("list", 131)
  U <- array(0, dim = c(1699, 131))
  
  # 新增：存储EVT拟合相关结果
  gpd_low_fits <- vector("list", 131)    # 存储下尾部GPD拟合对象
  gpd_high_fits <- vector("list", 131)   # 存储上尾部GPD拟合对象
  thresholds <- matrix(NA, nrow = 131, ncol = 2, 
                       dimnames = list(NULL, c("lower", "upper"))) # 存储阈值
  exceed_low_list <- vector("list", 131)  # 存储下尾部超出量
  exceed_high_list <- vector("list", 131) # 存储上尾部超出量
  
  for (j in 1:131) {
    # GARCH模型拟合部分保持不变
    spec1 <- ugarchspec(
      variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
      distribution.model = "sged")
    
    warn_msg <- NULL
    garch[[j]] <- withCallingHandlers(
      tryCatch({
        ugarchfit(spec = spec1, data = f[, j], out.sample = 0,
                  solver = "solnp", solver.control = list(trace = 0, maxit = 5000))
      }, error = function(e) {
        error_cols <<- c(error_cols, j)
        return(NULL)
      }),
      warning = function(w) {
        warn_msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )
    
    # 错误处理部分保持不变
    if (!is.null(warn_msg)) {
      warn_cols <- c(warn_cols, j)
      if (grepl("convergence problem|optim gave code", warn_msg)) {
        non_converge_cols <- c(non_converge_cols, j)
      }
    }
    
    if (is.null(garch[[j]])) {
      sigma[, j] <- rep(NA, 1699)
      mu_forecast[j] = NA
      sigma_forecast[j] = NA
      CDF[[j]] <- NULL
      inverse_CDF[[j]] <- NULL
      U[, j] = NA
      gpd_low_fits[[j]] <- NULL
      gpd_high_fits[[j]] <- NULL
      exceed_low_list[[j]] <- NA
      exceed_high_list[[j]] <- NA
      next
    }
    
    sigma_j <- tryCatch({
      as.numeric(sigma(garch[[j]]))
    }, error = function(e) {
      return(rep(NA, 1699))
    })
    
    if (length(sigma_j) != 1699 || any(is.na(sigma_j))) {
      sigma[, j] <- rep(NA, 1699)
      mu_forecast[j] = NA
      sigma_forecast[j] = NA
      CDF[[j]] <- NULL
      inverse_CDF[[j]] <- NULL
      U[, j] = NA
      gpd_low_fits[[j]] <- NULL
      gpd_high_fits[[j]] <- NULL
      exceed_low_list[[j]] <- NA
      exceed_high_list[[j]] <- NA
      next
    }
    
    sigma[, j] <- sigma_j
    fore <- ugarchforecast(garch[[j]], n.ahead = 1)
    mu_forecast[j] <- fitted(fore)[1]
    sigma_forecast[j] <- sigma(fore)[1]
    z[, j] <- residuals(garch[[j]], standardize = TRUE)
    
    # EVT拟合部分
    u_low <- quantile(z[, j], 0.05)
    u_high <- quantile(z[, j], 0.95)
    thresholds[j, ] <- c(u_low, u_high)
    
    # 计算并存储超出量
    exceed_low <- -(z[z[, j] < u_low, j] - u_low)
    exceed_high <- z[z[, j] > u_high, j] - u_high
    exceed_low_list[[j]] <- exceed_low  # 存储下尾部超出量
    exceed_high_list[[j]] <- exceed_high # 存储上尾部超出量
    
    # 存储GPD拟合对象
    gpd_low_fits[[j]] <- tryCatch(gpd(exceed_low, threshold = 0), error = function(e) NULL)
    gpd_high_fits[[j]] <- tryCatch(gpd(exceed_high, threshold = 0), error = function(e) NULL)
    
    if (!is.null(gpd_low_fits[[j]]) && !is.null(gpd_high_fits[[j]])) {
      xi_low <- gpd_low_fits[[j]]$par.ests["xi"]
      beta_low <- gpd_low_fits[[j]]$par.ests["beta"]
      xi_high <- gpd_high_fits[[j]]$par.ests["xi"]
      beta_high <- gpd_high_fits[[j]]$par.ests["beta"]
      
      ecdf_mid <- ecdf(z[z[, j] >= u_low & z[, j] <= u_high, j])
      
      # CDF和逆CDF函数定义保持不变
      CDF[[j]] <- function(x) {
        p <- numeric(length(x))
        for (i in seq_along(x)) {
          if (x[i] < u_low) {
            y <- -(x[i] - u_low)
            p[i] <- 0.05 * pgpd(y, xi = xi_low, beta = beta_low)
          } else if (x[i] > u_high) {
            y <- x[i] - u_high
            p[i] <- 0.95 + 0.05 * pgpd(y, xi = xi_high, beta = beta_high)
          } else {
            p[i] <- 0.05 + 0.90 * ecdf_mid(x[i])
          }
        }
        return(p)
      }
      
      inverse_CDF[[j]] <- function(p) {
        q <- numeric(length(p))
        for (i in seq_along(p)) {
          if (p[i] < 0.05) {
            q[i] <- u_low - qgpd(p[i] / 0.05, xi = xi_low, beta = beta_low)
          } else if (p[i] > 0.95) {
            q[i] <- u_high + qgpd((p[i] - 0.95) / 0.05, xi = xi_high, beta = beta_high)
          } else {
            q[i] <- quantile(z[z[, j] >= u_low & z[, j] <= u_high, j], 
                             probs = (p[i] - 0.05) / 0.90)
          }
        }
        return(q)
      }
      
      U[, j] <- tryCatch(CDF[[j]](z[, j]), error = function(e) rep(NA, length(z[, j])))
    } else {
      CDF[[j]] <- NULL
      inverse_CDF[[j]] <- NULL
      U[, j] <- NA
    }
  }
  
  
  U[, non_converge_cols] <- NA
  # 找出全为 NA 的列索引
  na_cols <- which(colSums(is.na(U)) == nrow(U))
  # 获取 mu_forecast 和 sigma_forecast 中的异常列
  na_mu <- which(is.na(mu_forecast) | is.infinite(mu_forecast))
  na_sigma <- which(is.na(sigma_forecast) | is.infinite(sigma_forecast))
  
  na_cols<- unique(c(
    which(colSums(is.na(U)) == nrow(U)),  # U中全NA的列
    na_mu,                                # mu_forecast异常列
    na_sigma                              # sigma_forecast异常列
  ))
  # 删除这些列
  if (length(na_cols) == 0) {
    U_1 <- U
  } else {
    U_1 <- U[, -na_cols, drop = FALSE]
  }
  
  
  RVINE<-RVineStructureSelect(
    data = U_1,              # U 是一个矩阵，行是样本，列是变量
    familyset = 1,      # 考虑的 copula 类型（1:10 覆盖常见类型）
    type = 0,              # 0 代表 R-Vine 结构
    selectioncrit = "AIC", # 模型选择准则，可以选 BIC / AIC / loglik
    indeptest = TRUE,      # 是否进行独立性检验
    rotations = TRUE       # 是否考虑旋转 copula（更灵活）
  )
  # 进行蒙特卡洛模拟，生成 10000 个样本
  set.seed(123)  # 可复现
  Sim <- RVineSim(1000, RVINE)  # 得到新的均匀分布样本
  
  sim_real <- matrix(NA, nrow = 1000, ncol = ncol(U_1))# 将模拟的样本从 Copula 空间反变换回原始尺度（需估计边缘分布）
  if (length(na_cols) == 0) {
    inverse_CDF_1 <- inverse_CDF
    mu_forecast_1 <- mu_forecast
    sigma_forecast_1 <- sigma_forecast
    garch_1 = garch
    sigma_1 = sigma
    log_diff_data_1 = log_diff_data_clean
    f_1 = f
  } else {
    inverse_CDF_1 <- inverse_CDF[-na_cols]
    mu_forecast_1 <- mu_forecast[, -na_cols, drop = FALSE]
    sigma_forecast_1 <- sigma_forecast[, -na_cols, drop = FALSE]
    garch_1 = garch[-na_cols]
    sigma_1 = sigma[, -na_cols, drop = FALSE]
    log_diff_data_1 = log_diff_data_clean[, -na_cols, drop = FALSE]
    f_1 = f[, -na_cols, drop = FALSE]
  }
  r=array(0,dim = c(1000, ncol(U_1)))
  for (m in 1:ncol(U_1)) {
    sim_real[, m] <- inverse_CDF_1[[m]](Sim[, m])
    r[,m] = mu_forecast_1[m] + coef(garch_1[[m]])[2] * tail(f_1[, m], 1) +  coef(garch_1[[m]])[3] * tail(sigma_1[, m], 1) + sigma_forecast_1[m] * sim_real[,m]
  }
  
  write.csv(r, file = paste0("/Users/mac/Desktop/期刊结果/r/r_", i, ".csv"), row.names = FALSE)
  write.csv(as.matrix(log_diff_data_1[1699+i,]), file = paste0("/Users/mac/Desktop/期刊结果/第一次剔除后日收益/收益_", i, ".csv"))
  

  print(i)
}
