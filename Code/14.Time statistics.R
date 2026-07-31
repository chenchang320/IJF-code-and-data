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

data = read.csv('/Users/mac/Desktop/GitHub-English version/data/Actual application dataset.csv')
data = data[c(-1,-2), -1]

rownames(data) <- data[[1]] 
data <- data[,-1]  

data_1 <- data[, !(as.numeric(data[1, ]) == 0)]

data_1[] <- lapply(data_1, as.numeric)

log_diff_data <- apply(data_1, 2, function(x) diff(log(x)))

log_diff_data <- as.data.frame(log_diff_data)

bad_cols <- sapply(log_diff_data, function(col) {
  any(!is.finite(col))  
})

log_diff_data_clean <- log_diff_data[, !bad_cols]
log_diff_data_clean = as.matrix(log_diff_data_clean)

VAR = matrix(0,974,4)
CVAR = matrix(0,974,4)
Mean=matrix(0,974,4)



time_stage1 <- numeric()  # stage1
time_stage2 <- numeric()  # stage2


for (i in 601:700) {

  t1_pre_rvine <- Sys.time()
  
  f = log_diff_data_clean[i:(1698+i), ]
  

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
  

  gpd_low_fits <- vector("list", 131)    
  gpd_high_fits <- vector("list", 131)   
  thresholds <- matrix(NA, nrow = 131, ncol = 2, 
                       dimnames = list(NULL, c("lower", "upper")))
  exceed_low_list <- vector("list", 131) 
  exceed_high_list <- vector("list", 131) 
  
  for (j in 1:131) {

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
    

    u_low <- quantile(z[, j], 0.05)
    u_high <- quantile(z[, j], 0.95)
    thresholds[j, ] <- c(u_low, u_high)
    

    exceed_low <- -(z[z[, j] < u_low, j] - u_low)
    exceed_high <- z[z[, j] > u_high, j] - u_high
    exceed_low_list[[j]] <- exceed_low  
    exceed_high_list[[j]] <- exceed_high
    

    gpd_low_fits[[j]] <- tryCatch(gpd(exceed_low, threshold = 0), error = function(e) NULL)
    gpd_high_fits[[j]] <- tryCatch(gpd(exceed_high, threshold = 0), error = function(e) NULL)
    
    if (!is.null(gpd_low_fits[[j]]) && !is.null(gpd_high_fits[[j]])) {
      xi_low <- gpd_low_fits[[j]]$par.ests["xi"]
      beta_low <- gpd_low_fits[[j]]$par.ests["beta"]
      xi_high <- gpd_high_fits[[j]]$par.ests["xi"]
      beta_high <- gpd_high_fits[[j]]$par.ests["beta"]
      
      ecdf_mid <- ecdf(z[z[, j] >= u_low & z[, j] <= u_high, j])
      

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

  na_cols <- which(colSums(is.na(U)) == nrow(U))

  na_mu <- which(is.na(mu_forecast) | is.infinite(mu_forecast))
  na_sigma <- which(is.na(sigma_forecast) | is.infinite(sigma_forecast))
  
  na_cols<- unique(c(
    which(colSums(is.na(U)) == nrow(U)),  
    na_mu,                               
    na_sigma                             
  ))

  if (length(na_cols) == 0) {
    U_1 <- U
  } else {
    U_1 <- U[, -na_cols, drop = FALSE]
  }
  
  t2_pre_rvine <- Sys.time()
  

  t1_rvine_to_r <- Sys.time()
  
  RVINE<-RVineStructureSelect(
    data = U_1,              
    familyset = 1,     
    type = 0,              
    selectioncrit = "AIC",
    indeptest = TRUE,      
    rotations = TRUE       
  )

  set.seed(123)  
  Sim <- RVineSim(1000, RVINE) 
  
  sim_real <- matrix(NA, nrow = 1000, ncol = ncol(U_1))
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
  
  t2_rvine_to_r <- Sys.time()
  

  duration1 <- as.numeric(difftime(t2_pre_rvine, t1_pre_rvine, units = "secs"))
  duration2 <- as.numeric(difftime(t2_rvine_to_r, t1_rvine_to_r, units = "secs"))
  time_stage1 <- c(time_stage1, duration1)
  time_stage2 <- c(time_stage2, duration2)
  

  cat("=====", i, "=====\n")
  cat("stage1：", round(duration1, 2), "seconds\n")
  cat("stage2：", round(duration2, 2), "seconds\n\n")
  
  print(i)
}
# Output the top ten average time
cat("====================\n")
cat("stage1：", round(mean(time_stage1[1:10]), 2), "seconds\n")
cat("stage1：", round(mean(time_stage2[1:10]), 2), "seconds\n")








greedy_reduce_condition_number <- function(r_matrix, cond_threshold = 100, min_assets = 5, verbose = TRUE) {
  
  current_indices <- 1:ncol(r_matrix)
  
  compute_condition <- function(mat) {
    X_centered <- scale(mat, center = TRUE, scale = FALSE)
    kappa <- kappa(X_centered)
    return(kappa)
  }
  
  current_cond <- compute_condition(r_matrix[, current_indices])
  
  if (verbose) {
    cat(sprintf("Initial asset number: %d, Initial condition number: %.2f\n", length(current_indices), current_cond))
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
      cat(sprintf("Remove assets %d，remain %d ，Current condition number: %.2f\n", worst_idx, length(current_indices), current_cond))
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
  
  
  scad_penalty <- function(beta, lambda, a = 3.7) {
    sapply(abs(beta), function(x) {
      if (x <= lambda) lambda * x
      else if (x <= a * lambda) (-x^2 + 2 * a * lambda * x - lambda^2) / (2 * (a - 1))
      else (lambda^2 * (a + 1)) / 2
    }) |> sum()
  }
  
  
  objective <- function(beta) {
    
    if (any(is.na(beta)) || any(is.nan(beta)) || any(is.infinite(beta))) return(1e10)
    
    port_returns <- as.vector(r %*% beta)
    
    
    if (any(is.na(port_returns)) || any(is.nan(port_returns)) || any(is.infinite(port_returns))) return(1e10)
    
    VaR_est <- quantile(port_returns, tau, type = 8)
    losses <- pmax(VaR_est - port_returns, 0)
    CVaR <- VaR_est + mean(losses) / tau
    
    mean_return <- mean(port_returns)
    penalty <- scad_penalty(beta, lambda, a)
    
    return(CVaR - gamma * mean_return + penalty)
  }
  
  
  
  equality_constraint <- function(beta) {
    sum(beta) - 1
  }
  
  
  init_beta <- rep(1 / p, p)
  
  
  res <- nloptr(
    x0 = init_beta,
    eval_f = objective,
    eval_g_eq = equality_constraint,
    lb = rep(0, p),  
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


time_stage3 <- numeric()  # stage3

VAR = matrix(0,974,4)
CVAR = matrix(0,974,4)
Mean=matrix(0,974,4)
for (num in 601:700) {
  t1_rvine_weight <- Sys.time()
  
  r = as.matrix(read.csv(paste0('/Users/mac/Desktop/GitHub-English version/data/r_forecast/r_', num, ".csv")))
  log_diff_data_final = read.csv(paste0("/Users/mac/Desktop/GitHub-English version/data/Actual differential return/return_", num, ".csv"))[,-1]
  result <- greedy_reduce_condition_number(r, cond_threshold = 100, min_assets = 80)
  r_reduced <- result$r_reduced
  selected_indices <- result$selected_indices
  final_cond <- result$final_condition
  log_diff_data_final <- log_diff_data_final[selected_indices, drop = FALSE]

  gamma_val <- 5
  min_nonzero_weights <- 25
  threshold <- 0.001
  weight <- NULL  
  
  repeat {

    g_2 <- optimize_portfolio_nonnegative(r_reduced, tau = 0.05, gamma = gamma_val, lambda = 1)
    weight <- g_2$weights
    

    nonzero_count <- sum(abs(weight) > threshold)
    
    if (nonzero_count >= min_nonzero_weights || gamma_val <= 1) {
      break
    }
    

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
  
  n_assets <- ncol(r_reduced)
  weights_ewp <- rep(1 / n_assets, n_assets)  
  simulated_returns_ewp <- r_reduced %*% weights_ewp  
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
  
  t2_rvine_weight <- Sys.time()
  
  duration3 <- as.numeric(difftime(t2_rvine_weight, t1_rvine_weight, units = "secs"))
  time_stage3 <- c(time_stage3, duration3)
  

  cat("=====", i, "=====\n")
  cat("stage3：", round(duration3, 2), "seconds\n")
  
  print(num)
  
}

cat("==========Output the top ten average time==========\n")
cat("stage3：", round(mean(time_stage3[1:10]), 2), "seconds\n")

