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
library(R.utils)
library(callr)


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


run_scad_repeat_timeout <- function(
    r_reduced,
    tau,
    timeout_seconds = 90,
    gamma_start = 5,
    min_nonzero_weights = 25,
    threshold = 0.001
) {
  
  result <- tryCatch(
    
    callr::r(
      func = function(
    r_reduced,
    tau,
    gamma_start,
    min_nonzero_weights,
    threshold,
    optimize_fun
      ) {
        
        gamma_val <- gamma_start
        weight <- NULL
        nonzero_count <- 0
        
        repeat {
          
          g <- optimize_fun(
            r = r_reduced,
            tau = tau,
            gamma = gamma_val,
            lambda = 1
          )
          
          weight <- g$weights
          
          if (
            is.null(weight) ||
            anyNA(weight) ||
            any(!is.finite(weight))
          ) {
            stop("Optimize without returning valid weights")
          }
          
          nonzero_count <- sum(
            abs(weight) > threshold
          )
          
          if (
            nonzero_count >= min_nonzero_weights ||
            gamma_val <= 1
          ) {
            break
          }
          
          gamma_val <- gamma_val - 1
        }
        
        list(
          success = TRUE,
          weights = weight,
          gamma = gamma_val,
          nonzero_count = nonzero_count
        )
      },
    
    args = list(
      r_reduced = r_reduced,
      tau = tau,
      gamma_start = gamma_start,
      min_nonzero_weights = min_nonzero_weights,
      threshold = threshold,
      optimize_fun = optimize_portfolio_nonnegative
    ),
    

    timeout = timeout_seconds,
    
    libpath = .libPaths()
    ),
    
    error = function(e) {
      
      list(
        success = FALSE,
        weights = NULL,
        gamma = NA_real_,
        nonzero_count = NA_integer_,
        message = conditionMessage(e)
      )
    }
  )
  
  result
}



Mean=matrix(0,974,6)
for (num in 1:974) {
  r = as.matrix(read.csv(paste0('/Users/mac/Desktop/GitHub-English version/Intermediate data/r_forecast/r_', num, ".csv")))
  log_diff_data_final = read.csv(paste0("/Users/mac/Desktop/GitHub-English version/Intermediate data/Actual differential return/return_", num, ".csv"))[,-1]
  result <- greedy_reduce_condition_number(r, cond_threshold = 100, min_assets = 80)
  r_reduced <- result$r_reduced
  selected_indices <- result$selected_indices
  final_cond <- result$final_condition
  log_diff_data_final <- log_diff_data_final[selected_indices, drop = FALSE]
  

  
  scad_failed <- FALSE
  
  
  # ---------------- tau = 0.01 ----------------
  result_scad_1 <- run_scad_repeat_timeout(
    r_reduced = r_reduced,
    tau = 0.01,
    timeout_seconds = 300,
    gamma_start = 5,
    min_nonzero_weights = 25,
    threshold = 0.001
  )
  
  if (!isTRUE(result_scad_1$success)) {
    
    cat(
      sprintf(
        "num = %d，tau = 0.01 excess or error ，Skip the current num：%s\n",
        num,
        result_scad_1$message
      )
    )
    
    scad_failed <- TRUE
    
  } else {
    
    weight_1 <- result_scad_1$weights
  }
  
  
  # ---------------- tau = 0.05 ----------------
  if (!scad_failed) {
    
    result_scad_2 <- run_scad_repeat_timeout(
      r_reduced = r_reduced,
      tau = 0.05,
      timeout_seconds = 300,
      gamma_start = 5,
      min_nonzero_weights = 25,
      threshold = 0.001
    )
    
    if (!isTRUE(result_scad_2$success)) {
      
      cat(
        sprintf(
          "num = %d，tau = 0.05 excess or error ，Skip the current num：%s\n",
          num,
          result_scad_2$message
        )
      )
      
      scad_failed <- TRUE
      
    } else {
      
      weight_2 <- result_scad_2$weights
    }
  }
  
  
  # ---------------- tau = 0.10 ----------------
  if (!scad_failed) {
    
    result_scad_3 <- run_scad_repeat_timeout(
      r_reduced = r_reduced,
      tau = 0.10,
      timeout_seconds = 300,
      gamma_start = 5,
      min_nonzero_weights = 25,
      threshold = 0.001
    )
    
    if (!isTRUE(result_scad_3$success)) {
      
      cat(
        sprintf(
          "num = %d，tau = 0.10 excess or error ，Skip the current num：%s\n",
          num,
          result_scad_3$message
        )
      )
      
      scad_failed <- TRUE
      
    } else {
      
      weight_3 <- result_scad_3$weights
    }
  }
  
  

  if (scad_failed) {
    next
  }
  
  
  #min-var
  mo=minvariancePortfolio(as.timeSeries(r_reduced), spec = portfolioSpec(), constraints = "LongOnly")
  qo=getWeights(mo)
  weights_minvar=matrix(qo,ncol(r_reduced),1)
  
  #mean-var
  Spec = portfolioSpec()
  setTargetReturn(Spec) = mean(colMeans(r_reduced))
  mo_1=efficientPortfolio(as.timeSeries(r_reduced), Spec , constraints = "LongOnly")
  qo_1=getWeights(mo_1)
  weights_meanvar=matrix(qo_1,ncol(r_reduced),1)
  
  #ewp
  n_assets <- ncol(r_reduced)
  weights_ewp <- rep(1 / n_assets, n_assets)  
  
  Mean[num,1]=(exp(log_diff_data_final) - 1)%*%weight_1
  Mean[num,2]=(exp(log_diff_data_final) - 1)%*%weight_2
  Mean[num,3]=(exp(log_diff_data_final) - 1)%*%weight_3
  Mean[num,4]=(exp(log_diff_data_final) - 1)%*%weights_minvar
  Mean[num,5]=(exp(log_diff_data_final) - 1)%*%weights_ewp
  Mean[num,6]=(exp(log_diff_data_final) - 1)%*%weights_meanvar
  print(num)
  
  write.csv(Mean, file = '/Users/mac/Desktop/GitHub-English version/Intermediate data/Daily return.csv', row.names = FALSE)
}






'============================================'
# Take the last "num" as an example
log_diff_name = read.csv(paste0("/Users/mac/Desktop/GitHub-English version/data/Actual differential return/return_", num, ".csv"))[,1]
log_diff_name <- log_diff_name[selected_indices, drop = FALSE]

group <- gsub("^X", "", log_diff_name)

# ==================================
values <- weight_1
values[abs(values) < 0.01] <- 0  



data_ggp <- data.frame(group = as.factor(group), values = values)


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


values <- weight_2
values[abs(values) < 0.01] <- 0  



data_ggp <- data.frame(group = as.factor(group), values = values)


library(ggplot2)
ggplot(data_ggp, aes(x = group, y = values)) +
  geom_bar(stat = "identity") +
  labs(title = "Optimal portfolio weights (tau = 0.05)", 
       x = "Stock Code", 
       y = "Weight") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.title.position = "plot")

# ==================================


values <- weight_3
values[abs(values) < 0.01] <- 0  



data_ggp <- data.frame(group = as.factor(group), values = values)


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
