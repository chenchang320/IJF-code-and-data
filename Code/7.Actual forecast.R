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



for (i in 1:974) {
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
  
  write.csv(r, file = paste0("/Users/mac/Desktop/GitHub-English version/data/r_forecast/r_", i, ".csv"), row.names = FALSE)
  write.csv(as.matrix(log_diff_data_1[1699+i,]), file = paste0("/Users/mac/Desktop/GitHub-English version/data/Actual differential return/return_", i, ".csv"))
  
  
  print(i)
}
