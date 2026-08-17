library(readxl)
library(forecast)
library(tseries)
library(MTS)
library(DistributionUtils)
library(dplyr)
library(rugarch)
library(fGarch)
library(timeDate)
library(timeSeries)
library(fBasics)
library(mistr)
library(quantmod)
data1=read_excel('/Users/mac/Desktop/GitHub-English version/Raw data/Ten industry dataset.xlsx')
select_data=data1[c(1:2674),c(-1)]
colnames(select_data)=c(1:10)


log.return=function(x){
  diff(log(x))
}

log_return=matrix(0,nrow = 2673,ncol = 10)
for (d in c(1:10)){
  e=select_data[[d]]
  h=as.numeric(e)
  log_return[,d]=log.return(h)
  
}  
f=log_return


data_outline=function(x){
  m=mean(x)
  d=max(x)
  xd=min(x)
  s=sd(x)
  kur=kurtosis(x)
  ske=skewness(x)
  R=max(x)-min(x)
  data.frame( Mean=m,max=d,min=xd,std_dev=s,Skewness=ske, Kurtosis=kur, R=R)
}

describe_statistic=data.frame()
zhengtai_test=list()
danweigen_test=list()
before_LB_test=list()
before_arch=list()
for (i in  1:10){
  describe_statistic=rbind(describe_statistic,data_outline(f[,i]))
  zhengtai_test[[i]]=jarque.bera.test(f[,i]) 
  danweigen_test[[i]]=adf.test(f[,i])
  before_LB_test[[i]]=Box.test(f[,i],10,type='Ljung-Box')
  before_arch[[i]]=archTest(residuals(auto.arima(f[,i])),lag=5)
}



garch = list()
sigma = array(0, dim = c(2673, 10))
z = array(0, dim = c(2673, 10))
fit = list()
CDF = list()
inverse_CDF = list()
U = array(0, dim = c(2673, 10))


logLik_values_1 = numeric(10)


spec1 = ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 1), arfima = FALSE),
  distribution.model = "ged"
)


for (i in 1:10) {
  garch[[i]] = ugarchfit(
    spec = spec1, 
    data = f[, i],
    out.sample = 0,
    solver = "solnp",
    solver.control = list(trace = 0)
  )
  

  sigma[, i] <- as.numeric(sigma(garch[[i]]))
  z[, i] <- residuals(garch[[i]], standardize = TRUE)
  

  logLik_values_1[i] <- likelihood(garch[[i]])
  

  fit[[i]] = GNG_fit(z[, i], start = c(break1 = -2, break2 = 1.5, mean = 0, sd = 1, shape1 = 0.1, shape2 = 0.1))
  CDF[[i]] = function(x) { mistr::p(distribution(fit[[i]]), x) }
  inverse_CDF[[i]] = function(x) { mistr::q(distribution(fit[[i]]), x) }
  U[, i] = CDF[[i]](z[, i])
}





garch = list()
sigma = array(0, dim = c(2673, 10))
z = array(0, dim = c(2673, 10))
fit = list()
CDF = list()
inverse_CDF = list()
U = array(0, dim = c(2673, 10))


logLik_values_0 = numeric(10)


spec1 = ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), arfima = FALSE),
  distribution.model = "ged"
)


for (i in 1:10) {
  garch[[i]] = ugarchfit(
    spec = spec1, 
    data = f[, i],
    out.sample = 0,
    solver = "solnp",
    solver.control = list(trace = 0)
  )
  

  sigma[, i] <- as.numeric(sigma(garch[[i]]))
  z[, i] <- residuals(garch[[i]], standardize = TRUE)
  

  logLik_values_0[i] <- likelihood(garch[[i]])
  

  fit[[i]] = GNG_fit(z[, i], start = c(break1 = -2, break2 = 1.5, mean = 0, sd = 1, shape1 = 0.1, shape2 = 0.1))
  CDF[[i]] = function(x) { mistr::p(distribution(fit[[i]]), x) }
  inverse_CDF[[i]] = function(x) { mistr::q(distribution(fit[[i]]), x) }
  U[, i] = CDF[[i]](z[, i])
}





LRT_values <- numeric(10)
p_values <- numeric(10)



for (i in 1:10) {

  LRT_values[i] <- -2 * (logLik_values_0[i] - logLik_values_1[i])
  

  p_values[i] <- 1 - pchisq(LRT_values[i], df = 2)
}


print(LRT_values)
print(p_values)

