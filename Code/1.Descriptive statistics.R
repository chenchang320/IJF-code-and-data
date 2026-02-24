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
data1=read_excel('/Users/mac/Desktop/git-hub/数据/10支上证行业指数.xlsx')
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


garch=list()
sigma=array(0,dim = c(2673,10))
z=array(0,dim = c(2673,10))
fit=list()
CDF=list()
inverse_CDF=list()
U=array(0,dim = c(2673,10))
spec1=ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),mean.model = list(armaOrder = c(1, 1), arfima = FALSE), distribution.model = "ged")
for (i in 1:10) {
  garch[[i]]=ugarchfit(spec = spec1, data=f[,i],out.sample=0,solver="solnp",solver.control=list(trace=0))
  sigma[,i]<-as.numeric(sigma(garch[[i]]))#提取条件标准差#
  z[,i]<-residuals(garch[[i]],standardize=TRUE)#提取标准化残差#
  fit[[i]]=GNG_fit(z[,i], start = c(break1 = -2, break2 =1.5, mean = 0, sd =1,shape1 = 0.1, shape2= 0.1))
  CDF[[i]]=function(x){mistr::p(distribution(fit[[i]]),x)}
  inverse_CDF[[i]]<-function(x){mistr::q(distribution(fit[[i]]),x)}
  U[,i]=CDF[[i]](z[,i])
  
}

after_LB_test=list()
after_arch=list()
for (j in 1:10) {
  after_LB_test[[j]]=Box.test(z[,j],5,type='Ljung-Box')
  after_arch[[j]]=archTest(z[,j],lag=20)
}
