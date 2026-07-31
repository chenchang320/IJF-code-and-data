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
data1=read_excel('/Users/mac/Desktop/GitHub-Englishi version/data/Ten industry dataset.xlsx')
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

