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
library(fPortfolio)
data1=read_excel('/Users/mac/Desktop/最初的沪深300中的158支个股股票数据.xlsx')
select_data=data1[c(3:2676),-1]
colnames(select_data)=c(1:158)
select_data=select_data[,c(-23,-74,-133)]
colnames(select_data)=c(1:155)
select_data=select_data[,c(-19,-66,-77,-78,-82,-85,-99,-100,-102,-109,-112,-129,-136,-144,-154)]
colnames(select_data)=c(1:140)
select_data=select_data[,c(-2,-7,-40,-62,-112)]
colnames(select_data)=c(1:135)
select_data=select_data[,c(-41,-78,-91)]
colnames(select_data)=c(1:132)
select_data=select_data[,-21]
colnames(select_data)=c(1:131)
select_data=select_data[,-91]
colnames(select_data)=c(1:130)

log.return=function(x){
  diff(log(x))
}

log_return=matrix(0,nrow = 2673,ncol = 130)
for (d in c(1:130)){
  e=select_data[[d]]
  h=as.numeric(e)
  log_return[,d]=log.return(h)
  
}  
f=log_return



z=array(0,dim = c(2673,130))
consigma=array(0,dim = c(2673,130))
cancha=array(0,dim = c(2673,130))
r=array(0,dim = c(2673,130))
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), mean.model = list(armaOrder = c(0,0)), distribution.model = "ged")

for (i in 1:130) {
  z[,i]=as.numeric(residuals(ugarchfit(spec = spec1, data=f[,i],out.sample=0,solver="solnp",solver.control=list(trace=0)),standardize=TRUE))
  consigma[,i]=as.numeric(sigma(ugarchfit(spec = spec1, data=f[,i],out.sample=0,solver="solnp",solver.control=list(trace=0))))
  cancha[,i]=z[,i]*consigma[,i]
  r[,i]=cancha[,i]+coef(ugarchfit(spec = spec1, data=f[,i],out.sample=0,solver="solnp",solver.control=list(trace=0)))[[1]]
}



VAR=matrix(0,974,18)
CVAR=matrix(0,974,18)
r_1=array(0,dim = c(1700,130,974))
for (j in 1:974) {
  r_1[,,j]=r[j:(j+1699),]

  beta_hat=matrix(0,1700,129)
  for (a in 2:130) {
    for (b in 1:1700) {
      beta_hat[b,a-1]=r_1[b,1,j]-r_1[b,a,j]
      
    }
    
  }
  g_1=rq.fit.scad(beta_hat,r_1[,1,j],tau = 0.01,alpha = 3.7,lambda = 0.6)
  g_2=rq.fit.scad(beta_hat,r_1[,1,j],tau = 0.05,alpha = 3.7,lambda = 0.6)
  g_3=rq.fit.scad(beta_hat,r_1[,1,j],tau = 0.1,alpha = 3.7,lambda = 0.6)
  beta_end_1=matrix(c(1-sum(g_1[["coefficients"]]),g_1[["coefficients"]]),130,1)
  beta_end_2=matrix(c(1-sum(g_2[["coefficients"]]),g_2[["coefficients"]]),130,1)
  beta_end_3=matrix(c(1-sum(g_3[["coefficients"]]),g_3[["coefficients"]]),130,1)
  y_end_1=r_1[,,j]%*%beta_end_1
  y_end_2=r_1[,,j]%*%beta_end_2
  y_end_3=r_1[,,j]%*%beta_end_3
  
  VaR_ninetynine_1<-sort(y_end_1)[10]
  ES_ninetynine_1<-sum(sort(y_end_1)[1:10])/10
  VaR_ninetyfive_1<-sort(y_end_1)[50]
  ES_ninetyfive_1<-sum(sort(y_end_1)[1:50])/50
  VaR_ninety_1<-sort(y_end_1)[100]
  ES_ninety_1<-sum(sort(y_end_1)[1:100])/100
  VaR_ninetynine_2<-sort(y_end_2)[10]
  ES_ninetynine_2<-sum(sort(y_end_2)[1:10])/10
  VaR_ninetyfive_2<-sort(y_end_2)[50]
  ES_ninetyfive_2<-sum(sort(y_end_2)[1:50])/50
  VaR_ninety_2<-sort(y_end_2)[100]
  ES_ninety_2<-sum(sort(y_end_2)[1:100])/100
  VaR_ninetynine_3<-sort(y_end_3)[10]
  ES_ninetynine_3<-sum(sort(y_end_3)[1:10])/10
  VaR_ninetyfive_3<-sort(y_end_3)[50]
  ES_ninetyfive_3<-sum(sort(y_end_3)[1:50])/50
  VaR_ninety_3<-sort(y_end_3)[100]
  ES_ninety_3<-sum(sort(y_end_3)[1:100])/100
  
  #EWP
  beta_ewp=1/130
  ewp_vector=matrix(beta_ewp,130,1)
  y_ewp=r%*%ewp_vector
  VaR_ewp_ninetynine<-sort(y_ewp)[10]
  ES_ewp_ninetynine<-sum(sort(y_ewp)[1:10])/10
  VaR_ewp_ninetyfive<-sort(y_ewp)[50]
  ES_ewp_ninetyfive<-sum(sort(y_ewp)[1:50])/50
  VaR_ewp_ninety<-sort(y_ewp)[100]
  ES_ewp_ninety<-sum(sort(y_ewp)[1:100])/100
  
  
  #min-var
  mo=minvariancePortfolio(as.timeSeries(r_1[,,j]), spec = portfolioSpec(), constraints = "LongOnly")
  qo=getWeights(mo)
  beta_minvar=matrix(qo,130,1)
  y_minvar=r_1[,,j]%*%beta_minvar
  VaR_minvar_ninetynine<-sort(y_minvar)[10]
  ES_minvar_ninetynine<-sum(sort(y_minvar)[1:10])/10
  VaR_minvar_ninetyfive<-sort(y_minvar)[50]
  ES_minvar_ninetyfive<-sum(sort(y_minvar)[1:50])/50
  VaR_minvar_ninety<-sort(y_minvar)[100]
  ES_minvar_niniety<-sum(sort(y_minvar)[1:100])/100
  #mean=var
  Spec = portfolioSpec()
  setTargetReturn(Spec) = mean(colMeans(r_1[,,j]))
  mo_1=efficientPortfolio(as.timeSeries(r_1[,,j]), Spec , constraints = "LongOnly")
  qo_1=getWeights(mo_1)
  beta_meanvar=matrix(qo_1,130,1)
  y_meanvar=r_1[,,j]%*%beta_meanvar
  VaR_meanvar_ninetynine<-sort(y_meanvar)[10]
  ES_meanvar_ninetynine<-sum(sort(y_meanvar)[1:10])/10
  VaR_meanvar_ninetyfive<-sort(y_meanvar)[50]
  ES_meanvar_ninetyfive<-sum(sort(y_meanvar)[1:50])/50
  VaR_meanvar_ninety<-sort(y_meanvar)[100]
  ES_meanvar_ninety<-sum(sort(y_meanvar)[1:100])/100
  
  VAR[j,1]=VaR_ninetynine_1
  VAR[j,2]=VaR_ninetynine_2
  VAR[j,3]=VaR_ninetynine_3
  VAR[j,4]=VaR_ewp_ninetynine
  VAR[j,5]=VaR_minvar_ninetynine
  VAR[j,6]=VaR_meanvar_ninetynine
  VAR[j,7]=VaR_ninetyfive_1
  VAR[j,8]=VaR_ninetyfive_2
  VAR[j,9]=VaR_ninetyfive_3
  VAR[j,10]=VaR_ewp_ninetyfive
  VAR[j,11]=VaR_minvar_ninetyfive
  VAR[j,12]=VaR_meanvar_ninetyfive
  VAR[j,13]=VaR_ninety_1
  VAR[j,14]=VaR_ninety_2
  VAR[j,15]=VaR_ninety_3
  VAR[j,16]=VaR_ewp_ninety
  VAR[j,17]=VaR_minvar_ninety
  VAR[j,18]=VaR_meanvar_ninety
  CVAR[j,1]=ES_ninetynine_1
  CVAR[j,2]=ES_ninetynine_2
  CVAR[j,3]=ES_ninetynine_3
  CVAR[j,4]=ES_ewp_ninetynine
  CVAR[j,5]=ES_minvar_ninetynine
  CVAR[j,6]=ES_meanvar_ninetynine
  CVAR[j,7]=ES_ninetyfive_1
  CVAR[j,8]=ES_ninetyfive_2
  CVAR[j,9]=ES_ninetyfive_3
  CVAR[j,10]=ES_ewp_ninetyfive
  CVAR[j,11]=ES_minvar_ninetyfive
  CVAR[j,12]=ES_meanvar_ninetyfive
  CVAR[j,13]=ES_ninety_1
  CVAR[j,14]=ES_ninety_2
  CVAR[j,15]=ES_ninety_3
  CVAR[j,16]=ES_ewp_ninety
  CVAR[j,17]=ES_minvar_niniety
  CVAR[j,18]=ES_meanvar_ninety
  print(j)
}
