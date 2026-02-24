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
z=array(0,dim = c(2673,10))
fit=list()
CDF=list()
inverse_CDF=list()
U=array(0,dim = c(2673,10))
spec1=ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),mean.model = list(armaOrder = c(1, 1), arfima = FALSE), distribution.model = "ged")
for (i in 1:10) {
  z[,i]=as.numeric(residuals(ugarchfit(spec = spec1, data=f[,i],out.sample=0,solver="solnp",solver.control=list(trace=0)),standardize=TRUE))
  fit[[i]]=GNG_fit(z[,i], start = c(break1 = -2, break2 =1.5, mean = 0, sd =1,shape1 = 0.1, shape2= 0.1))
  CDF[[i]]=function(x){mistr::p(distribution(fit[[i]]),x)}
  inverse_CDF[[i]]<-function(x){mistr::q(distribution(fit[[i]]),x)}
  U[,i]=CDF[[i]](z[,i])
  
}
RVINE<-RVineStructureSelect(U,c(1:6),rotations=T,se=T,type=0)
#str(RVINE)
#summary(RVINE)
contour(RVINE)
RVINEM <- RVineMatrix(Matrix = RVINE$Matrix, family = RVINE$family,par = RVINE$par, par2 = RVINE$par2,names = RVINE$names)
plot(RVINEM)

RVineTreePlot(RVINE,edge.labels='family')
