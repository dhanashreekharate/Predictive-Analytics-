data<-read.csv("HW1_data.csv",sep=",",header=TRUE)
MLC<-read.table("MLC.csv",sep=",",header=TRUE)
MLC<-read.table("MLC.csv",sep=",",header=TRUE)
x1<-MLC$Location;
x2<-MLC$Week;
y<-MLC$Efficiency
fn <- function(p) {yhat<-p[1]+p[2]*x1+p[4]*exp(p[3]*x2); sum((y-yhat)^2)} 


out<-nlm(fn,p=c(1,0,-.5,-.1),hessian=TRUE)

theta<-out$estimate#parameter estimates
theta
###we will use the following later, for finding SEs and CIs#######
MSE<-out$minimum/(length(y) - length(theta))  #estimate of the error variance
InfoMat<-out$hessian/2/MSE  #observed information matrix
CovTheta<-solve(InfoMat)
SE<-sqrt(diag(CovTheta))  #standard errors of parameter estimates
MSE
CovTheta
SE


#2a and 3a
data<-read.csv("HW1_data.csv",sep=",",header=TRUE)
x<- data$x
y <- data$ï..y
fun <- function(i){yhat <- i[1]+i[2]*x; sum((y-yhat)^2)}
output <- nlm(fun,p=c(1,0),hessian= TRUE)
theta1<-output$estimate #parameter estimates
MSE<-output$minimum/(length(y) - length(theta1))  #estimate of the error variance
InfoMat<-output$hessian/2/MSE  # info matrix
CovTheta<-solve(InfoMat)  #covariance matrix
SE<-sqrt(diag(CovTheta))  #standard errors of parameter estimates
MSE
CovTheta
SE

#2b
fun2 <- function(x,p){p[1]+p[2]*x}
output2<-nls(y~fun2(x,p),start=list(p=c(0,-0.5)),trace=TRUE)
theta2<-output2$estimate #parameter estimates
MSE1<-output2$minimum/(length(y) - length(theta2))  #estimate of the error variance
InfoMat1<-output2$hessian/2/MSE1  # info matrix
CovTheta1<-solve(InfoMat1)  #covariance matrix
SE1<-sqrt(diag(CovTheta1))  #standard errors of parameter estimates
MSE1
CovTheta1
SE1

#3b
Cov<-vcov(output2)
SE1<- sqrt(diag(Cov))

#3c
crude <- confint.default(output2,level=0.90)



library(boot)   #need to load the boot package
datafit<-function(Z,i,theta0,x_pred) {
Zboot<-Z[i,]
x1<-Zboot[[1]];y<-Zboot[[2]]
fn <- function(p) {yhat<-p[1]+p[2]*x1; sum((y-yhat)^2)} 
out<-nlm(fn,p=theta0)
theta<-out$estimate 
y_pred<- theta[1]+theta[2]*x_pred[1]} #predicted response
databoot<-boot(data, datafit, R=20000, theta0=c(1,-0.05), x_pred=c(16))
databoot

VarYhat<-var(databoot$t); 
VarYhat
SEYhat<-sqrt(VarYhat); 
SEYhat
plot(databoot)  
boot.ci(databoot,conf=c(.9,.95,.99),type=c("norm","basic"))


#6 same AIC 
fun2 <- function(x,p){p[1]+p[2]*x}
output2<-nls(y~fun2(x,p),start=list(p=c(0,-0.5)),trace=TRUE)
n <- nrow(data)
AIC<- -2*as.numeric(logLik(output2))/n+2*1/n
fun3 <- function(x,p){p[1]+p[2]*sqrt(x)}
output3 <-nls(y~fun3(x,p),start=list(p=c(0,-0.5)),trace=TRUE)
AIC_1 <- -2*as.numeric(logLik(output3))/n+2*1/n
out<-cv.glm(data, output3, function(y,phat) -mean(log(phat)*y+log(1-phat)*(1-y)), K=9)



#5
library(boot)   #need to load the boot package
data<-read.csv("HW1_data.csv",sep=",",header=TRUE)
datafit1<-function(Z,i,theta0,x_pred) {
  Zboot<-Z[i,]
  x1<-Zboot[[1]];y<-Zboot[[3]]
  fn <- function(p) {yhat<-p[1]+p[2]*x1; sum((y-yhat)^2)} 
  out<-nlm(fn,p=theta0)
  theta<-out$estimate 
  y_pred<- theta[1]+theta[2]*x_pred[1]} #predicted response
databoot<-boot(data, datafit, R=5000, theta0=c(1,2), x_pred=c(27))
databoot
VarYhat<-var(databoot$t); VarYhat
SEYhat<-sqrt(VarYhat); SEYhat
plot(databoot)  
boot.ci(databoot,conf=c(.95),type=c("norm","basic"))
Yhat0<-databoot$t0
Yhatboot<-databoot$t
g.hat <- 40.63  ;
MSE <- 4.68  #from slide 46
e<-rnorm(nrow(Yhatboot), mean=0, sd=sqrt(MSE))  
Yboot<-Yhatboot-e
Yquant<-quantile(Yboot,prob=c(0.95))
L<-2*Yhat0-Yquant[2]
U<-2*Yhat0-Yquant[1]
hist(Yboot,100)
c(L,U) #more complex PI
#simpler PI
SEY<-sqrt(var(Yhatboot)+MSE); SEY
c(g.hat-qnorm(0.95)*SEY, g.hat+qnorm(0.95)*SEY) #simpler PI



#7 model 2 is better
CVInd <- function(n,K) {  #n is sample size; K is number of parts; returns K-length list of indices for each part
  m<-floor(n/K)  #approximate size of each part
  r<-n-m*K  
  I<-sample(n,n)  #random reordering of the indices
  Ind<-list()  #will be list of indices for all K parts
  length(Ind)<-K
  for (k in 1:K) {
    if (k <= r) kpart <- ((m+1)*(k-1)+1):((m+1)*k)  
    else kpart<-((m+1)*r+m*(k-r-1)+1):((m+1)*r+m*(k-r))
    Ind[[k]] <- I[kpart]  #indices for kth part of data
  }
  Ind
}

Nrep<-20 #number of replicates of CV
K<-10  #K-fold CV on each replicate
n.models = 2 #number of different models to fit and compare
FitFun1 <- function(x,p){ p[1]+p[2]*x }
FitFun2 <- function(x,p) { p[1]+p[2]*sqrt(x) }
n=nrow(data)
y<-data$ï..y
yhat=matrix(0,n,n.models)
MSE<-matrix(0,Nrep,n.models)
for (j in 1:Nrep) {
  Ind<-CVInd(n,K)
  for (k in 1:K) {
    out1<-nls(ï..y~FitFun1(x,p),data=data[-Ind[[k]],],start=list(p=c(-.15,-.55)))
    yhat[Ind[[k]],1]<-as.numeric(predict(out1,data[Ind[[k]],]))
    #yhat_res[Ind[[k]],1]<- as.numeric(resid(out1,data[Ind[[k]],]))
    out2<-nls(ï..y~FitFun2(x,p),data=data[-Ind[[k]],],start=list(p=c(-.15,-.55)))
    yhat[Ind[[k]],2]<-as.numeric(predict(out2,data[Ind[[k]],]))
    #yhat_res_2[Ind[[k]],2]<-as.numeric(resid(out2,data[Ind[[k]],]))
  } #end of k loop
  MSE[j,]=apply(yhat,2,function(x) sum((y-x)^2))/n
} #end of j loop

MSE
MSEAve<- apply(MSE,2,mean); MSEAve #averaged mean square CV error
MSEsd <- apply(MSE,2,sd); MSEsd   #SD of mean square CV error
r2<-1-MSEAve/var(y); r2  #CV r^2

AIC<- -2*as.numeric(logLik(out1))/n+2*2/n
AIC_1<- -2*as.numeric(logLik(out2))/n+2*3/n








