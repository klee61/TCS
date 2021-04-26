#This code will replicate the study in Figure 1. 

rm(list = ls())
source("funcode.R"); library(mvtnorm); library(MASS); 

d=2; N=3e5; M=3e4; alpha=0.1
vector_tau = exp(seq(log(1),log(0.001),len=20))

f.test=matrix(,M,4); result=array(,dim=c(length(vector_tau),4,4)); 
FN=FP=matrix(,length(vector_tau),4)
p=opt.ind=rep(0,4)
test.sim=conv_testing=array(,dim=c(M,d,4));

train.sim=rmvnorm(N,mean = rep(0,d),sigma = diag(d))
omega_p=range(train.sim)
conv_training=convert.marginal(train.sim,train.sim)
low.limit = quantile(dmvnorm(train.sim,mean = rep(0,d),sigma = diag(d)),alpha)

for(j in 1:4){ 
  if(j==4){
    test.sim[,,j]=matrix(runif(M*d,min=omega_p[1],max=omega_p[2]),ncol=d)
  }else{test.sim[,,j]=rmvnorm(M,mean = rep(0,d),sigma = diag(d)*j)
  ii = which(rowSums(test.sim[,,j]<omega_p[1] | test.sim[,,j]>omega_p[2])>0)
  while(length(ii)>0){
    test.sim[ii,,j]=rmvnorm(length(ii),mean = rep(0,d),sigma = diag(d)*j)
    ii = which(rowSums(test.sim[,,j]<omega_p[1] | test.sim[,,j]>omega_p[2])>0)
  }}
  
  f.test[,j] = dmvnorm(test.sim[,,j],mean = rep(0,d),sigma = diag(d)) 
  conv_testing[,,j]<-as.matrix(convert.marginal(train.sim,test.sim[,,j]))
  out=cred.est1(test=conv_testing[,,j],train.sim=train.sim,train=conv_training,tau=vector_tau,alpha=alpha,p.alpha=0.95,f.test[,j],low.limit)
  result[,,j]=out$y.out
  
  p[j]=mean(f.test[,j]<low.limit)
  FP[,j]=out$y.out[,3]*p[j]; FN[,j]=out$y.out[,4]*(1-p[j])
  opt.ind[j]=which.min(FP[,j]+FN[,j])
  }

dev.new(); 
par(mfrow=c(2,2))
plot(vector_tau,FP[,1]+FN[,1],log='x',type='l',xlab=expression(tau),ylab='total loss')
points(vector_tau[opt.ind[1]],FP[opt.ind[1],1]+FN[opt.ind[1],1],pch=3)
lines(vector_tau,FP[,2]+FN[,2],lty=2)
points(vector_tau[opt.ind[2]],FP[opt.ind[2],2]+FN[opt.ind[2],2],pch=3)
lines(vector_tau,FP[,3]+FN[,3],lty=3)
points(vector_tau[opt.ind[3]],FP[opt.ind[3],3]+FN[opt.ind[3],3],pch=3)
lines(vector_tau,FP[,4]+FN[,4],lty=4)
points(vector_tau[opt.ind[4]],FP[opt.ind[4],4]+FN[opt.ind[4],4],pch=3)
legend("topright",lty=c(1:4),legend=c("c=1","c=2","c=3","uniform"))

plot(vector_tau,result[,2,1],log='x',xlab=expression(tau),ylab='coverage',col=result[,1,1]+1)

plot(vector_tau,FP[,1],log='x',type='l',ylim=c(0,0.035),xlab=expression(tau),ylab='false positive')
lines(vector_tau,FP[,2],lty=2)
lines(vector_tau,FP[,3],lty=3)
lines(vector_tau,FP[,4],lty=4)

plot(vector_tau,FN[,1],log='x',type='l',ylim=c(0,0.37),xlab=expression(tau),ylab='false negative')
lines(vector_tau,FN[,2],lty=2)
lines(vector_tau,FN[,3],lty=3)
lines(vector_tau,FN[,4],lty=4)

