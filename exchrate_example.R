#This code will construct a joint credible set of a posterior approximation and estimate operational coverages in Section 5.3.1.

rm(list=ls())
library(gk); library(EasyABC)
source("funcode.R"); # library(HDInterval);library(mclust)

nz=read.csv('NZUSrate.csv')
data=2e2*log(nz$rate[2:1085]/nz$rate[1:1084])

#summary stats
p.oct=seq(from=0,to=1,length.out=8)
sum_stat_obs=quantile(data,p.oct)

rate_model<-function(x){
  quantile(rgk(1084,A=x[1],B=(x[2]),g=x[3],k=x[4],c = 0.8),seq(from=0,to=1,length.out=8))  }

rate_prior=list(c("unif",-1,1),c("unif",0,1),c("unif",-5,5),c("unif",0,10))

ABC.train<-ABC_mcmc(method="Marjoram", model=rate_model,prior=rate_prior,summary_stat_target=sum_stat_obs,n_between_sampling=5,n_rec=3.1e5)
ABC.test<-ABC_mcmc(method="Marjoram", model=rate_model,prior=rate_prior,summary_stat_target=sum_stat_obs,n_between_sampling=5,n_rec=4e4)

#Use the below for ABC simulation using Wegmann method.
#ABC.train <-ABC_mcmc(method="Wegmann", model=rate_model,prior=rate_prior,summary_stat_target=sum_stat_obs,n_between_sampling=5, n_rec=3.01e5)
#ABC.test <-ABC_mcmc(method="Wegmann", model=rate_model,prior=rate_prior,summary_stat_target=sum_stat_obs,n_between_sampling=5, n_rec=3.1e4)


## approx credible set using Marjoram
train.sim=ABC.train$param[-c(1:1e4),]
test.sim=ABC.test$param[-c(1:1e4),]
conv_training <- convert.marginal(train.sim,train.sim)
conv_testing <- convert.marginal(train.sim,test.sim)
#vector.tau = seq(15,1,len=10)
vector.tau = seq(50,10,len=10)
out=cred.est2(test=conv_testing,train.sim=train.sim,train=conv_training,tau=vector.tau,alpha=0.1,p.alpha=0.95)


#Let's simulate samples from the exact posterior.
log_prior = function(theta) {
  return(sum(dunif(theta[1],-1,1,log=TRUE)+dunif(exp(theta[2]),0,1,log=TRUE)+dunif(theta[3],-5,5,log=TRUE)+dunif(theta[4],0,10,log=TRUE)))
}

b0=c(0,-0.9,-0.8,3.5); Sigma0=diag(c(0.1,1,5,5))
mcmcout=mcmc(data, N=1e5, logB=TRUE, get_log_prior=log_prior,theta0=b0, Sigma0=Sigma0)

exam.sim = mcmcout[seq(1e4,1e5,by=3,]; 
exam.sim[,2]=exp(exam.sim[,2])

conv_exam=convert.marginal(train.sim,exam.sim)
r=cov.prob(out$cre.set,conv_exam)
#exact operational calibration estimation 
r$coverage_prob

#Regression approach using BART
library(BART)
### Reg-tree estimate
nn=1e5
b.bart=cbind(runif(nn,-1,1), runif(nn,0,1), runif(nn,-5,5), runif(nn,0,10))
y.bart=matrix(,nn,8); d.bart=rep(0,nn)
for(i in 1:nn){
  y0=rgk(1e3,A=b.bart[i,1],B=b.bart[i,2],g=b.bart[i,3],k=b.bart[i,4],c = 0.8) 
  y.bart[i,]=quantile(y0,p.oct)
  d.bart[i]=mean((y.bart[i,]-sum_stat_obs)^2)
}

####Use the half of particles to build a reg-tree. 
conv_sim=convert.marginal(train.sim,b.bart[d.bart<quantile(d.bart,0.2),])
c.bart=cov.prob(out$cre.set,conv_sim)$y
bart.fit=pbart(as.matrix(y.bart[d.bart<quantile(d.bart,0.2),]),c.bart)
bart.pred = predict(bart.fit,matrix(sum_stat_obs,nrow=1))
bart.pred$prob.test.mean
sd(bart.pred$prob.test)

