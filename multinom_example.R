#This code consruct a DET and joint credible set for posterior probabilities in Section 5.2.2.

rm(list = ls())
source("funcode.R"); 
library(MCMCpack); library(mvtnorm); library(MASS); 
library(foreign); 
hsb <- read.dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta")

N=3e5; #training sample size 
M=3e4; #test sample size
alpha=0.1; #target quantile level

vector_tau = seq(0.5,0.01,len=10) #potential tau values
x=c(1,1,0,1,36) #predictor values of interest

#Get the train samples
post3 <- MCMCmnl(prog~female + ses + write, mcmc.method="RWM", B0=0,verbose=500, mcmc=N*10, thin=10, tune=0.5,data=hsb)
train.sim=matrix(,nrow(post3),3)
train.sim[,3]=1/(1+exp(post3[,c(1:5)*2-1]%*%x)+exp(post3[,c(1:5)*2]%*%x))
train.sim[,1]=exp(post3[,c(1:5)*2-1]%*%x)*train.sim[,3]
train.sim[,2]=exp(post3[,c(1:5)*2]%*%x)*train.sim[,3]

#Get the test samples  
post3 <- MCMCmnl(prog~female + ses + write,mcmc.method="RWM", B0=0,verbose=500, mcmc=M*10, thin=10, tune=0.5,data=hsb)
test.sim=matrix(,nrow(post3),3)
test.sim[,3]=1/(1+exp(post3[,c(1:5)*2-1]%*%x)+exp(post3[,c(1:5)*2]%*%x))
test.sim[,1]=exp(post3[,c(1:5)*2-1]%*%x)*test.sim[,3]
test.sim[,2]=exp(post3[,c(1:5)*2]%*%x)*test.sim[,3]

#Train and test samples are transformed.   
conv_training=convert.marginal(train.sim,train.sim)
conv_testing=as.matrix(convert.marginal(train.sim,test.sim))

#Find the DET and joint credible set using Algorithm 2.  
out=cred.est2(test=conv_testing,train.sim=train.sim,train=conv_training,tau=vector_tau,alpha=alpha,p.alpha=0.95)

#Scatter plot and credible set estimation
xx=c(1,2); f1=range(train.sim[,xx[1]]); f2=range(train.sim[,xx[2]])
dev.new(); plot(exam.sim[,xx],col='darkgrey',xlim=range(train.sim[,xx[1]]),ylim=range(train.sim[,xx[2]]),xlab='general',ylab='vocation'); 
for(i in 1:length(out$cre.set)){pt=out$cre.set[[i]][,xx];
pt[,1]=pt[,1]*diff(f1)+f1[1]; pt[,2]=pt[,2]*diff(f2)+f2[1]
lines(c(pt[1,1],pt[2,1],pt[2,1],pt[1,1],pt[1,1]),c(pt[1,2],pt[1,2],pt[2,2],pt[2,2],pt[1,2]),col='red') 
}

#Get M-samples to examine the credible set estimation
post3 <- MCMCmnl(prog~female + ses + write,mcmc.method="RWM", B0=0,verbose=500, mcmc=M*10, thin=10, tune=0.5,data=hsb)
exam.sim=matrix(,nrow(post3),3)
exam.sim[,3]=1/(1+exp(post3[,c(1:5)*2-1]%*%x)+exp(post3[,c(1:5)*2]%*%x))
exam.sim[,1]=exp(post3[,c(1:5)*2-1]%*%x)*exam.sim[,3]
exam.sim[,2]=exp(post3[,c(1:5)*2]%*%x)*exam.sim[,3]

conv_exam=convert.marginal(train.sim,exam.sim)

#Coverage estimation
r=cov.prob(out$cre.set,conv_exam)
r$coverage_prob
    
