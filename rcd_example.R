rm(list = ls())
source("rcdmcmc.R") 
SubS=3; T=3e5*SubS; M=3e4*SubS; d=9; alpha=0.1 

OP.train=rcd.mcmc(nd,L,U,T,SubS,w=0.25,scale.move=TRUE,beta=1)
train.sim=OP.train[,c(2:10)]; 
f.train=OP.train[,11]+OP.train[,12]
low.limit=quantile(f.train,alpha)

OP.test=rcd.mcmc(nd,L,U,M,SubS,w=0.25,scale.move=TRUE,beta=1)
test.sim=OP.test[,c(2:10)]; 
f.test=OP.test[,11]+OP.test[,12];

OP.exam=rcd.mcmc(nd,L,U,M,SubS,w=0.25,scale.move=TRUE,beta=1)
exam.sim=OP.exam[,c(2:10)]; 
f.exam=OP.exam[,11]+OP.exam[,12];

source("funcode.R"); 
conv_training=convert.marginal(train.sim,train.sim)
conv_testing=convert.marginal(train.sim,test.sim)
conv_exam=convert.marginal(train.sim,exam.sim)

vector.tau=c(seq(0.5,0.01,len=10),0.005)
out=cred.est1(test=conv_testing,train.sim=train.sim,train=conv_training,tau=vector.tau,alpha=alpha,p.alpha=0.95,f.test,low.limit)

r=cov.prob(out$cre.set,conv_exam); 
r$coverage_prob

