#This code constructs a DET and a credible set for the Banana shaped distribution in Section 5.2. 

rm(list = ls())
library(ks); 
library(HDInterval)
source("funcode.R"); 
source('banana.r');
set.seed(1234)

#number of train samples
N=3e5; 
#number of test samples
M=3e4;
#nominal level
alpha= 0.1 

#tau is a smoothing parameter and this is equivalent to tau in the paper.
#potential tau-values
vector_tau = seq(1,0.01,len=30)

#train data
train.sim=banana(N=N); 
#test data
test.sim = banana(N=M)
#transformation of train-data to a unit hyper0rectangular using min/max values of train data
conv_training=convert.marginal(train.sim,train.sim)
#transformation of test-data
conv_testing=convert.marginal(train.sim,test.sim)

#density evaluation at each test sample
f.test = d.banana(x=test.sim); 
#alpha-density level estimation 
low.limit = quantile(d.banana(x=train.sim),alpha)

#DET at a particular tau value, 50
DET = density_estimation(tau=50,sz=dim(train.sim),train.sim=train.sim,dataset=conv_training)
#The output cells are sorted in decreasing order of density

#Piecewise density approximation at each cell
DET$uni_density 
      
#Probability at each cell
DET$uni_prob          

#Boundary of each cell
DET$uni_space
# each list is a 2 x length of coordinates matrix; 
# lower (1st row) and upper (2nd row) boundaries 
# column indicates lower and upper boundaries of each coordinates 

#Histrogical partition decision
DET$partition_decision

#Train samples in each cell
DET$uni_point


#credible set at the level, (1-alpha) using Algorithm 1
out1=cred.est1(test=conv_testing,train.sim=train.sim,train=conv_training,tau=vector_tau,alpha=alpha,p.alpha=0.95,f.test,low.limit)

#credible set at the level, (1-alpha) using Algorithm 2
out2=cred.est2(test=conv_testing,train.sim=train.sim,train=conv_training,tau=vector_tau,alpha=alpha,p.alpha=0.95)

#validation samples and density
test.mc=banana(N=M); ftest = d.banana(x=test.mc); 
#validation sample transformation
conv_true=convert.marginal(train.sim,test.mc)
#proportion of validation samples in a credible set estimation from out1
r1=cov.prob(out1$cre.set,conv_true); 
r1$coverage_prob
  
#proportion of validation samples in a credible set estimation from out2
r2=cov.prob(out2$cre.set,conv_true); 
r2$coverage_prob


#Scatter plot with a joint credible set 
f1=train.sim[,1]; f2=train.sim[,2]
ii=which(ftest<low.limit)
dev.new(); plot(test.mc[ii,1:2],col='dark green',xlim=range(train.sim[,1]),ylim=range(train.sim[,2]),ylab='y',xlab='x')
points(test.mc[-ii,1:2],col='grey')
for(i in 1:length(out2$cre.set)){pt=out1$cre.set[[i]][,1:2]; 
  pt[,1]=min(f1)+diff(range(f1))*pt[,1]; pt[,2]=min(f2)+diff(range(f2))*pt[,2]
  lines(c(pt[1,1],pt[2,1],pt[2,1],pt[1,1],pt[1,1]),c(pt[1,2],pt[1,2],pt[2,2],pt[2,2],pt[1,2]),col='red') 
}


  