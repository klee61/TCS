-------------------------------------------------------------------------------
convert.marginal(train,test)
-------------------------------------------------------------------------------
Describtion: Construct a transformation function using min/max values of train samples and transform test samples according to the transformation function

Arguments:
train - Train samples used to construct a transformation function.
test - Test samples.

Values: Transformed test samples.

-------------------------------------------------------------------------------
cov.prob(cre.set,conv.test)
-------------------------------------------------------------------------------
Description: Proportion of test samples in the credible set.

Arguments: 
cre.set - Joint credible set.
conv.test - transformed test samples.

Values:
Proportion of test samples in the joint credible set.

-------------------------------------------------------------------------------
density_estimation(tau,sz=dim(train.sim),train.sim=train.sim,dataset=conv_training)
-------------------------------------------------------------------------------
Description: Fit a density estimation tree using a particular smoothing value, tau.

Arguments:
tau - tau value.
sz - Dimension of train samples.
train.sim - Train samples.
dataset - Transformed train samples.

Values: DET is fitted on a unit rectangular.
uni_prob - Probability at each hyper-rectangular.        
uni_space - Boundary of each hyper-rectangular.
	    Each list is a 2 x length of coordinates matrix; 
	    Lower and upper bounds on the 1st and 2nd rows respectively.
	    Each column shows bounds of each coordinate.
partition_decision - Histrogical partition decision.
uni_point - List of train samples in each hyper-rectangular. 

-------------------------------------------------------------------------------
cred.est1(test,train.sim,train,tau,alpha,p.alpha,f.test,low.limit)
-------------------------------------------------------------------------------
Description: Find the optimal DET and estimate a joint credible set estimation using Algorithm 1

Arguments: 
test - Transformed test samples. 
train.sim - Train samples (n x d).
train - Transformed train samples.
tau - A vector of potential tau values. It should be in a decreasing order.
alpha - (1-alpha) target coverage of a credible set.
p.alpha - Confidence level for a coverage testing. The default value is 0.95.
f.test - Density estimates for test samples.
low.limit - Density level.

Values: 
cre.set - Joint credible set. It is a list of hyper-rectangulars.
n.tree - Number of leaves.
prob - Proportion of test samples in cre.set.
tau.out - tau-value of the optimal DET.
uni_density - Density of each hyper-rectangulars.
uni_prob - Probabitliy of each hyper-rectangulars.
uni_space - Boundary of each hyper-rectangulars.
y.out - length(tau)x4 measurement matrix; 
	For each tau value, the coverage test result (1st column), 
        coverage estimation (2nd column), false positive (3rd column) 
        and false negative (4th column) are estimated. Here 0 means the coverage test passes.

-------------------------------------------------------------------------------
cred.est2(test,train.sim,train,tau,alpha,p.alpha)
-------------------------------------------------------------------------------
Description: Find the optimal DET and estimate a joint credible set estimation using Algorithm 2

Arguments: 
test - Transformed test samples. 
train.sim - Train samples (n x d).
train - Transformed train samples.
tau - A vector of potential tau values. It should be in a decreasing order.
alpha - (1-alpha) target coverage of a credible set.
p.alpha - Confidence level for a coverage testing. The default value is 0.95.

Values: 
cre.set - Joint credible set. It is a list of hyper-rectangulars.
n.tree - Number of leaves.
prob - Proportion of test samples in cre.set.
tau.out - tau-value of the optimal DET.
uni_density - Density of each hyper-rectangulars.
uni_prob - Probabitliy of each hyper-rectangulars.
uni_space - Boundary of each hyper-rectangulars.
y.out - length(tau)x4 measurement matrix; 
	For each tau value, the coverage test result (1st column), 
        coverage estimation (2nd column), false positive (3rd column) 
        and false negative (4th column) are estimated. Here 0 means the coverage test passes.
