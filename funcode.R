###################################################################################
# marginal pdf to marginal cdf

convert.marginal = function(data0,data1){
  d=ncol(data0); data_con = matrix(,nrow(data1),d)
  for(i in 1:d){ F = range(data0[,i])
  data_con[,i]=(data1[,i]-F[1])/diff(F)
  data_con[(data_con[,i]<0),i]=0; data_con[(data_con[,i]>1),i]=1; 
  }
  return(as.data.frame(data_con))
}


###################################################################################
# star-discrepancy
star.disc = function(uni_points,uni_space,tau=0.01,tot.n,m=10){
  
  n=length(uni_points); dimension=ncol(uni_points[[1]]); star=ul=rep(0,n)
  for(i in 1:n){ convert_point=matrix(,nrow(uni_points[[i]]),dimension); 
  for(d in 1:dimension){ 
    convert_point[,d] = (uni_points[[i]][,d]-uni_space[[i]][1,d])/uni_space[[i]][2,d]
  }
  convert_area = c(1:(m-1))/(m-1)
  prob_area = convert_area^dimension
  #star discrepancy calculation
  c=matrix(1,nrow(convert_point),(m-1))
  for( d in 1:dimension){
    c=c*sapply(convert_area,function(x) as.numeric(convert_point[,d]<x))
  }
  star[i] = max(abs(colMeans(c)-convert_area^dimension))
  ul[i] = tau*sqrt(tot.n)/nrow(uni_points[[i]])
  }
  return(list(star=star,ul=ul))
}

###################################################################################
#A coverage probability of test data in the credible set is estimated.
#two inputs; credible set estimation from cred.est and test data of interest.
#two ouputs; coverage probability, indicator
cov.prob = function(credible_set,test_data){
  #y is a vector of indicators; 1 means a test point is in the credible set and otherwise, it is 0.
  y=rep(0,nrow(test_data))
  for(i in 1:length(credible_set)){ ind=rep(0,nrow(test_data))
  for(d in 1:ncol(test_data)){ 
    ind=ind+(test_data[,d]<=credible_set[[i]][1,d])+(test_data[,d]>credible_set[[i]][2,d])  
  }
  y=y+as.numeric(ind==0)
  }
  coverage_prob=mean(y)
  return(list(coverage_prob=coverage_prob,y=y))
}


##############################################################################

density_estimation = function(m=5,tau=1.5,sz,train.sim,dataset=NULL,data_points=NULL,data_space=NULL,data_prob_region=NULL){
  
  if(length(data_points)==0){
    data_points=list(dataset); data_space=list(replicate(ncol(dataset),c(0,1))); 
    data_prob_region=1 }
  
  N = sz[1]
  dimension = sz[2]
  partition_decision=NULL
  uni_point=NULL
  uni_space=NULL
  uni_prob=NULL
  index_point=0
  
  index_of_partition=0 
  to_test_point=NULL
  to_test_space=NULL
  to_test_prob_region=NULL
  
  # starting loop over data points set
  # parameters: 
  ## space: space of interest;
  
  ## dimension: dimension;
  ## x_set: points of in this area
  while(length(data_points)>0){
    index_point=index_point+1 # index of data point sets
    
    # points set
    x_set=as.data.frame(data_points[[index_point]])
    space=data_space[[index_point]]
    prob_region=data_prob_region[index_point]
    

    if(nrow(x_set)>2 & nrow(unique(x_set))>1){
      #check discrepancy and partition based on max gap
      #define the domain of space in multi dimensions
      gap=partition=null.space=matrix(,m-2,dimension)
      #generate matrix for partition position and gap
      for(d in 1:dimension){
        split_point = seq(space[1,d],space[2,d],len=m)[2:(m-1)]
        ind=sapply(split_point ,function(x) x_set[,d]<x)
        null.space[,d]=colSums(ind)
        gap[,d]=abs(colMeans(ind)-c(1:(m-2))/(m-1))
        partition[,d]=split_point
      }

      convert_point=NULL
      for(d in 1:dimension){ convert_point = cbind(convert_point,(x_set[,d]-space[1,d])/(space[2,d]-space[1,d])) }
      
      #scaled areas
      subarea=rbind(partition,space[2,]) 
      convert_area = c(1:20)/20
      prob_area = convert_area^dimension
      #star discrepancy calculation
      s_discrepancy=NULL
      c=matrix(1,nrow(convert_point),20)
      for( d in 1:dimension){
        c=c*sapply(convert_area,function(x) as.numeric(convert_point[,d]<x))
      }
      
      if(sum(!is.na(gap))>0){split.cond=(max(abs(colMeans(c)-convert_area^dimension))<=tau*sqrt(N)/nrow(x_set))
      }else{split.cond=TRUE} 
      
      # evaluate star discrepancy
      if (split.cond){
        uni_point=c(uni_point,list(x_set))
        uni_space=c(uni_space,list(space))
        uni_prob=c(uni_prob,prob_region)
      }else{
        if(sum(colSums(null.space==0,na.rm=TRUE)>m/2)>0){
          max_gap=matrix(,1,2)
          max_gap[1,2]=which.max(colSums(null.space==0,na.rm=TRUE))[1]
          max_gap[1,1]=max(which(null.space[,max_gap[1,2]]==0))
        }else{
          #generate position of max gap
          max_gap=which(gap==max(gap,na.rm=TRUE),arr.ind = TRUE)
        }
        # generate partition point
        #par_decision includes (dimension,position,gap)
        par_decision=c(max_gap[1,2],partition[max_gap[1,1],max_gap[1,2]],gap[max_gap[1,1],max_gap[1,2]])
        
        #print(par_decision)===============
        partition_decision=rbind(partition_decision,par_decision)
        # two regions are generated and the domain is divided into two parts
        space_low=space_high=space
        #region in lower half
        space_low[2,par_decision[1]] = par_decision[2]
        #points in lower half
        indL = which(x_set[,par_decision[1]]<par_decision[2])
        point_low = x_set[indL,]
        prob_region_low=prob_region*nrow(point_low)/nrow(x_set)
        #region in higher half
        space_high[1,par_decision[1]]=par_decision[2]
        #points in higher half
        point_high=x_set[which(x_set[,par_decision[1]]>=par_decision[2]),]
        prob_region_high=prob_region-prob_region_low

        # store the points and space
        to_test_point=c(to_test_point, list(point_low,point_high))
        to_test_space=c(to_test_space,list(space_low,space_high))
        to_test_prob_region=c(to_test_prob_region,prob_region_low,prob_region_high)
      }
      
    }else{uni_point=c(uni_point,list(x_set))
    uni_space=c(uni_space,list(space))
    uni_prob=c(uni_prob,prob_region)
    }
    
    # when finishing all data point sets,\ set all parameter to default
    # renew data point sets
    if(index_point==length(data_points)){
      data_points=to_test_point
      data_space=to_test_space
      data_prob_region=to_test_prob_region
      to_test_point=to_test_space=to_test_prob_region=NULL
      index_point=0 
    }
    index_of_partition= index_of_partition+1
  }
  
  #sort hyper-rectangulars in decreasing order of density 
  p=rep(1,length(uni_prob)); 
  for(j in 1:sz[2]){ F=range(train.sim[,j]); 
  for(k in 1:length(uni_prob)){
    p[k]=p[k]/diff(uni_space[[k]][,j])/diff(F)
  }}
  den=uni_prob*p
  sort.prob=sort(den,decreasing=TRUE,index.return=TRUE)
  
  return(list(uni_point=uni_point[sort.prob$ix],uni_space=uni_space[sort.prob$ix],uni_prob=uni_prob[sort.prob$ix],uni_density=den[sort.prob$ix], partition_decision=partition_decision))
}

###################################################################################
#Credible set estimation using Algorithm 2

cred.est2 = function(test,train.sim,train,tau,alpha,p.alpha){
  
  y.out=matrix(,length(tau),2)
  i=1; N=nrow(test); p.alpha=(1-p.alpha)/2
  result = density_estimation(tau=tau[i], sz=dim(train),train.sim=train.sim,dataset=train)
  uni_space = result$uni_space; uni_prob = result$uni_prob; uni_point = result$uni_point
  ind=c(1:which.min(abs(cumsum(uni_prob)-(1-alpha))))
  cred.set=uni_space[ind]; 
  r=cov.prob(cred.set,test)
  del=abs(r$coverage_prob-(1-alpha))
  pp=pnorm(r$coverage_prob,(1-alpha),sqrt((1-alpha)*alpha/N))
  
  y.out[1,2]=r$coverage_prob	
  y.out[1,1]=sum((pp<p.alpha)|(pp>(1-p.alpha))) 
  
  for(i in 2:length(tau)){  
    result = density_estimation(tau=tau[i], sz=dim(train),train.sim=train.sim,dataset=train)
    uni_space = result$uni_space; uni_prob = result$uni_prob; uni_point = result$uni_point
    ind=c(1:which.min(abs(cumsum(uni_prob)-(1-alpha))))
    cred.set=uni_space[ind]; 
    r=cov.prob(cred.set,test)
    pp=pnorm(r$coverage_prob,(1-alpha),sqrt((1-alpha)*alpha/N))
    
    y.out[i,2]=r$coverage_prob
    y.out[i,1]=sum((pp<p.alpha)|(pp>(1-p.alpha)))
  }
  
  if(sum(y.out[,1]==0)==0){tau0=tau[which.min(abs(y.out[,2]-(1-alpha)))] 
  }else{tau0=min(tau[y.out[,1]==0])}  
  result = density_estimation(tau=tau0, sz=dim(train),train.sim=train.sim,dataset=train)
  uni_space = result$uni_space; uni_prob = result$uni_prob; uni_point = result$uni_point
  ind=c(1:which.min(abs(cumsum(uni_prob)-(1-alpha))))
  cred.set=uni_space[ind]; 
  r=cov.prob(cred.set,test) 
  
  return(list(cre.set=cred.set,n.tree=length(uni_space),tau.out=tau0,prob=r$coverage_prob,uni_space=uni_space,uni_prob=uni_prob,uni_density=result$uni_density,y.out=y.out))
}

###################################################################################

#Credible set estimation using Algorithm 1
cred.est1 = function(test,train.sim,train,tau,alpha,p.alpha,f.test,f.low){
  
  y.out=matrix(,length(tau),4)
  N=nrow(test); p.alpha=(1-p.alpha)/2
  result = density_estimation(tau=tau[1], sz=dim(train),train.sim=train.sim,dataset=train)
  uni_space = result$uni_space; uni_prob = result$uni_prob; uni_point = result$uni_point
  ind=c(1:which.min(abs(cumsum(uni_prob)-(1-alpha))))
  cred.set=uni_space[ind]; 
  r=cov.prob(cred.set,test)
  del=abs(r$coverage_prob-(1-alpha))
  pp=pnorm(r$coverage_prob,(1-alpha),sqrt((1-alpha)*alpha/N))
  
  y.out[1,2]=r$coverage_prob	
  y.out[1,1]=sum((pp<p.alpha)|(pp>(1-p.alpha)))
  r=cov.prob(cred.set,test[f.test<f.low,])
  y.out[1,3]=r$coverage_prob; 
  r=cov.prob(cred.set,test[f.test>f.low,])
  y.out[1,4]=1-r$coverage_prob; 
  
  for(i in 2:length(tau)){  
    result = density_estimation(tau=tau[i],sz=dim(train),train.sim=train.sim,data_points=uni_point,data_space=uni_space,data_prob_region=uni_prob)
    uni_space = result$uni_space; uni_prob = result$uni_prob; uni_point = result$uni_point
    ind=c(1:which.min(abs(cumsum(uni_prob)-(1-alpha))))
    cred.set=uni_space[ind]; 
    r=cov.prob(cred.set,test)
    del=abs(r$coverage_prob-(1-alpha))
    pp=pnorm(r$coverage_prob,(1-alpha),sqrt((1-alpha)*alpha/N))
    
    y.out[i,2]=r$coverage_prob
    y.out[i,1]=sum((pp<p.alpha)|(pp>(1-p.alpha)))
    r=cov.prob(cred.set,test[f.test<f.low,])
    y.out[i,3]=r$coverage_prob; 
    r=cov.prob(cred.set,test[f.test>f.low,])
    y.out[i,4]=1-r$coverage_prob;   
  }
  
  if(sum(y.out[,1]==0)>0){ i0=which(y.out[,1]==0)
  tau0=tau[i0[which.min(y.out[i0,3])]]
  }else{tau0=tau[which.min(mean(f.test<f.low)*y.out[,3]+mean(f.test>=f.low)*y.out[,4])]
  print("The coverage is poor and the total loss minimizer is used.")
  }
  
  result = density_estimation(tau=tau0, sz=dim(train),train.sim=train.sim,dataset=train)
  uni_space = result$uni_space; uni_prob = result$uni_prob; uni_point = result$uni_point
  ind=c(1:which.min(abs(cumsum(uni_prob)-(1-alpha))))
  cred.set=uni_space[ind]; 
  r=cov.prob(cred.set,test)
  
  return(list(cre.set=cred.set,n.tree=length(uni_space),tau.out=tau0,prob=r$coverage_prob,uni_space=uni_space,uni_prob=uni_prob,uni_density=result$uni_density,y.out=y.out))
}
