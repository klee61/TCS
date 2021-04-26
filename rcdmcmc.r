
#load the data
rcc.dat<-read.table("http://www.stats.ox.ac.uk/~nicholls/BayesMethods/SHRCC2013.txt",sep=",",header=TRUE)
str(rcc.dat)
attach(rcc.dat)

Earliest=1000; Latest=1;
y.BP=Latest:Earliest
mu=spline(CAL.BP, X14C.age, xmin = 0.9*Latest, xmax = 1.1*Earliest, xout=y.BP)$y
err=spline(CAL.BP,Error, xmin = 0.9*Latest, xmax = 1.1*Earliest, xout=y.BP)$y

#Data - come from a Moari rivermouth settlement in NZ
#NZ 7758 (1, 1) 580 47
#NZ 7761 (2, 1) 600 50
#NZ 7757 (3, 1) 537 44
#NZ 7756 (4, 1) 670 47
#NZ 7755 (5, 1) 646 47
#WK 2589 (5, 2) 630 35
#NZ 7771 (6, 1) 660 46

#radiocarbon dates
y=c(580,600,537,670,646,630,660); nd=length(y)
#measurement errors
d=c(47,50,44,47,47,35,46)
#layers
s=c(1,2,3,4,5,5,6)
#the laboratory date identifiers
nm=c("NZ 7758","NZ 7761","NZ 7757","NZ 7756","NZ 7755","WK 2589","NZ 7771")

#the log likelihood for calendar dates theta
llk<-function(theta) {
  sig=sqrt(d^2+err[theta]^2)
  llk=sum(dnorm(y,mean=mu[theta],sd=sig,log=T))
  return(llk)
}
#the log prior for the shrinkage prior
lpr<-function(psi) {
  S=psi[2]-psi[1]
  lpr=-nd*log(S)-log(U-L-S)
  return(lpr)
}

#these are very conservative bounds on the settlement dates
#current best estimates suggest NZ was only settled around 700 BP 
#so a lower bound of 1000 is quite conservative. The upper bound
#comes from information independent of the dates
L=500; U=1000; x=L:U;

#MCMC sampler very simple
rcd.mcmc<-function(nd,L,U,T,Subs,shrink=TRUE,w,scale.move=TRUE,beta) {
  #initialise the start state for the MCMC
  theta=sort(runif(nd,min=L,max=U))
  psi=c((L+theta[1])/2,(theta[nd]+U)/2)
  llo=llk(theta)*beta
  lpo=lpr(psi)
  OP=matrix(NA,T/SubS,1+2+nd+2)
  colnames(OP)<-c('span','psi1','psi2','thet1','thet2','thet3','thet4','thet5','thet6','thet7','llkd','lpri')
  
  #run the MCMC T steps
  for (t in 1:T) {
    #go through each of the parameters in turn
    for (i in 1:nd) {
      thetap=theta; 
      #the proposal is uniform in allowed range
      if (shrink) {
        thetap[i]=runif(1,min=psi[1],max=psi[2]);
      } else {
        thetap[i]=runif(1,min=L,max=U)
      }
      lln=llk(thetap)*beta; #no change to prior
      logMHR=lln-llo
      if (log(runif(1))<logMHR) {
        theta=thetap
        llo=lln; 
      }
    }
    
    #if we are using the shrinkage prior (psi etc) update the extra parameters
    if (shrink) {
      psip=psi; psip[1]=runif(1,min=L,max=min(theta));
      lpn=lpr(psip)
      logMHR=lpn-lpo
      if (log(runif(1))<logMHR) {
        psi=psip
        lpo=lpn
      }
      psip=psi; psip[2]=runif(1,min=max(theta),max=U);
      lpn=lpr(psip)
      logMHR=lpn-lpo
      if (log(runif(1))<logMHR) {
        psi=psip
        lpo=lpn
      }
      
      if (scale.move) {
        
        delta=runif(1,w,1/w)
        m=mean(c(psi,theta))
        psip=m+(psi-m)*delta
        thetap=m+(theta-m)*delta
        if ( (psip[2]<U) & (psip[1]>L) ) {
          lpn=lpr(psip)
          lln=llk(thetap)*beta
          lqq=(nd-1)*log(delta)
          logMHR=lln-llo+lpn-lpo+lqq
          if (log(runif(1))<logMHR) {
            psi=psip
            theta=thetap
            lpo=lpn
            llo=lln
          }
        }
      }
    }
    #subsample - take a sample every SubS steps - store it in OutPut
    if (!t%%SubS) {
      span=diff(psi)
      OP[t/SubS,]=c(span,psi,theta,llo/beta,lpo)
    }
  }
  return(OP)
}

