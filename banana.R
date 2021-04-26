banana = function(A=0.5, B=0, C1=3, C2=3, N, keep=10, init=10) {
  R = init*keep + N*keep
  x1 = x2 = 0
  bimat = matrix(double(2*N), ncol=2)
  for (r in 1:R) {
    x1 = rnorm(1,mean=(B*x2+C1) / (A*(x2^2)+1), sd=sqrt(1/(A*(x2^2)+1)))
    x2 = rnorm(1,mean=(B*x2+C2) / (A*(x1^2)+1), sd=sqrt(1/(A*(x1^2)+1)))
    if (r>init*keep && r%%keep==0) {
      mkeep = r/keep
      bimat[mkeep-init,] = c(x1,x2)
    }
  }
  return(bimat)
}

d.banana = function(x,A=0.5, B=0, C1=3, C2=3) {
  f = dnorm(x[,1],mean=(B*x[,2]+C1) / (A*(x[,2]^2)+1), sd=sqrt(1/(A*(x[,2]^2)+1)))
  f = f*dnorm(x[,2],mean=(B*x[,2]+C2) / (A*(x[,1]^2)+1), sd=sqrt(1/(A*(x[,1]^2)+1)))
  return(f)
}
