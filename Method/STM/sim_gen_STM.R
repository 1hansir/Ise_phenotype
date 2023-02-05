require(MASS)
sim.gen.STM = function(n,p,beta0,sigmaZ,Nperturb=0)
{
  Z = mvrnorm(n,rep(0,p),sigmaZ^2*(0.2+0.8*diag(1,p)))
  u = runif(n)
  T = exp((log(u/(1-u)) - Z%*%beta0)/3)*4 #h^-1 (g^-1(u) - beta'Z)
  C = runif(n,0,12)
  delta = (T<=C)
  
  return(list(delta = delta, C = C, 
              Z = Z))
}