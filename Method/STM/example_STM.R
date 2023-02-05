source("sim_gen_STM.R")

# simulation parameters
# n = 500
# p = 10
# beta0=c(0.7,0.7,0.7,-0.5,-0.5,-0.5,0.3,0.3,0.3,0)
n = 1500
p = 3
beta0=c(0.7,-0.5,0.3)
sigmaZ = 1

Nrep = 10000

# Data generation
dat = sim.gen.STM(n,p,beta0,sigmaZ)

# You implement of two algorithms
# According to my code, the implicit profiling algorithm 
# should be twice faster with half number of iterations. 

source("utility.R")
source("implicit_profiling.R")

# Bandwidth and kernel
h = sd(dat$C)/(sum(dat$delta))^0.25
KC = dnorm(as.matrix(dist(dat$C/h,diag=T,upper=T)))/h

ip.fit = implicit.profile(dat$delta,dat$Z,KC, dat$C, h.loop = F)

# Create the prediction matrix
source("pred_STM.R")
pred.mat = pred.STM(ip.fit, dat$Z, tgrid = seq(0,10,0.1))
