
# Master function to compute measures from MATA paper
#------------------------------------------------------------------
# inc.prob: matrix of predicted incidence probabilities
# fu: vector of follow up times
# Xobs: vector of observation times from annotation
# Dobs: vector of event status from annotation
# time.grid: vector of grid for predicted event probability (shared)
# targetFPR: the target FPR to select cutoff
MASTA.meas = function(inc.prob, fu, Xobs, Dobs, time.grid, targetFPR)
{
  if(missing(time.grid))
  {
    time.grid = as.numeric(colnames(inc.prob))
  }
  DT = pred.DT(fu, time.grid, inc.prob)
  thres.list = quantile(DT$piC, probs = seq(0,1,0.01))
  Xhat = pred.Xhat(fu, DT, thres.list)
  Dhat = outer(DT$piC, thres.list, ">=")
  if(missing(targetFPR))
  {
    APE = pred.APE(Xobs, Xhat, thres.list)
    u = APE$best.cut
  }else{
    FPR = match.FPR(Dhat,Dobs, thres.list, targetFPR)
    u = FPR$match
  }
  D.best = as.numeric(DT$piC >= u)
  X.best = Xhat[,which.min(abs(thres.list-u))]
  return(pred.u(D.best, X.best, Dobs, Xobs))
}

pred.DT = function(fu, time.grid, inc.prob)
{
  n = length(fu)
  K = length(time.grid)
  out = data.frame(piC = rep(NA, n),
                   That = rep(NA, n))
  time.diff = c(diff(time.grid),0)
  
  for (i in 1:n) 
  {
    bf.C = time.grid <= fu[i]
    C.pos = sum(bf.C)
    if(C.pos == 0)
    {
      out$piC[i] = 0 
      out$That[i] = fu[i]
      next
    } 

      # print(c(i,C.pos))
    out$piC[i] = inc.prob[i,C.pos]
    out$That[i] = sum((1-inc.prob[i,bf.C])*time.diff[bf.C])
    if(C.pos < K)
    {
      out$That[i] =  out$That[i] - (1-inc.prob[i,C.pos])*(time.grid[C.pos+1]-fu[i])
    }
  }
  
  return(out)
}

pred.Xhat = function(fu, DT, cut.off)
{
  Du = outer(DT$piC, cut.off, ">=")
  out = fu*(1-Du) + Du*DT$That
  colnames(out) = cut.off
    
  return(out)
}

pred.APE = function(Xobs, Xhat, cut.off)
{
  ape.tab = data.frame(cut.off = cut.off, 
                       APE = apply(abs(Xobs-Xhat), 2, mean))
  best.cut = cut.off[which.min(ape.tab$APE)]
  best.ape = min(ape.tab$APE)
  
  return(list(best = best.ape, 
              best.cut = best.cut, 
              ape = ape.tab))
}

match.FPR = function(Dhat, Dobs, cut.off, targetFPR)
{
  fpr.tab = data.frame(cut.off = cut.off, 
                       FPR = apply(Dhat*(1-Dobs), 2, mean)/(1-mean(Dobs)))
  match.fpr = cut.off[which.min(abs(fpr.tab$FPR - targetFPR))]
  
  return(list(match = match.fpr, 
              fpr = fpr.tab))
}

pred.u = function(Du,Xu, Dobs,Xobs, symmetry = TRUE)
{
  TPR = sum(Du*Dobs)/sum(Dobs)
  FPR = sum(Du*(1-Dobs))/sum(1-Dobs)
  CCR = mean(Du == Dobs)
  APE = mean(abs(Xu-Xobs))
  
  tri.Xobs = outer(Xobs,Xobs, "<=")
  if(!symmetry)
  {
    tri.Xobs = lower.tri(tri.Xobs)*tri.Xobs
  }
  tri.XX = tri.Xobs*outer(Xu,Xu, "<=")
  Cu = sum(tri.XX)/sum(tri.Xobs)
  
  tri.Xobs = outer(Xobs,Xobs, "<=")*Dobs
  if(!symmetry)
  {
    tri.Xobs = lower.tri(tri.Xobs)*tri.Xobs
  }
  tri.XX = tri.Xobs*outer(Xu,Xu, "<=")*Du
  Cup = sum(tri.XX)/sum(tri.Xobs)
  
  out = c(TPR,FPR, CCR,  Cu, Cup, APE)
  names(out) = c("TPRu", "FPRu", "CCRu", "Cu", "Cu+", "APEu")
  return(out)
}
