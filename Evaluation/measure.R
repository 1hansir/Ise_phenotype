require(pROC)
# event = valid.rec
# fu = last.Month
# time.grid = time.MASTA
# inc.prob = pred.MASTA

#
# Some measures I invented
#==================================================================

# Time specific AUC
#------------------------------------------------------------------
# event: vector of event times: event times from full labeled validation dataset
#For censored patients, max(fu.time) + 100
# fu: vector of follow up times: follow-up times from full labeled dataset
# time.grid: vector of grid for predicted event probability (shared): 0:500, e.g. (or largest follow-up time); spacing of 0.1, 1, 5, etc.
# inc.prob: matrix of predicted incidence probabilities : ROWs, validation patients, COLUMNS: preds for each time grid point
# For new patients, predict their cumulative prob. dist. ( f.hat(time.grid.point,X) )
# Write prediction function using object
auc.tspec = function(event, fu, time.grid, inc.prob)
{
  nt = length(time.grid)
  auc.t = rep(0, nt)

  for(i in 1:nt)
  {
    Y.dich = as.numeric(event <= time.grid[i])
    atrisk =  time.grid[i] <= fu

    if(sum(Y.dich[atrisk]) > 0 & sum(Y.dich[atrisk]) < length(Y.dich[atrisk])){
      roc.MASTA = roc(Y.dich[atrisk], inc.prob[atrisk,i],
                     levels = c(0,1), direction = "<")
      auc.t[i] = roc.MASTA$auc
    }else{auc.t[i]=NA}
  }
  return(auc.t)
}

# auc.tspec = function(event, fu, time.grid, inc.prob)
# {
#   nt = length(time.grid)
#   auc.t = rep(0, nt)
#   
#   for(i in 1:nt)
#   {
#     Y.dich = as.numeric(event <= time.grid[i])
#     atrisk =  time.grid[i] <= fu
#     
#     roc.MASTA = roc(Y.dich[atrisk], inc.prob[atrisk,i], 
#                     levels = c(0,1), direction = "<")
#     auc.t[i] = roc.MASTA$auc
#   }
#   return(auc.t)
# }
# 
# 


# Time specific PPV at given cutoff
#------------------------------------------------------------------
# event: vector of event times
# fu: vector of follow up times
# time.grid: vector of grid for predicted event probability (shared)
# inc.prob: matrix of predicted incidence probabilities
ppv.tspec = function(event, fu, time.grid, inc.prob, spec = 0.95)
{
  nt = length(time.grid)
  ppv.t = rep(0, nt)
  
  for(i in 1:nt)
  {
    Y.dich = as.numeric(event <= time.grid[i])
    atrisk =  time.grid[i] <= fu
    
    roc.MASTA = roc(Y.dich[atrisk], inc.prob[atrisk,i], 
                    levels = c(0,1), direction = "<")
    if(all(roc.MASTA$specificities<spec))
    {
      warning("No match for spericificties")
    }
    spec.thres = roc.MASTA$thresholds[min(which(roc.MASTA$specificities>=spec))]
    ppv.t[i] = mean(Y.dich[(inc.prob[,i] >= spec.thres)&atrisk])
  }
  return(ppv.t)
}

# Individual level area between curves (observed counting process vs predicted)
#-------------------------------------------------------------------
# event: vector of event times
# fu: vector of follow up times
# time.grid: vector of grid for predicted event probability (shared)
# inc.prob: matrix of predicted incidence probabilities
AbC.prob = function(event, fu, time.grid, inc.prob)
{
  n = length(event)
  ngrid = length(time.grid)
  AbC = rep(NA, n)
  
  for (i in 1:n) 
  {
    bf.event = time.grid <= event[i]
    aft.event = (!bf.event) & (time.grid <= fu[i])
    bf.pos = sum(bf.event)
    AbC[i] = sum(c(inc.prob[i,bf.event]*diff(c(0,time.grid[bf.event])), 
        inc.prob[i,bf.pos]*(event[i] - time.grid[bf.pos])))
    if(any(aft.event))
    {
        sum(c((1-inc.prob[i,bf.pos])*(time.grid[bf.pos+1] - event[i]), 
        (1-inc.prob[i,aft.event])*diff(c(time.grid[aft.event], fu[i]))))
    }
  }

  return(AbC)
}

#
# Measures from Yuri's paper
#==================================================================

# Pooled AUC, weighted by window length (as if no weights if all
#  windows are length 1)
#-------------------------------------------------------------------
# event: vector of event times
# fu: vector of follow up times
# time.grid: vector of grid for predicted event probability (shared)
# inc.prob: matrix of predicted incidence probabilities
auc.pooled = function(event, fu, time.grid, inc.prob,
                      weight.time = TRUE)
{
  n = length(event)
  
  Y.dich = as.numeric(outer(event,time.grid, "<="))
  if(weight.time)
  {
    wgt = as.numeric(outer(fu,time.grid, ">="))*rep(diff(c(0,time.grid)),each = n)
  }else{
    wgt = as.numeric(outer(fu,time.grid, ">="))
  }
  roc.fit = roc.wgt(Y.dich, as.numeric(unlist(inc.prob)), wgt)
  
  return(roc.fit$auc)
}

# Pooled F1 score, weighted by window length (as if no weights if all
#  windows are length 1)
#-------------------------------------------------------------------
# event: vector of event times
# fu: vector of follow up times
# time.grid: vector of grid for predicted event probability (shared)
# inc.prob: matrix of predicted incidence probabilities
F1.pooled = function(event, fu, time.grid, inc.prob, thres, spec = .95)
{ 
  n = length(event)
  Y.dich = as.numeric(outer(event,time.grid, "<="))
  wgt = as.numeric(outer(fu,time.grid, ">="))*rep(diff(c(0,time.grid)),each = n)
  wgt = wgt/sum(wgt)
  
  if(missing(thres))
  {
    roc.fit = roc.wgt(Y.dich, as.numeric(unlist(inc.prob)), wgt)
    thres = min(roc.fit$thres[roc.fit$specificity >= spec])
  }
  pred = as.numeric(inc.prob >= thres)
  TP = sum(pred*Y.dich*wgt)
  TN = sum((1-pred)*(1-Y.dich)*wgt)
  P = sum(Y.dich*wgt)
  N = sum((1-Y.dich)*wgt)
  
  return(2/(P/TP+N/TN))
}

# Individual level area between curves (observed counting process vs 
#  dichotomized predicted probability)
#-------------------------------------------------------------------
# event: vector of event times
# fu: vector of follow up times
# time.grid: vector of grid for predicted event probability (shared)
# inc.prob: matrix of predicted incidence probabilities
AbC.dich = function(event, fu, time.grid, inc.prob, thres)
{
  n = length(event)
  ngrid = length(time.grid)
  if(missing(thres))
  {
    Y.dich = as.numeric(outer(event,time.grid, "<="))
    wgt = as.numeric(outer(fu,time.grid, ">="))*rep(diff(c(0,time.grid)),each = n)
    wgt = wgt/sum(wgt)
    
    prev = sum(Y.dich*wgt)
    
    pred = as.numeric(inc.prob)
    order.pred = order(pred, decreasing = TRUE)
    
    pred.prev = 0
    for(i in 1:length(order.pred))
    {
      thres = pred[order.pred[i]]
      pred.prev = pred.prev + wgt[order.pred[i]]
      if(pred.prev > prev)
        break
    }
  }
  
  inc.dich = (inc.prob >= thres)
  
  AbC = rep(NA, n)
  
  for (i in 1:n) 
  {
    bf.event = time.grid <= event[i]
    aft.event = (!bf.event) & (time.grid <= fu[i])
    bf.pos = sum(bf.event)
    AbC[i] = sum(c(inc.dich[i,bf.event]*diff(c(0,time.grid[bf.event])), 
                   inc.dich[i,bf.pos]*(event[i] - time.grid[bf.pos])))
    if(any(aft.event))
    {
      sum(c((1-inc.dich[i,bf.pos])*(time.grid[bf.pos+1] - event[i]), 
            (1-inc.dich[i,aft.event])*diff(c(time.grid[aft.event], fu[i]))))
    }
  }
  
  return(AbC)
}

# Null prediction
#-------------------------------------------------------------------
# event: vector of event times
# fu: vector of follow up times
# time.grid: vector of grid for predicted event probability (shared)
# inc.prob: matrix of predicted incidence probabilities
inc.null = function(event, fu, time.grid)
{
  n = length(event)
  Y.dich = outer(event,time.grid, "<=")
  wgt = outer(fu,time.grid, ">=")*matrix(rep(diff(c(0,time.grid)),each = n),n)
  
  return(matrix(
    rep(apply(Y.dich*wgt, 2, sum)/apply(wgt, 2, sum),
        each = n), n))
}