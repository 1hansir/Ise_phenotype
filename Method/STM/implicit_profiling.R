implicit.profile = function(data,delta,Z, KC, C,
                            init = rep(0,ncol(Z)), tol=1e-7,
                     maxit = 100, min.factor = 0.75,
                     ls.factor = 0.75, max.move = 1,
                     h.loop = T, glink = expit, 
                     dglink = dexpit)
{

  p = ncol(data$Z)
  init = glm(delta~log(C)+Z,data=data,family=binomial)$coef[3:(p+2)]
  if(any(is.na(init)))
    {
      data$Z = data$Z[,!is.na(init)]
      #Z.sd = Z.sd[!is.na(tmpinit)]
      p = ncol(data$Z)
      init = glm(delta~log(C)+Z,data=data,family=binomial)$coef[3:(p+2)]
    }


  start = Sys.time()
  n = nrow(Z)
  KCd = drop(KC%*%delta)
  hC = rep(0,n)
  oldscore = NULL

  for(k in 1:maxit) 
  {
    lp = drop(Z%*%init)
    hC.flag = rep(TRUE,n)
    gij = wZbar = matrix(0,n,n)
    hHess = rep(0,n)
    for(kk in 1:ifelse(h.loop, maxit, 1))
    {
      # print(hC[hC.flag])
      lp.flag = outer(hC[hC.flag],lp,"+")
      gij[hC.flag,] = glink(lp.flag)
      tmp = KC[hC.flag,]*gij[hC.flag,]
      wZbar[hC.flag,] = KC[hC.flag,]*dglink(lp.flag)
      if(sum(hC.flag)>=2)
      {
        hscore = apply(tmp,1,sum)-KCd[hC.flag]
        hHess[hC.flag] = apply(wZbar[hC.flag,],1,sum)
      }else
      {
        hscore = sum(tmp)-KCd[hC.flag]
        hHess[hC.flag] = sum(wZbar[hC.flag,])
      }
      
      dhC = hscore/hHess[hC.flag]
      dhC = sign(dhC)*pmin(abs(dhC),max.move)
      kk.flag = abs(hscore) > tol
      if(!any(kk.flag))
        break
      hC[hC.flag][kk.flag] = hC[hC.flag][kk.flag] - dhC[kk.flag]
      hC.flag[hC.flag] = kk.flag
    }
    if(kk >= maxit)
      stop("Numerical error when computing h0(Ci)")
    Zbar =  (wZbar%*%Z) / hHess 
    
    gi = glink(hC+lp)
    bscore = drop(t(Z)%*% (delta - gi))
    if(!is.null(oldscore))
      if(((sum(oldscore^2)*min.factor) <= sum(bscore^2)))
      {
        init = init+dinit
        dinit = dinit*ls.factor
        if(max(abs(dinit))<tol)
        {
          if(max(abs(oldscore)) > 1e-6)
            warning(paste("Algorithm stops in line-search. Target tol: ",
                        tol, ". Current tol: ", max(abs(oldscore)),
                        ". ", sep = ''))
          break
        }
        init = init - dinit
        #print(init)
        next
      }
    oldscore = bscore
    bHess = t(dglink(hC+lp)*Z) %*% (Zbar-Z)
    dinit = solve(bHess,bscore)
    if(all(abs(bscore)<tol))
      break
    #print(rbind(init,bscore,dinit))
    #print(init)
    init = init - dinit
  }
  if(k >=maxit)
    stop("Numerical error when computing beta_delta")

  C.order = order(C)
  hfun = stepfun(C[C.order],c(min(hC),hC[C.order]))
  
  run.time = Sys.time() - start
  return(list(beta = init, 
              hfun= hfun,
              run.time = run.time, 
              step = k))
}

