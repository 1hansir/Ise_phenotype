rm(list = objects())

# .libPaths(c("~/R-3.6.1/library",.libPaths()))

library(MASS)
library(survival)
library(doParallel)
library(doRNG)
library(glmnet)

source("source/utility.R")
source("source/beta_delta.R")
source("source/beta_delta_perturb.R")
source("source/Sk_sym.R")
source("source/Sk_sym_perturb.R")
source("source/W_hat_adaPCA_v2.R")
source("source/cv_concordance.R")


np = 10
cl = makeCluster(getOption("cl.cores", np))
registerDoParallel(cl)
# clusterEvalQ(cl, .libPaths(c("~/R-3.6.1/library",.libPaths())))

out.dir = "dataresult/obesity-MAP/"
link = expit
dlink = dexpit
Nperturb = 1000
  
  load("data/obesity-MAP-popPCA-oldage.rda")
  data$IPW = data$ipw
  n = nrow(data)
  
  Z.sd = apply(data$Z, 2, sd)
  data$Z = scale(data$Z)
  p = ncol(data$Z)

  set.seed(531)
  registerDoRNG(seed = 525)

  tmp.save = paste0(out.dir,"SSL_pca_tmp_IPW_oldage.rda")
  if(file.exists(tmp.save))
  {
    load.list <- load(tmp.save)
  }else{
    load.list <- NULL
  }
  
  
  #    beta_delta and its perturbation
  #------------------------------------------
  
  if(!any("betadelta" == load.list))
  {
    
    print("Initial estimation...")
    start.time = Sys.time()
    fold.pos = 1:m
    h = sd(data$C[fold.pos])/(sum(data$delta[fold.pos]))^0.25
    KC = dnorm(as.matrix(dist(data$C[fold.pos]/h,diag=T,upper=T)))/h
    std.wgt = data$IPW[fold.pos]/mean(data$IPW[fold.pos])
    
    tmpinit = glm(delta~log(C)+Z,data=data,family=binomial,
                  weights = IPW, subset = fold.pos)$coef[3:(p+2)]
    if(any(is.na(tmpinit)))
    {
      data$Z = data$Z[,!is.na(tmpinit)]
      Z.sd = Z.sd[!is.na(tmpinit)]
      p = ncol(data$Z)
      npc = p-1 
      tmpinit = glm(delta~log(C)+Z,data=data,family=binomial)$coef[3:(p+2)]
    }
        
    betadelta = init.beta.perturb(data$delta[fold.pos],data$Z[fold.pos,],KC,
                          init = tmpinit,
                          V = std.wgt,
                          maxit = 1000,
                          min.factor = 1, tol=1e-8)
    
    run.time = Sys.time() - start.time
    print(paste("Finished in ", format(run.time, digits = 0, nsmall = 2)))
    
    save(betadelta, h, 
         file = tmp.save)
  }
  
  
  #    Sk
  #------------------------------------------
  if(!any("Sk" == load.list))
  {
    print("Sk...")
    start.time = Sys.time()
    Sk = rep(0,p*2)
    beta.std = betadelta/sqrt(sum(betadelta^2))
    lp = drop(data$Z %*% beta.std)
    sdlp = sd(lp)
    h1 = sdlp/(sum(data$delta1))^0.3
    Sk[1:p] = Sk_sym(matrix(lp),data$Z,data$X1,data$delta1,
                     data$Cstar,dnorm,h1)
    h2 = sdlp/(sum(data$delta2))^0.3
    Sk[p+1:p] = Sk_sym(matrix(lp),data$Z,data$X2,data$delta2,
                       data$Cstar,dnorm,h2)
    
    run.time = Sys.time() - start.time
    print(paste("Finished in ", format(run.time, digits = 0, nsmall = 2)))
    
    save(betadelta,  h, 
         Sk, h1, h2, 
         file = tmp.save)
  }
  #    perturbation
  #------------------------------------------
  if(!any("Skb" == load.list))
  {
    print("Perturbation...")
    start.time = Sys.time()
    fold.pos = 1:m
    KC = dnorm(as.matrix(dist(data$C[fold.pos]/h,diag=T,upper=T)))/h
    tmpout<- foreach(iperturb = 1:Nperturb,.packages = c("MASS","glmnet")) %dopar%
      {
        V = rbeta(n,0.5,1.5)*4
        div = TRUE
        while(div)
        {
          div = FALSE
          std.weight = (V[fold.pos]*data$IPW[fold.pos])/mean(V[fold.pos]*data$IPW[fold.pos])
          fit = try(betak <-
                      init.beta.perturb(data$delta[fold.pos],data$Z[fold.pos,],KC,
                                        std.weight, 
                                        init = betadelta,
                                        maxit = 1000,
                                        min.factor = 1, tol=1e-9))
          
          if(class(fit) == "try-error")
          {
            print(paste("Error in iteration ", iperturb))
            div = TRUE
            V = rbeta(n,0.5,1.5)*4
            break
          }
        }
        
        betak.std = betak/sqrt(sum(betak^2))
        lp = data$Z %*% betak.std
        
        Skb = c(Sk_sym_perturb(lp,data$Z,data$X1,data$delta1,
                               data$Cstar,dnorm,h1,
                               matrix(V)),
                Sk_sym_perturb(lp,data$Z,data$X2,data$delta2,
                               data$Cstar,dnorm,h2,
                               matrix(V)))
        return(list(betak = betak, Skb=Skb))
      }
    
    Skb =  do.call(rbind,sapply(tmpout, '[',"Skb"))
    betak =  do.call(cbind,sapply(tmpout, '[',"betak"))
    
    rm(tmpout)
    run.time = Sys.time() - start.time
    print(paste("Finished in ", format(run.time, digits = 0, nsmall = 2)))
    
    save(betadelta, h,  
         Sk, h1, h2, 
         Skb, betak, 
         file = tmp.save)
  }
  
  #    Estimating W matrix
  #------------------------------------------
  if(!any("W.hat" == load.list))
  {
    print("What...")
    start.time = Sys.time()
    W.hat = W_hat_adaPCA_ridge(betak,betadelta, Skb)$W.hat
    
    run.time = Sys.time() - start.time
    print(paste("Finished in ", format(run.time, digits = 0, nsmall = 2)))
    
    save(betadelta, h,  
         Sk, h1, h2, 
         Skb, betak, 
         W.hat, 
         file = tmp.save)
  }
  
  #------------------------------------------
  # CV for adaptive soft-thresholding initial
  #------------------------------------------ 
  
  Nfold = 5
  Nlam = 100
  min.lam = 1e-4
  
  
  if(!any("lp.init.cv" == load.list))
  {
    print("CV for initial estimation...")
    start.time = Sys.time()
    init.cv = matrix(0,p,Nfold)
    foldid = rep(0,n)
    pos1 = 1:m
    fold.pos = 1:m
    h = sd(data$C[fold.pos])/(sum(data$delta[fold.pos]))^0.25
    KC = dnorm(as.matrix(dist(data$C[fold.pos]/h,diag=T,upper=T)))/h
    
    div = TRUE
    while (div) 
    {
      div = FALSE
      foldid[pos1] = rep_len(1:Nfold,length(pos1))[sample(length(pos1))]
      foldid[-pos1] = rep_len(1:Nfold,n-length(pos1))[sample(n-length(pos1))]
      lp.init.cv = matrix(0,n,Nfold)
      for(ifold in 1:Nfold)
      {
        # i.label = which(foldid[1:m] != ifold)
        i.train = which(foldid[1:m] != ifold)
        
        fit = try(init.cv[,ifold] <- init.beta.perturb(data$delta[i.train],data$Z[i.train,],
                                                       KC[i.train,i.train],
                                                       init = betadelta,
                                                       V = data$IPW[i.train], 
                                                       maxit = 1000,
                                                       min.factor = 1, tol=1e-8)
        )
        
        if(class(fit) == "try-error")
        {
          print(paste("Error in fold ", ifold))
          div = TRUE
          break
        }
        lp.init.cv[,ifold] = drop(data$Z %*% init.cv[,ifold])
      }
    }
    
    run.time = Sys.time() - start.time
    print(paste("Finished in ", format(run.time, digits = 0, nsmall = 2)))
    
    save(betadelta, h,  
         Sk, h1, h2, 
         Skb, betak, 
         W.hat, 
         lp.init.cv, foldid, init.cv, 
         file = tmp.save)
  }
  # 
  # ada.factor = 1/(abs(betadelta)*apply(data$Z, 2, sd))
  # lambda = c(0,exp(seq(log(min.lam),log(min(max(betadelta^2),m^(-.75)))
  #                      ,length.out = Nlam))[-Nlam])
  # 
  # #beta.cv = matrix(0,p,Nfold)
  # cv.concord = matrix(0,Nfold,Nlam)
  # for(ifold in 1:Nfold)
  # {
  #   i.train = which(foldid != ifold)
  #   i.test = which(foldid == ifold)
  #   beta.cv = init.cv[,ifold] 
  #   beta.thres.cv = sign(beta.cv)*matrix(pmax(0,abs(beta.cv)-outer(ada.factor,lambda,"*")),p)
  #   lp.test = data$Z[i.test,] %*% beta.thres.cv
  #   tmp.concord = (cv.concordance(lp.test, 
  #                                 data$X1[i.test],data$delta1[i.test],
  #                                 data$Cstar[i.test]) + 
  #                    cv.concordance(lp.test, 
  #                                   data$X2[i.test],data$delta2[i.test],
  #                                   data$Cstar[i.test]))
  #   cv.concord[ifold,] = tmp.concord[,1]/tmp.concord[,2]
  # }
  # 
  # best.ilam = which.max(apply(cv.concord,2,mean))
  # lambda.best = lambda[best.ilam]
  # thres.best = lambda.best * ada.factor
  # betadelta.thres = sign(betadelta) * pmax(0,abs(betadelta) - thres.best)
  
  
  # CV for adaptive soft-thresholding
  #------------------------------------------
  
  
  betaG = betadelta - drop(W.hat %*% Sk)
  ada.factor = 1/abs(betaG)
  # ada.factor = rep(1,p)
  max.lambda = m^(-0.75)
  lambda = c(0,exp(seq(log(min.lam),log(max.lambda),length.out = Nlam-1)))
  
  # beta.std = betadelta/sqrt(sum(betadelta^2))
  # lp = drop(data$Z %*% beta.std)
  # sdlp = sd(lp)
  # h1 = sdlp/(sum(data$delta1))^0.3
  # h2 = sdlp/(sum(data$delta2))^0.3
  
  cv.concord = matrix(0,Nfold,Nlam)
  for(ifold in 1:Nfold)
  {
    i.train = which(foldid != ifold)
    i.test = which(foldid == ifold)
    Sk.cv = rep(0,p*2)
    Sk.cv[1:p] = Sk_sym(lp.init.cv[i.train,ifold],data$Z[i.train,],
                        data$X1[i.train],data$delta1[i.train],
                        data$Cstar[i.train],dnorm,h1)
    Sk.cv[p+1:p] = Sk_sym(lp.init.cv[i.train,ifold],data$Z[i.train,],
                          data$X2[i.train],data$delta2[i.train],
                          data$Cstar[i.train],dnorm,h2)
    beta.cv = init.cv[,ifold] - drop (W.hat %*% Sk.cv)
    beta.thres.cv = sign(beta.cv)*matrix(pmax(0,abs(beta.cv)-outer(ada.factor,lambda,"*")),p)
    lp.test = data$Z[i.test,] %*% beta.thres.cv
    tmp.concord = (cv.concordance(lp.test, 
                                  data$X1[i.test],data$delta1[i.test],
                                  data$Cstar[i.test]) + 
                     cv.concordance(lp.test, 
                                    data$X2[i.test],data$delta2[i.test],
                                    data$Cstar[i.test]))
    cv.concord[ifold,] = tmp.concord[,1]/tmp.concord[,2]
  }
  
  best.ilam = which.max(apply(cv.concord,2,mean))
  lambda.best = lambda[best.ilam]
  
  
  # One-step estimator and its permutation
  #------------------------------------------
  
  beta.init = betadelta
  betas.init = betak
  betaG = (
    betadelta - drop(W.hat %*% Sk)
  )
  betaGs = (
    betak - W.hat %*% t(Skb)
  )
  thres.best = 1/abs(betaG) * lambda.best
  betaG.thres = sign(betaG) * pmax(0,abs(betaG) - thres.best)
  betaGs.thres = sign(betaGs) * pmax(0,abs(betaGs) - thres.best)
  betaG.thres.std = (betaG.thres/
                       sqrt(sum(betaG.thres^2))*
                       sqrt(sum(betaG^2)))
  
  betaGs.thres.std = (betaGs.thres/
                        rep(sqrt(apply(betaGs.thres^2,2,sum))/
                              sqrt(apply(betaGs^2,2,sum)),each = p))
  
  thresGs = 1/abs(betaGs) * lambda.best
  betaGs.thres2 = sign(betaGs) * pmax(0,abs(betaGs) - thresGs)
  betaGs.thres2.std = (betaGs.thres2/
                         rep(sqrt(apply(betaGs.thres2^2,2,sum))/
                               sqrt(apply(betaGs^2,2,sum)),each = p))

    save(betaG, betaGs, beta.init,betas.init,
         betaG.thres,betaGs.thres,betaGs.thres2,
         betaG.thres.std, betaGs.thres.std,betaGs.thres2.std,
         Z.sd, 
         file = paste0(out.dir,"SSL_MAP_pca_IPW_oldage.rda"))



stopCluster(cl)

rm(list = objects())
