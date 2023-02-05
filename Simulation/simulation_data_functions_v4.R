#### load packges and functions ####
library(MASS)
library(pROC)
library(dplyr)
source(paste0(mdir,"/Method/MAP_S/MAPflex.R"))
source(paste0(mdir,"/Method/MAP_S/library_v2.R"))
source(paste0(mdir,"/Method/MAP_S/MAP_0.1.0/MAP/R/FUN_NLP_PheWAS.R"))
# source(("addMAP_1mo_timeWindow.R"))
#==================================================================================
# n: number of observations
# p: number of features

##### Step1_generateX #####
# @family: the type of X-features
# c(rep("binary",4), rep("count",4), rep("numeric",2)) indicates 4 binary, 4 counting and 2 numeric covariates
# @rho: correlation matrix of MVN = (1-rho)I + rho
# @coef: coefs for transforming MVN to Binary or Numeric
#----------------------------------------------------------
generate_X = function(n, p = 10, family = c(rep("binary",4), rep("count",4), rep("numeric",2)),
                      rho = 0.1, coef = c(0.3,2)){
  # generate MVN
  Sigma = (1-rho)*diag(1,nrow = p) + rho*matrix(1, p, p)
  X = mvrnorm(n = n, rep(0,10), Sigma)
  # convert to different types
  for (j in 1:p) {
    if(family[j] == 'binary'){
      X[,j] = ifelse(X[,j]>coef[1],1,0)
    }
    
    if(family[j] == 'count'){
      X[,j] = sapply(X[,j],FUN = function(x){
        rpois(1,lambda = pnorm(coef[2]*x))})
    }
    
    if(family[j] == 'numeric'){
      X[,j] = pnorm(X[,j])
    }
  }
  
  X.scaled = scale(X,scale = FALSE)
  return(X.scaled)
}

##### Step2_generateT #####
# @X: feature matrix
# @family: hazard model of event time, "linear", "quadratic"
# @b0,beta: coefficients of the linear components in the model
# @B: Matrix of the quadratic components in the quadratic model
# @coef: coefs for cox model
# --------------------------------------------------------------
generate_T = function(n, X, b0 = 0, beta = c(runif(4,0.1,0.2),runif(4,-0.2,-0.1),runif(2,-0.3,-0.2)),
                      B = diag(c(0.1,0.1,-0.2,-0.2,-0.1,-0.1,-0.2,-0.2,0.1,-0.2)),
                      family = "linear", coef = 0.005){
  if(family == 'linear'){
    T = apply(X, 1, FUN = function(X_i){
      rexp(1, exp(t(X_i)%*%beta+b0)*coef*12)
    })
  }
 
  if(family == 'quadratic'){
    quadratic_term = apply(X,1,function(X_i){
      (b0 + t(X_i)%*%beta + t(X_i)%*%B%*%X_i)
    })
    ## justify for C_max 
    # Corrected a bug in event time generation
    piecewise_hazard_generate = function(i){
      h = function(u,i){
        exp(quadratic_term[i]/(u+1))*coef
      }
      Ti = 0
      
      for(u in 0:34)
      {
        tmpT = rexp(1,h(u,i))
        if(tmpT < 12)
        {
          Ti = Ti + tmpT
          break
        }else{
          Ti = Ti + 12
        }
      }
      return(Ti/12)
    }
    T = sapply(1:n, piecewise_hazard_generate)
  }
  return(T)
}

##### Step3_generateC #####
# @T: Event time vector
# @family: hazard model of event time, "linear", "quadratic"
# @rate: desired censor rate(not needed in this version)
# @t_max,c_min: the time period when a censor might happen 
# --------------------------------------------------------------
generate_C = function(n, T, family, rate = 0.5){
  if(family == 'linear'){
    t_max = 33
    c_min = 8
    C = runif(n, c_min, t_max+1) 
  }
  if(family == 'quadratic'){
    t_max = 33
    c_min = 8
    C = runif(n, c_min, t_max+1) 
  }
  return(C)
}

##### Step4_generateS #####
# @T: Event time vector
# @C: Censoring time vector
# @type: 'numerical': generate single numerical surrogate
#        'count': generate three counting surrogates
# --------------------------------------------------------------
generate_S = function(n, T, C, type){
  C = floor(C)
  if(type == "numerical"){
    num_S = function(T,C){
      sapply(1:C, function(t){
        if(t<T){return(rbeta(1,1,2))}
        # can use different parameters for beta distribution
        if(t>=T){return(rbeta(1,1+(t-T)/2,1))}
      })
    }
    S = sapply(1:n, function(i){num_S(T[i],C[i])})
    }
  # can use different parameters for poison distribution
  if(type == "count"){
    hos_S = function(T,C){
      sapply(1:C, function(t){
        if(t<T){return(rpois(1,0.5))}
        if(t>=T & t<=T+1){return(rpois(1,2))}
        if(t>T+1){return(rpois(1,1))}
      })
    }
    S.1 = sapply(1:n, function(i){hos_S(T[i],C[i])})
    
    dio_S = function(T,C,S){
      sapply(1:C, function(t){
        if(t<T)
        {return(rpois(1,0.05*S[t]))}
        if(t>=T & t<=T+1)
        {return(rpois(1,3*S[t]))}
        if(t>T+1)
        {return(rpois(1,2*S[t]))}
      })
    }
    S.2 = sapply(1:n, function(i){dio_S(T[i],C[i],S.1[[i]])})
    
    cui_S = function(T,C,S){
      sapply(1:C, function(t){
        if(t<T)
        {return(rpois(1,0.2*S[t]))}
        if(t>=T & t<=T+1)
        {return(rpois(1,6*S[t]))}
        if(t>T+1)
        {return(rpois(1,4*S[t]))}
      })
    }
    S.3 = sapply(1:n, function(i){cui_S(T[i],C[i],S.1[[i]])})
    
    S = sapply(1:n, function(i){rbind(S.1[[i]],S.2[[i]],S.3[[i]])})
  }
  return(S)
}

##### Step5_generateP #####
# return the true risk of the model(each patient with a risk vector at every month)
# for AUC(risk vs. status)
# @X: feature matrix
# @S: Surrogates
# @family: family argument used in generate_T
# @b0,beta: coefficients of the linear components in the model
# @B: Matrix of the quadratic components in the quadratic model
# --------------------------------------------------------------
generate_P = function(n, X, S, b0 = 0, beta = c(runif(4,0.1,0.2),runif(4,-0.2,-0.1),runif(2,-0.3,-0.2)),
                      B = diag(c(0.1,0.1,-0.2,-0.2,-0.1,-0.1,-0.2,-0.2,0.1,-0.2)),
                      family = "linear",type = 'count', coef = 0.005){
  if(type == 'numerical'){fu=sapply(S, length)}
  else{fu = sapply(S, ncol)}
  
  S.id  = rep(1:length(S), fu)
  S.mth = unlist(sapply(fu, function(x) 1:x))
  Riskmth = sapply(fu, rep, x= NA)
  for (i in 1:length(S.id)){
    if(family == 'linear'){
      linear_term = apply(X,1,function(X_i){
        (b0 + t(X_i)%*%beta)
      })
      cdf = function(i,t){
          pexp(t, exp(linear_term[i])*coef*12)
      }
      Riskmth[[S.id[i]]][S.mth[i]] = cdf(S.id[i],S.mth[i])
    }
    
    if(family == 'quadratic'){
      quadratic_term = apply(X,1,function(X_i){
        (b0 + t(X_i)%*%beta + t(X_i)%*%B%*%X_i)
      })
      cdf = function(i,t){
          a = as.numeric(35)
          for (u in 0:34) {
            a[u+1] = coef*
              max(0,min((t*12-u*12),12))*exp(
                quadratic_term[i]/(u+1) # Corrected a bug in cumulative hazard
              )
          }
          return(1-exp(-sum(a)))
      }
      Riskmth[[S.id[i]]][S.mth[i]] = cdf(S.id[i],S.mth[i])
    }
  }
  
  return(P = Riskmth)
}
##### Step5_generateMAP #####
# generate MAP if the Surrogates are counting variables
# @S: Surrogates
# --------------------------------------------------------------
generate_MAP = function(S){
  fu = sapply(S, ncol)
  Scum = sapply(S, function(x) apply(x,1,cumsum))
  S.merge = do.call("rbind", Scum)
  S.id = rep(1:length(S), fu)
  S.mth = unlist(sapply(fu, function(x) 1:x))
  MAPfit.merge = MAP_PheWAS_main(dat.icd = data.frame(ID=1:nrow(S.merge), 
                                                      feature = S.merge[,2]), 
                                 dat.nlp = data.frame(ID=1:nrow(S.merge), 
                                                      feature =  S.merge[,3]), 
                                 dat.note = data.frame(ID=1:nrow(S.merge), 
                                                       ICDday =  S.merge[,1]),
                                 nm.phe = "feature",
                                 nm.ID = "ID",
                                 nm.utl = "ICDday",
                                 p.icd = 0.001, n.icd = 10, 
                                 p.nlp = 0.001, n.nlp = 10,
                                 yes.con=FALSE, yes.nlp=FALSE)
  MAPmth = sapply(fu, rep, x= NA)
  for (i in 1:length(S.id))
  {
    MAPmth[[S.id[i]]][S.mth[i]] = MAPfit.merge$MAP[i]
  }
  return(MAPmth)
}

##### Step6_generateData #####
# @T.family: 'family' argument for generate_T
# @T.b0, T.beta, T.B: coefficients argument for generate_T
# @S.type: 'type' argument for generate_S
# @C.rate: 'rate' argument for generate_C(not needed in this version)
# @filename: data matrix for simulation, '.csv' file will be saved in working dir
# @dataname: summary of data, '.Rds' file will be saved in working dir
# -----------------------------------------------------------------------
gen_data = function(n, p=10, T.family, S.type, T.b0, T.beta, T.B, T.coef,
                    C.rate, filename, dataname){
  X = generate_X(n,p,rho = 0.1)
  T = generate_T(n, X=X, family = T.family, 
                 b0 = T.b0, beta = T.beta, B = T.B,
                 coef = T.coef)
  C = generate_C(n, T, family =T.family, rate = C.rate)
  S = generate_S(n, T, C, type = S.type)
  if(S.type == 'numerical'){MAP = S}
  if(S.type == 'count'){MAP = generate_MAP(S)}
  I = ifelse(T<=C,1,0)
  #P = generate_P(n, X, S,  b0 = T.b0, beta = T.beta, B = T.B,
  #                       family = T.family, coef = T.coef, type = S.type)
  saveRDS(list(X=X,T=T,C=C,S=S,I=I,
               #P=P,
               MAP=MAP),file = dataname)
  
  gen_single_patient = function(i, S.type){
    Censor = floor(C[i])
    T_star = floor(T[i])
    ID = rep(i,Censor+1)
    if(T_star < Censor){Y = c(rep(0,T_star),rep(1,Censor-T_star+1))}
    if(T_star >= Censor){Y = rep(0,Censor+1)}
    t = 0:Censor
    if(S.type == 'numerical'){S_t = c(0,S[[i]])}
    if(S.type == 'count'){S_t = rbind(rep(0,3),t(S[[i]]))}
    MAP_t = c(0,MAP[[i]])
    X_feature = matrix(X[i,],1,p)[rep(1,Censor+1),]
    single_patient = data.frame(ID = ID, Y = Y, T = t, X = X_feature, S = S_t, MAP = MAP_t)
    single_patient = as.matrix(single_patient)
    return(single_patient)
  }
  
  
  data = data.frame()
  for (i in 1:n) {
    data = rbind(data,gen_single_patient(i, S.type = S.type))
    
    # monitor
    if (i %% 1000 == 0){
        print(paste0(i/n*100,'%'))
    }
    
  }
  write.csv(data, file = filename)
}

##### Step7_split #####
# @data.name: the data file to be splitted to labeled and unlabeled
# -----------------------------------------------------------------------
label_split = function(data.name,label.num=1000,total=37736){
  data = read.csv(paste0(mdir,'/Simulation/','/SimDat/',data.name,'/',data.name,'.csv'))
  #print(data, ID %in% 1:label.num)
  write.csv(filter(data,ID %in% 1:label.num),paste0(mdir,'/Simulation/','/SimDat/',data.name,'/',data.name,'_labeled.csv'))
  write.csv(filter(data,ID %in% (label.num+1):total),paste0(mdir,'/Simulation/','/SimDat/',data.name,'/',data.name,'_unlabeled.csv'))
}



getauc = function(summary.file){
  Summary = readRDS(summary.file)
  time.grid = 11:35
  T = Summary$T
  C = Summary$C
  P = Summary$P
  risk = t(sapply(P, function(x) x[time.grid]))
  nt = length(time.grid)
  auc.t = rep(0, nt)
  return(auc.tspec(T,C,time.grid,risk))
}


risk.auc = function(n.repeat=100,data.file,summary.file,train.file,test.file,result.file){
  auc_TS = matrix(NA,length(time.grid),n.repeat)
  Dat = read.csv(data.file)
  Train = read.csv(train.file)
  Test = read.csv(test.file)
  Summary = readRDS(summary.file)
  P = Summary$P
  T = Summary$T
  C = Summary$C
  X = Summary$X
  Y = Summary$I
  for (i in 1:n.repeat) {
    train = Train[,i]
    test = Test[,i]
    risk = sapply(P, function(x) x[time.grid])
    auc_TS[,i] = auc.tspec(T[test],C[test],time.grid,risk)
  }
  saveRDS(auc_TS,result.file)
}

