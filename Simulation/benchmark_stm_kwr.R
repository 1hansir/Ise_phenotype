source(paste0(mdir,"/Evaluation/measure.R"))


## benchmark method ##
# ----------------------------#
# kernel weighted ts regression #
kwts.regression = function(X,C,S,Y,train,test,unlabel,time.grid){
  N = length(Y)
  dt.1 = data.frame(X,S,Y)
  h = sd(C)/sum(Y)^0.25
  K = function(t){1/(h*sqrt(2*pi))*exp(-(C-t)^2/(2*h^2))}
  
  # regression for imputed outcome
  kwip.regression = function(t){
    model = glm(Y~., data = dt.1, family = 'binomial',weights = K(t),subset =train)
    imputed.values = predict(model,data.frame(X,S),type = 'response')
    return(c(Y[1:1000],imputed.values[1001:N]))
  }
  
  # ts.weighted logistic regression
  prediction = matrix(NA,length(test),length(time.grid))
  for (t in time.grid) {
    dt.2 = data.frame(X,Y=kwip.regression(t))
    model = glm(Y~., data = dt.2,  weights = K(t),subset = c(train,unlabel))
    prediction[,t-time.grid[1]+1] = 1/(1+exp(-predict(model,data.frame(X[test,]))))
  }
  
  return(prediction)
}


simulation_kw = function(n.repeat=100,summary.file,train.file,test.file,unlabeled.file,result.file,time.grid){
  auc_TS = matrix(NA,length(time.grid),n.repeat)
  Train = read.csv(train.file)
  Test = read.csv(test.file)
  Unlabeled = read.csv(unlabeled.file)
  Summary = readRDS(summary.file)
  S = sapply(Summary$MAP, function(x) rev(x)[1])
  T = Summary$T
  C = Summary$C
  X = Summary$X
  Y = Summary$I
  for (i in 1:n.repeat) {
    train = Train[,i]
    test = Test[,i]
    unlabel = Unlabeled[,1]
    prediction = kwts.regression(X,C,S,Y,train,test,unlabel,time.grid)
    auc_TS[,i] = auc.tspec(T[test],C[test],time.grid,prediction)
  }
  saveRDS(auc_TS,result.file)
}  





## semiparametric transformation #
# -------------------------------------#
source(paste0(mdir,'/SemiSupervised_RiskPrediction/Method/STM/pred_STM.R'))
source(paste0(mdir,'/SemiSupervised_RiskPrediction/Method/STM/implicit_profiling.R'))
source(paste0(mdir,'/SemiSupervised_RiskPrediction/Method/STM/utility.R'))

stm.regression = function(X,C,S,Y,train,test,time.grid){
  h = sd(C[train])/(sum(Y[train]))^0.25
  # PCA
  Z = prcomp(X,center = TRUE,scale. = TRUE)$x[,1:8]
  # STM
  ip.fit = implicit.profile(Y[train],Z[train,],KC = dnorm(as.matrix(dist(C[train]/h,diag=T,upper=T)))/h, 
                            C[train], h.loop = F,
                            glink = expexp, dglink = dexpexp)
  pred.mat = pred.STM(ip.fit, Z[test,], tgrid = time.grid,glink = expexp)
  return(pred.mat)
}

simulation_stm = function(n.repeat=100,summary.file,train.file,test.file,result.file,time.grid){

  auc_TS = matrix(NA,length(time.grid),n.repeat)
  Train = read.csv(train.file)
  Test = read.csv(test.file)
  Summary = readRDS(summary.file)
  S = sapply(Summary$MAP, function(x) rev(x)[1])
  T = Summary$T
  C = Summary$C
  X = Summary$X
  Y = Summary$I
  for (i in 1:n.repeat) {
    train = Train[,i]
    test = Test[,i]
    prediction = stm.regression(X,C,S,Y,train,test,time.grid)
    auc_TS[,i] = auc.tspec(T[test],C[test],time.grid,prediction)
  }
  saveRDS(auc_TS,result.file)
}



