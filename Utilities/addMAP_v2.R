rm(list=ls())
source('~/Dropbox/SemiSupervised_RiskPrediction/Utilities/MAPflex.R')
source("~/Dropbox/SemiSupervised_RiskPrediction/Utilities/library_v2.R")
source("~/Dropbox/SemiSupervised_RiskPrediction/Utilities/MAP_0.1.0/MAP/R/FUN_NLP_PheWAS.R")
source("~/Dropbox/SemiSupervised_RiskPrediction/Utilities/addMAP_1mo_timeWindow.R")
library(data.table)
library(MAP)


getMAP_month = function(EHR.dat,ICD.nm){
  
  code.nm = grep("PheCode",names(EHR.dat))
  EHR.dat$utl = apply(EHR.dat[,c(code.nm)],1,sum)
    
  EHR.sets = split(EHR.dat,EHR.dat$ID)
  
  # Generate MAP
  fu = sapply(EHR.sets, function(x){nrow(x)}); summary(fu)
  
  Scum = sapply(EHR.sets, function(x) apply(x[,c(ICD.nm,"utl")],2,cumsum))
  S.merge = do.call("rbind", Scum)
  S.id = rep(1:length(EHR.sets), fu)
  S.mth = unlist(sapply(fu, function(x) 1:x))
  
  MAPfit.merge = MAP_PheWAS_main(dat.icd = data.frame(ID=1:nrow(S.merge), 
                                                      feature = S.merge[,1]), 
                                 dat.nlp = data.frame(ID=1:nrow(S.merge), 
                                                      feature =  S.merge[,1]), 
                                 dat.note = data.frame(ID=1:nrow(S.merge), 
                                                       ICDday =  S.merge[,2]),
                                 nm.phe = "feature",
                                 nm.ID = "ID",
                                 nm.utl = "ICDday",
                                 p.icd = 0.001, n.icd = 10, 
                                 p.nlp = 0.001, n.nlp = 10,
                                 yes.con=FALSE, yes.nlp=FALSE)
 
  ###############
  event.time = sapply(EHR.sets, function(x){x$T[grep(1,x$Y)[1]]})

  # Evaluate pooled monthly AUC
  cs.mth = rep(0,length(S.id))
  
  for (i in 1:length(S.id))
  {
    cs.mth[i] = as.numeric(event.time[S.id[i]] <= S.mth[i])
  }
  MAP.mth.res = MAPfit.merge$MAP[!is.na(cs.mth)]; CS.MTH = cs.mth[!is.na(cs.mth)]
  #AUC.mth.pool = auc(cs.mth[!is.na(cs.mth)], drop(MAPfit.merge$MAP)[!is.na(cs.mth)]) # pooled AUC
  if(length(unique(EHR.dat$Y)) == 1){ AUC.mth.pool = NULL} else { AUC.mth.pool = auc(CS.MTH, MAP.mth.res) } # pooled AUC

  
 ################  
 
  # #Evaluate monthly AUC 
  MAPmth = sapply(fu, rep, x= NA)
  for (i in 1:length(S.id))
  {
    MAPmth[[S.id[i]]][S.mth[i]] = MAPfit.merge$MAP[i]
  }

  # auc.mth = rep(NA, max(fu))
  # for (mth in 1:max(fu))
  # {
  #   sel = S.mth==mth
  #   if(length(unique(cs.mth[sel]))==2)
  #   {
  #     auc.mth[mth] = auc(cs.mth[sel], 
  #                        drop(MAPfit.merge$MAP)[sel])
  #   }
  # }
  
######################
  
  #Evaluate final status AUC
  cs = sapply(EHR.sets, function(x) x$Y[nrow(x)])
  if(length(unique(EHR.dat$Y)) == 1){ AUC.final.status = NULL} else {   AUC.final.status = auc(cs, sapply(MAPmth, function(x) rev(x)[1])) }# current status AUC
  
  #, MAP.AUC.mth = auc.mth
  return(list(MAP = MAPfit.merge$MAP, 
              MAP.AUC.pool = AUC.mth.pool, MAP.AUC.final = AUC.final.status))
}



# phe.nm = "T2D"
# all.dat = list.files(paste0("~/Dropbox/RiskPrediction_Deep/Data/",phe.nm),full.names=T,pattern = "labeled.csv")
# EHR.dat = fread(all.dat[1],data.table=F)
# 
# 
# labeled_MAP = getMAP_month(EHR.dat)
# 
# labeled_MAP$MAP.AUC.pool
# labeled_MAP$MAP.AUC.final
# 







