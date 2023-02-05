rm(list=ls())
source('~/Dropbox/SemiSupervised_RiskPrediction/Utilities/MAPflex.R')
source("~/Dropbox/SemiSupervised_RiskPrediction/Utilities/library_v2.R")
source("~/Dropbox/SemiSupervised_RiskPrediction/Utilities/MAP_0.1.0/MAP/R/FUN_NLP_PheWAS.R")
source("~/Dropbox/SemiSupervised_RiskPrediction/Utilities/addMAP_1mo_timeWindow.R")
library(data.table)
library(MAP)

phe.nm = "T2D"

all.dat = list.files(paste0("~/Dropbox/RiskPrediction_Deep/Data/",phe.nm),full.names=T,pattern = "labeled.csv")
dat.labeled = fread(all.dat[1],data.table=F)
labeled.sets = split(dat.labeled,dat.labeled$ID)

all.dat_RP = list.files(paste0("~/Dropbox/RiskPrediction_Deep/Data/",phe.nm),full.names=T,pattern = "labeled_RP.csv")
dat.labeled_RP = fread(all.dat_RP[1],data.table=F)

# Generate MAP
fu = sapply(labeled.sets, function(x){nrow(x)}); summary(fu)

Scum = sapply(labeled.sets, function(x) apply(x[,c("PheCode:250.2","utl")],2,cumsum))
S.merge = do.call("rbind", Scum)
S.id = rep(1:length(labeled.sets), fu)
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


MAPmth = sapply(fu, rep, x= NA)
for (i in 1:length(S.id))
{
  MAPmth[[S.id[i]]][S.mth[i]] = MAPfit.merge$MAP[i]
}

dat.labeled$MAP = dat.labeled_RP$MAP = MAPfit.merge$MAP
dat.labeled$S = NULL

fwrite(dat.labeled,file=all.dat[1]); fwrite(dat.labeled_RP,file=all.dat_RP[1])


################Unlabeled: will create a unified function to take any single dataset as input, either labeled or unlabeled


phe.nm = "T2D"

dat.unlabeled = fread(all.dat[2],data.table=F)
unlabeled.sets = split(dat.unlabeled,dat.unlabeled$ID)

dat.unlabeled_RP = fread(all.dat_RP[2],data.table=F)

# Generate MAP
fu = sapply(unlabeled.sets, function(x){nrow(x)}); summary(fu)

Scum = sapply(unlabeled.sets, function(x) apply(x[,c("PheCode:250.2","utl")],2,cumsum))
S.merge = do.call("rbind", Scum)
S.id = rep(1:length(unlabeled.sets), fu)
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


MAPmth = sapply(fu, rep, x= NA)
for (i in 1:length(S.id))
{
  MAPmth[[S.id[i]]][S.mth[i]] = MAPfit.merge$MAP[i]
}

dat.unlabeled$MAP = dat.unlabeled_RP$MAP = MAPfit.merge$MAP
dat.unlabeled$S = NULL

fwrite(dat.unlabeled,file=all.dat[2]); fwrite(dat.unlabeled_RP,file=all.dat_RP[2])



event.time = sapply(labeled.sets, function(x){x$T[grep(1,x$Y)[1]]})

# Evaluate AUC
cs.mth = rep(0,length(S.id))
for (i in 1:length(S.id))
{
  cs.mth[i] = as.numeric(event.time[S.id[i]] <= S.mth[i])
}
auc(cs.mth[!is.na(cs.mth)], drop(MAPfit.merge$MAP)[!is.na(cs.mth)]) # pooled AUC
auc(cs.mth[!is.na(cs.mth)], dat.labeled$MAP[!is.na(cs.mth)]) # pooled AUC

cs = as.numeric(simdat$T <= simdat$C)
auc(cs, sapply(MAPmth, function(x) rev(x)[1])) # current status AUC
auc.mth = rep(NA, max(fu))
for (mth in 1:max(fu))
{
  sel = S.mth==mth
  if(length(unique(cs.mth[sel]))==2)
  {
    auc.mth[mth] = auc(cs.mth[sel], 
                       drop(MAPfit.merge$MAP)[sel])
  }
}





dt_2.2_labeled$S = dt_2.2_labeled[,]
# Generate MAP
fu = sapply(dt_2.2_labeled[,c("S.1","S.2","S.3")], ncol)
summary(fu)

Scum = sapply(dt_2.2_labeled[,c("S.1","S.2","S.3")], function(x) apply(x,1,cumsum))
S.merge = do.call("rbind", Scum)
S.id = rep(1:length(dt_2.2_labeled$S), fu)
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
