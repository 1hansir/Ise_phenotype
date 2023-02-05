
#Risk Prediction: MAP for silver standard labels at each one year time window
#Datasets: HF, T1DM and T2DM
#Note: T1DM is a juvenile onset disease, so may be less compelling from a clinical point of view for risk prediction. 
#As in, we know that genetics play a very large, if not the largest, role in determining patient risk of disease (and also [<perhaps>] age)


rm(list=ls())
source('~/Dropbox/SemiSupervised_RiskPrediction/Utilities/MAPflex.R')
source("~/Dropbox/SemiSupervised_RiskPrediction/Utilities/library_v2.R")
source("~/Dropbox/SemiSupervised_RiskPrediction/Utilities/MAP_0.1.0/MAP/R/FUN_NLP_PheWAS.R")
source("~/Dropbox/SemiSupervised_RiskPrediction/Utilities/addMAP_1mo_timeWindow.R")
library(data.table)
library(MAP)


sim.dat = list.files("~/Dropbox/SemiSupervised_RiskPrediction/Simulation/data_samples",full.names=T, pattern = ".csv")

# Read in saved simulation dataset 2, the one with counts as surrogate variables S.1, S.2 and S.3. Replace indices with correct values
dt2.labeled = fread(sim.dat[3],data.table=F)
dt2.unlabeled = fread(sim.dat[4],data.table=F)


############MAP testing and troubleshooting, for labeled dataset 2.2: you cam skip this part
dt_2.2_labeled = fread("~/Dropbox/SemiSupervised_RiskPrediction/Simulation/data_samples/data.2.2.csv",data.table=F)
#EHR.dat = dt_2.2_labeled; mo.window = 1; ICD_NLP.nm = c("S.2","S.3")
TEST_2.2 = add_MAP_SimDat(dt_2.2_labeled,1,c("S.2","S.3"))
############
dt_1.2_labeled = fread("~/Dropbox/SemiSupervised_RiskPrediction/Simulation/data_samples/data.1.2.csv",data.table=F)
#EHR.dat = dt_2.2_labeled; mo.window = 1; ICD_NLP.nm = c("S.2","S.3")
TEST_1.2 = add_MAP_SimDat(dt_1.2_labeled,1,c("S.2","S.3"))
#summary(TEST_2.2$MAP); summary(dt_2.2_labeled$S.2); summary(dt_2.2_labeled$S.3); summary(log(dt_2.2_labeled$S.2+1)); summary(log(dt_2.2_labeled$S.3+1))
fwrite(TEST_2.2,file = "~/Dropbox/SemiSupervised_RiskPrediction/Simulation/dt1.1-2.2/data.2.2.labeled_MAP.csv")
TEST_2.2 = fread("~/Dropbox/SemiSupervised_RiskPrediction/Simulation/data_samples/data.2.2_MAP.csv", data.table=F)

dat.sets = split(TEST_1.2,TEST_1.2$ID)

MAP.dat = do.call(rbind,lapply(dat.sets, function(x){x[nrow(x),c("Y","MAP")] }))

auc(MAP.dat[,1],MAP.dat[,2])



lab.sets = split(dt_2.2_labeled,dt_2.2_labeled$ID)
lab.dat = do.call(rbind,lapply(lab.sets,function(x){x[nrow(x),]}))




MAP.dat = list()

patient_num = as.matrix(0:(nrow(lab.dat)-1))
colnames(patient_num) = "patient_num"
MAP.dat$ID = patient_num

ICD_NLP.vec = lab.dat[,c("S.2","S.3")]

#ICD_NLP = as.matrix(ICD_NLP.vec); ICD_NLP[is.na(ICD_NLP)] = 0
ICD_NLP = ICD_NLP.vec; ICD_NLP[is.na(ICD_NLP)] = 0
#colnames(ICD_NLP) = paste(phe.nm,"_ICD_NLP",sep="")
colnames(ICD_NLP) = ICD_NLP.nm
MAP.dat$mat = ICD_NLP

utl = as.matrix(lab.dat$S.1)
colnames(utl) = "utl"
MAP.dat$note = utl

#Compute MAP probabilities, based on ICD_NLP codes only
# if (length(unique(ICD_NLP[,1])) <= 2 | length(unique(ICD_NLP[,2])) <= 2 | 
#     mean(ICD_NLP[,1] > 0) < 0.001 | mean(ICD_NLP[,2] > 0) < 0.001) {
#   MAP.probs[,i]= rep(0,nrow(MAP.probs))
# } else {
MAP.prob = tryCatch(MAP_flex_core(data = MAP.dat,  yes.con = FALSE, full.output = TRUE))
MAP.probs = MAP.prob$scores.all[,2]
# }

auc(lab.dat$Y,MAP.probs)


############

#Run lines below to append MAP probabilities to data, to labeled and unlabeled samples separately
fwrite(add_MAP_SimDat(dt2.labeled,1,c("S.2","S.3")),file = "~/Dropbox/SemiSupervised_RiskPrediction/Simulation/dt2.labeled_MAP.csv")
fwrite(add_MAP_SimDat(dt2.unlabeled,1,c("S.2","S.3")),file = "~/Dropbox/SemiSupervised_RiskPrediction/Simulation/dt2.unlabeled_MAP.csv")


