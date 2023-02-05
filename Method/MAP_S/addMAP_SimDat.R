
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


sim.dat = list.files("~/Dropbox/SemiSupervised_RiskPrediction/Simulation",full.names=T, pattern = ".csv")

# Read in saved simulation dataset 2, the one with counts as surrogate variables S.1, S.2 and S.3. Replace indices with correct values
dt2.labeled = fread(sim.dat[3],data.table=F)
dt2.unlabeled = fread(sim.dat[4],data.table=F)


############MAP testing and troubleshooting, for labeled dataset 2.2: you cam skip this part
dt_2.2_labeled = fread("~/Dropbox/SemiSupervised_RiskPrediction/Simulation/dt1.1-2.2/data.2.2.labeled.csv",data.table=F)
#EHR.dat = dt_2.2_labeled; mo.window = 1; ICD_NLP.nm = c("S.2","S.3")
TEST_2.2 = add_MAP_SimDat(dt_2.2_labeled,1,c("S.2","S.3"))
#summary(TEST_2.2$MAP); summary(dt_2.2_labeled$S.2); summary(dt_2.2_labeled$S.3); summary(log(dt_2.2_labeled$S.2+1)); summary(log(dt_2.2_labeled$S.3+1))
fwrite(TEST_2.2,file = "~/Dropbox/SemiSupervised_RiskPrediction/Simulation/dt1.1-2.2/data.2.2.labeled_MAP.csv")
############

#Run lines below to append MAP probabilities to data, to labeled and unlabeled samples separately
fwrite(add_MAP_SimDat(dt2.labeled,1,c("S.2","S.3")),file = "~/Dropbox/SemiSupervised_RiskPrediction/Simulation/dt2.labeled_MAP.csv")
fwrite(add_MAP_SimDat(dt2.unlabeled,1,c("S.2","S.3")),file = "~/Dropbox/SemiSupervised_RiskPrediction/Simulation/dt2.unlabeled_MAP.csv")


