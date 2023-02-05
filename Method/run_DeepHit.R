
###################Commands to enter in O2 prior to launching R.
#####Setup for necessary modules and virtual environment with installed packages

# srun --pty -p interactive -t 0-10:00 --mem 75G bash
# module purge
# module load gcc/6.2.0 python/3.7.4 R/3.6.1
# unset PYTHONPATH
# source ~/Survival_env/bin/activate
# R

args = commandArgs(trailingOnly = TRUE)

library(data.table)
# library(distr6)
library(dplyr)
library(reticulate)
library(survival)
library(survivalmodels)

# use_virtualenv("~/Survival_env")

# TODO: change the home dir
source("./set_mdir.R")

#phe.nm = "SimDat.1"; num.labels = 200; dat.type = "cum"
# phe.nm = args[1]
#num.labels will be equal to 200, by default, as this is the only version we are considering...(n.train in get_Folds.R)

#"cum": cumulative counts; "stacked": stacked counts

dat.type = 'cum'

phe.nm = 'T2D'

#phe.nm = "SimDat.1.1"; Observation_years = 5; ntrain = 300; ntest = 170


###################Load labeled data for phenotype of interest

print('--------------------------------')
print(paste0('Training Process of: ', phe.nm))
print('--------------------------------')

# TODO: change the path to utilize T2D data;
dat.nlong.lab = fread(paste0(mdir,"/Real/T2D/NonLong/",phe.nm,"_labeled_",dat.type,"Counts.csv"),data.table=F)

feat.keep = colnames(dat.nlong.lab)[grep("X",colnames(dat.nlong.lab))]

#Load IDs of patients in training and testing sets

# TODO: change the path to utilize T2D data;
train_folds = fread(paste0(mdir,"/Real/T2D/train_patients.csv"), data.table=F)
test_folds = fread(paste0(mdir,"/Real/T2D/test_patients.csv"),data.table = F)


############ DeepHit implementation


for (f in 1:50){       # f in 1:50

print(paste0("replicate:",f))
train.inds = match(train_folds[,f],dat.nlong.lab[,"ID"])
test.inds = match(test_folds[,f],dat.nlong.lab[,"ID"])

# train_data = data.frame(dat.nlong.lab[train.inds,-c(3:4)])
# test_data = data.frame(dat.nlong.lab[test.inds,-c(3:4)])
train_data = data.frame(dat.nlong.lab[train.inds, c("Times","Y",feat.keep)])
test_data = data.frame(dat.nlong.lab[test.inds,c("Times","Y",feat.keep)])

#Defining cuts as max(train_data$Times)+1 allows for obtaining a prediction at each month, from 0 to the maximum event/censoring time
#in the training set
set.seed(181694)
print("start training")
DeepHit_Model = deephit(Surv(Times, Y,type="right") ~ . , data = train_data, cuts = max(train_data$Times)+1, epochs = 25,
                        best_weights = TRUE, verbose = TRUE, frac = 0.2)

#Obtain cumulative incidence probabilities, instead of survival probabilities
print("start evaluating")
DeepHit_predictions = 1 - predict(DeepHit_Model,newdata = test_data, distr6= F, type="survival")

# TODO: change the path to utilize T2D data;
fwrite(data.frame(DeepHit_predictions),file=paste0(mdir,"/Real/Results/",phe.nm,"-DeepHit_predictions_labels-Rep_",f,".csv"))

}


