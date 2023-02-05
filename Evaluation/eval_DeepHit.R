
# TODO: change the home dir
source("./set_mdir.R")

source(paste0(mdir,"/Evaluation/measure.R"))
source(paste0(mdir,"/Evaluation/measure_MASTA.R"))
args = commandArgs(trailingOnly = TRUE)

library(data.table)
library(pROC)

#phe.nm = "T2D"; num.labels = 120; dat.type = "cum"
# phe.nm = args[1]

#"cum": cumulative counts; "stacked": stacked counts
dat.type = 'cum'   # stands for cumulative

phe.nm = 'T2D'


############DeepHit evaluation##############

# TODO: change the path to utilize T2D data;
dat.nlong.lab = fread(paste0(mdir,"/Real/T2D/NonLong/",phe.nm,"_labeled_",dat.type,"Counts.csv"),data.table=F)

# TODO: change the path to utilize T2D data;
all.dat = list.files(paste0(mdir,"/Real/T2D"), full.names=T, pattern = "labeled.csv")
dat.labeled = fread(all.dat[1],data.table=F)
max_censoring_t = quantile(dat.nlong.lab$F.times,0.95, na.rm=TRUE)
print(max_censoring_t)
dat.sets = split(dat.labeled,dat.labeled$ID)
# print(dat.sets$`1`)

# TODO: change the path to utilize T2D data;
test_folds = fread(paste0(mdir,"/Real/T2D/test_patients.csv"),data.table = F)

#Time-specific AUCs: one AUC per month after month 24: 25-324
# TODO: Adapt the length of the matrix's second dimension(equals to max_time - baseline_time e.g. 30 - 5 in simulation data)
AUC.tspec.deephit = matrix(0,max_censoring_t,50)

#Current status AUCs
Cstat.deephit = matrix(0,50,1)

set.seed(181694)

for (f in (1:50)){

# TODO: change the path to utilize T2D data;
DeepHit_predictions = read.csv(file=paste0(mdir,"/Real/Results/",phe.nm,"-DeepHit_predictions_labels-Rep_",f,".csv"))

preds.T = as.integer(colnames(DeepHit_predictions))

test.inds = test_folds[,f]
test_data = data.frame(dat.nlong.lab[match(test.inds,dat.nlong.lab$ID),])

dat_test = dat.sets[match(test.inds,dat.nlong.lab$ID)]
# print(match(test.inds,dat.nlong.lab$ID))
# Include event times following baseline period. This is especially important, given that patients may already have ICD codes at the first baseline time point.
event = do.call(rbind,lapply(dat_test ,function(x){post.baseline.T = x$T; post.baseline.T[grep(1,x$Y)[1]]}))  # first event_onset time after 5 years
# print(event)
# TODO: adapt the time_grid
time.grid = 1:max_censoring_t

event[is.na(event)] = max(test_data$Times) + 1
# print(dim(event))

f.time = fu = test_data$F.times # Time_final

# TODO: adapt the time range

#print(if(max_t==30))
inc.prob = DeepHit_predictions[,1:max_censoring_t]


#############Time-specific AUCs
AUC.tspec.deephit[,f] = auc.tspec(event, f.time, time.grid , inc.prob)


Y.final = test_data$Y

preds.final = apply(as.matrix(seq_along(test.inds)),1,function(i){DeepHit_predictions[i,test_data$F.times[i]]})


Cstat.deephit[f,1] = auc(Y.final,preds.final,levels = c(0,1), direction = "<")   # this setence produce output, why?
# Cstat.deephit[f,1] = auc(Y.final,preds.final)

}

# print(AUC.tspec.deephit)

# TODO: adapt the time range
# TODO: change the path to utilize T2D data;
fwrite(cbind(1:max_censoring_t,AUC.tspec.deephit),file=paste0(mdir,
                                             "/Evaluation/Results/DeepHit/T2D/",phe.nm,"_tspec_AUC-ALL.deephit_labels.csv"))

# TODO: adapt the time range
# TODO: change the path to utilize T2D data;
fwrite(cbind(1:max_censoring_t,apply(AUC.tspec.deephit,1,mean),apply(AUC.tspec.deephit,1,sd)),file=paste0(mdir,
                                                "/Evaluation/Results/DeepHit/T2D/",phe.nm,"_tspec_AUC.deephit_labels.csv"))

# TODO: change the path to utilize T2D data;
fwrite(Cstat.deephit,file=paste0(mdir,
                                "/Evaluation/Results/DeepHit/T2D/",phe.nm,"_Cstat-ALL.deephit_labels.csv"))
# TODO: change the path to utilize T2D data;
write.csv(c(mean(Cstat.deephit),sd(Cstat.deephit)),file=paste0(mdir,
                                                            "/Evaluation/Results/DeepHit/T2D/",phe.nm,"_Cstat.deephit_labels.csv"))



