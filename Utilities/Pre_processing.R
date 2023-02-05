
library(MAP)
library(ff)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

source("./set_mdir.R")
total_data_num = as.numeric(args[1])
# 37736
labeled_data_num = as.numeric(args[2])
# 1000
num.train = args[3]

num.test = args[4]

data.file = args[5] # == phe.nm

min_t = as.numeric(args[6])

#----------------------------------------
# data generation #
#----------------------------------------
source(paste0(mdir,'/Simulation/simulation_data_functions_v4.R'))
source(paste0(mdir,"/Utilities/get_Folds.R"))

Time_offset = function(dat.sets){

  dat_offset = do.call(rbind,lapply(dat.sets,function(patient.dat){
    First_date = patient.dat$T[1]
    patient.dat$T = patient.dat$T - First_date
    # print(patient.dat)
    # print('-------------------------')
    colnames(patient.dat) = c('V1','X','ID','T','Y','X.1','X.2','X.3','X.4','X.5','X.6','X.7','X.8','X.9','X.10','X.11','X.12','S','MAP')
    patient.dat = data.frame(patient.dat)

  }))

  return(dat_offset)
 }


label_split = function(data.name,label.num=1000,total=37736){
  data = read.csv(paste0(mdir,'/Real/T2D/',data.name,'.csv'))
  # print(data, 'ID' %in% 1:label.num)

  write.csv(data[data$ID %in% 1:label.num,],paste0(mdir,'/Real/T2D/',data.name,'_labeled_mid.csv'))
  write.csv(data[!(data$ID %in% 1:label.num),],paste0(mdir,'/Real/T2D/',data.name,'_unlabeled_mid.csv'))

  all.dat = list.files(paste0(mdir,"/Real/T2D/"),full.names=T, pattern = "labeled_mid.csv")

  dat.labeled = fread(all.dat[1],data.table=F);
  dat.unlabeled = fread(all.dat[2],data.table=F);

  labeled.sets = split(dat.labeled,dat.labeled$ID); unlabeled.sets = split(dat.unlabeled,dat.unlabeled$ID)

  fwrite(Time_offset(labeled.sets),paste0(mdir,'/Real/T2D/',data.name,'_labeled.csv'))
  fwrite(Time_offset(unlabeled.sets),paste0(mdir,'/Real/T2D/',data.name,'_unlabeled.csv'))
}

label_split(data.file,label.num = labeled_data_num, total = total_data_num)
getFolds_RP_v2(phe.nm = data.file, Observation_years = min_t, ntrain = num.train, ntest = as.numeric(num.test) )
