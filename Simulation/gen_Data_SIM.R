# TODO: change the home dir
mdir = "/home/yihan/anaconda3/envs/risk-prediction"

library(MAP)
library(ff)
library(data.table)

args = commandArgs(trailingOnly = TRUE)


total_data_num = as.numeric(args[1])
# 37736
labeled_data_num = as.numeric(args[2])
# 1000
num.train = args[3]

#"cum": cumulative counts; "stacked": stacked counts
num.test = args[4]

#----------------------------------------
# data generation #
#----------------------------------------
source(paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/simulation_data_functions_v4.R'))
source(paste0(mdir,"/SemiSupervised_RiskPrediction/Utilities/get_Folds_SIM.R"))

# interaction matrix for data.2.1 and data.2.2
A<-matrix(rep(c(1,-2),50),10)
pmean <- function(x,y) (x+y)/2
A[] <- pmean(A, matrix(A, nrow(A), byrow=TRUE))
s<-matrix(rep(1,100),10)-2*A+diag(rep(c(0.5,-4.5),5))


# data.1.1
if (!(file.exists(paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/SimDat/SimDat.1.1/','SimDat.1.1.csv')))){

    gen_data(total_data_num, p=10, T.family = 'linear', S.type = 'numerical', 
         T.b0 = -3, T.beta = 0.5*c(1,-1,-2,-2, 1,-1,-2,-2, 1,-2), T.coef = 0.05,
         C.rate = 0.5, filename = paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/SimDat/SimDat.1.1/','SimDat.1.1.csv'),dataname =paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/SimDat/SimDat.1.1/','SimDat.1.1.Rds'))
}
print('1.1 ready')

# data.1.2
if (!(file.exists(paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/SimDat/SimDat.1.2/','SimDat.1.2.csv')))){

    gen_data(total_data_num, p=10, T.family = 'linear', S.type = 'count', 
         T.b0 = -3, T.beta = 0.5*c(1,-1,-2,-2, 1,-1,-2,-2, 1,-2), T.coef = 0.05,
         C.rate = 0.5, filename = paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/SimDat/SimDat.1.2/','SimDat.1.2.csv'),dataname =paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/SimDat/SimDat.1.2/','SimDat.1.2.Rds'))
}
print('1.2 ready')

# data.2.1
if (!(file.exists(paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/SimDat/SimDat.2.1/','SimDat.2.1.csv')))){

    gen_data(total_data_num, p=10, T.family = 'quadratic', S.type = 'numerical', 
         T.b0 = -30, T.beta = 0.5*c(1,-1,-2,-2, 1,-1,-2,-2, 1,-2), T.B = -1*s, T.coef = 0.05,
         C.rate = 0.5, filename = paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/SimDat/SimDat.2.1/','SimDat.2.1.csv'),dataname =paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/SimDat/SimDat.2.1/','SimDat.2.1.Rds'))
}
print('2.1 ready')

# data.2.2
if (!(file.exists(paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/SimDat/SimDat.2.2/','SimDat.2.2.csv')))){

    gen_data(total_data_num, p=10, T.family = 'quadratic', S.type = 'count', 
         T.b0 = -30, T.beta = 0.5*c(1,-1,-2,-2, 1,-1,-2,-2, 1,-2), T.B = -1*s, T.coef = 0.05,
         C.rate = 0.5, filename = paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/SimDat/SimDat.2.2/','SimDat.2.2.csv'),dataname =paste0(mdir,'/SemiSupervised_RiskPrediction/Simulation/SimDat/SimDat.2.2/','SimDat.2.2.Rds'))
}
print('2.2 ready')
# split labeled and unlabeled


for (data.file in c('SimDat.1.1','SimDat.1.2','SimDat.2.1','SimDat.2.2')) {
  label_split(data.file,label.num = labeled_data_num, total = total_data_num)
  getFolds_RP_v2(phe.nm = data.file, Observation_years = 5, ntrain = num.train, ntest = as.numeric(num.test) )
}

