
library(MAP)
library(ff)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

source("./set_mdir.R")


# 1000
num.train = args[1]

#----------------------------------------
# data generation #
#----------------------------------------

source(paste0(mdir,"/Utilities/get_Folds.R"))

getFolds_RP_v2(phe.nm = 'T2D', Observation_years = -1, ntrain = as.numeric(num.train))
