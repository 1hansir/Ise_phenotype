
rm(list=ls())
# print('1')
#args = commandArgs(trailingOnly = TRUE)

library(data.table)
library(dplyr)


# SimDat.1
# TODO: change the home dir
source("./set_mdir.R")

###################Load labeled data for phenotype of interest

# source(paste0(mdir,"/Utilities/get_Folds.R"))
phe.nm = 'T2D'

cumDat = function(dat.sets){
  
  
  dat.cum = do.call(rbind,lapply(dat.sets,function(long.dat){
    #Retain only data up to month 324, 25 years after 2-year baseline observation period


    # TODO: adapt the maximum time(30 for simulation data)
    #x = long.dat[long.dat$T <= 20,]  # Maybe for simulation data, we dont need this line
    x = long.dat
    
    #Event time of interest: (minimum of) positive diagnosis or censoring time
    e.time = ifelse(length(grep(1,x$Y)) > 0, x$T[grep(1,x$Y)[1]],x$T[nrow(x)])
    surv.info = t(as.matrix(c(e.time,x[nrow(x),c("Y","ID","T")])))  # (event_onset_time, Y_final, ID, T_final <-- Last row
    colnames(surv.info) = c("Times","Y","ID","F.times")
    

    # Define baseline data vector

    end.base = which(x$T >= 0)[1]         # ==1 here
    #print(end.base)
    if (is.na(end.base)){
      print(1)
      cov.info = rep(0,ncol(x)-5)   # 这个5并不是min_t!!!
    } else {
      X = x[1:(end.base-1),]
      #if (length(X) == 18){
      #  print(X)
      # }
      if (is.null(dim(X))){        # if the X has only one row 
        # print(0)
        cov.info = apply(X[-c(1:3,ncol(X)-1,ncol(X))],2,sum) 
      } else {
        cov.info = apply(X[1:nrow(X),-c(1:3,ncol(X)-1,ncol(X))],2,sum)
      }
    }

    # print(cov.info)
    DAT = unlist(c(surv.info,cov.info)); 

    # (event_onset time, final_Y, ID, censoring time, Y_sum(of the first 5 years), T_sum(of the first 5 years, which is always 10)), sum of first 5 years' X_i
  }))
  
  return(dat.cum)
}



all.dat = list.files(paste0(mdir,"/Real/T2D/"),full.names=T, pattern = "labeled.csv")

#Make sure the selected files are the ones containing EXCLUSIVELY labeled and unlabeled subjects
dat.labeled = fread(all.dat[1],data.table=F); dat.unlabeled = fread(all.dat[2],data.table=F)


labeled.sets = split(dat.labeled,dat.labeled$ID); unlabeled.sets = split(dat.unlabeled,dat.unlabeled$ID)

#####################
#Compress all longitudinal EHR data into vectors of cumulative code counts up to month 24 (this gives us one vector of baseline data for each patient),
# follow-up time F.times as last observed T, time of interest Times as minimum of observed event time and censoring time,
# and censoring status as last observed value of Y, last current status value


dat.cum.lab = cumDat(labeled.sets); dat.cum.unlab = cumDat(unlabeled.sets); dat.cum.all = cumDat(c(labeled.sets,unlabeled.sets))
colnames(dat.cum.lab)[1:4] = colnames(dat.cum.unlab)[1:4] = colnames(dat.cum.all)[1:4] = c("Times","Y","ID","F.times")


fwrite(as.data.frame(dat.cum.lab), file = paste0(mdir,"/Real/T2D/NonLong/",phe.nm,"_labeled_cumCounts.csv"))
fwrite(as.data.frame(dat.cum.unlab), file = paste0(mdir,"/Real/T2D/NonLong/",phe.nm,"_unlabeled_cumCounts.csv"))
fwrite(as.data.frame(dat.cum.all), file = paste0(mdir,"/Real/T2D/NonLong/",phe.nm,"_all_cumCounts.csv"))

print(paste0(phe.nm,' convertion complete'))


