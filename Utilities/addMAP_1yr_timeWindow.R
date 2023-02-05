
#Risk Prediction: MAP for silver standard labels at each one year time window
#Datasets: HF, T1DM and T2DM
#Note: T1DM is a juvenile onset disease, so may be less compelling from a clinical point of view for risk prediction. 
#As in, we know that genetics play a very large, if not the largest, role in determining patient risk of disease (and also [<perhaps>] age)


rm(list=ls())
source('~/Dropbox/Phenotyping_Deep/MAPflex.R')
source("~/Dropbox/Phenotyping_Deep/Isabella/R.funs/Utilities/library_v2.R")
source("~/Dropbox/Phenotyping_Deep/Isabella/MAP_0.1.0/MAP/R/FUN_NLP_PheWAS.R")
source("~/Dropbox/RiskPrediction_Deep/Code/addMAP_1mo_timeWindow.R")
library(data.table)
library(ff)
library(MAP)
library(SemiEstimate)


T2D.dat = list.files("~/Dropbox/RiskPrediction_Deep/Data/T2D",full.names=T)

#Includes filter positive AND filter negative patients
T2D.fpos = fread(T2D.dat[6],data.table=F)
npos = length(unique(T2D.fpos$ID))
#Filter negative patients
T2D.fneg.0 = fread(T2D.dat[5],data.table=F); 

code.8 = T2D.fneg.0[T2D.fneg.0$`PheCode:250.2` == 8,]
code.8$ID




#Check that patient IDs from filter positive and negative sets do not overlap
common.IDs = intersect(T2D.fpos$ID,T2D.fneg.0$ID)
#common.IDs:104062 101140 101859 100558 102331 106865 100108 102491 105079 104756

#Compare patients w/ same IDs in filter positive and negative sets
#Patient 104062
T2D.fpos[T2D.fpos$ID == common.IDs[1],"PheCode:250.2"]
#[1] 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 1 1 0 0 0 0 3 0 1 3 1 0 0
T2D.fneg.0[T2D.fneg.0$ID == common.IDs[1],"PheCode:250.2"]
#[1] 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 1 0 1 0 0 0 0 3 0 1 2 1 0

#Patient 101140
T2D.fpos[T2D.fpos$ID == common.IDs[2],"PheCode:250.2"]
#[1] 0 0 0 0 0 0 0 0 0 3 1 0 0 0 0 0 2 0 0 0 0 0 0 1
T2D.fneg.0[T2D.fneg.0$ID == common.IDs[2],"PheCode:250.2"]
#[1] 0 3 1 2 0 1

#Patient 101859
T2D.fpos[T2D.fpos$ID == common.IDs[3],"PheCode:250.2"]
#[1] 1 2 0 0 1 0 0 0 0 0 0 1 0 1 0 0 0 0 1 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0
T2D.fneg.0[T2D.fneg.0$ID == common.IDs[3],"PheCode:250.2"]
#[1] 1 2 1 0 0 0 1 0 1 0 1 2 0 0 0 0 0 0

T2D.fneg.0$`PheCode:250.2`


#Demographic + GRS data frame with pseudo patient IDs
T2D.demo.dat.0 = fread("~/Dropbox/SemiSupervised_RiskPrediction/Data Analysis/MGB_demo_GRS.csv", data.table=F)
#Pseudo patient ID to true patient ID mapping
ID.map = fread("~/Dropbox/RiskPrediction_Deep/Data/T2D/Utilities/patient_fake_id_map.csv", data.table=F)
true.IDs = match(ID.map$pat.did,T2D.demo.dat.0$PATIENT_NUM)
summary(true.IDs)
TRUE.IDs = ID.map$pat.true[true.IDs]

#Create demographic + GRS data frame with true IDs instead of pseudo IDs
T2D.demo.dat = T2D.demo.dat.0; T2D.demo.dat$PATIENT_NUM = TRUE.IDs
head(T2D.demo.dat.0); head(T2D.demo.dat)
dim(T2D.demo.dat.0); dim(T2D.demo.dat)



T2D.fneg.sep = split(T2D.fneg.0,T2D.fneg.0$ID); 
for (i in seq_along(T2D.fneg.sep)){
  T2D.fneg.sep[[i]]$ID = i + max(T2D.fpos$ID) 
}
T2D.fneg.1 = do.call(rbind,T2D.fneg.sep[1:1000])
pos.feats = setdiff(colnames(T2D.fpos),colnames(T2D.fneg.1))
T2D.fneg.sup = matrix(0,nrow(T2D.fneg.1),length(pos.feats)-3)
colnames(T2D.fneg.sup) = pos.feats[-c(1:2,108)]
T2D.names = intersect(colnames(T2D.fpos),colnames(T2D.fneg.1))
T2D.fneg = cbind(T2D.fneg.1[,T2D.names],T2D.fneg.sup)
T2D.lab = rbind(T2D.fpos[,colnames(T2D.fneg)],T2D.fneg)
#Unlabeled patients
T2D.unlab = fread(T2D.dat[7],data.table=F)[,colnames(T2D.fneg)]


HF.dat = list.files("~/Dropbox/RiskPrediction_Deep/Data/HF",full.names=T)

#Includes filter positive AND filter negative patients
HF.fpos = fread(HF.dat[3],data.table=F)
HF.fneg.0 = fread(HF.dat[4],data.table=F); 
HF.fneg.sep = split(HF.fneg.0,HF.fneg.0$ID); 
for (i in seq_along(HF.fneg.sep)){
  HF.fneg.sep[[i]]$ID = i + max(HF.fpos$ID) 
}
HF.fneg = do.call(rbind,HF.fneg.sep[1:1950])
HF.names = intersect(colnames(HF.fpos),colnames(HF.fneg))
HF.lab = rbind(HF.fpos[,HF.names],HF.fneg[,HF.names])
#Unlabeled patients
HF.unlab = fread(HF.dat[7],data.table=F)[,HF.names]


fwrite(add_MAP(HF.lab,1,"PheCode:428.1","HF"),file = HF.dat[1])
#fwrite(add_MAP(HF.unlab,1,"PheCode:428.1","HF"),file = HF.dat[6])

fwrite(add_MAP(T2D.lab,1,"PheCode:250.2","T2D"),file = T2D.dat[1])
#fwrite(add_MAP(T2D.unlab,1,"PheCode:250.2","T2D"),file = T2D.dat[5])


T2D.lab = fread(T2D.dat[1],data.table=F)
T2D.lab$MAP[is.na(T2D.lab$MAP)] = 0
fwrite(T2D.unlab,file=T2D.dat[1])

T2D.unlab = fread(T2D.dat[5],data.table=F)
T2D.unlab$MAP[is.na(T2D.unlab$MAP)] = 0
fwrite(T2D.unlab,file=T2D.dat[5])


HF.lab = fread(HF.dat[1],data.table=F)
HF.lab$MAP[is.na(HF.lab$MAP)] = 0
fwrite(HF.unlab,file=HF.dat[1])

HF.unlab = fread(HF.dat[5],data.table=F)
HF.unlab$MAP[is.na(HF.unlab$MAP)] = 0
fwrite(HF.unlab,file=HF.dat[5])



# T2D.dat.0 = fread(all.dat[4],data.table=F)
# T2D.dat.1 = fread(all.dat[10],data.table=F); #T2D.dat.1$Y = 9
# shared.cols.T2D = intersect(colnames(T2D.dat.0),colnames(T2D.dat.1))
# T2D.dat = rbind(T2D.dat.0[,-c(1,2,166)],T2D.dat.1[,-c(1,165)])
# 
# HF.dat = fread(all.dat[12],data.table=F)
# HF.dat.0 = read.csv(all.dat[5])
# HF.dat.2 = read.csv(all.dat[7])
# HF.dat_new = HF.dat
# HF.dat_new$MAP = HF.dat.MAPs$MAP
# HF.dat_MAPs = EHR.DAT
# write.csv(HF.dat_MAPs,file = "/Users/isabella/Dropbox/RiskPrediction_Deep/Data/Dates_De_HF_PROD4_keser_V2.csv")
# #T1DM.dat = read.csv(all.dat[2])
# T2DM.dat = read.csv(all.dat[2])
# T2DM.dat.2 = read.csv(all.dat[3])
# 
# EHR.dat = T2D.dat.1; ICD.nm = "PheCode:250.2"; phe.nm = "T2DM"
# # EHR.dat = T1DM.dat; ICD.nm = "PheCode.250.1"; phe.nm = "T1DM"
# EHR.dat = HF.dat; ICD.nm = "PheCode:428.1"; phe.nm = "HF"
# EHR.dat = HF.dat.2; ICD.nm = c("PheCode.428.1","PheCode.428.2","PheCode.428.3","PheCode.428.4"); phe.nm = "HF"




# T2DM.sets = split(T2DM.dat,T2DM.dat$new_num); 
# HF.sets = split(HF.dat,HF.dat$new_num)
# for(i in seq_along(T2DM.sets)) {print(c(nrow(T2DM.sets[[i]]),unique(T2DM.sets[[i]]$H)))}
# for (i in seq_along(HF.sets)) {print(c(nrow(HF.sets[[i]]),unique(HF.sets[[i]]$H)))}
# 
 


#Create MAP probabilities at each time point for all patients.
create_MAP = function(EHR.dat,ICD.nm,phe.nm){
  #Separate longitudinal data by ID
  code.nm = grep("PheCode",names(EHR.dat))
  EHR.dat$utl = apply(EHR.dat[,c(code.nm)],1,sum)
  EHR.sets.0= split(EHR.dat,EHR.dat$ID)
  #;x$T = x$T/12 ; x
  TIMES = lapply(EHR.sets.0,function(x){times = x$T/12; ceiling(times - min(times))} )
  outlier.times = which(unlist(lapply(TIMES,max))  > 40)
  if(length(outlier.times) > 0) { EHR.sets.1 = EHR.sets.0[-c(outlier.times)]} else { EHR.sets.1 = EHR.sets.0 }
  EHR.sets = lapply(EHR.sets.1,function(x){x$T = x$T - min(x$T); x})
  #lapply(EHR.sets.0,function(x){times = x$T/12; if(max(ceiling(times - min(times))) > 35) {NULL}})
  #EHR.times = unique(unlist(TIMES))
  EHR.times = unique(unlist(TIMES[-c(outlier.times)]))
  #EHR.times = sort(unique(ceiling(EHR.dat$T/12))); yrs.total = ceiling(max(EHR.times)/12)
  
  MAP.probs = matrix(0,length(EHR.sets),1); EHR.dat$MAP = rep(0,nrow(EHR.dat))
  #MAP.probs = matrix(0,length(EHR.sets),length(EHR.times)); EHR.dat$MAP = rep(0,nrow(EHR.dat))
  
#for (i in seq_along(EHR.times)){

  MAP.dat = list()
  
  patient_num = as.matrix(0:(length(EHR.sets)-1))
  colnames(patient_num) = "patient_num"
  MAP.dat$ID = patient_num
  
  ICD.vec = unlist(lapply(EHR.sets,function(x){sum(x[,ICD.nm])}))
  ICD = as.matrix(ICD.vec)
  colnames(ICD) = paste(phe.nm,"_ICD",sep="")
  MAP.dat$mat = ICD
  
  
  utl = as.matrix(unlist(lapply(EHR.sets,function(x){sum(x[,"utl"])})))
  colnames(utl) = "utl"
  MAP.dat$note = utl
  
  #Compute MAP probabilities, based on ICD codes only
  if (length(unique(ICD)) <= 2 | mean(ICD > 0) < 0.01) {
    MAP.probs= rep(0,nrow(MAP.probs))
  } else {
    MAP.prob= tryCatch(MAP_flex_core(data = MAP.dat,  yes.con = FALSE, full.output = TRUE))
    MAP.probs = MAP.prob$scores.all[,2]
  }
  
  
  #j: patient; i: year
  for (j in seq_along(EHR.sets)){
    dat = EHR.sets[[j]]
    EHR.sets[[j]][, "MAP" ] = MAP.probs[j]
  }

  
  EHR.DAT = do.call(rbind,EHR.sets)
  # Times = unlist(lapply(EHR.sets,function(x){x[nrow(x),"T"]}))
  # Labels = unlist(lapply(EHR.sets,function(x){x[nrow(x),"Y"]}))
  
 return(EHR.DAT)
}




EHR.dat = HF.dat.labeled; ICD.nm = "PheCode:428.1"; phe.nm = "HF"

#Create MAP probabilities at each time point for all patients.
# create_MAP = function(EHR.dat,ICD.nm,phe.nm){
#   #Separate longitudinal data by ID
#   code.nm = grep("PheCode",names(EHR.dat))
#   EHR.dat$utl = apply(EHR.dat[,c(code.nm)],1,sum)
#   EHR.sets.0= split(EHR.dat,EHR.dat$ID)
#   #;x$T = x$T/12 ; x
#   #TIMES = lapply(EHR.sets.0,function(x){times = x$T/12; ceiling(times - min(times))} )
#   # TIMES = 
#   # outlier.times = which(unlist(lapply(TIMES,max))  > 40)
#   # EHR.sets.1 = EHR.sets.0[-c(outlier.times)]
#   # EHR.sets.1 = EHR.sets.0
#   #EHR.sets = lapply(EHR.sets.1,function(x){x$T = x$T - min(x$T); x})
#   EHR.sets= EHR.sets.0
#   #lapply(EHR.sets.0,function(x){times = x$T/12; if(max(ceiling(times - min(times))) > 35) {NULL}})
#   #EHR.times = unique(unlist(TIMES))
#   #EHR.times = unique(unlist(TIMES[-c(outlier.times)]))
#   #EHR.times = sort(unique(ceiling(EHR.dat$T/12))); yrs.total = ceiling(max(EHR.times)/12)
# 
#   MAP.probs = matrix(0,max(EHR.dat$T)+1,length(EHR.sets)+1); 
#   MAP.probs[,1] = 0:max(EHR.dat$T)
#   #EHR.dat$MAP = rep(0,nrow(EHR.dat))
#   #MAP.probs = matrix(0,length(EHR.sets),length(EHR.times)); EHR.dat$MAP = rep(0,nrow(EHR.dat))
# 
#   #for (i in seq_along(EHR.times)){
#   for (i in seq_len(max(EHR.dat$T))){
#     #t = EHR.times[i] + 1
#     print(i);
#     #print(t)
#     MAP.dat = list()
# 
#     patient_num = as.matrix(0:(length(EHR.sets)-1))
#     colnames(patient_num) = "patient_num"
#     MAP.dat$ID = patient_num
# 
#     ICD.vec = unlist(lapply(EHR.sets,function(x){sum(x[x$T <= i,ICD.nm])}))
#     ICD = as.matrix(ICD.vec)
#     colnames(ICD) = paste(phe.nm,"_ICD",sep="")
#     MAP.dat$mat = ICD
# 
# 
#     utl = as.matrix(unlist(lapply(EHR.sets,function(x){sum(x[x$T <= i,"utl"])})))
#     colnames(utl) = "utl"
#     MAP.dat$note = utl
# 
#     #Compute MAP probabilities, based on ICD codes only
#     if (length(unique(ICD)) <= 2 | mean(ICD > 0) < 0.01) {
#       MAP.probs[i+1,-1] = rep(0,nrow(MAP.probs))
#     } else {
#       MAP.prob= tryCatch(MAP_flex_core(data = MAP.dat,  yes.con = FALSE, full.output = TRUE))
#       MAP.probs[i+1,-1] = t(MAP.prob$scores.all[,2])
#     }
# 
# 
#   #   #j: patient; i: year
#   #   for (j in seq_along(EHR.sets)){
#   #     dat = EHR.sets[[j]]
#   #     EHR.sets[[j]][dat$T/12 > (i-1) & dat$T/12 <= i, "MAP" ] = MAP.probs[j,i]
#   #   }
#   # 
#    }
#   # 
#   # EHR.DAT = do.call(rbind,EHR.sets)
#   # for (l in seq_along(EHR.sets)){
#   #   EHR.DAT = rbind(EHR.DAT,EHR.sets[[l]])
#   # }
#   # EHR.DAT$MAP[is.na(EHR.DAT$MAP)] = 0
#   # 
#   # return(EHR.DAT)
#   return(MAP.probs)
# }
# 
# fwrite(MAP.probs,file=paste0("~/Dropbox/RiskPrediction_Deep/Longitudinal_Data/MAP_longitudinal/",phe.nm,"_MAP.probs.csv"))
# fwrite(MAP.probs,file=paste0("~/Dropbox/RiskPrediction_Deep/Longitudinal_Data/MAP_longitudinal/",phe.nm,"_MAP.probs.csv"))
# 

# T2DM.dat.MAPs = create_MAP(T2DM.dat[,-166],"PheCode.250.2","T2DM")
# summary(T2DM.dat.MAPs[,"MAP"])
# 
# T1DM.dat.MAPs = create_MAP(T1DM.dat,"PheCode.250.1","T1DM")
# summary(T1DM.dat.MAPs[,"MAP"])
# 
# HF.dat.MAPs = create_MAP(HF.dat,"PheCode.428.1","HF")
# HF.dat.MAPs$ID = HF.dat.MAPs$new_num
# HF.dat.MAPs$new_num = HF.dat.MAPs$X = NULL
# summary(HF.dat.MAPs[,"MAP"])

T2D.unlabeled = fread(all.dat[11],data.table=F)
View(T2D.unlabeled[T2D.unlabeled$ID == 66,])

T2D.u.sets = split(T2D.unlabeled,T2D.unlabeled$ID)
unlist(lapply(T2D.u.sets,function(x){if (nrow(x) <= 1) return (unique(x$ID))}))

# 
# HF.dat = fread(all.dat[7],data.table=F)
# T2D.dat = fread(all.dat[4],data.table=F)
# HF.dat.labeled = create_MAP(fread(all.dat[7],data.table=F),"PheCode:428.1","HF")
# T2D.dat.labeled = create_MAP(fread(all.dat[4],data.table=F),"PheCode:250.2","T2DM")
# HF.dat.unlabeled = create_MAP(fread(all.dat[14],data.table=F),"PheCode:428.1","HF")
# T2D.dat.unlabeled = create_MAP(fread(all.dat[11],data.table=F),"PheCode:250.2","T2DM")
# 
# 
# #[,-grep("PheCode:428.1",colnames(HF.dat.unlabeled))]
# 
# fwrite(HF.dat.labeled,
#        file = "/Users/isabella/Dropbox/RiskPrediction_Deep/Data/Dates_De_HF_PROD4_keser_V2.csv")
# fwrite(T2D.dat.labeled[,-1],
#        file = "/Users/isabella/Dropbox/RiskPrediction_Deep/Data/Dates_De_De_identied_De_identied_T2D_total__codified_V2.csv")
# fwrite(HF.dat.unlabeled,
#        file = "/Users/isabella/Dropbox/RiskPrediction_Deep/Data/HF_biobank_3000_unlabeled_de_identified_V2.csv")
# fwrite(T2D.dat.unlabeled,
#        file = "/Users/isabella/Dropbox/RiskPrediction_Deep/Data/De_T2D_unlabeled__codified_ONEICD_NON_overlap_V2.csv")
# 
# 
# 
# #write.csv(T1DM.dat.MAPs,file = "~/Dropbox/RiskPrediction_Deep/Data/T1DM.dat_withMAPs.csv", row.names=F)
# write.csv(T2DM.dat.MAPs,file = "/Users/isabella/Dropbox/RiskPrediction_Deep/Data/Dates_De_De_identied_De_identied_T2D_total__codified_V2.csv", row.names = F)
# write.csv(HF.dat_MAPs,file = "/Users/isabella/Dropbox/RiskPrediction_Deep/Data/Dates_De_HF_PROD4_keser_V2.csv",row.names=F)
# write.csv(HF.dat.MAPs,file = "~/Dropbox/RiskPrediction_Deep/Data/HF.dat_withMAPs.csv", row.names = F)
# 
# read.csv(all.dat[7])
# 
# T2DM.dat.MAPs = EHR.DAT
# fwrite(T2DM.dat.MAPs,file = "/Users/isabella/Dropbox/RiskPrediction_Deep/Data/De_T2D_unlabeled__codified_ONEICD_NON_overlap_V3.csv")
# 
# T2DM.dat.MAPs_new = fread(file = "/Users/isabella/Dropbox/RiskPrediction_Deep/Data/De_T2D_unlabeled__codified_ONEICD_NON_overlap_V2.csv",data.table = F)
# T2DM.dat.MAPs_new$MAP[is.na(T2DM.dat.MAPs_new$MAP)] = 0
# fwrite(T2DM.dat.MAPs_new,file = "/Users/isabella/Dropbox/RiskPrediction_Deep/Data/De_T2D_unlabeled__codified_ONEICD_NON_overlap_V2.csv")
# 
# 
# HF.dat.MAPs = EHR.DAT
# fwrite(HF.dat.MAPs,file = "/Users/isabella/Dropbox/RiskPrediction_Deep/Data/HF_biobank_3000_unlabeled_de_identified_V2.csv")
# 
# HF.dat.MAPs_new = fread(file = "/Users/isabella/Dropbox/RiskPrediction_Deep/Data/HF_biobank_3000_unlabeled_de_identified_V2.csv",data.table=F)
# 
# 
# 
# 
# T1DM.dat = read.csv(file = "~/Dropbox/RiskPrediction_Deep/Data/Dates_De_De_identied_De_identied_T1D_total__codified.csv")
# T2DM.dat = fread(file = "~/Dropbox/RiskPrediction_Deep/Data/De_T2D_unlabeled__codified_ONEICD_NON_overlap_V2.csv",data.table=F)
# HF.dat = read.csv(file = "~/Dropbox/RiskPrediction_Deep/Data/HF.dat_withMAPs.csv")
# 

#Create MAP probabilities at each time point for all patients.
#create_MAP = function(EHR.dat,ICD.nm,phe.nm){
  #Separate longitudinal data by ID
  code.nm = grep("PheCode",names(EHR.dat))
  EHR.dat$utl = apply(EHR.dat[,code.nm],1,sum)
  EHR.sets = split(EHR.dat,EHR.dat$new_num)
  EHR.times = unique(EHR.dat$T); yrs.total = ceiling(max(EHR.times)/12)
  
  MAP.probs = matrix(0,length(EHR.sets),yrs.total)
  
  for (i in seq_len(yrs.total)){
    print(i)
    MAP.dat = list()
    
    patient_num = as.matrix(0:(length(EHR.sets)-1))
    colnames(patient_num) = "patient_num"
    MAP.dat$ID = patient_num
    
    ICD.vec = unlist(lapply(EHR.sets,function(x){sum(x[x$T/12 <= i,ICD.nm])}))
    ICD = as.matrix(ICD.vec)
    colnames(ICD) = paste(phe.nm,"_ICD",sep="")
    MAP.dat$mat = ICD
    
    
    utl = as.matrix(unlist(lapply(EHR.sets,function(x){sum(x[x$T/12 <= i,"utl"])})))
    colnames(utl) = "utl"
    MAP.dat$note = utl
    
    #Compute MAP probabilities, based on ICD codes only
    if (length(unique(ICD)) <= 2 | mean(ICD > 0) < 0.01) {
      MAP.probs[,i] = rep(0,nrow(MAP.probs))
    } else {
      MAP.prob= tryCatch(MAP_flex_core(data = MAP.dat,  yes.con = FALSE, full.output = TRUE))
      MAP.probs[,i] = MAP.prob$scores.all[,2]
    }
    
    # k = 0
    # 
    # while (i - k > 0 & max(MAP.probs[,i-k]) == 0 & max(ICD) > 1) { 
    #   if (sum(dat$T/12 > (i-k-1) & dat$T/12 <= (i-k)) == 0){k = k+1 } 
    #   print(k)
    # }
    
    for (j in seq_along(EHR.sets)){
      dat = EHR.sets[[j]]
      EHR.sets[[j]][,paste0("MAP.Y",i)] = rep(0,nrow(EHR.sets[[j]]))
      if (i  == 1 ){ ind = dat$T/12 <= i} else {ind = dat$T/12 <= i & dat[,paste0("MAP.Y",i-1)] == 0}
      EHR.sets[[j]][ ind, paste0("MAP.Y",i) ] = MAP.probs[j,i]
    }
    
  }
  
  EHR.DAT = NULL

  return(EHR.DAT)
#}
