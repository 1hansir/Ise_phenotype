
#Risk Prediction: MAP for silver standard labels at each one year time window
#Datasets: HF, T1DM and T2DM
#Note: T1DM is a juvenile onset disease, so may be less compelling from a clinical point of view for risk prediction. 
#As in, we know that genetics play a very large, if not the largest, role in determining patient risk of disease (and also [<perhaps>] age)



library(data.table)
library(ff)
library(MAP)




compute_MAP = function(EHR.sets,i,ICD.nm,phe.nm){
  
  MAP.dat = list()
  
  patient_num = as.matrix(0:(length(EHR.sets)-1))
  colnames(patient_num) = "patient_num"
  MAP.dat$ID = patient_num
  
  ICD.vec = unlist(lapply(EHR.sets,function(x){sum(x[x$T <= i,ICD.nm])}))
  ICD = as.matrix(ICD.vec); ICD[is.na(ICD)] = 0
  colnames(ICD) = paste(phe.nm,"_ICD",sep="")
  MAP.dat$mat = ICD
  
  utl = as.matrix(unlist(lapply(EHR.sets,function(x){sum(x[x$T <= i,"utl"])})))
  colnames(utl) = "utl"
  MAP.dat$note = utl
  
  #Compute MAP probabilities, based on ICD codes only
  if (length(unique(ICD)) <= 2 | mean(ICD > 0) < 0.01) {
    MAP.probs= rep(0,nrow(MAP.probs))
  } else {
    MAP.prob= tryCatch(MAP_flex_core(data = MAP.dat,  yes.con = FALSE, full.output = TRUE))
    MAP.probs = MAP.prob$scores.all[,2]
  }
  
  return(MAP.probs)
}

#Need to modify EHR data so that there is one row for each month
#Create MAP probabilities at each time point for all patients.
add_MAP = function(EHR.dat,mo.window,ICD.nm){
  #Separate longitudinal data by ID
  code.nm = grep("PheCode",names(EHR.dat))
  EHR.dat$utl = apply(EHR.dat[,c(code.nm)],1,sum)
  EHR.sets.0= split(EHR.dat,EHR.dat$ID)
  TIMES = lapply(EHR.sets.0,function(x){times = x$T/mo.window; ceiling(times - min(times))} )
  EHR.sets = lapply(EHR.sets.0,function(x){x$T = x$T - min(x$T); x})
  All.times = unlist(lapply(EHR.sets,function(x){max(x$T)}))
  MAP.probs = matrix(0,length(EHR.sets),max(All.times)+1); EHR.dat$MAP = rep(0,nrow(EHR.dat))
  
  for (i in seq_len(max(All.times)+1)){
  
    MAP.dat = list()
    
    patient_num = as.matrix(0:(length(EHR.sets)-1))
    colnames(patient_num) = "patient_num"
    MAP.dat$ID = patient_num
    
    ICD.vec = unlist(lapply(EHR.sets,function(x){sum(x[x$T <= i,ICD.nm])}))
    ICD = as.matrix(ICD.vec); ICD[is.na(ICD)] = 0
    colnames(ICD) = "Main_ICD"
    MAP.dat$mat = ICD
    
    utl = as.matrix(unlist(lapply(EHR.sets,function(x){sum(x[x$T <= i,"utl"])})))
    colnames(utl) = "utl"
    MAP.dat$note = utl
    
    #Compute MAP probabilities, based on ICD codes only
    if (length(unique(ICD)) <= 2 | mean(ICD > 0) < 0.001) {
      MAP.probs[,i]= rep(0,nrow(MAP.probs))
    } else {
      MAP.prob = tryCatch(MAP_flex_core(data = MAP.dat,  yes.con = FALSE, full.output = TRUE))
      MAP.probs[,i] = MAP.prob$scores.all[,2]
    }
   
    #j: patient; i: year
    for (j in seq_along(EHR.sets)){
      dat = EHR.sets[[j]]
      EHR.sets[[j]][dat$T/mo.window > (i-1) & dat$T/mo.window <= i, "MAP" ] = MAP.probs[j,i]
    }
  } 
  
  EHR.DAT = do.call(rbind,EHR.sets)
  EHR.DAT$MAP[is.na(EHR.DAT$MAP)] = 0
  return(EHR.DAT)
}


add_MAP_SimDat = function(EHR.dat,mo.window,ICD_NLP.nm){
  #Separate longitudinal data by ID
  EHR.dat$utl = EHR.dat$S.1
  EHR.sets.0= split(EHR.dat,EHR.dat$ID)
  TIMES = lapply(EHR.sets.0,function(x){times = x$T/mo.window; ceiling(times - min(times))} )
  EHR.sets = lapply(EHR.sets.0,function(x){x$T = x$T - min(x$T); x})
  All.times = unlist(lapply(EHR.sets,function(x){max(x$T)}))
  MAP.probs = matrix(0,length(EHR.sets),max(All.times)+1); EHR.dat$MAP = rep(0,nrow(EHR.dat))
  
  for (i in seq_len(max(All.times)+1)){
    
    MAP.dat = list()
    
    patient_num = as.matrix(0:(length(EHR.sets)-1))
    colnames(patient_num) = "patient_num"
    MAP.dat$ID = patient_num
    
    ICD_NLP.vec = matrix(0,length(EHR.sets),length(ICD_NLP.nm))
    for (s in seq_along(ICD_NLP.nm)){
      ICD_NLP.vec[,s] = unlist(lapply(EHR.sets,function(x){sum(x[x$T <= i,ICD_NLP.nm[s]])}))
    }
    
    #ICD_NLP = as.matrix(ICD_NLP.vec); ICD_NLP[is.na(ICD_NLP)] = 0
    ICD_NLP = ICD_NLP.vec; ICD_NLP[is.na(ICD_NLP)] = 0
    #colnames(ICD_NLP) = paste(phe.nm,"_ICD_NLP",sep="")
    colnames(ICD_NLP) = ICD_NLP.nm
    MAP.dat$mat = ICD_NLP
    
    utl = as.matrix(unlist(lapply(EHR.sets,function(x){sum(x[x$T <= i,"utl"])})))
    colnames(utl) = "utl"
    MAP.dat$note = utl
    
    #Compute MAP probabilities, based on ICD_NLP codes only
    if (length(unique(ICD_NLP[,1])) <= 2 | length(unique(ICD_NLP[,2])) <= 2 | 
        mean(ICD_NLP[,1] > 0) < 0.001 | mean(ICD_NLP[,2] > 0) < 0.001) {
      MAP.probs[,i]= rep(0,nrow(MAP.probs))
    } else {
      MAP.prob = tryCatch(MAP_flex_core(data = MAP.dat,  yes.con = FALSE, full.output = TRUE))
      MAP.probs[,i] = MAP.prob$scores.all[,2]
    }
    
    
    #j: patient; i: year
    for (j in seq_along(EHR.sets)){
      dat = EHR.sets[[j]]
      EHR.sets[[j]][dat$T/mo.window > (i-1) & dat$T/mo.window <= i, "MAP" ] = MAP.probs[j,i]
    }
  } 
  
  EHR.DAT = do.call(rbind,EHR.sets)
  EHR.DAT$MAP[is.na(EHR.DAT$MAP)] = 0
  return(EHR.DAT)
}

