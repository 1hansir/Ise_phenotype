source(paste0(mdir,'/Evaluation/measure.R'))

# Transform the csv files in folder to list
#-----------------------------------------------
# fdir: path of the folder
# 
# Assume the csv files are in format
#   file name: patient id
#   col 1: status at period 
#         0 - before onset/censoring
#         1 - after onset before censoring
#        NA - after censoring
#   col 2: predicted risk
#   col 3: period index
# Output: 
#   info: a data frame with id, follow up time
#         and event time (NA if no event)
#   pred: a list for predicted risk for periods
#         1 through end of followup

pred.csv.to.list=function(fdir)
{
  csvlist = list.files(fdir)
  npred = length(csvlist)
  
  idlist = gsub("\\.csv",'',csvlist)
  
  pred = vector("list",npred)
  names(pred) = idlist
  
  info = data.frame(id = idlist)
  info$event = Inf
  info$fu = NA
  
  for(i in 1:npred)
  {
    tab = read.csv(paste(fdir, csvlist[i],sep='/'),
                   header = FALSE)
    info$fu[i] = tab$V3[max(which(!is.na(tab$V1)))]    # 由于中间有空缺，因此使用绝对坐标，而不再使用之前的累加方法
    #info$fu[i] = sum(!is.na(tab$V1))-1
    info$event[i] = min(which(tab$V1==1))
    
    pred[[i]] = tab$V2[1:info$fu[i]+1]
  }
  
  return(list(info = info, 
              pred = pred))
}


## get main auc ###
## input: result sub-folders directory for main method
## output: length(time.grid)*100 matrix, rows are time.grids, columns are replications
#--------------------------
auc.main.method = function(dir,time.grid){
  auc_TS = matrix(NA,length(time.grid),100)
  get_auc = function(i){
    list = pred.csv.to.list(paste0(dir,'Rep_',i,'/'))
    cs = list$info$event < list$info$fu
    pred = sapply(list$pred, function(x) rev(x)[1])
    inc.prob = t(sapply(list$pred, function(x) x[time.grid]))
    return(auc.tspec(list$info$event, list$info$fu, time.grid, inc.prob))
  }
  for (i in 1:100) {              # temporarily, we just run for 5 repulicate
    auc_TS[,i] = tryCatch({get_auc(i)},
      error=function(cond) {return(NA)})
    }
  # print(auc_TS)
  return(auc_TS)

}

# create CI dataframe #
# input: length(time.grid)*100 matrix, rows are time.grids, columns are replications
# output: dataframe of ts-AUC for plot,  $mean,$5-th percentage, $95-th percentage, $time.grid
#------------------------------------
ts.auc.CI = function(data,time.grid){
  data.frame(time = time.grid,
             mean = apply(data, 1, function(x) mean(x, na.rm=TRUE)),
             low = apply(data, 1, function(x) quantile(x,0.05, na.rm=TRUE)),
             up = apply(data, 1, function(x) quantile(x,0.95, na.rm=TRUE)))
}
