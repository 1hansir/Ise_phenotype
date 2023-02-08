source("./set_mdir.R")
# TODO: change the home dir

# load packages
library(MAP)
library(ff)
library(data.table)
library(ggplot2)
library(ggpubr)

#----------------------------------------
# data generation #
#----------------------------------------
source(paste0(mdir,'/Real/benchmark_stm_kwr.R'))

source(paste0(mdir,'/Real/pred_csv_to_list.R'))

dir.benchmark = paste0(mdir,'/Evaluation/Results/DeepHit/T2D/')


args = commandArgs(trailingOnly = TRUE)


max_time = as.numeric(args[1])
#data.file = args[2]    # == phe.nm
data.file = 'T2D'
max_time = 25
#----------------------------------------
# benchmark methods #
#----------------------------------------

# stm.ts.auc

# TODO: adapt(or abandon) the 'for' loop and set data.file




simulation_stm(50,paste0(mdir,'/Real/T2D/',data.file,'.Rds'),
        paste0(mdir,'/Real/T2D/train_patients.csv'),
        paste0(mdir,'/Real/T2D/test_patients.csv'),
        result.file=paste0(mdir,'/Real/Benchmark/stm.auc.',data.file,'.Rds'),time.grid = 1:max_time)



print('stm method complete')

# kw.ts.auc
# TODO: adapt(or abandon) the 'for' loop and set data.file



#simulation_kw(50,paste0(mdir,'/Real/T2D/',data.file,'.Rds'),
#        paste0(mdir,'/Real/T2D/train_patients_',num_train,'.csv'),
#        paste0(mdir,'/Real/T2D/test_patients.csv'),
#        paste0(mdir,'/Real/T2D/unlabeled_patients.csv'),
#        result.file=paste0(mdir,'/Real/Benchmark/kw.auc.',data.file,'.Rds'),time.grid = (min_t+1):max_t)



print('kw method complete')

#----------------------------------------
# plot out #
#----------------------------------------


plot_methods_dh = function(phe.nm){

    # TODO: change the path to utilize T2D data;
    main_dir = paste0(mdir, "/Real/T2D/Main_result/",phe.nm,"_results_RETTAIN_patient_longitudinal_predictions/")   #where we store the result from main method
    # dir3 = paste0("C:/Users/sirius/Dropbox/2022summer/intern/Data/FinalData0716/SimDat.",data.name,"/result3/")

    # TODO: adapt the time_grid
    auc.main2 = ts.auc.CI(auc.main.method(paste0(main_dir),time.grid = 1:max_time),1:max_time)

    # TODO: adapt the time_grid
   # auc.kw = ts.auc.CI(readRDS(paste0(mdir,'/Real/Benchmark/kw.auc.',phe.nm,'.Rds')),(min_t+1):max_t)

    # TODO: adapt the time_grid
    auc.stm = ts.auc.CI(readRDS(paste0(mdir,'/Real/Benchmark/stm.auc.',phe.nm,'.Rds')),1:max_time)

    # TODO: adapt the time_grid
    auc.deephit = ts.auc.CI(read.csv(paste0(dir.benchmark,phe.nm,'_tspec_AUC-ALL.deephit_labels.csv'))[,2:51],1:max_time)

    ggdf =  data.frame(time = c( auc.main2$time,
                              # auc.kw$time,
                              auc.stm$time,
                               auc.deephit$time),
                     auc = c( auc.main2$mean,
                             # auc.kw$mean,
                              auc.stm$mean,
                              auc.deephit$mean),
                     low = c( auc.main2$low,
                             # auc.kw$low,
                             auc.stm$low,
                              auc.deephit$low),
                     up = c( auc.main2$up,
                            # auc.kw$up,
                            auc.stm$up,
                             auc.deephit$up),
                     method = c(rep('Main Method',length(auc.main2$time)),
                              #  rep('Kernel Weighted Time-specific Regression',length(auc.kw$time)),
                               rep('Semiparametric Transformation Model',length(auc.stm$time)),
                                rep('Deep Hit',length(auc.deephit$time)))
    )
  # ggplot
  ggplot(ggdf, aes(x = time , y = auc, color = method))+    # x = time - 5(represents postbaseline)
    theme_bw() +
    geom_smooth(aes(ymin = low, ymax = up, fill = method),
                stat = "identity")+  scale_alpha(range = c(0.1, 0.4)) + 
    ylab("AUC") + xlab("Year") +
    xlim(c(1,max_time)) + ylim(c(0,1)) +
    # TODO: change the path to utilize T2D data
    ggtitle(paste0(phe.nm)) + theme(plot.title = element_text(hjust = 0.5),legend.position="right", legend.direction="vertical")
}

# TODO: change the phe.nm

g_ =plot_methods_dh(data.file)

#g.1.1=plot_methods_dh('SimDat.1.1')
#g.1.2=plot_methods_dh('SimDat.1.2')
#g.2.1=plot_methods_dh('SimDat.2.1')
#g.2.2=plot_methods_dh('SimDat.2.2')


pdf(paste0(mdir,"/Real/Plots/Time specific AUC plots.pdf"),
    width = 13, height = 8)
ggarrange(g_ ,
          common.legend = TRUE, legend = "right")
dev.off()

# Cstatis
print('Cstat for Main-Method: ')
c.stat_main_v1 = function(data.name){
  main_dir = paste0(mdir, "/Real/T2D/Main_result/Binary/")
  # main method setting AUC
  n.repeat = 4
  C = numeric(n.repeat)
  for (i in 1:n.repeat) {
    dat = read.csv(paste0(main_dir,data.name,'_results_RETTAIN/Rep_',i,'/Cstat_results.csv'))
    C[i] = auc(dat$Y,dat$Predition,levels=c(0,1),direction='<')
  }
 return(c(mean(C),sd(C)))
}

c.stat_main_v2 = function(data.name){      # 从整体Summary得到Cstatis
  main_dir = paste0(mdir, "/Real/T2D/Main_result/")
  # main method setting AUC
  n.repeat = 4

  C = numeric(n.repeat)
  for (i in 0:n.repeat-1) {
    dat = read.csv(paste0(main_dir,data.name,"_results_RETTAIN_Attention_metric_binaryMLP.txt"))
    # print(unlist(strsplit(dat[,7] , split = ":"))[2])
    C[i] = as.numeric(unlist(strsplit(dat[i,7] , split = ":"))[2])
    #print(C[i])
  }

 # print(C)
 return(c(mean(C),sd(C)))
}

c.stat_main_v1('T2D')
c.stat_main_v2('T2D')