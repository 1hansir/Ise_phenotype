#!/bin/bash

#module purge
#module load gcc/6.2.0 python/3.7.4 R/3.6.1
#unset PYTHONPATH
#source ~/Survival_env/bin/activate
#pip3 install numpy
#pip3 install torch
#pip3 install pycox

cd C:\Users\NORTH\Dropbox\SemiSupervised_RiskPrediction\Yihan Updates


RScript Utilities/Pre_processing.R 12919 600 300 160 T2D 3

#--------------------------Train and Evaluate DeepHit
RScript Utilities/convert_data_NonLongitudinal.R T2D 3 20


RScript Method/run_DeepHit.R 300 T2D


RScript Evaluation/eval_DeepHit.R 300 T2D 3 20

#--------------------------Train and Evaluate method

python Method/Method_AE/train_AE2.py --phe_nm=T2D --epochs=30 --maximum_visit=20 --colums_max=17     # 17= 5+12(X.12)

#--------------------------Plot the result

RScript Real/plot_main_benchmark_tune.R 300 T2D 3 20