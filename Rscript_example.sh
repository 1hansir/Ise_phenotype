#!/bin/bash


# Switch to target directory
cd C:\Users\NORTH\Dropbox\REAL\SemiSupervised_RiskPrediction

#--------------------------Split data to train set and test set

RScript Utilities/get_fold_new.R 30000 1000 300 170 T2D 5

#--------------------------Train and Evaluate DeepHit
RScript Utilities/convert_data_NonLongitudinal.R T2D 5 30


RScript Method/run_DeepHit.R 300 T2D


RScript Evaluation/eval_DeepHit.R 300 T2D 5 30

#--------------------------Train and Evaluate method

python Method/Method_AE/train_AE2.py --phe_nm=T2D --epochs=4

#--------------------------Plot the result

RScript Real/plot_main_benchmark_tune.R 300 T2D 5 30



