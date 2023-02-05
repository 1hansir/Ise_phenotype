#!/bin/bash

#module purge
#module load gcc/6.2.0 python/3.7.4 R/3.6.1
#unset PYTHONPATH
#source ~/Survival_env/bin/activate
#pip3 install numpy
#pip3 install torch
#pip3 install pycox

# Switch to target directory
cd C:\Users\NORTH\Dropbox\SemiSupervised_RiskPrediction\Yihan Updates

#--------------------------Split data to train set and test set

RScript Utilities/create_longdata.R

RScript Utilities/Pre_processing_v2.R ${num_train}
# 1.total_number 2.labeled_num 3.num_train 4. num_test 5.phe.nm 6. observation time
# e.g. RScript Utilities/get_fold_new.R 30000 1000 300 170 T2D 5

#--------------------------Train and Evaluate DeepHit
RScript Utilities/convert_data_NonLongitudinal.R
# 1.phe.nm 2. observation time(baseline time)(5 for simulation) 3. maximum time
# e.g. RScript Utilities/convert_data_NonLongitudinal.R T2D 5 30

RScript Method/run_DeepHit.R
# 1.num_train   2. phe.nm
# e.g. RScript Method/run_DeepHit.R 300 T2D

RScript Evaluation/eval_DeepHit.R
#  1.num_labels  2. phe.nm 3. observation_time 4. maximum time
# e.g. RScript Evaluation/eval_DeepHit.R 300 T2D 5 30

#--------------------------Train and Evaluate method
python Method/Method_AE/train_AE2.py --maximum_visits=26
# e.g python Method/Method_AE/train_AE2.py --phe_nm=T2D --epochs=4

#--------------------------Plot the result
RScript Real/plot_main_benchmark_tune.R ${max_time}
#  1.num_labels 2. phe.nm 3. observation_time 4. maximum time
# e.g. RScript Real/plot_main_benchmark_tune.R 300 T2D 5 30



