# -*- coding: utf-8 -*-
"""
TRAIN LAUNCHER
"""
import pandas as pd

from decimal import Decimal
from utilize_AE import get_data_from_csv
from GRU_AE2 import Attention_train
import os
import argparse
import random
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
from tensorflow.python.client import device_lib
import shutil

def get_available_gpus():
  local_device_protos = device_lib.list_local_devices()
  return [x.name for x in local_device_protos if x.device_type == 'GPU']

print(get_available_gpus())
os.environ['CUDA_VISIBLE_DEVICES'] = "0"
# 开启GPU加速

# TODO: change home directory
mdir = 'C:/Users/NORTH/Dropbox/SemiSupervised_RiskPrediction/Yihan Updates'    # your home file

# dirr_data: directory of data = dirr_save(will be save to the same directory)
# columns_min: first column of feature in the data(from index 0)   5 in simdata
# columns_max: last column of feature in the data    15 in simdata
# observation_months:         if int(data_Time[rowi])> dic_patient_firstDate[str(Patient_num[rowi])] + Observation_months:
#             dic_map_results.setdefault(str(Patient_num[rowi]), []).append(float(MAP_results[rowi]))
# year_window: total_month = int(int(max_date-min_date)/year_window)+1    set to 1
# visit_maximum: the corresponding size of the final saved data
# window_h: used to generate kernal: set to 6


def parse_arguments(parser):
    """Read user arguments"""
    parser.add_argument('--train_directory', type=str,
                        default=mdir+"/Real/T2D/",
                        help='Directory of train data ')
    parser.add_argument('--train_filename', type=str, default="_train.csv",
                        help='Filename of the train data')
    parser.add_argument('--train_filename_RP', type=str, default="_train.csv",
                        help='Filename of the train RP data')
    parser.add_argument('--test_directory', type=str,
                        default=mdir+"/Real/T2D/",
                        help='Directory of test data ')
    parser.add_argument('--test_filename', type=str, default="_test.csv",
                        help='Filename of the test data')
    parser.add_argument('--test_filename_RP', type=str, default="_test.csv",
                        help='Filename of the test RP data')
    # parser.add_argument('--embedding_filename', type=str, default="embedding.csv",
    #                     help='Filename of the embedding data in the train_directory ')
    parser.add_argument('--unlabel_directory', type=str,
                        default=mdir+"/Real/T2D/",
                        help='Directory of unlabeled data ')
    parser.add_argument('--unlabel_filename', type=str, default="_unlabeled.csv",
                        help='Filename of the unlabeled data')
    parser.add_argument('--unlabel_filename_RP', type=str, default="_unlabeled",
                        help='Filename of the unlabeled RP data')
    parser.add_argument('--file_risk_factors', type=str, default="T2D_risk_factors.csv",
                        help='Filename for list of risk factors')
    parser.add_argument('--save_directory', type=str,
                        default=mdir+"/Real/T2D/Main_result",
                        help='Directory to save the results')
    parser.add_argument('--results_filename', type=str, default="results_RETTAIN",
                        help='Filename to save the result data')
    parser.add_argument('--colums_min', type=int, default= 5,
                        help='data beginning column index  ')
    parser.add_argument('--colums_max', type=int, default=17,
                        help='data end column index')
    parser.add_argument('--epochs', type=int, default = 40,                   # 1.1 40
                        help='training epoches')
    # parser.add_argument('--embedding_dim', type=int, default=80,
    #                     help='dimension of embedding')
    parser.add_argument('--width_kernel', type=int, default=3,              # 1.1\1.2 3
                        help='width of kernel in months')
    parser.add_argument('--Observation_times', type=int, default=5,
                        help='Observation_times to predict')
    # parser.add_argument('--number_labels', type=int, default=400,
    #                    help='how many labeled patients to use to  train')
    parser.add_argument('--kernel_flag', type=int, default=1,
                        help='if use kernels to reweights: 1 to use; 0 not to use')
    parser.add_argument('--weight_unlabel', type=float, default=0.1,
                        help='weight of MAP to do semi-supervised')
    parser.add_argument('--time_window', type=int, default=3,
                        help='time_window: the time window to aggregate visits')
    parser.add_argument('--maximum_visits', type=int, default=26,     # 26 for T2D
                        help='maximum month-level visits for each patient; otherwise will be truncated')
    parser.add_argument('--Rep_num', type=int, default=1,
                        help='Replication of algorithm implementation')
    parser.add_argument('--flag_imputation', type=int, default=1,
                        help='Indicate whether or not to incorporate imputation step into model training')
    parser.add_argument('--phe_nm', type=str, default='T2D',
                        help='The phenotypes name to be predicted')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    print('---Model runs')
    PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)           # this will print default value when print args.help()
    ARGS = parse_arguments(PARSER)
    train_directory = ARGS.train_directory
    test_directory = ARGS.test_directory
    unlabel_directory = ARGS.unlabel_directory
    save_directory = ARGS.save_directory  # make sure this save dirr exist

    time_window = ARGS.time_window
    # watch out for this time_window
    width_kernel = ARGS.width_kernel
    #    embedding_dim = ARGS.embedding_dim
    Observation_times = ARGS.Observation_times
    kernel_flag = ARGS.kernel_flag  # if using kernels
    weight_unlabel = ARGS.weight_unlabel

    file_risk_factors = ARGS.file_risk_factors
    Rep_num = ARGS.Rep_num
    flag_imputation = ARGS.flag_imputation

    colums_min = ARGS.colums_min  #### begnning afer the  "Time" "Patient" "Label"
    colums_max = ARGS.colums_max
    maximum_visits = ARGS.maximum_visits
    epochs = ARGS.epochs  # how many epoches: 40
    phe_nm = ARGS.phe_nm

    # phe_nm = ARGS.phe_nm
    # Rep_num = ARGS.Rep_num
    phe_nm = 'T2D'
    print('phe_nm: ' ,phe_nm)
    results_filename = '/' + phe_nm + '_' + ARGS.results_filename  # results_RETTAIN


    filename_save_longitudinal = str(results_filename).split(".csv")[0]
    dirr_t = save_directory + filename_save_longitudinal + "_patient_longitudinal_predictions/"
    if os.path.exists(dirr_t):
        shutil.rmtree(dirr_t, ignore_errors=True)  # 如果存在这个文件夹，直接删除
        print('Previous Results deleted')

    for Rep_num in range(1,51):

        train_filename = phe_nm + "_train/REP_" + str(Rep_num) + ".csv"
        train_filename_RP = phe_nm + "_train/REP_"+ str(Rep_num) +".csv"

        test_filename = phe_nm + "_test/REP_"+ str(Rep_num) +".csv"
        test_filename_RP = phe_nm + "_test/REP_"+ str(Rep_num) +".csv"
        # embedding_filename = ARGS.embedding_filename  #in the train_directory
        # ARGS.train_directory  # ARGS.unlabel_directory
        # unlabel_filename =  ARGS.train_filename[0:3]+"_unlabeled__codified_ONEICD.csv" # #ARGS.train_filename[0:3]+"_unlabeled__codified.csv" #ARGS.unlabel_filename

        unlabel_filename= phe_nm + "_unlabeled.csv"
        unlabel_filename_RP= phe_nm + "_unlabeled.csv"



        print("---train_directory: ", train_directory)
        print("---train_filename: ", train_filename)
    #    print("---embedding_filename: ", embedding_filename)
        print("---test_directory: ", test_directory)
        print("---test_filename: ", test_filename)
        print("---unlabel_filename: ", unlabel_filename)


        dic_colums = {"Time": "T", "Patient": "ID", "Label": "Y", "MAP": "MAP"}
        ###use different column names for Time, here is "Patient_num;  patient column here is: "Patient_num"
        if not os.path.exists(save_directory):
            os.mkdir(save_directory)

        ###################Labeled patients################

        train_patients = (pd.read_csv(train_directory + "train_patients.csv"))["X" + str(Rep_num)]
        train_patients = list(train_patients); train_patients = [int(train_patients) for train_patients in train_patients]   #是否是train_patient?
        test_patients = (pd.read_csv(test_directory + "test_patients.csv"))["X" + str(Rep_num)]
        test_patients = list(test_patients); test_patients = [int(test_patients) for test_patients in test_patients]
        unlabeled_patients = (pd.read_csv(unlabel_directory + "unlabeled_patients.csv"))['unlabeled_patients']
        unlabeled_patients = list(unlabeled_patients); unlabeled_patients = [int(unlabeled_patients) for unlabeled_patients in unlabeled_patients]


        print("\n--------getting data for {}--------\n".format(phe_nm), "Rep: ", Rep_num)

        print("\n--------getting training data--------\n", "Rep: ", Rep_num)
        get_data_from_csv(train_directory, train_filename, train_filename_RP,
                          train_directory, colums_min=colums_min, colums_max=colums_max,
                          file_risk_factors=file_risk_factors,
                          visit_maximum=maximum_visits, weight_unlabel=weight_unlabel,
                          dic_items=dic_colums, year_window=1, train_mode="train", window_h=width_kernel,
                          patient_targeted=train_patients, Observation_months=Observation_times,
                          kernal_flag=kernel_flag)
        print("\n--------getting testing data--------\n", "Rep: ", Rep_num)
        get_data_from_csv(test_directory, test_filename, test_filename_RP,
                          test_directory, colums_min=colums_min, weight_unlabel=weight_unlabel,
                          colums_max=colums_max, file_risk_factors=file_risk_factors, visit_maximum=maximum_visits,
                          dic_items=dic_colums, year_window=1, train_mode="test", window_h=width_kernel,
                          patient_targeted=test_patients, Observation_months=Observation_times,
                          kernal_flag=kernel_flag)

        if Rep_num == 1:
        # Train and Test file must be generated for every replicate, but unlabeled data only need to be generate once.

            print("\n--------getting unlabeled data--------\n")
            get_data_from_csv(unlabel_directory, unlabel_filename,unlabel_filename_RP,
                                      unlabel_directory, colums_min=colums_min, weight_unlabel = weight_unlabel,
                                      colums_max=colums_max, file_risk_factors = file_risk_factors,visit_maximum=maximum_visits,
                                      dic_items=dic_colums,year_window=1,train_mode="unlabeled",window_h=width_kernel,
                                      patient_targeted=unlabeled_patients[1:],Observation_months=Observation_times,kernal_flag = kernel_flag)


        print("--------training--------",phe_nm,'Rep:',Rep_num)

        train_file = str(train_filename).split(".csv")[0] + str(train_filename).split(".csv")[1]
        test_file = str(test_filename).split(".csv")[0] + str(test_filename).split(".csv")[1]
        unlabel_file = str(unlabel_filename).split(".csv")[0] + str(unlabel_filename).split(".csv")[1]


        Attention_train(dirr_train=train_directory + "PKL_files/",
                      filename_train=train_file + "_train_kernel_"+str(Decimal(kernel_flag))+"_weight_unlabel_"+str(round(Decimal(weight_unlabel),5))+".pkl",
                      dirr_test=test_directory + "PKL_files/",
                      filename_test=test_file + "_test_kernel_"+str(Decimal(kernel_flag))+"_weight_unlabel_"+str(round(Decimal(weight_unlabel),5))+".pkl",
                      dirr_unlabel=unlabel_directory + "PKL_files/",
                      filename_unlabel=unlabel_file + "_unlabeled_kernel_"+str(Decimal(kernel_flag))+"_weight_unlabel_"+str(round(Decimal(weight_unlabel),5))+".pkl",
                      dirr_save=save_directory,
                      filename_save=results_filename,
                      epochs=epochs, colums_min=colums_min,
                      colums_max=colums_max,max_visits=maximum_visits,
                      year_window=time_window,weight_unlabel=weight_unlabel,
                      flag_imputation = flag_imputation,epochs_imputation=10,Rep = Rep_num)  # original 20

        print("\n--------Finish Training for {}!--------".format(phe_nm), ' Rep:', Rep_num,"\n")


