from decimal import Decimal
import pandas as pd
import numpy as np
import pickle
from sklearn.utils import shuffle
import math

'''
Three main steps:
1. generate data for all and risk_prediction only
2. generate kernel_weight for every patients
3. 

'''

# train_directory = "~/Dropbox/RiskPrediction_Deep/Data/T2D/"
# train_filename = "Dates_De_De_identied_De_identied_T2D_total__codified_V2.csv"
# embedding_filename = "embedding_codified_T2D_codes_PCA80.csv"
# test_directory = "~/Dropbox/RiskPrediction_Deep/Data/T2D/"
# test_filename  = "Dates_De_De_identied_De_identied_T2D_total__codified_V2.csv"
# file_risk_factors  = "T2D_risk_factors.csv"
# unlabel_directory = "~/Dropbox/RiskPrediction_Deep/Data/T2D/"
# unlabel_filename = "De_T2D_unlabeled__codified_ONEICD_NON_overlap_V2.csv"
# save_directory  = "~/Dropbox/RiskPrediction_Deep/Results/T2D/"
# results_filename  = "T2D_kernel_1_numlabels_200-weight_unlabel_0.5-Rep_1.csv"
# Rep_num = 1
# colums_min = 3
# colums_max = 13
# year_window = 3
# embedding_dim = 80
# width_kernel = 6
# maximum_visits = visit_maximum = 300
# number_labels = 120
# Observation_months = 24
# kernal_flag = 1
# weight_unlabel = 0.5
# epochs = 80
#
#
#
# train_patients = (pd.read_csv(train_directory + "train_patients_" + str(number_labels) + ".csv"))["X" + str(Rep_num)]
# train_patients = list(train_patients)
# test_patients = (pd.read_csv(test_directory + "test_patients.csv"))["X" + str(Rep_num)]
# test_patients = list(test_patients)
# #Check that training and testing sets don't overlap
# set(test_patients).intersection(train_patients)
# unlabeled_patients = (pd.read_csv(unlabel_directory + "unlabeled_patients.csv"))['unlabeled_patients']
# unlabeled_patients = list(unlabeled_patients)
#
# patient_targeted = test_patients
# patient_targeted
# train_mode = "test"
# dirr_data = dirr_save = test_directory
# filename = test_filename
# filename = "T2D_labeled_RP.csv"
# filename_RP = "T2D_labeled_RP.csv"
# window_h = 6


# dirr_data: directory of data = dirr_save(will be save to the same directory)
# RP
# columns_min: first column of feature in the data(from index 0)   5 in simdata
# columns_max: last column of feature in the data    15 in simdata
# observation_months:         if int(data_Time[rowi])> dic_patient_firstDate[str(Patient_num[rowi])] + Observation_months:
#             dic_map_results.setdefault(str(Patient_num[rowi]), []).append(float(MAP_results[rowi]))
# year_window: total_month = int(int(max_date-min_date)/year_window)+1    set to 1
# visit_maximum: the corresponding size of the final saved data
# window_h: used to generate kernal: set to 6

# file_risk_factors is useless
def get_data_from_csv(dirr_data, filename, filename_RP,dirr_save,colums_min,colums_max,visit_maximum,file_risk_factors,
                              dic_items,year_window=1,train_mode="test",window_h=6,weight_unlabel = 0.4,
                      patient_targeted=["dafda"],Observation_months=2, kernal_flag=1):
    print ("------get data from csv:  ",dirr_data,": ",filename)
    number_maximum=50                 # maximum number of code counts: what is this for?
    dirr=dirr_data
    dirr_save=dirr_save
    df = pd.read_csv(dirr + filename)
    df_RP = pd.read_csv(dirr + filename_RP)

    # risk_factors = (pd.read_csv(dirr + file_risk_factors))["Risk_Var"]
    # feat_risk_keep = []
    # for i in range(len(risk_factors)):
    #     if("RXNORM" not in risk_factors[i]):
    #         feat_risk_keep.append(risk_factors[i])
    # feat_risk_inds = []
    # all_nms = list(df.columns)
    # for i in range(len(feat_risk_keep)):
    #     if (feat_risk_keep[i] in all_nms):
    #         feat_risk_inds.append(all_nms.index(feat_risk_keep[i]))

    # print("feature_column:", df.columns[colums_min:colums_max])
    data_array_all=np.array(df[df.columns[colums_min:colums_max]])
    data_array_RP_all = np.array(df_RP[df.columns[colums_min:colums_max]])
    #print ("df.columns: ",df.columns)
    data_Time=list(df["T"])
    Patient_NUM = list(df["ID"])
    Patient_num= [int(Patient_NUM) for Patient_NUM in Patient_NUM]
    Y_label=list(df["Y"])
    MAP_results = list(df["MAP"])


    dic_patient_firstDate = {}
    dic_patient_label={}
    dic_patient_lastDate = {}
    dic_map_results = {}

    dic_patient_label[str(Patient_num[0])] = int(Y_label[0])
    dic_patient_firstDate[str(Patient_num[0])] = int(data_Time[0])
    dic_patient_lastDate[str(Patient_num[0])] = int(data_Time[0])
    # initialization


    for rowi in range(1, len(data_Time)):
        dic_patient_label[str(Patient_num[rowi])] = int(Y_label[rowi])

        # print("Row:", rowi, ". New patient: ", not Patient_num[rowi - 1] == Patient_num[rowi],
        # ". Post baseline observation period: ", int(data_Time[rowi])> dic_patient_firstDate[str(Patient_num[rowi])]+Observation_months,
        # ". Patient ID: ", Patient_num[rowi])
        if not Patient_num[rowi - 1] == Patient_num[rowi]:
            dic_patient_firstDate[str(Patient_num[rowi])] = int(data_Time[rowi])
            dic_patient_lastDate[str(Patient_num[rowi-1])] = int(data_Time[rowi-1])
        # smart!
        
        if int(data_Time[rowi])> dic_patient_firstDate[str(Patient_num[rowi])] + Observation_months:
            dic_map_results.setdefault(str(Patient_num[rowi]), []).append(float(MAP_results[rowi]))
        # dict.setdefault(key, value): check if key exist, if not, create one and value it. 
    ############

    #################getting map results#########using the last five predictions as the pseudo-label

    key_feature1 = list(df.columns[1])         ###if need key features
    key_feature2 = list(df.columns[2])         ###if need key features

    print ("------data_array_all.shape: ",np.array(data_array_all).shape)
    print ("------data_array_RP_all.shape: ",np.array(data_array_RP_all).shape)
    min_date=min(np.array(data_Time,dtype=np.int))
    max_date = max(np.array(data_Time, dtype=np.int))
    total_month = int(int(max_date-min_date)/year_window)+1
    total_codes= colums_max-colums_min

    # df_embedding = pd.read_csv(dirr_data+file_embedding)
    #
    # data_embedding_all=[]
    # for colmi in df.columns[colums_min:colums_max]:
    #     data_embedding_all.append(df_embedding[colmi])
    # data_embedding_all = np.array(data_embedding_all, dtype=np.float)
    # print("---data_embedding_all: ", data_embedding_all.shape)
    # data_embedding_all = data_embedding_all.transpose()
    patients_total = {}
    patients_total_RP = {}
    dates_patient_total={}
    patient_total_valid={}

    for rowi in range(len(data_Time)):
        patient = str(Patient_num[rowi])
        ##having data after observation period, and in data set for relevant training phase (i.e. train, test, or unlabeled)
        if (int(patient) in patient_targeted) and str(patient) in dic_map_results :
            if (int(patient) in patient_targeted or 'ALL' in patient_targeted) and str(patient) in dic_map_results :  ##having data after observation period,this is unnecessary
                if not patient in patients_total:
                    patient_col = np.ones(shape=(total_month, 1)) * int(Patient_num[rowi])
                    label_col = np.zeros(shape=(total_month, 1))
                    time_col = np.array(list(range(total_month))) * 0
                    time_col = np.expand_dims(time_col, axis=-1)
                    patient_time_col = np.concatenate((patient_col, time_col), axis=-1)
                    patient_time_col = np.concatenate((patient_time_col, label_col), axis=-1)

                    
                    patients_total[patient] = np.concatenate((patient_time_col, np.zeros(shape=(total_month,
                                                                                    total_codes+3))), -1)  ###########total_codes+2: phecodes
                    patients_total_RP[patient] = np.concatenate((patient_time_col, np.zeros(shape=(total_month,
                                                                              total_codes+3))), -1)  ###########total_codes+2: phecodes
                    #  —— ID | T | label | total_codes+3+
                    #  total_month
                    #
                    
                if patient in patients_total and patient in dic_patient_firstDate and patient in dic_patient_lastDate:
                    index_row = int(int(int(data_Time[rowi]) - min_date) / year_window)

                    if int(data_Time[rowi]) < dic_patient_firstDate[patient] + Observation_months:      ######observation data
                        patient_total_valid[patient]=1
                        patients_total[patient][index_row, 6:6+colums_max-colums_min] += data_array_all[rowi]

                        patients_total_RP[patient][index_row, 6:6+colums_max-colums_min] += data_array_RP_all[rowi]

                        #########make sure the code counts no more than number_maximum to dominate the embedding
                        patients_total[patient][index_row, 6:6 + colums_max - colums_min] = \
                            np.minimum(patients_total[patient][index_row, 6:6 + colums_max - colums_min], number_maximum)


                        patients_total_RP[patient][index_row, 6:6 + colums_max - colums_min] = \
                            np.minimum(patients_total_RP[patient][index_row, 6:6 + colums_max - colums_min],
                                       number_maximum)

                        #################
                        # the censor time as the Keyfeature1
                        patients_total[patient][index_row, 5] = 1
                        patients_total[patient][index_row, 3] = dic_patient_lastDate[patient] - dic_patient_firstDate[patient] -Observation_months
                        #############the MAP results as the keyfeature2: using the last 20 windows MAP results
                        patients_total[patient][index_row, 4] = np.mean(dic_map_results[patient][-1:])  ###the map results: the last 3 monthely visits
                        #patients_total[patient][index_row, 4] = np.mean
                        patients_total[patient][index_row, 2] = int(dic_patient_label[patient])  # int(Y_label[rowi])
                        patients_total[patient][index_row, 1] = int(data_Time[rowi]-(dic_patient_firstDate[patient] +Observation_months))
                        patients_total[patient][index_row, 0] = int(Patient_num[rowi])
                        
                        #############
                        patients_total_RP[patient][index_row, 5] = 1
                        patients_total_RP[patient][index_row, 3] = dic_patient_lastDate[patient] - dic_patient_firstDate[patient] -Observation_months
                        #############the MAP results as the keyfeature2: using the last 20 windows MAP results
                        patients_total_RP[patient][index_row, 4] = np.mean(dic_map_results[patient][-1:])  ###the last one?
                        #patients_total[patient][index_row, 4] = np.mean
                        patients_total_RP[patient][index_row, 2] = int(dic_patient_label[patient]) # int(Y_label[rowi])
                        patients_total_RP[patient][index_row, 1] = int(data_Time[rowi]-(dic_patient_firstDate[patient] + Observation_months))
                        patients_total_RP[patient][index_row, 0] = int(Patient_num[rowi])

                        dates_patient_total.setdefault(patient, []).append(int(data_Time[rowi]))

                        # the first 6 cols
                        # ID | T - (firstData + Observation_months) | label_final | last_Date - First_Date - Observation_t(censoring time) | MAP result(mean) | 1(what's 1 for?)
                        # first 3 rows is as before

                    else:       ###### data after the observation: to provide labels for different censor times
                        ######observed data for imputation##
                        patients_total[patient][index_row, 6:6 + colums_max - colums_min] += data_array_all[rowi]
                        patients_total[patient][index_row, 6:6 + colums_max - colums_min] = \
                            np.minimum(patients_total[patient][index_row, 6:6 + colums_max - colums_min],
                                       number_maximum)
                        patients_total[patient][index_row, 5] = 1
                        patients_total[patient][index_row, 3] = dic_patient_lastDate[patient] - dic_patient_firstDate[
                            patient] - Observation_months
                        #############the MAP results as the keyfeature2: using the last 20 windows MAP results
                        patients_total[patient][index_row, 4] = np.mean(
                            dic_map_results[patient][-1:])  ###the map results: the last 3 monthely visits(or 1?)
                        # patients_total[patient][index_row, 4] = np.mean
                        patients_total[patient][index_row, 2] = int(dic_patient_label[patient]) #int(Y_label[rowi])
                        patients_total[patient][index_row, 1] = int(
                            data_Time[rowi] - (dic_patient_firstDate[patient] + Observation_months))
                        patients_total[patient][index_row, 0] = int(Patient_num[rowi])
                        
                        #####   Important: RP stands for risk_pridiction, which does not contains the data after the observation time
                        patients_total_RP[patient][index_row, 5] = 0 #### set to be zero since future risk factors cannot be observed
                        patients_total_RP[patient][index_row, 3] = dic_patient_lastDate[patient] - dic_patient_firstDate[
                            patient] - Observation_months
                        #############the MAP results as the keyfeature2: using the last 20 windows MAP results
                        patients_total_RP[patient][index_row, 4] = np.mean(
                            dic_map_results[patient][-1:])  ###the map results: the last 3 monthely visits
                        # patients_total[patient][index_row, 4] = np.mean
                        patients_total_RP[patient][index_row, 2] = int(dic_patient_label[patient]) #int(Y_label[rowi])
                        patients_total_RP[patient][index_row, 1] = int(
                            data_Time[rowi] - (dic_patient_firstDate[patient] + Observation_months))
                        patients_total_RP[patient][index_row, 0] = int(Patient_num[rowi])

                        dates_patient_total.setdefault(patient, []).append(int(data_Time[rowi]))

        #if patient == str(663):
           # print(patients_total[str(663)][:, 6:6 + colums_max - colums_min])

    # print(patients_total[str(663)][:, 6:6 + colums_max - colums_min])



    #################getting kernal weights for each patient

    if not window_h>0:
        print("--------without kernel weights for each patient---------\n")
    else:
        print("--------getting kernel weights for each patient---------\n")
    dic_patient_span={}
    for patienti in patients_total:
        if patienti in dic_patient_lastDate and patienti in dic_patient_firstDate:
            time_span=dic_patient_lastDate[patienti]-dic_patient_firstDate[patienti]-Observation_months   # time_span = censoring_time
            if time_span>0:
                dic_patient_span[patienti]=time_span

    time_spans_all=np.array(list(dic_patient_span.values()))   # an array contains time_span for every patients
    dic_patient_weight={}

    for patienti in dic_patient_span:
        if window_h > 0:
            time_span=dic_patient_span[patienti]
            weights=(1.0/(window_h*math.sqrt(2.0*math.pi)))*np.exp(-np.square((time_span-time_spans_all)/window_h)/2.0)
            dic_patient_weight[patienti]=np.average(weights)
            # patient_weight is the average of it to every other patients?

        else:
            dic_patient_weight[patienti] = 1

    #########################################getting kernal weights for each patient
    weights_sum = np.sum(np.array(list(dic_patient_weight.values())))
    weights_total=len(dic_patient_weight)
    weight_adjust=weights_total*1.0/weights_sum   # what's this?
    ########to make sure the total weights to be data size


    # the key in dic_patient_span satisfies (if (int(patient) in patient_targeted or 'ALL' in patient_targeted) /
    # and str(patient) in dic_map_results <- [map results after observation time] and time-span > 0)
    numbers_total = dic_patient_span.keys()
    # Must be the same when changed to patient_total_vaild

    data_total = []
    data_risk_total = []
    data_total_counts=[]
    patient_num_total=[]
    label_total = []
    date_total=[]
    weight_total = []
    key_feature1_total=[]
    key_feature2_total=[]
    squence_maximu=0
    feature_d=colums_max-colums_min
    # feature_d = total_code

    '''
    key_feature1_total: span_time
    key_feature2_total: MAP_result
    '''

    for patient_i in numbers_total:
        # print('patient:', patients_total[patient_i])
        data_i = patients_total[patient_i]
        data_RP_i = patients_total_RP[patient_i]
        #print ("----data_i.shape: ",data_i.shape)
        rows, cols = data_i.shape # patient_num, date, Y, not included
        data_temp = []
        data_risk_temp = []
        label_total_temp = []
        patient_num_total_temp = []
        date_total_temp = []
        weight_temp = []
        key_feature1_total_temp=[]
        key_feature2_total_temp=[]
        # feature_d=colums_max-colums_min
        # risk_feat_use =
        # feature_risk = index(risk_factors)
        visit_i_valid_last=0
        for visit_i in range(rows):
            if np.sum(data_i[visit_i, 5]) == 1 and len(data_temp)<visit_maximum:  ########all data are included
                visit_i_valid_last=visit_i
                data_visit_i=data_i[visit_i]
                data_visit_i[5]=data_visit_i[5]+int(dic_patient_span[patient_i])
                # print(data_visit_i[-feature_d:])
                data_temp.append(data_visit_i[-feature_d:])
                ############
                data_visit_RP_i=data_RP_i[visit_i]
                data_visit_RP_i[5]=data_visit_RP_i[5]+int(dic_patient_span[patient_i])
                data_risk_temp.append(data_visit_RP_i[-feature_d:])
                
                ############
                key_feature2_total_temp.append(data_i[visit_i, 4])          # MAP result
                key_feature1_total_temp.append(data_i[visit_i, 3])          # span time
                label_total_temp.append(data_i[visit_i, 2])                 # label_final
                date_total_temp.append(data_i[visit_i, 1])                  # T
                patient_num_total_temp.append(data_i[visit_i, 0])           # patient_ID

                if (train_mode=="train" or train_mode=="unlabeled") and patient_i in dic_patient_weight : ###########only train use weights
                    weight_temp.append(weight_adjust*dic_patient_weight[patient_i])
                else:
                    weight_temp.append(1)       #######weights only for MLP binary prediction

        # What are these 2 blocks for?
        if len(data_temp)>squence_maximu:
            squence_maximu=len(data_temp)

        # Padding,to keep the size correspond
        for add_i in range(visit_maximum - len(data_temp)):
            #data_temp.append(np.zeros(shape=(embeddings_d)))
            data_temp.append(np.zeros(shape=(feature_d)))
            data_risk_temp.append(np.zeros(shape=feature_d))
            key_feature2_total_temp.append(data_i[visit_i_valid_last, 4])
            key_feature1_total_temp.append(data_i[visit_i_valid_last, 3])
            label_total_temp.append(data_i[visit_i_valid_last, 2])
            date_total_temp.append(-1)          #########
            patient_num_total_temp.append(data_i[visit_i_valid_last, 0])
            if (train_mode=="train" or train_mode=="unlabeled")  and patient_i in dic_patient_weight:  ###########only train use weights
                weight_temp.append(weight_adjust*dic_patient_weight[patient_i])   #dic_patient_weight[patient_i]
            else:
                weight_temp.append(1)
        #print ("0000data_temp.shape: ",np.array(data_temp).shape)
        #print("0000data_temp.shape: ", np.array(data_temp).shape)
        # print(np.array(data_temp).reshape(-1, feature_d))
        data_total.append(np.array(data_temp).reshape(-1, feature_d))
        data_risk_total.append(np.array(data_risk_temp).reshape(-1, feature_d))
        key_feature2_total.append(key_feature2_total_temp)
        key_feature1_total.append(key_feature1_total_temp)
        label_total.append(label_total_temp)
        date_total.append(date_total_temp)
        patient_num_total.append(patient_num_total_temp)
        weight_total.append(weight_temp)
        # print("------data_total[0]", data_total[0])



    label_total = np.array(label_total)
    label_total = np.expand_dims(label_total, axis=-1)
    weight_total = np.array(weight_total)
    patient_num_total=np.array(patient_num_total)
    data_total = np.array(data_total)
    data_risk_total = np.array(data_risk_total)
    key_feature1_total = np.mean(np.array(key_feature1_total), axis=-1) #/ (np.sum(data_total, axis=(1, 2))+1.0)
    key_feature2_total = np.mean(np.array(key_feature2_total), axis=-1) #/ (np.sum(data_total, axis=(1, 2))+1.0)
    # the average over time-dimension

    key_feature1_total = np.expand_dims(np.array(key_feature1_total),-1)
    key_feature2_total = np.expand_dims(np.array(key_feature2_total), -1)
    key_feature12_total=np.concatenate([key_feature1_total,key_feature2_total],axis=-1)

    # print("----------------------------------------------year_window: ", year_window)
    print("------train_mode:  ", train_mode, "\n------patient_total_valid: ",
          len(patient_total_valid))
    print("------data_total: ", np.array(data_total).shape)
    print("------data_risk_total: ", np.array(data_risk_total).shape)
    print("------label_total.shape: ", label_total.shape)
    print("------weight_total.shape: ", weight_total.shape)
    print("------patient_num_total.shape: ", patient_num_total.shape)
    print("------date_total.shape: ", np.array(date_total).shape)
    print("------key_feature12_total.shape: ", key_feature12_total.shape)
    print("------key_feature1_total.shape: ", key_feature1_total.shape)
    print("------key_feature1_total[0](span_time): ",key_feature1_total[0])

    # print("------data_total[0]",np.array(data_total)[0])
    # print("------data_risk_total[0]: ", np.array(data_risk_total)[0])
    # print("np.array(label_total)[0]: ", label_total[0])
    # print("np.array(weight_total)[0]: ", weight_total[0])

    filename_save =  filename + "_"+str(train_mode)+"_kernel_"+str(Decimal(kernal_flag))+"_weight_unlabel_"+str(round(Decimal(weight_unlabel),5))+".pkl"
    save_name = str(filename_save).split(".csv")[0] + str(filename_save).split(".csv")[1]
    with open(dirr_save + "PKL_files/" + save_name, 'wb') as fid:
        pickle.dump((np.array(data_total),np.array(data_risk_total),np.array(label_total),
                     np.array(patient_num_total),np.array(date_total),
                     np.array(weight_total),
                     key_feature1_total,key_feature2_total,key_feature12_total),fid, protocol=4)
    

        # pickle.dump((np.array(data_total),np.array(label_total),
        #              np.array(patient_num_total),np.array(date_total),
        #              np.array(weight_total),np.array(data_embedding_all),
        #              key_feature1_total,key_feature2_total,key_feature12_total),fid, protocol=4)

