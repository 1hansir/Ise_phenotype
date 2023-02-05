from numpy import asarray
from numpy import unique
from numpy import argmax
from tensorflow.keras.datasets.mnist import load_data
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import Conv2D
from tensorflow.keras.layers import MaxPool2D
from tensorflow.keras.layers import Flatten
from tensorflow.keras.layers import Dropout
from sklearn.utils import shuffle
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score, average_precision_score
import numpy as np
import pandas as pd
import csv
import random
import logging,os
import pickle
from collections import deque
from tensorflow.keras import layers,models
import tensorflow as tf
import tensorflow as tf
from tensorflow.keras.layers import Dense, Flatten, Conv2D
from tensorflow.keras import Model
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_auc_score, f1_score,average_precision_score,accuracy_score,precision_score, precision_recall_curve
from tensorflow.keras.models import load_model



tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)


def Attention_train(dirr_train,filename_train,dirr_test, filename_test, dirr_unlabel,filename_unlabel,
                    dirr_save, filename_save,
                  colums_min,colums_max,  max_visits,
            epochs=50,batch_size=256,epoch_show=1,feature_d=161,
                    year_window=1,MLP_only=False,weight_unlabel=0.4,
                    weight_super=1.0,flag_imputation=True,epochs_imputation=15,Rep = 1):

    #flag_imputation = True
    if flag_imputation>0 or flag_imputation==True:
        flag_imputation=True
    else:
        flag_imputation = False
    feature_d=colums_max-colums_min
    feature_d_imput = colums_max - colums_min
    train_epoches=epochs
    batch_size=batch_size
    epoch_show=epoch_show
    max_lengh=max_visits
    codes_total=colums_max-colums_min
    smooth_weight=0.0001
    weights_binary_MLP=1.0
    weights_keyfeature=0.001
    optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)
    optimizer_ft = tf.keras.optimizers.Adam(learning_rate=0.001)
    
    with open(dirr_train + filename_train, 'rb') as fid:
        (data_train, data_risk_train, y_train, patient_num_train,
         date_train, weights_train,
         key_feature1_train,key_feature2_train,key_feature12_train) = pickle.load(fid)
    with open(dirr_test + filename_test, 'rb') as fid:
        (data_test, data_risk_test, y_test, patient_num_test, date_test,
         weights_test,
         key_feature1_test,key_feature2_test,key_feature12_test) = pickle.load(fid)
    with open(dirr_unlabel + filename_unlabel, 'rb') as fid:
        (data_train_unlabel, data_risk_train_unlabel, y_train_unlabel, patient_num_train_unlabel,
         date_train_unlabel, weights_train_unlabel,
         key_feature1_unsuper,key_feature2_unsuper,key_feature12_unsuper) = pickle.load(fid)

    print("--------------------------data_train.shape: ", np.array(data_test).shape)
    print("--------------------------data_risk_train.shape: ", np.array(data_test).shape)
    print("--------------------------y_train.shape: ", np.array(y_train).shape)
    print("--------------------------date_train.shape: ", np.array(date_train).shape)
    print("--------------------------patient_num_train.shape: ", np.array(patient_num_train).shape)
    print("--------------------------weights_train.shape: ", np.array(weights_train).shape)
    print("--------------------------data_test.shape: ", np.array(data_test).shape)
    print("--------------------------data_risk_test.shape: ", np.array(data_test).shape)
    print("--------------------------y_test.shape: ", np.array(y_test).shape)
    print("--------------------------patient_num_test.shape: ", np.array(patient_num_test).shape)
    print("--------------------------weights_test.shape: ", np.array(weights_test).shape)
    print("--------------------------data_train_unlabel.shape: ", np.array(data_train_unlabel).shape)
    print("--------------------------data_risk_train_unlabel.shape: ", np.array(data_risk_train_unlabel).shape)

    data_unsuper=[]
    data_risk_unsuper=[]
    weights_unsuper = []
    data_train_new=[]
    data_risk_train_new=[]
    y_train_new=[]
    patient_num_train_new=[]
    date_train_new=[]
    weights_train_new=[]
    train_num=batch_size*20
    key_features_train=[]
    key_features_unsuper=[]
    for i in range(train_num):
        sample_i=i%(len(data_train))
        key_features_train.append(key_feature12_train[sample_i])
        data_train_new.append(data_train[sample_i])
        data_risk_train_new.append(data_risk_train[sample_i])
        y_train_new.append(y_train[sample_i])
        patient_num_train_new.append(patient_num_train[sample_i])
        date_train_new.append(date_train[sample_i])
        weights_train_new.append(weights_train[sample_i])
        sample_i_unsuper=random.randint(0, len(data_train_unlabel) - 1)
        data_unsuper.append(data_train_unlabel[sample_i_unsuper])
        data_risk_unsuper.append(data_risk_train_unlabel[sample_i_unsuper])
        weights_unsuper.append(weights_train_unlabel[sample_i_unsuper])
        key_features_unsuper.append(key_feature12_unsuper[sample_i_unsuper])

    data_train=np.array(data_train_new)
    data_risk_train=np.array(data_risk_train_new)
    y_train = np.array(y_train_new)
    patient_num_train = np.array(patient_num_train_new)
    date_train = np.array(date_train_new)
    key_features_train = np.array(key_features_train)
    weights_train = np.array(weights_train_new)
    data_unsuper = np.array(data_unsuper)
    weights_unsuper = np.array(weights_unsuper)
    key_features_unsuper = np.array(key_features_unsuper)


    print ("--------------------------data_train.shape: ",np.array(data_train).shape)
    print ("--------------------------data_risk_train.shape: ",np.array(data_risk_train).shape)

    train_loss = tf.keras.metrics.Mean(name='train_loss')
    train_loss_imputation = tf.keras.metrics.Mean(name='train_loss_imputation')
    train_loss_map = tf.keras.metrics.Mean(name='train_loss_map')
    train_loss_classfication = tf.keras.metrics.Mean(name='train_smooth_loss_classfication')

    train_keyfeature_loss = tf.keras.metrics.Mean(name='train_keyfeature_loss')
    train_metric = tf.keras.metrics.AUC(name='train_auc', )
    train_loss_mono = tf.keras.metrics.Mean(name='train_loss_mono' )
    train_loss_smooth = tf.keras.metrics.Mean(name='train_loss_smooth')
    train_loss_constrast = tf.keras.metrics.Mean(name='train_loss_constrast')
    valid_loss = tf.keras.metrics.Mean(name='valid_loss')
    valid_smooth_loss = tf.keras.metrics.Mean(name='valid_smooth_loss')
    valid_metric = tf.keras.metrics.AUC(name='test_auc')
    ds_train = tf.data.Dataset.from_tensor_slices((data_train, data_risk_train,y_train,
                                                   patient_num_train, date_train, weights_train,
                                                   data_unsuper, data_risk_unsuper,weights_unsuper,
                                                   key_features_train,key_features_unsuper)) \
        .shuffle(buffer_size=500).batch(batch_size) \
        .prefetch(tf.data.experimental.AUTOTUNE).cache()
    ds_test = tf.data.Dataset.from_tensor_slices((data_test, data_risk_test, y_test,
                                                  patient_num_test, date_test, weights_test,
                                                  data_test, data_risk_test, weights_test,
                                                  key_feature12_test,key_feature12_test)) \
        .shuffle(buffer_size=500).batch(batch_size) \
        .prefetch(tf.data.experimental.AUTOTUNE).cache()

    def Model_GRUs():
#       embedding_d=embedding_dim
        inputs1 = layers.Input(shape=(max_lengh,feature_d,))
        inputs_unsuper = layers.Input(shape=(max_lengh, feature_d,))
        inputs_risk = layers.Input(shape=(max_lengh,feature_d,))
        inputs_risk_unsuper = layers.Input(shape=(max_lengh, feature_d,))

        inputs_keyfeature = layers.Input(shape=(2,))
        inputs_keyfeature_unsuper = layers.Input(shape=( 2,))
        selection_weight_label=layers.Input(shape=(max_lengh,))
        selection_weight_unlabel=layers.Input(shape=(max_lengh,))



        #############impute the labels using the MLP netowrk overe the glbal embeddin
        ##########MAP prediction as the second key feature
        FCN0_impute = layers.Dense(20, activation=tf.nn.relu, name="inputation/FCN0")
        GRU_imput = tf.keras.layers.GRU(units=80, return_sequences=False, activation=tf.nn.relu, name="inputation/GRU1",recurrent_dropout=0.1)
        #GRU_imput2 = tf.keras.layers.GRU(units=30, return_sequences=True, activation=tf.nn.relu)
        FCN1_impute = layers.Dense(20, activation=tf.nn.relu, name="inputation/FCN1")
        FCN2_impute = layers.Dense(2, activation=tf.nn.tanh, name="inputation/FCN2")
        Predictor_impute = layers.Dense(1, activation=None, name="inputation/Predictor")

        out_imput_label=FCN0_impute(inputs1)
        out_imput_label=GRU_imput(out_imput_label)
        out_imput_unlabel = FCN0_impute(inputs_unsuper)
        out_imput_unlabel=GRU_imput(out_imput_unlabel)

        out_imput_label = tf.nn.dropout(out_imput_label, rate=0.1)
        out_imput_label = FCN1_impute(out_imput_label)
        out_imput_label = tf.nn.dropout(out_imput_label, rate=0.1)
        out_imput_unlabel = tf.nn.dropout(out_imput_unlabel, rate=0.1)
        out_imput_unlabel = FCN1_impute(out_imput_unlabel)
        out_imput_unlabel = tf.nn.dropout(out_imput_unlabel, rate=0.1)

        out_imput_label =FCN2_impute(out_imput_label)
        out_imput_unlabel =FCN2_impute(out_imput_unlabel)

        out_imput_label = tf.concat([inputs_keyfeature[:, 1:], out_imput_label], axis=-1)
        out_imput_unlabel = tf.concat([inputs_keyfeature_unsuper[:, 1:], out_imput_unlabel], axis=-1)
        out_imput_label=Predictor_impute(out_imput_label)
        out_imput_unlabel = Predictor_impute(out_imput_unlabel)
        print("---out_imput_unlabel: ", out_imput_unlabel)
        #############################imputation end##########



        #############using LSTM to learn the representation of the observation periods and do prediction
        FCN_pre = layers.Dense(9, activation=tf.nn.relu)
        GRU_EN = tf.keras.layers.GRU(units=15, return_sequences=False, activation=tf.nn.relu,go_backwards=True)
        GRU_DE1 = tf.keras.layers.GRU(units=25, return_sequences=True, activation=tf.nn.relu)
        GRU_Bidirectional = tf.keras.layers.Bidirectional(tf.keras.layers.GRU(units=20, return_sequences=True,
                 activation=tf.nn.relu),merge_mode='ave')

        GRU_DE2 = tf.keras.layers.GRU(units=30, return_sequences=True, activation=tf.nn.relu,recurrent_dropout=0.1)
        GRU_FCN1 = layers.Dense(7, activation=tf.nn.relu)
        Predictior = layers.Dense(1, activation=None,kernel_regularizer=tf.keras.regularizers.L1L2(0.01))
        Predictior_unlabel = layers.Dense(1, activation=None,kernel_regularizer=tf.keras.regularizers.L1L2(0.01))
############encoding the input sequence

        feature_label=FCN_pre(inputs_risk)  #GRU1(inputs1)
        feature_unlabel=FCN_pre(inputs_risk_unsuper)  # GRU1(inputs_unsuper)
        feature_label = GRU_EN(feature_label)  # GRU1(inputs1)
        feature_unlabel = GRU_EN(feature_unlabel)  # GRU1(inputs_unsuper)

############decoding the  sequence
        feature_label_in=tf.expand_dims(feature_label, 1)
        feature_unlabel_in = tf.expand_dims(feature_unlabel, 1)
        for rowi in range(1, max_lengh):
            temp1 = tf.expand_dims(feature_label, 1)
            feature_label_in = tf.concat([feature_label_in, temp1], axis=1)
            temp2 = tf.expand_dims(feature_unlabel, 1)
            feature_unlabel_in = tf.concat([feature_unlabel_in, temp2], axis=1)
        print (" ---feature_unlabel_in: ", feature_unlabel_in)
        print (" ---feature_label_in: ", feature_label_in)

        out_label=GRU_DE1(feature_label_in)
        #out_label=GRU_DE2(out_label)
        out_label = tf.nn.dropout(out_label, rate=0.1)
        fcn_label = GRU_FCN1(out_label)
        out_unlabel=GRU_DE1(feature_unlabel_in)
        #out_unlabel=GRU_DE2(out_unlabel)
        out_unlabel=tf.nn.dropout(out_unlabel,rate=0.1)
        fcn_unlabel=GRU_FCN1(out_unlabel)

        loss_smooth_label=tf.reduce_sum(tf.reduce_mean(tf.abs(fcn_label[:,0:-1,:]-fcn_label[:,1:,:]),axis=-1))
        loss_smooth_unlabel=tf.reduce_sum(tf.reduce_mean(tf.abs(fcn_unlabel[:,0:-1,:]-fcn_unlabel[:,1:,:]),axis=-1))
        loss_smooth=(loss_smooth_label+loss_smooth_unlabel)*0.5

        #out_label = tf.nn.dropout(fcn_label, rate=0.1)
        out_label=Predictior(out_label)
        #out_unlabel = tf.nn.dropout(fcn_unlabel, rate=0.1)
        out_unlabel=Predictior_unlabel(out_unlabel)

        loss_mono_Y=tf.reduce_sum(tf.maximum(out_label[:,0:-1,:]-out_label[:,1:,:],0))
        loss_mono_Y_unlabel=tf.reduce_sum(tf.maximum(out_unlabel[:,0:-1,:]-out_unlabel[:,1:,:],0))
        loss_mono=(loss_mono_Y+loss_mono_Y_unlabel)*0.5

        out_label = tf.reduce_sum(out_label * tf.expand_dims(selection_weight_label, axis=-1), axis=1)
        out_unlabel = tf.reduce_sum(out_unlabel * tf.expand_dims(selection_weight_unlabel, axis=-1), axis=1)

        fcn_label = tf.reduce_sum(fcn_label * tf.expand_dims(selection_weight_label, axis=-1), axis=1)
        print("---Prediction_label: ", out_label)

        model = models.Model(inputs = [inputs1,inputs_unsuper,
                                       inputs_keyfeature,inputs_keyfeature_unsuper,
                                       selection_weight_label,selection_weight_unlabel,
                                       inputs_risk,inputs_risk_unsuper],
                             outputs = [out_label,out_unlabel,out_imput_label,
                                        out_imput_unlabel,loss_mono,loss_smooth,fcn_label]) #  outputs_sque  outputs  outputs_fused
        return model

    def train_step(model, data_train,labels,weights,
                   data_unsuper,weights_unsuper,smooth_weight_in,
                                keyfeature,keyfeature_unsuper,weight_super,data_risk, data_risk_unsuper):
        weight_impute=weight_super
        threshold_MLP_same=9
        threshold_MLP_differ=10
        with tf.GradientTape(persistent=True) as tape:
            weights = tf.cast(weights, tf.float32)
            weights_unsuper = tf.cast(weights_unsuper, tf.float32)
            #print("weights ori shape: ", weights.shape)
           # print("weights ori : ", weights)
            weights = tf.reduce_mean(weights, axis=1)  #### here the weights only for MLP
            weights_unsuper = tf.reduce_mean(weights_unsuper, axis=1)
            # print (" weights: ",weights[0:20])
            # print(" weights_unsuper: ", weights_unsuper[0:20])
            # print(" weights mean: ", np.mean(weights))
            # print(" weights_unsuper mean: ", np.mean(weights_unsuper))
            ############to randomly change censor time

            batch_num=len(data_train.numpy())
            selection_weight_label=np.zeros((batch_num,max_lengh))
            selection_weight_unlabel = np.zeros((batch_num, max_lengh))
            keyfeature_temp=np.array(keyfeature)
            keyfeature_unsuper_temp = np.array(keyfeature_unsuper)

            for i in range(batch_num):
                selection_weight_label[i, min(int(keyfeature_temp[i,0]),max_lengh-1)]=1
                selection_weight_unlabel[i,  min(int(keyfeature_unsuper_temp[i, 0]),max_lengh-1)] = 1
            selection_weight_label=np.array(selection_weight_label)
            selection_weight_unlabel = np.array(selection_weight_unlabel)

            MLP_prediction, MLP_prediction_unlabel,MLP_prediction_imputed,\
            MLP_prediction_unlabel_imputed,loss_mono,loss_smooth,fcn_labeled= \
                model([data_train,data_unsuper,keyfeature,keyfeature_unsuper,
                       selection_weight_label,selection_weight_unlabel,
                       data_risk, data_risk_unsuper], training=True)

            ########## label loss
            labels_MLP = tf.reduce_sum(labels, axis=1)
            labels_MLP = tf.cast(tf.greater_equal(labels_MLP, 1.0), tf.float32)
            #print(labels_MLP)
            indices = tf.range(start=0, limit=tf.shape(data_train)[0], dtype=tf.int32)
            shuffled_indices = tf.random.shuffle(indices)
            shuffled_fcn_labeled = tf.gather(fcn_labeled, shuffled_indices)
            shuffled_labels_MLP = tf.gather(labels_MLP, shuffled_indices)

            fcn_batch = tf.reduce_sum(tf.abs(fcn_labeled - shuffled_fcn_labeled), -1)
            Y_batch= tf.reduce_sum(abs(labels_MLP - shuffled_labels_MLP), -1)
            distance_same = (1 - Y_batch) * tf.maximum(fcn_batch - threshold_MLP_same, 0.0)
            distance_differ = Y_batch * tf.maximum(threshold_MLP_differ - fcn_batch, 0.0)
            distance_constastive = tf.reduce_mean(distance_same + distance_differ)

            MLP_prediction = tf.nn.sigmoid(MLP_prediction)
            #print(MLP_prediction)
            #print ("-----labels_MLP sum: ",np.sum(labels_MLP.numpy()))
            loss_binary_MLP = tf.keras.losses.binary_crossentropy(labels_MLP, MLP_prediction) *weights
            ############# label patients with MAP loss
            MLP_prediction_imputed = tf.nn.sigmoid(MLP_prediction_imputed)
            loss_binary_imputation = tf.keras.losses.binary_crossentropy(labels_MLP, MLP_prediction_imputed) #*weights
            #############
            loss_mono=loss_mono/batch_num
            loss_smooth=loss_smooth/batch_num
            ##########MAP loss for unlabeled patients
            # print ("keyfeature_unsuper[:,1]: ", np.array(keyfeature_unsuper[:,1]).shape)
            # print("-----unlabeled labels sum: ", np.sum(keyfeature_unsuper[:,1]))
            MLP_prediction_unlabel_imputed = tf.nn.sigmoid(MLP_prediction_unlabel_imputed)
            MLP_prediction_unlabel=tf.nn.sigmoid(MLP_prediction_unlabel)

            if flag_imputation==True:
                loss_binary_MLP_unlabel = tf.keras.losses.binary_crossentropy(MLP_prediction_unlabel_imputed,
                                                                              MLP_prediction_unlabel)  * weights_unsuper
            else:
                loss_binary_MLP_unlabel = tf.keras.losses.binary_crossentropy(keyfeature_unsuper[:,1],
                                                                              MLP_prediction_unlabel) * weights_unsuper
            loss_imputation=loss_binary_imputation
            loss =  loss_binary_MLP*1.0+loss_mono*0.1+loss_binary_MLP_unlabel*0.8#+loss_binary_MLP_unlabel*0.1 # loss_binary_MLP_unlabel*0.2+  +distance_constastive*0.1

        all_variables = model.trainable_variables
        variable_inputation = [v for v in all_variables if  'inputation' in v.name ]  #if  'inputation' in v.name
        variable_risk = [v for v in all_variables if not 'inputation' in v.name]
        gradients_imputation = tape.gradient(loss_imputation, variable_inputation)
        gradients = tape.gradient(loss, variable_risk)

        if flag_imputation==True:
            if weight_super<0.1:   ########only do imputation
                optimizer_ft.apply_gradients(grads_and_vars=zip(gradients_imputation, variable_inputation))
            else:
                optimizer.apply_gradients(grads_and_vars=zip(gradients, variable_risk))
        else:
            optimizer.apply_gradients(grads_and_vars=zip(gradients, variable_risk))

        train_loss.update_state(loss_binary_MLP)
        train_loss_map.update_state(loss_binary_MLP_unlabel)
        train_loss_classfication.update_state(loss)


        train_keyfeature_loss.update_state(loss)
        train_loss_imputation.update_state(loss_imputation)
        train_loss_smooth.update_state(loss_smooth)
        train_loss_constrast.update_state(distance_constastive)
        #train_loss_mono.update_state(loss_mono)
        return loss.numpy(),MLP_prediction.numpy(),labels_MLP.numpy(),loss_mono.numpy()

    def valid_step(model,data_test, labels,weights,data_unsuper,weights_unsuper,
                   smooth_weight,keyfeature, keyfeature_unsuper,data_risk, data_risk_unsuper):
    # def valid_step(model,data_test, labels,weights,data_embedding_test,
    #                data_embedding_test_unsuper,data_unsuper,weights_unsuper,smooth_weight
    #                ,keyfeature, keyfeature_unsuper):

        batch_num = len(data_test.numpy())
        selection_weight_label = np.zeros((batch_num, max_lengh))
        selection_weight_unlabel = np.zeros((batch_num, max_lengh))

        keyfeature_temp = np.array(keyfeature)
        keyfeature_unsuper_temp = np.array(keyfeature_unsuper)
        for i in range(batch_num):
            selection_weight_label[i, min(int(keyfeature_temp[i, 0]),max_lengh-1)] = 1
            selection_weight_unlabel[i, min(int(keyfeature_unsuper_temp[i, 0]),max_lengh-1)] = 1
        selection_weight_label = np.array(selection_weight_label)
        selection_weight_unlabel = np.array(selection_weight_unlabel)

        MLP_prediction, MLP_prediction_unlabel,MLP_prediction_MAP, \
        MLP_prediction_unlabel_MAP,loss_mono,loss_smooth,fcn_label=\
        model([data_test,data_unsuper,keyfeature,
               keyfeature_unsuper,selection_weight_label,selection_weight_unlabel,
               data_risk, data_risk_unsuper])
        # model([data_test,data_embedding_test,data_unsuper,data_embedding_test_unsuper,
        #                keyfeature, keyfeature_unsuper,keyfeature, keyfeature_unsuper])
        MLP_prediction  = tf.nn.sigmoid(MLP_prediction)
        labels_MLP = tf.reduce_sum(labels, axis=1)
        labels_MLP = tf.cast(tf.greater_equal(labels_MLP, 1.0), tf.float32)
        #print(labels_MLP.numpy().shape)
        #print(MLP_prediction.numpy().shape)
        loss_binary_MLP = tf.keras.losses.binary_crossentropy(labels_MLP, MLP_prediction)
        loss = loss_binary_MLP
        valid_loss.update_state(loss)
        return loss.numpy(),MLP_prediction.numpy(),labels_MLP.numpy()

    def train_model(model, ds_train, ds_valid, epochs):
        print ("--------begin training.....")
        epoch_num=-1
        epoch_show=1
        AUC_test_total=[]
        PPV_test_total=[]
        AUC_MLP_test_total=[]
        threshold = 0.5
        smooth_weight_temp=smooth_weight
        while(epoch_num<epochs):
            epoch_num+=1
            if epoch_num<3:
                smooth_weight_temp=smooth_weight/2.0
            else:
                smooth_weight_temp=smooth_weight /2.0+smooth_weight *(epoch_num-1 )/epochs   #  *(epoch_num-1)/epochs
            predictions_train_total=[]
            labels_train_total = []
            weights_train_total = []
            predictions_test_total = []
            labels_test_total = []
            weights_test_total = []
            date_test_total = []
            patient_num_test_total = []
            MLP_prediction_total=[]
            labels_MLP_total=[]
            MLP_prediction_total_test = []
            labels_MLP_total_test = []
            i_number=0
            weight_super=1.0
            loss_mono_total=[]

            if flag_imputation==True:
                if epoch_num< epochs_imputation:
                    weight_super=0.00000001
                    print ('--------------------------imputaion learning----------------- epoch_num: ',epoch_num, " total epoch: ", epochs)
                else:
                    print ('--------------------------imputaion fixed; semi-supervised training----------------- epoch_num: ',epoch_num, " total epoch: ", epochs)
                    weight_super=1.0
            else:
                print ("  --------------------------no imputation and MAP as the psedu-label for unlabeled  -- epoch_num: ",epoch_num, " total epoch: ", epochs)
            if epoch_num>0:
                for data_train, data_risk_train,y_train, patient_num_train, \
                    date_train, weights_train,data_unsuper,data_risk_unsuper,weights_unsuper,\
                        keyfeature_train,keyfeature_unsuper in ds_train:
                    keyfeature_train=np.array(keyfeature_train)

                    #print("keyfeature_train ori: ", keyfeature_train.shape)
                    loss_out,MLP_prediction,labels_MLP,loss_mono=train_step(model,data_train, y_train ,
                                    weights_train,data_unsuper,weights_unsuper,smooth_weight_temp,
                                    keyfeature_train,keyfeature_unsuper,weight_super, data_risk_train,data_risk_unsuper)
                    # loss_out,MLP_prediction,labels_MLP,loss_mono=train_step(model,data_train, y_train ,
                    #                 weights_train,data_embedding_all_temp,
                    # data_embedding_unsuper_temp,data_unsuper,weights_unsuper,smooth_weight_temp,
                    #                 keyfeature_train,keyfeature_unsuper,weight_super)
                    loss_mono_total.append(loss_mono)
                    if i_number == 0:
                        labels_train_total = np.array(y_train)
                        weights_train_total = np.array(weights_train)
                        MLP_prediction_total= np.array(MLP_prediction)
                        labels_MLP_total = np.array(labels_MLP)
                    else:
                        labels_train_total = np.concatenate((labels_train_total, np.array(y_train)), axis=0)
                        weights_train_total = np.concatenate((weights_train_total, np.array(weights_train)), axis=0)
                        MLP_prediction_total = np.concatenate((MLP_prediction_total, np.array(MLP_prediction)), axis=0)
                        labels_MLP_total = np.concatenate((labels_MLP_total, np.array(labels_MLP)), axis=0)
                    i_number += 1
                i_number = 0
                if epoch_num%epoch_show==0 or epoch_num>epochs-3:
                    for data_test, data_risk_test,y_test, patient_num_test, date_test, \
                        weights_test,data_unsuper,data_risk_unsuper,weights_unsuper,\
                        keyfeature_test,keyfeature_test_unsuper in ds_valid:
                    # for data_test, y_test, patient_num_test, date_test, \
                    #     weights_test,data_embedding_test_temp,data_embedding_test_temp_temp,\
                    #     data_unsuper,weights_unsuper,\
                    #     keyfeature_test,keyfeature_test_unsuper in ds_valid:

                        loss_out,MLP_prediction,labels_MLP=valid_step(model,
                                    data_test,y_test, weights_test,data_unsuper,weights_unsuper,
                              smooth_weight_temp,keyfeature_test, keyfeature_test_unsuper,data_risk_test,data_risk_unsuper)
                    # loss_out,MLP_prediction,labels_MLP=valid_step(model,
                    #             data_test,y_test, weights_test,
                    #     data_embedding_test_temp,data_embedding_test_temp_temp,data_unsuper,weights_unsuper,
                    #       smooth_weight_temp,keyfeature_test, keyfeature_test_unsuper)
                        if i_number == 0:
                            labels_test_total = np.array(y_test)
                            weights_test_total = np.array(weights_test)
                            date_test_total = np.array(date_test)
                            patient_num_test_total = np.array(patient_num_test)
                            MLP_prediction_total_test = np.array(MLP_prediction)
                            labels_MLP_total_test = np.array(labels_MLP)
                        else:
                            labels_test_total = np.concatenate((labels_test_total, np.array(y_test)), axis=0)
                            weights_test_total = np.concatenate((weights_test_total, np.array(weights_test)), axis=0)
                            date_test_total = np.concatenate((date_test_total, np.array(date_test)), axis=0)
                            patient_num_test_total = np.concatenate((patient_num_test_total, np.array(patient_num_test)), axis=0)
                            MLP_prediction_total_test = np.concatenate((MLP_prediction_total_test, np.array(MLP_prediction)), axis=0)
                            labels_MLP_total_test = np.concatenate((labels_MLP_total_test, np.array(labels_MLP)),
                                                                    axis=0)

                        i_number += 1
                        print("prediction train[0:10]: ", MLP_prediction_total[0:10,:].tolist())
                        print("prediction test [0:10]: ", MLP_prediction_total_test[0:10,:].tolist())

                if epoch_num % epoch_show == 0 or epoch_num > epochs - 3 :

                    AUC_MLP_train=roc_auc_score(y_true=labels_MLP_total,
                                                 y_score=MLP_prediction_total,
                                                 average='macro')
                    AUC_MLP_test = roc_auc_score(y_true=labels_MLP_total_test,
                                                  y_score=MLP_prediction_total_test,
                                                  average='macro')

                    AUC_MLP_test_total.append(AUC_MLP_test)

                    MLP_prediction_total_binary = np.array(np.greater(MLP_prediction_total, threshold), dtype=np.int32)
                    MLP_prediction_total_test_binary = np.array(np.greater(MLP_prediction_total_test, threshold), dtype=np.int32)
                    tn, fp, fn, tp = confusion_matrix(labels_MLP_total_test, MLP_prediction_total_test_binary).ravel()
                    specificity_test_MLP = tn / (tn + fp)
                    sensitivity_test_MLP = tp / (tp + fn)
                    PPV_test_MLP = tp / (tp + fp)
                    NPV_test_MLP = tn / (tn + fn)
                    # Fall out or false positive rate
                    FPR = fp / (fp + tn)
                    # False negative rate
                    FNR = fn / (tp + fn)
                    acc_train_MLP=accuracy_score(labels_MLP_total,MLP_prediction_total_binary)
                    acc_test_MLP=accuracy_score(labels_MLP_total_test,MLP_prediction_total_test_binary)
                    logs = '--- epoch {}  loss_labeled: {}, loss_unlabeled: {},  loss_mono: {},  loss_imputation: {}, ' \
                           'AUC_train: {},acc_train:{}, AUC_test:{}  acc_test:{}  PPV_test:{},' \
                           'specif_test_:{},' \
                           'loss_smooth:{},' \
                           'loss_contrast:{}'
                    tf.print(tf.strings.format(logs, (epoch_num,train_loss.result(),train_loss_map.result(),np.mean(loss_mono_total),
                                                      train_loss_imputation.result(),AUC_MLP_train,
                                                          acc_train_MLP,AUC_MLP_test,acc_test_MLP,PPV_test_MLP,
                                         specificity_test_MLP, train_loss_smooth.result(),train_loss_constrast.result())))
                    if epoch_num == epochs - 2 :
                        if True:
                            print("---------saving--- ", dirr_save, ": ", filename_save)
                            #print ("MLP_prediction_total_test: ",MLP_prediction_total_test)
                            print ("len(MLP_prediction_total_test): ",np.array(MLP_prediction_total_test).shape)
                            print("len(labels_MLP_total_test): ", np.array(labels_MLP_total_test).shape)
                            dataframe = pd.DataFrame({'Predition':np.array(MLP_prediction_total_test).reshape(len(MLP_prediction_total_test),),
                                                      'Y': np.array((labels_MLP_total_test)).reshape(len(labels_MLP_total_test),)})

                            dirr_temp = dirr_save + "/Binary" + filename_save + "/Rep_"+str(Rep)+"/"
                            if not os.path.exists(dirr_temp):
                                os.makedirs(dirr_temp)  # 创建多级目录
                            dataframe.to_csv(dirr_temp+"Cstat_results.csv", index=True, sep=',')
                            print("---------saving ends successfully--- ", dirr_save, ": ", filename_save)


                    if epoch_num == epochs-2 :
                        if True:
                            f = open(dirr_save + filename_save + "_Attention_metric_binaryMLP.txt", 'a')
                            f.write("AUC_train_MLP :%4f , ACC_train_MLP:%4f, ,FPR:%4f, FNR:%4f, NPV_MLP:%4f "
                                    " ,AUC_test_MLP:%4f, ACC_test_MLP:%4f "
                                    ", Speci_MLP:%4f ,Sensitivity_MLP:%4f, PPV_MLP:%4f"
                                     % (
                                AUC_MLP_train, acc_train_MLP,FPR,FNR, NPV_test_MLP,AUC_MLP_test,
                                acc_test_MLP,specificity_test_MLP,
                                sensitivity_test_MLP,
                                PPV_test_MLP))
                            f.write("\r")
                            f.close()

            train_loss.reset_states()
            train_loss_map.reset_states()
            valid_loss.reset_states()
            train_metric.reset_states()
            valid_metric.reset_states()
            train_loss_mono.reset_states()
            train_loss_imputation.reset_states()
            train_loss_smooth.reset_states()
            train_loss_constrast.reset_states()


        flag_longitudinal=True
                    ####################
        print ("--------------------------------computing patient longitudinal prediction results---------------")
        filename_save_longitudinal=str(filename_save).split(".csv")[0]
        dirr_save_patient=dirr_save+filename_save_longitudinal+"_patient_longitudinal_predictions/Rep_"+str(Rep)+"/"
        print("--------------savging dirr : ", dirr_save_patient)

        if not os.path.exists(dirr_save_patient):
            os.makedirs(dirr_save_patient)   # 创建多级目录
        if True:
            if True:
                if True:
                    if flag_longitudinal==True:
                        for data_test, data_risk_test,y_test, patient_num_test, date_test, \
                            weights_test, data_unsuper, data_risk_unsuper,weights_unsuper, \
                            keyfeature_test, keyfeature_test_unsuper in ds_valid:
                            ############predicting future 25 years
                            censor_T_max=max_visits
                            y_test_value=np.array(y_test.numpy())
                            print("y_test_value shape before squeeze: ", y_test_value.shape)
                            y_test_value=np.squeeze(y_test_value,axis=-1)
                            date_test_value = np.array(date_test.numpy())

                            for censor_t in range(0,censor_T_max):
                                if censor_t%20==0:
                                    print ("censor_T_max: ",censor_T_max, "months;  current censor_t: ",censor_t)
                                keyfeature_test_i = np.array(keyfeature_test)
                                keyfeature_test_temp=np.ones(shape=np.array(keyfeature_test[:,0]).shape)*censor_t
                                keyfeature_test_i[:,0]=keyfeature_test_temp
                                loss_out, MLP_prediction, labels_MLP = valid_step(model,
                                                                              data_test, y_test, weights_test,
                                                                              data_unsuper, weights_unsuper,
                                                                              smooth_weight_temp, keyfeature_test_i,
                                                                              keyfeature_test_unsuper,data_risk_test,data_risk_unsuper)

                                for rowi in range(len(MLP_prediction)):
                                    visit_Times=list(date_test_value[rowi,:])
                                    visit_labels = list(y_test_value[rowi,:])
                                    if censor_t in visit_Times:
                                        index=visit_Times.index(censor_t)
                                        label_rowi=visit_labels[index]
                                    else:
                                        label_rowi="NA"

                                    savename=str(int(patient_num_test[rowi,0]))+".csv"
                                    # print ("savefile: ",dirr_save_patient)

                                    with open(dirr_save_patient + savename, 'a+', newline='') as f:
                                        writer = csv.writer(f)
                                        writer.writerow([label_rowi, float(MLP_prediction[rowi][0]),censor_t])


        print (" training and test end.....")

    model = Model_GRUs()
    savename_model=dirr_save + filename_save+".h5"
    flag_load=False
    if os.path.exists(savename_model) and flag_load==True:
        print ("---------------------------------------------loadding saved model....................")
        model = load_model(savename_model)
    train_model(model, ds_train, ds_test, epochs=train_epoches)
    if not os.path.exists(savename_model):
        print ("---------------------------------------------saving model....................")
        model.save(savename_model)
