#!/bin/bash

#Directory to hold all of the training data, trained weights and prediction data
#Even at very coarse resoluton (e.g. 3.75 degree)
#there will be over a Terabyte of data stored here 
ROOT_DIR=/scratch/user/troyarcomano/speedy-ml-test #'/path/to/storage'


###Related to Training the models####

#Whether you are using an pre-existing trained hybrid model
TRAINED=1

#Number of cores to get training data
NUM_DATA_PROCS=6

#Training Date Start
IYYYY=1990
IMM=01
IDD=01
IHH=00

#Training date end
FYYYY=2000
FMM=01
FDD=01
FHH=00

###Related to Predictions

#Prediction Start Date
PSYYYY=2000
PSMM=01
PSDD=02
PSHH=00

#Number of predictions
NUM_PREDS=1

#
SYNC_LENGTH=72

#Prediction Length In hours
predict_len=8760*10



