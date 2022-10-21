#!/bin/bash
#=======================================================================
# get_trained_couple_data.sh
#   This script downloads the necessary machine learning weights and ERA5 initial conditions 
#   needed for prediction 
#   This can take a long time on the order of a few hours
#=======================================================================
set -e

echo "Starting get_trained_couple_data.sh"

# Directory settings
cd ../

HYBRID=`pwd`
SCRIPTS=$HYBRID/scripts

# Source experiment configuration and time increment function
source $HYBRID/config.sh

DATADIR=$ROOT_DIR/ML_SPEEDY_WEIGHTS/

cd $SCRIPTS

if [ ! -d "$DATADIR" ]; then
  echo "making $DATADIR "
  mkdir $DATADIR
fi

#Got to data directory and start downloading the weights
cd $DATADIR

wget -O coupled_model_weights.tar.gz https://zenodo.org/record/7222831/files/coupled_model_weights.tar.gz?download=1

#Untar downloaded data and then delete tar file
tar -xf coupled_model_weights.tar.gz

rm coupled_model_weights.tar.gz


### Get initial conditions 
DATADIR=$ROOT_DIR/ERA_5/2007

if [ ! -d "$DATADIR" ]; then
  echo "making $DATADIR "
  mkdir -p $DATADIR
fi

cd $DATADIR

#Use wget to download weights from a remote server
wget -O regridded_era5.tar.gz https://zenodo.org/record/7222831/files/regridded_era5.tar.gz?download=1

#Untar downloaded data and then delete tar file
tar -xf regridded_era5.tar.gz

rm regridded_era5.tar.gz

echo "Finished downloading and extracting weight files and ICs files"




