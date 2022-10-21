#!/bin/bash
#=======================================================================
# get_training_preduction_data.sh
#   This script downloads the necessary ERA5 training and prediction data
#   needed for training and verifying the model
#   This can take a long time for 40 years of data even with multiple cores this
#   script can take days
#=======================================================================
set -e

# Directory settings
cd ../

HYBRID=`pwd`
SCRIPTS=$HYBRID/scripts

# Source experiment configuration and time increment function
source $HYBRID/config.sh

DATADIR=$ROOT_DIR/ERA_5/

cd $SCRIPTS

if [ ! -d "$DATADIR" ]; then
  echo "making $DATADIR "
  mkdir $DATADIR
fi

#Calculate final year of prediction based off config.sh
FINAL_YEAR="$(($PSYYYY+$predict_len/8760))"
if [ $FINAL_YEAR -gt 2021 ]; then
  $FINAL_YEAR=2021
fi

python download_pressure_grid.py $IYYYY $FINAL_YEAR -path=$DATADIR --np=$NUM_DATA_PROCS --vars temperature u_component_of_wind v_component_of_wind specific_humidity --plev 20 30 50 70 100 125 200 225 250 300 350 400 450 500 550 600 650 700 750 775 800 825 850 875 900 925 950 975 --grid 3.0 3.0

python download_single_level.py $IYYYY $FINAL_YEAR -path=$DATADIR --np=$NUM_DATA_PROCS --vars sea_surface_temperature surface_pressure total_precipitation --grid 3.0 3.0



