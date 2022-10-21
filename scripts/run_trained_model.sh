#!/bin/bash
#=======================================================================
# run_trained_model.sh
#   This script runs a trained coupled model. First it reads a controller filer
#   makes the necessary changes to the Fortran code then used make to create the 
#   executable file. It creates a tmp directory in the ROOT_DIR and runs it there
#=======================================================================
set -e

# Directory settings
cd ../

HYBRID=`pwd`
SCRIPTS=$HYBRID/scripts

# Source experiment configuration and time increment function
source $HYBRID/config.sh

DATADIR=$ROOT_DIR/ML_SPEEDY_WEIGHTS

TMPDIR=$ROOT_DIR/TMP

MODELDIR=$HYBRID/src

echo "Deleting old output directory"
rm -rf $TMPDIR

if [ ! -d "$TMPDIR" ]; then
  echo "making $TMPDIR"
  mkdir $TMPDIR
fi

echo "Copying model over"
cp $MODELDIR/* $TMPDIR

cd $TMPDIR

echo "make clean"
make clean

sed -i -r -e "s/n_ens = [0-9]+/n_ens = $n_ens/" $TMPDIR/mod_reservoir.f90
sed -i -r -e "s/relax_alpha = [0-9]+.[0-9]+d0/relax_alpha = $rtpp/" $TMPDIR/mod_reservoir.f90
sed -i -r -e "s/sigma_obs=[0-9]+.[0-9]+d[0-9]+/sigma_obs=$hor_loc/" $TMPDIR/mod_reservoir.f90
sed -i -r -e "s/sigma_obsv=[0-9]+.[0-9]+d[0-9]+/sigma_obsv=$ver_loc/" $TMPDIR/mod_reservoir.f90
sed -i -r -e "s/cov_infl_mul = -*[0-9]+.[0-9]+d0/cov_infl_mul = $cov_infl/" $TMPDIR/mod_reservoir.f90

make imp.exe COMPILER=gcc2021 



