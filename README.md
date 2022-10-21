# SPEEDY-ML

Fortran code for a hybrid model that combines an atmospheric general circulation model (SPEEDY) and a reservoir computing-based machine learning algorithm.

## Installation Instructions
Download the SPEEDY-ML model with:
<pre><code> $ git clone git@github.com:Arcomano1234/SPEEDY-ML.git
</code></pre>

## Prerequisites
SPEEDY-ML requires a MPI built NetCDF library, [ARPACK-NG](https://github.com/opencollab/arpack-ng), and Intel oneapi. 
This code has been tested with both intel 2020 and gfortran 2021 with openmpi. 

## Quick Setup
1. Fill in the necessary information in config.sh

2. If you want to train your model use scripts/get_training_prediction_data.sh. If you want to use an already trained model skip to step 5.

3. Run scripts/regrid_data.sh

4. After all of the necessary training data is downloaded and regridded to the SPEEDY horizontal and vertical grid run scripts/train_model.sh . Depending on the number of processors this can take anywhere from 40 minutes to a day.

5. Run scripts/run_trained_model.sh . This can be used with a small number of CPUs with the only computational requiring being 32 GB of memory alotted to the program.

6. After the predictions are done, a number of python scripts are provided to analyze the forecasts or climate simulations 

## Training Data
To download training data run scripts/get_training_prediction_data.sh . This script downloads ERA5 data using parallel threads with each thread downloading a year of data (so max threads would be number of years of data). 

## Trained Model 
Trained model weights for a hybrid atmospheric model and machine learning only ocean model are provided at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7222831.svg)](https://doi.org/10.5281/zenodo.7222831). The get_trained_coupled_data.sh will actually download these weights and move them into the appropriate directories.



