### Prerequisites ###
#Make account to download ERA 5 data
#
#Follow instructions at https://cds.climate.copernicus.eu/api-how-to
#
#Make sure CDS API url and key are installed to $HOME/.cdsapirc

#See (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview) for all variable, time dimensions, and grid

###EXAMPLES###
#python download_single_level.py 2019 2020 -path=/skydata2/troyarcomano/ERA_5/ --np=2 --vars sea_surface_temperature surface_pressure total_precipitation --grid 3.0 3.0
#downloads SST, surface pressure and 1 hour total preceip for 2019 and 2020 at a 3 degree by 3 degree grid 

import cdsapi
from calendar import monthrange
from multiprocessing import Pool, cpu_count
from functools import partial
from os import path
import logging
import argparse
import numpy as np

c = cdsapi.Client()

def month_downloader(client,month,year,variables,grid,directory):
    """
    Function that downloads hourly pressure level ERA5 reanalysis data to a netCDF file

    Input: client          : cdsapi client
           month           : integer, month to download
           year            : integer, year to download
           variables       : list of strings, variables wished to download
           grid            : list of reals, resolution of the data [x,y] in degrees
           directory       : string, directory to write the data to

    Output: None but does write file to directory
    """

    days = monthrange(year,month)[1]
    daysarray = np.arange(1,days+1)

    nparraytostr = lambda daysarray:'%02d'%daysarray
    dayslist = list(map(nparraytostr,daysarray))

    #reanalysis-era5-pressure-levels-preliminary-back-extension
    #reanalysis-era5-pressure-levels
    client.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type':'reanalysis',
            'variable':variables,
            'year':f'{year}',
            'month':f'{month:01}',
            'day':dayslist,
            'time':[
                 '00:00','01:00','02:00',
                 '03:00','04:00','05:00',
                 '06:00','07:00','08:00',
                 '09:00','10:00','11:00',
                 '12:00','13:00','14:00',
                 '15:00','16:00','17:00',
                 '18:00','19:00','20:00',
                 '21:00','22:00','23:00'
             ],
             'format':'netcdf',
             'grid'  : grid
        },
        f'{directory}era_5_m{month:01}_y{year}_single_level_data.nc')

def year_loop(capi,args,year):
    months = np.arange(1,13)
    for i in months:
        month_downloader(capi,i,year,args.vars,args.grid,args.path)

###Main part of the code
parser = argparse.ArgumentParser(description="Download Hourly ERA5 Pressure Level Variables")
parser.add_argument("startyear",     type=int, help="Starting year to download")
parser.add_argument("endyear",       type=int, help="Ending year to download")
parser.add_argument("-path",         type=str, default="",metavar='path',help="Path to write the files to")
parser.add_argument("--np",          type=int, default=1, metavar='num_procs',help='Number of processors to use to download. Np should be <= number of years to download (default=1)')
parser.add_argument("--vars",        type=str, metavar='var', nargs='+' , default='surface_pressure',help='Variable or variables to download (default is temperature')
parser.add_argument("--grid",        type=str, nargs=2 , default=[1.0,1.0], metavar=('gridx', 'gridy'), help='Resolution of ERA5 data in the x direction in degrees (default=1.0)')

args = parser.parse_args()

years = np.arange(args.startyear,args.endyear+1)

func = partial(year_loop,c,args)
pool = Pool(args.np)
pool.map(func,years)
pool.close()
 
