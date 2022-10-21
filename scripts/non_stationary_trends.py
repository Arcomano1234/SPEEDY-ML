import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from netCDF4 import Dataset
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.geoaxes import GeoAxes
import xarray as xr
import glob
from datetime import datetime, timedelta
from numba import jit
import calendar
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy.ndimage.filters import uniform_filter1d

@jit(nopython=True,fastmath=True)
def lin_interp(var,ps,target_pressures):
    speedy_sigma = np.array([0.025, 0.095, 0.20, 0.34, 0.51, 0.685, 0.835, 0.95])

    #var = np.asarray(var)

    ps = np.exp(ps) * 1000.0

    ygrid = np.shape(var)[2]
    xgrid = np.shape(var)[3]

    time_length = np.shape(ps)[0]

    var_pressure = np.zeros((time_length,len(speedy_sigma),ygrid,xgrid))

    for t in range(time_length):
        for i in range(len(speedy_sigma)):
            var_pressure[t,i,:,:] = speedy_sigma[i] * ps[t,:,:]


    regridded_data = np.zeros((time_length,len(target_pressures),ygrid,xgrid))

    for t in range(time_length):
        for i in range(ygrid):
            for j in range(xgrid):
                regridded_data[t,:,i,j] = np.interp(target_pressures,var_pressure[t,:,i,j],var[t,:,i,j])

    return regridded_data

def make_ds_time_dim(ds,timestep,startdate):
    begin_year_str = startdate.strftime("%Y-%m-%d")

    attrs = {"units": f"hours since {begin_year_str} "}

    ds = ds.assign_coords({"Timestep": ("Timestep", ds.Timestep.values*timestep, attrs)})
    ds = xr.decode_cf(ds)

    return ds

startdate = datetime(2000,1,1,0)
enddate = datetime(2011,1,1,0)

ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_2kbias_10_year_then_platue_speedy_atmo_onlytrial_12_29_2006_00.nc')

ds_hybrid = make_ds_time_dim(ds_hybrid,6,startdate)

level = 7

da = ds_hybrid.assign_coords(year_month=ds_hybrid.Timestep.dt.strftime("%Y-%m"))
result = da.groupby("Timestep.month") - da.groupby("Timestep.month").mean("Timestep")

temp = result.sel(Sigma_Level=level)['Temperature'][:].values


average_temp = np.mean(temp,axis=(1,2))
print(np.shape(average_temp))

time = np.arange(1,np.shape(average_temp)[0]+1)
time = time/1460
plt.plot(time,uniform_filter1d(average_temp,size=1460))

ticks = np.arange(1,11)
plt.xticks(ticks)

plt.show()

