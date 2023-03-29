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
from datetime import datetime as dt
from datetime import timedelta 
from numba import jit
import calendar
from mpl_toolkits.axes_grid1 import AxesGrid

@jit()
def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def get_truth_data_ps(startdate,enddate,x_slice,y_slice):
    path = '/scratch/user/troyarcomano/ERA_5/'

    hours = enddate - startdate
    hours = hours.days * 4 + 1

    current_year = startdate.year
    end_year = enddate.year

    start_index = 0
    i = 0
    while current_year <=end_year:
      file_path = f"{path}/{current_year}/era_5_y{current_year}_regridded_mpi.nc"
      hours_in_year = dt(current_year,12,31,23) - dt(current_year,1,1,0)
      time_slice = slice(0,hours_in_year.days*24 + 1,6)
      temp_ds = xr.open_dataset(file_path)
      data_ps = temp_ds['logp'].sel(Timestep=time_slice,Lon=x_slice,Lat=y_slice)
      if start_index == 0:
         truth_ps = np.zeros((hours,np.shape(data_ps)[1],np.shape(data_ps)[2]))
      time_length = np.shape(data_ps)[0]
      truth_ps[start_index:start_index+time_length,:,:] = data_ps.values

      print('start_index,start_index+time_length',start_index,start_index+time_length)


      current_year += 1
      i += 1
      start_index = start_index + time_length
  
    return truth_ps  

def get_obs_era5_timeseries(startdate,enddate,timestep,lat_slice,lon_slice,var):
    start_year = startdate.year
    end_year = enddate.year

    currentdate = startdate
    counter = 0
    while currentdate.year <= enddate.year:
        print(currentdate.year)
        try:
           ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{currentdate.year}/era_5_y{currentdate.year}_regridded_mpi.nc')
        except:
           ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{currentdate.year}/era_5_y{currentdate.year}_regridded_mpi_fixed_var.nc')

        begin_year = dt(currentdate.year,1,1,0)
        begin_year_str = begin_year.strftime("%Y-%m-%d")
        attrs = {"units": f"hours since {begin_year_str} "}
        ds_era = ds_era.assign_coords({"Timestep": ("Timestep", ds_era.Timestep.values, attrs)})
        ds_era = xr.decode_cf(ds_era)

        ds_era = ds_era[var].sel(Lat=lat_slice,Lon=lon_slice)
        if start_year == currentdate.year:
           ds_merged = ds_era
        else:
           ds_merged = xr.merge([ds_merged,ds_era])

        currentdate = currentdate + timedelta(hours=ds_era.sizes['Timestep'])

    time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"),timestep)
    return  ds_merged.sel(Timestep=time_slice)

@jit()
def latituded_weighted_rmse_3d(data,lats):

    weights = np.cos(np.deg2rad(lats))

    weights2d = np.zeros(np.shape(data))

    weights2d = np.tile(weights,(96,1))
    weights2d = np.transpose(weights2d)
    masked = np.ma.MaskedArray(data, mask=np.isnan(data))
    weighted_average = np.ma.average(masked,weights=weights2d)
    return weighted_average

#@jit(nopython=True,fastmath=True)
def time_series_atmo_total_weight(logp,lats):
    p = np.exp(logp) * 1000.0 * 100.0
    
    length = np.shape(logp)[0]
    total_weight = np.zeros((length))
    for i in range(length):
        total_weight[i] = latituded_weighted_rmse_3d(p[i,:,:],lats)

    return total_weight

def time_series_atmo_total_weight_scalar(var,lats):
    length = np.shape(var)[0]
    total_weight = np.zeros((length))
    for i in range(length):
        total_weight[i] = latituded_weighted_rmse_3d(var[i,:,:],lats)

    return total_weight

plotdir = '/home/troyarcomano/FortranReservoir/hybridspeedy/plots/'

#ds_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.1_beta_model_1.0_prior_1.0_overlap1_era5_6hrtimestep_19years_climate_simulation_t_ends_2000_tisrtrial_12_30_1999_00.nc')
ds_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_sigma0.5_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_multi_gaussian_noise_newest_version_32_processors_root_ssttrial_12_29_2006_00.nc')

lats = ds_hybrid_sim['Lat'][:]

region_slice = slice(-90,90)
lon_slice = slice(0,360)

hybrid_logp = ds_hybrid_sim['logp'].sel(Lon=lon_slice,Lat=region_slice)
ds_hybrid_sim.close()
hybrid_weight = time_series_atmo_total_weight(hybrid_logp.values,lats.values) # (np.exp(hybrid_logp) * 1000.0 * 5.100644719 * 10**14)/9.8
hybrid_weight = (hybrid_weight * 5.100644719 * 10**14)/9.8


ds_speedy_sim = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/climo_speedy/speedy_era_start_climo_sst_decade_sim12_31_1999_00.nc')
speedy_logp = ds_speedy_sim['logp'].sel(Lon=lon_slice,Lat=region_slice)
ds_speedy_sim.close()
speedy_weight = time_series_atmo_total_weight(speedy_logp.values,lats.values) # (np.exp(hybrid_logp) * 1000.0 * 5.100644719 * 10**14)/9.8
speedy_weight = (speedy_weight * 5.100644719 * 10**14)/9.8

ds_climo = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/regridded_era_climatology2001_2011.nc')
climo_logp = ds_climo['logp'].sel(Lon=lon_slice,Lat=region_slice)
ds_climo.close()
climo_weight = time_series_atmo_total_weight(climo_logp.values,lats.values) # (np.exp(hybrid_logp) * 1000.0 * 5.100644719 * 10**14)/9.8
climo_weight = (climo_weight * 5.100644719 * 10**14)/9.8
climo_weight = np.average(climo_weight)
total_climo_weight = climo_weight

startdate = dt(2000,1,1,0)
enddate = dt(2010,12,31,1)
ps_time = get_truth_data_ps(startdate,enddate,lon_slice,region_slice)
era_weight = time_series_atmo_total_weight(ps_time,lats.values)
era_weight = (era_weight * 5.100644719 * 10**14)/9.8

print(np.shape(era_weight))
print(np.shape(speedy_weight))
time_speedy = np.arange(1,np.shape(speedy_weight)[0]+1)
time_speedy = time_speedy / (4 * 365.25)

time_hybrid = np.arange(1,np.shape(hybrid_weight)[0]+1)
time_hybrid = time_hybrid / (4 * 365.25)

time_era = np.arange(1,np.shape(era_weight)[0]+1)
time_era = time_era / (4 * 365.25)

z = np.polyfit(time_hybrid,hybrid_weight,1)
p = np.poly1d(z)
print('z total mass',z,time_era[0:10])
trend_line = p(time_hybrid)

plt.rc('font', family='serif')

plt.plot(time_hybrid[0:-49],moving_average(hybrid_weight,50),color='#377eb8',ls='-',linewidth=2.0,label='Hybrid')
plt.plot(time_hybrid,trend_line,color='k',linewidth=2.0,ls='--',label='Hybrid Trend')
plt.axhline(y=np.average(hybrid_weight),color='#377eb8',ls='--',linewidth=2.0)

plt.plot(time_speedy[0:-49],moving_average(speedy_weight,50),color='#4daf4a',ls='-',linewidth=2.0,label='SPEEDY')
plt.axhline(y=np.average(speedy_weight),color='#4daf4a',ls='--',linewidth=2.0)

plt.plot(time_era[0:-149],moving_average(era_weight[0:-100],50),color='#e41a1c',ls='-',linewidth=2.0,label='ERA 5')
plt.axhline(y=climo_weight,color='#e41a1c',ls='--',linewidth=2.0)

plt.legend()
plt.ylabel('Total Atmospheric Mass (kg)',fontsize=16)
plt.xlabel('Years Into Simulation',fontsize=16)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.ylim(np.min(hybrid_weight),np.max(hybrid_weight))
plt.xlim(0,np.max(time_hybrid[0:-149]))

plt.show()

####Moisture conservation######
var = 'Specific-Humidity'
hybrid_sh = ds_hybrid_sim[var].sel(Lon=lon_slice,Lat=region_slice)
print(np.shape(hybrid_sh))
hybrid_sh = np.sum(hybrid_sh.values,axis=(1))
ds_hybrid_sim.close()
print(np.shape(hybrid_sh))
hybrid_weight = np.average(hybrid_sh,axis=(1,2))#time_series_atmo_total_weight_scalar(hybrid_sh,lats.values) # (np.exp(hybrid_logp) * 1000.0 * 5.100644719 * 10**14)/9.8
hybrid_weight = hybrid_weight * total_climo_weight * 0.001 / 9.81


ds_speedy_sim = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/climo_speedy/speedy_era_start_climo_sst_decade_sim12_31_1999_00.nc')
speedy_sh = ds_speedy_sim['Specific_Humidity'].sel(Lon=lon_slice,Lat=region_slice)
speedy_sh = np.sum(speedy_sh.values,axis=(1))
ds_speedy_sim.close()
speedy_weight = np.average(speedy_sh,axis=(1,2))#time_series_atmo_total_weight_scalar(speedy_sh,lats.values) # (np.exp(hybrid_logp) * 1000.0 * 5.100644719 * 10**14)/9.8
speedy_weight = speedy_weight * total_climo_weight * 0.001 / 9.81

ds_climo = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/regridded_era_climatology2001_2011.nc')
climo_sh = ds_climo['Specific_Humidity'].sel(Lon=lon_slice,Lat=region_slice)
climo_sh = np.sum(climo_sh.values,axis=(1))
ds_climo.close()
climo_weight = np.average(climo_sh,axis=(1,2))#time_series_atmo_total_weight_scalar(climo_sh,lats.values) # (np.exp(hybrid_logp) * 1000.0 * 5.100644719 * 10**14)/9.8
climo_weight = np.average(climo_weight)
climo_weight = climo_weight * total_climo_weight * 0.001 / 9.81


startdate = dt(2000,1,1,0)
enddate = dt(2010,12,31,1)
sh_time = get_obs_era5_timeseries(startdate,enddate,6,region_slice,lon_slice,'Specific_Humidity')
sh_time = sh_time['Specific_Humidity'].values #*
sh_time = np.sum(sh_time*1000,axis=(1))
era_weight = np.average(sh_time,axis=(1,2))#time_series_atmo_total_weight_scalar(sh_time,lats.values)
era_weight = era_weight * total_climo_weight * 0.001 / 9.81

print(np.shape(era_weight))
print(np.shape(speedy_weight))
time_speedy = np.arange(1,np.shape(speedy_weight)[0]+1)
time_speedy = time_speedy / (4 * 365.25)

time_hybrid = np.arange(1,np.shape(hybrid_weight)[0]+1)
time_hybrid = time_hybrid / (4 * 365.25)

time_era = np.arange(1,np.shape(era_weight)[0]+1)
time_era = time_era / (4 * 365.25)

z = np.polyfit(time_hybrid,hybrid_weight,1)
p = np.poly1d(z)
print('z water',z)
trend_line = p(time_hybrid)

plt.plot(time_hybrid[0:-49],moving_average(hybrid_weight,50),color='#377eb8',ls='-',linewidth=2.0,label='Hybrid')
plt.plot(time_hybrid,trend_line,color='k',linewidth=1.0,ls='--',label='Hybrid Trend')
plt.axhline(y=np.average(hybrid_weight),color='#377eb8',ls='--',linewidth=2.0)

#plt.plot(time_speedy[0:-49],moving_average(speedy_weight,50),color='#4daf4a',ls='-',linewidth=2.0,label='SPEEDY')
#plt.axhline(y=np.average(speedy_weight),color='#4daf4a',ls='--',linewidth=2.0)

plt.plot(time_era[0:-149],moving_average(era_weight[0:-100],50),color='#e41a1c',ls='-',linewidth=2.0,label='ERA 5')
plt.axhline(y=climo_weight,color='#e41a1c',ls='--',linewidth=2.0)

plt.legend()
plt.ylabel('Total Atmospheric Water Mass (kg) ',fontsize=16)
plt.xlabel('Years Into Simulation',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim(np.min(hybrid_weight)*0.98,np.max(hybrid_weight)*1.02)
plt.xlim(0,np.max(time_hybrid[0:-149]))
plt.show()
