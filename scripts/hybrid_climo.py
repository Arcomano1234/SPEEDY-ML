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

class StretchOutNormalize(plt.Normalize):
    def __init__(self, vmin=None, vmax=None, low=None, up=None, clip=False):
        self.low = low
        self.up = up
        plt.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.low, self.up, self.vmax], [0, 0.5-1e-9, 0.5+1e-9, 1]
        return np.ma.masked_array(np.interp(value, x, y))

@jit()
def rms(true,prediction):
    return np.sqrt(np.nanmean((prediction-true)**2))

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

plotdir = '/home/troyarcomano/FortranReservoir/vert_loc_hybridspeedy_leakage/plots/'

#ds_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.1_beta_model_1.0_prior_1.0_overlap1_era5_6hrtimestep_19years_climate_simulation_t_ends_2000_tisrtrial_12_30_1999_00.nc')
#ds_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.01_dp_beta_model_1_prior_0.0_overlap1_era5_6hrtimestep_tisr_full_years_0.00noise_speedy_states_climate_simtrial_12_31_1999_00.nc')
#ds_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_full_test_climate_all_tisr_longertrial_12_31_1999_00.nc')
#ds_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era5000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_2_vertlap_2_full_test_climate_all_tisr_longertrial_12_28_2009_00.nc')
#ds_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_slab_ocean_model_7d_0.1rho_10noise_beta0.001_20yearstrial_01_09_2010_00.nc')
#ds_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era3000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap2_vertlevels_4_vertlap_1_slab_ocean_model_falsetrial_01_02_2000_00.nc')
#ds_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era5000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap2_vertlevels_4_vertlap_3_slab_ocean_model_falsestrial_12_31_1999_00.nc')
#ds_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era5000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_4_vertlap_6_slab_ocean_model_falses_precip_truetrial_12_31_1999_00.nc')
ds_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_ml_only_ocean_beta_0.0001_sigma0.6_6diff_speedytrial_12_29_2006_00.nc')

print(ds_hybrid_sim)

lats = ds_hybrid_sim.Lat
lons = ds_hybrid_sim.Lon

target_pressures = [25,95,200,350,500,680,850,950]

temp_climo = np.zeros((1460,8,48,96))
v_wind_climo = np.zeros((1460,8,48,96))
u_wind_climo = np.zeros((1460,8,48,96))
sh_climo = np.zeros((1460,8,48,96))
logp_climo = np.zeros((1460,48,96))

start_year = 2000
start_date = datetime(2000,1,3,0)
#start_date = datetime(1999,12,29,0)

region_slice = slice(-90,90)
lon_slice = slice(0,360)
level = slice(0,8)


for i in range(10):
    year = start_year + i + 1
    start_climo_data = datetime(year,1,1,0)
    start_index = start_climo_data - start_date
    print(start_index.days)
    start_index = start_index.days * 4
    if(calendar.isleap(start_year+i+1)):
      time_slice = slice(start_index,start_index+240)
      print(start_index,start_index+240)
      temp_climo[0:240,:,:,:] += lin_interp(ds_hybrid_sim['Temperature'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      v_wind_climo[0:240,:,:,:] += lin_interp(ds_hybrid_sim['V-wind'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      u_wind_climo[0:240,:,:,:] += lin_interp(ds_hybrid_sim['U-wind'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      sh_climo[0:240,:,:,:] += lin_interp(ds_hybrid_sim['Specific-Humidity'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      logp_climo[0:240,:,:] +=  ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice)

      time_slice = slice(start_index+244,start_index+1464)
      temp_climo[240:1460,:,:,:] += lin_interp(ds_hybrid_sim['Temperature'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      v_wind_climo[240:1460,:,:,:] += lin_interp(ds_hybrid_sim['V-wind'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      u_wind_climo[240:1460,:,:,:] += lin_interp(ds_hybrid_sim['U-wind'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      sh_climo[240:1460,:,:,:] += lin_interp(ds_hybrid_sim['Specific-Humidity'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      logp_climo[240:1460,:,:] += ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice) 
    else: 
      time_slice = slice(start_index,start_index+1460)
      print(year,start_index,start_index+1460) 
      temp_climo += lin_interp(ds_hybrid_sim['Temperature'].sel(Lon=lon_slice,Lat=region_slice,Sigma_Level=level,Timestep=time_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      v_wind_climo += lin_interp(ds_hybrid_sim['V-wind'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      u_wind_climo += lin_interp(ds_hybrid_sim['U-wind'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      sh_climo += lin_interp(ds_hybrid_sim['Specific-Humidity'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      logp_climo += ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice)

ds_hybrid_sim.close()

temp_climo_hybrid = temp_climo/10.0
v_wind_climo_hybrid = v_wind_climo/10.0
u_wind_climo_hybrid = u_wind_climo/10.0
sh_climo_hybrid = sh_climo/10.0
logp_climo_hybrid = logp_climo/10.0
surfacep_climo_hybrid = np.exp(logp_climo_hybrid) * 1000.0


ds_speedy_sim = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/climo_speedy/speedy_era_start_climo_sst_decade_sim12_31_1999_00.nc')#'/scratch/user/troyarcomano/Predictions/Hybrid/speedy_era_start01_07_2000_01.nc')

print(ds_speedy_sim)

temp_climo = np.zeros((1460,8,48,96))
v_wind_climo = np.zeros((1460,8,48,96))
u_wind_climo = np.zeros((1460,8,48,96))
sh_climo = np.zeros((1460,8,48,96))
logp_climo = np.zeros((1460,48,96))

for i in range(10):
    year = start_year + i + 1
    start_climo_data = datetime(year,1,1,0)
    start_index = start_climo_data - start_date
    print(start_index.days)
    start_index = start_index.days * 4
    if(calendar.isleap(start_year+i+1)):
      time_slice = slice(start_index,start_index+240)
      print(start_index,start_index+240)
      temp_climo[0:240,:,:,:] += lin_interp(ds_speedy_sim['Temperature'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      v_wind_climo[0:240,:,:,:] += lin_interp(ds_speedy_sim['V-wind'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      u_wind_climo[0:240,:,:,:] += lin_interp(ds_speedy_sim['U-wind'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      sh_climo[0:240,:,:,:] += lin_interp(ds_speedy_sim['Specific_Humidity'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      logp_climo[0:240,:,:] +=  ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice)

      time_slice = slice(start_index+244,start_index+1464)
      temp_climo[240:1460,:,:,:] += lin_interp(ds_speedy_sim['Temperature'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      v_wind_climo[240:1460,:,:,:] += lin_interp(ds_speedy_sim['V-wind'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      u_wind_climo[240:1460,:,:,:] += lin_interp(ds_speedy_sim['U-wind'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      sh_climo[240:1460,:,:,:] += lin_interp(ds_speedy_sim['Specific_Humidity'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      logp_climo[240:1460,:,:] += ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice)
    else:
      time_slice = slice(start_index,start_index+1460)
      print(year,start_index,start_index+1460)
      temp_climo += lin_interp(ds_speedy_sim['Temperature'].sel(Lon=lon_slice,Lat=region_slice,Sigma_Level=level,Timestep=time_slice).values,ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      v_wind_climo += lin_interp(ds_speedy_sim['V-wind'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      u_wind_climo += lin_interp(ds_speedy_sim['U-wind'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      sh_climo += lin_interp(ds_speedy_sim['Specific_Humidity'].sel(Timestep=time_slice,Sigma_Level=level,Lon=lon_slice,Lat=region_slice).values,ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
      logp_climo += ds_speedy_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice)

ds_speedy_sim.close()

temp_climo_speedy = temp_climo/10.0
v_wind_climo_speedy = v_wind_climo/10.0
u_wind_climo_speedy = u_wind_climo/10.0
sh_climo_speedy = sh_climo/10.0
logp_climo_speedy = logp_climo/10.0
surfacep_climo_speedy = np.exp(logp_climo_speedy) * 1000.0

ds_climo = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/regridded_era_climatology2001_2011.nc')

time_slice = slice(0,8760,6)
temp_climo_era = lin_interp(ds_climo['Temperature'].sel(Lon=lon_slice,Lat=region_slice,Sigma_Level=level,Timestep=time_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
uwind_climo_era = lin_interp(ds_climo['U-wind'].sel(Lon=lon_slice,Lat=region_slice,Sigma_Level=level,Timestep=time_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
vwind_climo_era = lin_interp(ds_climo['V-wind'].sel(Lon=lon_slice,Lat=region_slice,Sigma_Level=level,Timestep=time_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
sh_climo_era = lin_interp(ds_climo['Specific_Humidity'].sel(Lon=lon_slice,Lat=region_slice,Sigma_Level=level,Timestep=time_slice).values,ds_hybrid_sim['logp'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice).values,target_pressures)
logp_climo_era = ds_climo['logp'].sel(Lon=lon_slice,Lat=region_slice,Timestep=time_slice)
surfacep_climo_era = np.exp(logp_climo_era) * 1000.0


temp_diff_hybrid = temp_climo_hybrid - temp_climo_era
temp_diff_speedy = temp_climo_speedy - temp_climo_era

uwind_diff_hybrid = u_wind_climo_hybrid - uwind_climo_era
uwind_diff_speedy = u_wind_climo_speedy - uwind_climo_era

surfacep_diff_hybrid = surfacep_climo_hybrid - surfacep_climo_era
surfacep_diff_speedy = surfacep_climo_speedy - surfacep_climo_era

sh_diff_hybrid = sh_climo_hybrid - sh_climo_era
sh_diff_speedy = sh_climo_speedy - sh_climo_era

##Months 
jan = np.arange(0,124)
feb = np.arange(124,236)
march = np.arange(236,360)
april = np.arange(360,480)
may = np.arange(480,604)
june = np.arange(604,724)
july = np.arange(724,848)
august = np.arange(848,972)
septemember = np.arange(972,1092)
october = np.arange(1092,1216)
november = np.arange(1216,1336)
december = np.arange(1336,1460)

months_jja = list(june)+list(july)+list(august)
months_mam = list(march)+list(april)+list(may)
months_djf = list(december)+list(jan)+list(feb)
months_son = list(septemember)+list(october)+list(november)
#months = list(december)+list(jan)+list(feb) #months_jja#list(december)+list(jan)+list(feb)

months = months_jja + months_mam + months_djf + months_son

##Temp
temp_error_diffs = abs(temp_diff_speedy) - abs(temp_diff_hybrid)

average_temp_error_diff = np.average(temp_error_diffs[months,:,:,:],axis=(0,3))

average_temp_diff_speedy = np.average(temp_diff_speedy[months,:,:,:],axis=(0,3))

average_temp_diff_hybrid = np.average(temp_diff_hybrid[months,:,:,:],axis=(0,3))

average_temp_era =  np.average(temp_climo_era[months,:,:,:],axis=(0,3))


#Uwind
uwind_error_diffs = abs(uwind_diff_speedy) - abs(uwind_diff_hybrid)

average_uwind_error_diff = np.average(uwind_error_diffs[months,:,:,:],axis=(0,3))

average_uwind_diff_speedy = np.average(uwind_diff_speedy[months,:,:,:],axis=(0,3))

average_uwind_diff_hybrid = np.average(uwind_diff_hybrid[months,:,:,:],axis=(0,3))

average_uwind_era =  np.average(uwind_climo_era[months,:,:,:],axis=(0,3))


#SH 
sh_error_diffs = abs(sh_diff_speedy) - abs(sh_diff_hybrid)

average_sh_error_diff = np.average(sh_error_diffs[months,:,:,:],axis=(0,3))

average_sh_diff_speedy = np.average(sh_diff_speedy[months,:,:,:],axis=(0,3))

average_sh_diff_hybrid = np.average(sh_diff_hybrid[months,:,:,:],axis=(0,3))

average_sh_era =  np.average(sh_climo_era[months,:,:,:],axis=(0,3))

##surface p 
surfacep_error_diffs = abs(surfacep_diff_speedy) - abs(surfacep_diff_hybrid)

average_surfacep_error_diff = np.average(surfacep_error_diffs[months,:,:],axis=0)

average_surfacep_diff_speedy = np.average(surfacep_diff_speedy[months,:,:],axis=0)

average_surfacep_diff_hybrid = np.average(surfacep_diff_hybrid[months,:,:],axis=0)

##surface p
surfacep_error_diffs_jja = abs(surfacep_diff_speedy) - abs(surfacep_diff_hybrid)

average_surfacep_error_diff_jja = np.average(surfacep_error_diffs[months_jja,:,:],axis=0)

average_surfacep_diff_speedy_jja = np.average(surfacep_diff_speedy[months_jja,:,:],axis=0)

average_surfacep_diff_hybrid_jja = np.average(surfacep_diff_hybrid[months_jja,:,:],axis=0)

#Average zonal pressure
average_zonal_pressure_era = np.average(surfacep_climo_era[months,:,:],axis=(0,2))


better_diff = np.average(temp_climo_speedy[months,:,:,:],axis=(0,3)) - np.average(temp_climo_era[months,:,:,:],axis=(0,3))
print('SPEEDY max,min, better average zonal temperature bias',np.max(better_diff),np.min(better_diff),np.average(better_diff))
print('SPEEDY max,min,average zonal temperature bias',np.max(average_temp_diff_speedy),np.min(average_temp_diff_speedy),np.average(average_temp_diff_speedy))
print('SPEEDY Temperature RMSE',rms(np.average(temp_climo_era[months,2:-1,:,:],axis=(0,3)),np.average(temp_climo_speedy[months,2:-1,:,:],axis=(0,3))))

better_diff = np.average(u_wind_climo_speedy[months,:,:,:],axis=(0,3)) - np.average(uwind_climo_era[months,:,:,:],axis=(0,3))
print('SPEEDY max,min,average zonal wind bias',np.max(better_diff),np.min(better_diff),np.average(better_diff))
print('SPEEDY zonal wind RMSE',rms(np.average(uwind_climo_era[months,2:-1,:,:],axis=(0,3)),np.average(u_wind_climo_speedy[months,2:-1,:,:],axis=(0,3))))

better_diff = np.average(temp_climo_hybrid[months,2:-1,:,:],axis=(0,3)) - np.average(temp_climo_era[months,2:-1,:,:],axis=(0,3))

print('Hybrid max,min,average zonal temperature bias',np.max(better_diff),np.min(better_diff),np.average(better_diff))
print('Hybrid temperature RMSE',rms(np.average(temp_climo_era[months,2:-1,:,:],axis=(0,3)),np.average(temp_climo_hybrid[months,2:-1,:,:],axis=(0,3))))
better_diff = np.average(u_wind_climo_hybrid[months,2:-1,:,:],axis=(0,3)) - np.average(uwind_climo_era[months,2:-1,:,:],axis=(0,3))
print('Hybrid max,min,average zonal wind bias',np.max(better_diff),np.min(better_diff),np.average(better_diff))
print('Hybrid zonal wind RMSE',rms(np.average(uwind_climo_era[months,2:-1,:,:],axis=(0,3)),np.average(u_wind_climo_hybrid[months,2:-1,:,:],axis=(0,3))))

diff_max = 21

cyclic_data, cyclic_lons = add_cyclic_point(average_surfacep_diff_speedy, coord=lons)
lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
projection = ccrs.PlateCarree()
axes_class = (GeoAxes,dict(map_projection=projection))
plt.rc('font', family='serif')

fig = plt.figure()
axgr = AxesGrid(fig, 111, axes_class=axes_class,
                nrows_ncols=(2, 2),
                axes_pad=0.9,
                cbar_location='right',
                cbar_mode='each',
                cbar_pad=0.2,
                cbar_size='3%',
                label_mode='')  # note the empty label_mode


ax1 = axgr[0]#plt.subplot(221,projection=ccrs.PlateCarree())
ax1.coastlines()

ax1.set_xticks(np.linspace(-180, 180, 5), crs=projection)
#ax1.set_yticks(np.linspace(-90, 90, 5), crs=projection)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
#lat_formatter = LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
#ax1.yaxis.set_major_formatter(lat_formatter)

plot = ax1.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=-1*diff_max,vmax=diff_max,cmap='seismic')

cbar = axgr.cbar_axes[0].colorbar(plot, extend='both')
#cbar.set_label_text(f'hPa',fontsize=16)
cbar.set_label(f'hPa',fontsize=16)
cbar.ax.tick_params(labelsize=16)

ax1.set_title(f"SPEEDY Surface Pressure Bias DJF",fontsize=18,fontweight="bold")

###########

diff_max = 21

cyclic_data, cyclic_lons = add_cyclic_point(average_surfacep_diff_hybrid, coord=lons)
lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

ax2 = axgr[1] #plt.subplot(222,projection=ccrs.PlateCarree())
ax2.coastlines()
ax2.set_xticks(np.linspace(-180, 180, 5), crs=projection)
#ax2.set_yticks(np.linspace(-90, 90, 5), crs=projection)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
#lat_formatter = LatitudeFormatter()
ax2.xaxis.set_major_formatter(lon_formatter)
#ax2.yaxis.set_major_formatter(lat_formatter)

plot_2 = ax2.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=-1*diff_max,vmax=diff_max,cmap='seismic')
cbar2 = axgr.cbar_axes[1].colorbar(plot_2, extend='both')
#cbar2.set_label_text(f'hPa',fontsize=16)
cbar2.set_label(f'hPa',fontsize=16)
cbar2.ax.tick_params(labelsize=16)

ax2.set_title(f"Hybrid Surface Pressure Bias DJF",fontsize=18,fontweight="bold")

diff_max = 21

cyclic_data, cyclic_lons = add_cyclic_point(average_surfacep_diff_speedy_jja, coord=lons)
lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

ax3 = axgr[2] #plt.subplot(223,projection=ccrs.PlateCarree())
ax3.coastlines()
ax3.set_xticks(np.linspace(-180, 180, 5), crs=projection)
#ax3.set_yticks(np.linspace(-90, 90, 5), crs=projection)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
#lat_formatter = LatitudeFormatter()
ax3.xaxis.set_major_formatter(lon_formatter)
#ax3.yaxis.set_major_formatter(lat_formatter)

plot_3 = ax3.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=-1*diff_max,vmax=diff_max,cmap='seismic')

cbar = axgr.cbar_axes[2].colorbar(plot_3,extend='both')
#cbar.set_label_text(f'hPa',fontsize=16)
cbar.set_label(f'hPa',fontsize=16)
cbar.ax.tick_params(labelsize=16)

ax3.set_title(f"SPEEDY Surface Pressure Bias JJA",fontsize=18,fontweight="bold")

###########
diff_max = 21

cyclic_data, cyclic_lons = add_cyclic_point(average_surfacep_diff_hybrid_jja, coord=lons)
lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

ax4 = axgr[3]#plt.subplot(224,projection=ccrs.PlateCarree())
ax4.coastlines()
ax4.set_xticks(np.linspace(-180, 180, 5), crs=projection)
#ax4.set_yticks(np.linspace(-90, 90, 5), crs=projection)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
#lat_formatter = LatitudeFormatter()
ax4.xaxis.set_major_formatter(lon_formatter)
#ax4.yaxis.set_major_formatter(lat_formatter)

plot_4 = ax4.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=-1*diff_max,vmax=diff_max,cmap='seismic')
cbar4 = axgr.cbar_axes[3].colorbar(plot_4,extend='both')
#cbar4.set_label_text(f'hPa',fontsize=16)
cbar4.set_label(f'hPa',fontsize=16)
cbar4.ax.tick_params(labelsize=16)

ax4.set_title(f"Hybrid Surface Pressure Bias JJA",fontsize=18,fontweight="bold")
'''
###########
diff_max = 50
cyclic_data, cyclic_lons = add_cyclic_point(average_surfacep_error_diff, coord=lons)
lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

ax3 = plt.subplot(313,projection=ccrs.PlateCarree())
ax3.coastlines()

plot_3 = ax3.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=-1*diff_max,vmax=diff_max,cmap='seismic')
cbar3 = fig.colorbar(plot_3, ax=ax3, extend='both')
cbar3.set_label(f'hPa',fontsize=16)
cbar3.ax.tick_params(labelsize=16)
ax3.set_title(f"Difference (SPEEDY - Hybrid)",fontsize=18,fontweight="bold")
'''
plt.show()
#plt.tight_layout()
#plt.savefig(f'{plotdir}hybrid_surfacep_climo_error.png')
#plt.savefig(f'{plotdir}hybrid_surfacep_djf_climo_error.pdf')
#plt.close("all")

pressure_levels = [25,95,200,350,500,680,850,950] 
pressure_levels = np.flip(pressure_levels)

ticktimes = [25,100,200,300,400,500,600,700,800,900]
plt.rc('font', family='serif')
x_tick_labels = [u'60\N{DEGREE SIGN}S', u'30\N{DEGREE SIGN}S',u'EQ', u'30\N{DEGREE SIGN}N',u'60\N{DEGREE SIGN}N']

fig, axes = plt.subplots(3,2,figsize=(18,10),sharex=True,sharey='row',constrained_layout=True)

norm=StretchOutNormalize(vmin=-8, vmax=8, low=-0.5, up=0.5)
####
#SPEEDY TEMP BIAS
####
ax = axes[0,0]#plt.subplot(321)
ax.invert_yaxis()  # Reverse the time order to do oldest first

clevs = np.arange(-8,9,0.5)

cf = ax.contourf(lats,pressure_levels,np.flipud(average_temp_diff_speedy),cmap='seismic',norm=norm,levels=clevs,extend='both')
cs = ax.contour(lats,pressure_levels,np.flipud(average_temp_era) - 273.15,levels=np.arange(-80,40,10),linestyles='--',colors='k')

ax.clabel(cs,fmt='%2.1f', colors='k', fontsize=14)

ax.fill_between(lats,y1=950, y2=average_zonal_pressure_era,color='grey',zorder=10)
ax.set_yticks(ticktimes)
ax.set_yticklabels(ticktimes,fontsize=16)
ax.set_ylabel('hPa',fontsize=18)
ax.set_title('SPEEDY Temperature Bias',fontsize=16,fontweight='bold')
ax.grid()
ax.minorticks_off()

for tic in ax.xaxis.get_major_ticks():
    tic.tick1On = tic.tick2On = False
    tic.label1On = tic.label2On = False

ax.set_ylim(bottom=950, top=100)



####
#HYBRID TEMP BIAS
####
ax = axes[0,1]#plt.subplot(322)
ax.invert_yaxis()  # Reverse the time order to do oldest first

cf_temp_hybrid = ax.contourf(lats,pressure_levels,np.flipud(average_temp_diff_hybrid),cmap='seismic',norm=norm,levels=clevs,extend='both')
cs = ax.contour(lats,pressure_levels,np.flipud(average_temp_era) - 273.15,levels=np.arange(-80,40,10),linestyles='--',colors='k')
ax.clabel(cs,fmt='%2.1f', colors='k', fontsize=14)

ax.fill_between(lats,y1=950, y2=average_zonal_pressure_era,color='grey',zorder=10)

ax.set_title('Hybrid Temperature Bias',fontsize=16,fontweight='bold')
ax.grid()
ax.minorticks_off()

ax.set_ylim(bottom=950, top=100)

for tic in ax.yaxis.get_major_ticks():
    tic.tick1On = tic.tick2On = False
    tic.label1On = tic.label2On = False

cbar = fig.colorbar(cf_temp_hybrid, ax=axes[0,:], extend='both')
cbar.set_label(f'Kelvin',fontsize=16)
cbar.ax.tick_params(labelsize=16)

#####
#SPEEDY UWIND BIAS
####
ax = axes[1,0]#plt.subplot(323)
ax.invert_yaxis()  # Reverse the time order to do oldest first

cf = ax.contourf(lats,pressure_levels,np.flipud(average_uwind_diff_speedy),cmap='seismic',norm=norm,levels=clevs,extend='both')
cs = ax.contour(lats,pressure_levels,np.flipud(average_uwind_era),levels=np.arange(-60,80,10),linestyles='--',colors='k')
ax.clabel(cs,fmt='%2.1f', colors='k', fontsize=14)

ax.fill_between(lats,y1=950, y2=average_zonal_pressure_era,color='grey',zorder=10)

ax.set_yticks(ticktimes)
ax.set_yticklabels(ticktimes,fontsize=16)
ax.set_ylabel('hPa',fontsize=18)
ax.set_title('SPEEDY Zonal Wind Bias',fontsize=16,fontweight='bold')
ax.grid()
ax.minorticks_off()

ax.set_ylim(bottom=950, top=100)

for tic in ax.xaxis.get_major_ticks():
    tic.tick1On = tic.tick2On = False
    tic.label1On = tic.label2On = False


####
#HYBRID UWIND BIAS
####
ax = axes[1,1]#plt.subplot(324)
ax.invert_yaxis()  # Reverse the time order to do oldest first

cf_hybrid_uwind = ax.contourf(lats,pressure_levels,np.flipud(average_uwind_diff_hybrid),cmap='seismic',norm=norm,levels=clevs,extend='both')
#cs = plt.contour(cf,levels=cf.levels[::2], colors='k')
cs = ax.contour(lats,pressure_levels,np.flipud(average_uwind_era),levels=np.arange(-60,80,10),linestyles='--',colors='k')
ax.clabel(cs,fmt='%2.1f', colors='k', fontsize=14)

ax.fill_between(lats,y1=950, y2=average_zonal_pressure_era,color='grey',zorder=10)

ax.set_title('Hybrid Zonal Wind Bias',fontsize=16,fontweight='bold')
ax.grid()
ax.minorticks_off()

ax.set_ylim(bottom=950, top=100)

for tic in ax.yaxis.get_major_ticks():
    tic.tick1On = tic.tick2On = False
    tic.label1On = tic.label2On = False

for tic in ax.xaxis.get_major_ticks():
    tic.tick1On = tic.tick2On = False
    tic.label1On = tic.label2On = False

cbar = fig.colorbar(cf_hybrid_uwind, ax=axes[1,:], extend='both')
cbar.set_label(f'm/s',fontsize=16)
cbar.ax.tick_params(labelsize=16)


#####
#SPEEDY SH BIAS
####

norm=StretchOutNormalize(vmin=-2, vmax=2, low=0, up=0)

ax = axes[2,0]#plt.subplot(325)
ax.invert_yaxis()  # Reverse the time order to do oldest first

cf = ax.contourf(lats,pressure_levels,np.flipud(average_sh_diff_speedy),cmap='seismic',norm=norm,levels=np.arange(-2,2.2,0.2),extend='both')
#cs = plt.contour(cf,levels=cf.levels[::2], colors='k')
cs = ax.contour(lats,pressure_levels,np.flipud(average_sh_era),levels=np.arange(0,20,2),linestyles='--',colors='k')
ax.clabel(cs,fmt='%2.1f', colors='k', fontsize=14)

ax.fill_between(lats,y1=950, y2=average_zonal_pressure_era,color='grey',zorder=10)

ax.set_xticks([-60, -30, 0, 30, 60])
ax.set_xticklabels(x_tick_labels,fontsize=16)
ax.set_xlabel('Latitude',fontsize=16)
ax.set_yticks(ticktimes)
ax.set_yticklabels(ticktimes,fontsize=16)
ax.set_ylabel('hPa',fontsize=18)
ax.set_title('SPEEDY Zonal Specific Humidity Bias',fontsize=16,fontweight='bold')
ax.grid()
ax.minorticks_off()

ax.set_ylim(bottom=950, top=400)


#####
#HYBRID SH BIAS
#####
ax = axes[2,1]#plt.subplot(326)
ax.invert_yaxis()  # Reverse the time order to do oldest first

cf = ax.contourf(lats,pressure_levels,np.flipud(average_sh_diff_hybrid),cmap='seismic',norm=norm,levels=np.arange(-2,2.2,0.2),extend='both')
#cs = plt.contour(cf,levels=cf.levels[::2], colors='k')
cs = ax.contour(lats,pressure_levels,np.flipud(average_sh_era),levels=np.arange(0,20,2),linestyles='--',colors='k')
ax.clabel(cs,fmt='%2.1f', colors='k', fontsize=14)

ax.fill_between(lats,y1=950, y2=average_zonal_pressure_era,color='grey',zorder=10)

ax.set_xticks([-60, -30, 0, 30, 60])
ax.set_xticklabels(x_tick_labels,fontsize=16)
ax.set_xlabel('Latitude',fontsize=16)
ax.set_title('Hybrid Zonal Specific Humidity Bias',fontsize=16,fontweight='bold')
ax.grid()
ax.minorticks_off()

for tic in ax.yaxis.get_major_ticks():
    tic.tick1On = tic.tick2On = False
    tic.label1On = tic.label2On = False

ax.set_ylim(bottom=950, top=400)

cbar = fig.colorbar(cf, ax=axes[2,:], extend='both')
cbar.set_label(f'g/kg',fontsize=16)
cbar.ax.tick_params(labelsize=16)

plt.savefig(f'{plotdir}hybrid_annual_three_panel_sst_coupled_precip_enso.png')
plt.savefig(f'{plotdir}hybrid_annual_three_panel_sst_coupled_precip_enso.pdf')
plt.close("all")
#plt.show()
