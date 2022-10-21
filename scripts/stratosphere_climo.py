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

def era_5_climo(start_year,end_year,lon_slice,region_slice,var,height):
    root_path = '/scratch/user/troyarcomano/ERA/'
    hours_in_year = 24*365
    xgrid = 96
    ygrid = 48
    average_var = np.zeros((hours_in_year,ygrid,xgrid))

    for current_year in range(start_year,end_year + 1):
        try:
           ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{current_year}/era_5_y{current_year}_regridded_mpi.nc')
        except:
           ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{current_year}/era_5_y{current_year}_regridded_mpi_fixed_var.nc')

        if current_year == start_year:
           shape = np.shape(ds_era[var].sel(Lon=lon_slice,Lat=region_slice,Sigma_Level=height))
           xgrid = shape[2]
           ygrid = shape[1]
           average_var = np.zeros((hours_in_year,ygrid,xgrid))

        if calendar.isleap(current_year):
           time_slice = slice(0,240*6)
           average_var[0:240*6,:,:] += ds_era[var].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice,Sigma_Level=height)

           time_slice = slice(244*6,1464*6)
           average_var[240*6:1460*6,:,:] += ds_era[var].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice,Sigma_Level=height)
        else:
           time_slice = slice(0,hours_in_year)
           average_var += ds_era[var].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice,Sigma_Level=height)

    return average_var/(end_year-start_year+1)

def get_obs_era5_timeseries(startdate,enddate,timestep,lat_slice,lon_slice,sigma):
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

        begin_year = datetime(currentdate.year,1,1,0)
        begin_year_str = begin_year.strftime("%Y-%m-%d")
        attrs = {"units": f"hours since {begin_year_str} "}
        ds_era = ds_era.assign_coords({"Timestep": ("Timestep", ds_era.Timestep.values, attrs)})
        ds_era = xr.decode_cf(ds_era)

        ds_era = ds_era['U-wind'].sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma)

        if start_year == currentdate.year:
           ds_merged = ds_era
        else:
           ds_merged = xr.merge([ds_merged,ds_era])

        currentdate = currentdate + timedelta(hours=ds_era.sizes['Timestep'])

    time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"),timestep)
    return ds_merged.sel(Timestep=time_slice)

def zonal_wind_mean_plot(ds_era,ds_hybrid,ds_speedy):
    ds_zmean_era = ds_era.mean(dim='Lon')
    ds_zmean_era = ds_zmean_era.mean(dim='Lat') 

    mean_era = ds_zmean_era.groupby(ds_era.Timestep.dt.dayofyear).mean("Timestep")['U-wind']
    std_era = ds_zmean_era.groupby(ds_era.Timestep.dt.dayofyear).std("Timestep")['U-wind']

    ds_zmean_hybrid = ds_hybrid.mean(dim='Lon')
    ds_zmean_hybrid = ds_zmean_hybrid.mean(dim='Lat')

    mean_hybrid = ds_zmean_hybrid.groupby(ds_hybrid.Timestep.dt.dayofyear).mean("Timestep")['U-wind']
    std_hybrid = ds_zmean_hybrid.groupby(ds_hybrid.Timestep.dt.dayofyear).std("Timestep")['U-wind']

    ds_zmean_speedy = ds_speedy.mean(dim='Lon')
    ds_zmean_speedy = ds_zmean_speedy.mean(dim='Lat')

    mean_speedy = ds_zmean_speedy.groupby(ds_speedy.Timestep.dt.dayofyear).mean("Timestep")['U-wind']
    std_speedy = ds_zmean_speedy.groupby(ds_speedy.Timestep.dt.dayofyear).std("Timestep")['U-wind']

    timestep = mean_era.dayofyear.values

    fig, axes = plt.subplots(1,3,figsize=(18,10),sharex=True,sharey=True,constrained_layout=True)

    mean_era = np.roll(mean_era,-180)
    mean_hybrid = np.roll(mean_hybrid,-180)
    mean_speedy = np.roll(mean_speedy,-180)
    
    std_era = np.roll(std_era,-180)
    std_hybrid = np.roll(std_hybrid,-180)
    std_speedy = np.roll(std_speedy,-180)


    x_labels_vals = [15,45,76,107,138,168,199,229,257,287,317,348]
    x_tick_labels = ['July','Aug','Sept','Oct','Nov','Dec','Jan','Feb','March','April','May','June']
    ####
    #ERA Strato
    ####
    ax = axes[0]

    ax.plot(timestep,mean_era)
    ax.fill_between(timestep,mean_era+std_era*2,mean_era-std_era*2, facecolor='grey', alpha=0.4) 
   
    ax.set_title('ERA5',fontsize=20,fontweight='bold')
    ax.set_xlim([np.min(timestep),np.max(timestep)]) 
    ax.set_ylim([-20,60])
    #ax.set_xlabel('Day of Year')
    ax.set_ylabel('20 hPa Zonal Mean Wind')

    ax.set_xticks(x_labels_vals)
    ax.set_xticklabels(x_tick_labels, rotation=30)

  
    ####
    #Hybrid Strato
    ####
    ax = axes[1]

    ax.plot(timestep,mean_hybrid)
    ax.fill_between(timestep,mean_hybrid+std_hybrid*2,mean_hybrid-std_hybrid*2, facecolor='grey', alpha=0.4)

    ax.set_title('Hybrid',fontsize=20,fontweight='bold')
    ax.set_xlim([np.min(timestep),np.max(timestep)])
    ax.set_ylim([-20,60])
    #ax.set_xlabel('Day of Year')
    #ax.set_ylabel('Hybrid 20 hPa Zonal Mean Wind (2000-2011)')
    ax.set_xticks(x_labels_vals)
    ax.set_xticklabels(x_tick_labels, rotation=30)

    ####
    #SPEEDY Strato
    ####
    ax = axes[2]

    ax.plot(timestep,mean_speedy)
    ax.fill_between(timestep,mean_speedy+std_speedy*2,mean_speedy+std_speedy-std_speedy*2, facecolor='grey', alpha=0.4)

    ax.set_title('SPEEDY',fontsize=20,fontweight='bold')
    ax.set_xlim([np.min(timestep),np.max(timestep)])
    ax.set_ylim([-20,60])
    #ax.set_xlabel('Day of Year')
    #ax.set_ylabel('SPEEDY 20 hPa Zonal Mean Wind (2000-2011)')
    ax.set_xticks(x_labels_vals)
    ax.set_xticklabels(x_tick_labels, rotation=30)
    plt.tight_layout()
    plt.show()

def make_ds_time_dim(ds,timestep,startdate):
    begin_year_str = startdate.strftime("%Y-%m-%d")

    attrs = {"units": f"hours since {begin_year_str} "}

    ds = ds.assign_coords({"Timestep": ("Timestep", ds.Timestep.values*timestep, attrs)})
    ds = xr.decode_cf(ds)
 
    return ds

def qbo_plot(ds_era,ds_hybrid,ds_speedy):
    lat_slice = slice(-35,35)
    lon_slice = slice(0,365) 

    sigma = 0

    ds_era = ds_era.sel(Lat=lat_slice,Lon=lon_slice)
    ds_era_zmean = ds_era.mean(dim='Lon') 

    ds_hybrid = ds_hybrid.sel(Lat=lat_slice,Lon=lon_slice)
    ds_hybrid_zmean = ds_hybrid.mean(dim='Lon')

    ds_speedy = ds_speedy.sel(Lat=lat_slice,Lon=lon_slice)
    ds_speedy_zmean = ds_speedy.mean(dim='Lon')

    x_tick_labels = [u'30\N{DEGREE SIGN}S', u'15\N{DEGREE SIGN}S',"EQ", u'15\N{DEGREE SIGN}N',u'30\N{DEGREE SIGN}N'] 

    fig, axes = plt.subplots(1,3,figsize=(18,10),sharex=True,sharey=True,constrained_layout=True)

    ax = axes[0]
    ax.contourf(ds_era_zmean.Lat.values,ds_era_zmean.Timestep.values,ds_era_zmean['U-wind'].values,levels=np.arange(-30,35,5),cmap='seismic',extend="both")
    ax.set_title('ERA5',fontsize=20,fontweight='bold')

    ax = axes[1]
    ax.contourf(ds_hybrid_zmean.Lat.values,ds_hybrid_zmean.Timestep.values,ds_hybrid_zmean['U-wind'].values,levels=np.arange(-30,35,5),cmap='seismic',extend="both")
    ax.set_title('Hybrid',fontsize=20,fontweight='bold')

    ax = axes[2]
    cf = ax.contourf(ds_speedy_zmean.Lat.values,ds_speedy_zmean.Timestep.values,ds_speedy_zmean['U-wind'].values,levels=np.arange(-30,35,5),cmap='seismic',extend="both")
    ax.set_title('SPEEDY',fontsize=20,fontweight='bold')

    plt.subplots_adjust(top=0.93, bottom=0.16, left=0.04, right=0.99, hspace=0.2, wspace=0.135)
    cax = plt.axes([0.125, 0.08, 0.80, 0.025])
    cbar = fig.colorbar(cf,cax,orientation='horizontal',fraction=0.1)
    cbar.set_label('m $s^{-1}$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    plt.show()

def qbo_plot_darpa_talk(ds_era,ds_hybrid,ds_hybrid_old,ds_speedy):
    lat_slice = slice(-35,35)
    lon_slice = slice(0,365)

    sigma = 0

    ds_era = ds_era.sel(Lat=lat_slice,Lon=lon_slice)
    ds_era_zmean = ds_era.mean(dim='Lon')

    ds_hybrid = ds_hybrid.sel(Lat=lat_slice,Lon=lon_slice)
    ds_hybrid_zmean = ds_hybrid.mean(dim='Lon')

    ds_hybrid_old = ds_hybrid_old.sel(Lat=lat_slice,Lon=lon_slice)
    ds_hybrid_zmean_old = ds_hybrid_old.mean(dim='Lon')

    ds_speedy = ds_speedy.sel(Lat=lat_slice,Lon=lon_slice)
    ds_speedy_zmean = ds_speedy.mean(dim='Lon')

    x_tick_labels = [u'30\N{DEGREE SIGN}S', u'15\N{DEGREE SIGN}S',"EQ", u'15\N{DEGREE SIGN}N',u'30\N{DEGREE SIGN}N']
    x_labels_vals = [-30,-15,0,15,30]
    
    fig, axes = plt.subplots(1,4,figsize=(18,10),sharex=True,sharey=True,constrained_layout=True)

    ax = axes[3]
    ax.contourf(ds_era_zmean.Lat.values,ds_era_zmean.Timestep.values,ds_era_zmean['U-wind'].values,levels=np.arange(-30,35,5),cmap='seismic',extend="both")
    ax.set_title('ERA5',fontsize=20,fontweight='bold')
    ax.set_xticks(x_labels_vals)
    ax.set_xticklabels(x_tick_labels)

    ax = axes[2]
    ax.contourf(ds_hybrid_zmean.Lat.values,ds_hybrid_zmean.Timestep.values,ds_hybrid_zmean['U-wind'].values,levels=np.arange(-30,35,5),cmap='seismic',extend="both")
    ax.set_title('Hybrid 3D Loc',fontsize=20,fontweight='bold')
    ax.set_xticks(x_labels_vals)
    ax.set_xticklabels(x_tick_labels)
 
    ax = axes[0]
    cf = ax.contourf(ds_speedy_zmean.Lat.values,ds_speedy_zmean.Timestep.values,ds_speedy_zmean['U-wind'].values,levels=np.arange(-30,35,5),cmap='seismic',extend="both")
    ax.set_title('SPEEDY',fontsize=20,fontweight='bold')
    ax.set_xticks(x_labels_vals)
    ax.set_xticklabels(x_tick_labels)

    ax = axes[1]
    ax.contourf(ds_hybrid_zmean_old.Lat.values,ds_hybrid_zmean_old.Timestep.values,ds_hybrid_zmean_old['U-wind'].values,levels=np.arange(-30,35,5),cmap='seismic',extend="both")
    ax.set_title('Hybrid 2D Loc',fontsize=20,fontweight='bold')
    ax.set_xticks(x_labels_vals)
    ax.set_xticklabels(x_tick_labels)

    plt.subplots_adjust(top=0.93, bottom=0.16, left=0.04, right=0.99, hspace=0.2, wspace=0.135)
    cax = plt.axes([0.125, 0.08, 0.80, 0.025])
    cbar = fig.colorbar(cf,cax,orientation='horizontal',fraction=0.1)
    cbar.set_label('m $s^{-1}$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    plt.show()

lat_slice = slice(-45,45)
lon_slice = slice(0,365)
sigma = 0

startdate = datetime(2010,1,1,0)
enddate = datetime(2018,12,31,0)

ds_era5 = get_obs_era5_timeseries(startdate,enddate,6,lat_slice,lon_slice,sigma)

#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era9000_15_30_30_beta_res0.1_beta_model_1_prior_0.0_overlap1_era5_6hrtimestep_sst_true_climo_inputtrial_01_02_2000_00.nc')
#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era5000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_2_vertlap_2_full_test_climate_all_tisr_longertrial_12_28_2009_00.nc')
ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_ml_only_ocean_beta_0.0001_sigma0.6_6diff_speedytrial_12_29_2006_00.nc')
#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era4000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_4_vertlap_2_full_test_climate_all_tisr_longertrial_12_28_2009_00.nc')

time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"))

ds_hybrid = make_ds_time_dim(ds_hybrid,6,startdate)


ds_speedy_sim = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/climo_speedy/speedy_era_start_climo_sst_decade_sim12_31_1999_00.nc')

#startdate_speedy = datetime(2000,1,1,0)
#enddate_speedy = datetime(2012,1,1,23)
startdate_speedy = datetime(2010,1,1,0)
enddate_speedy = datetime(2018,12,31,23)

time_slice_speedy = slice(startdate_speedy.strftime("%Y-%m-%d"),enddate_speedy.strftime("%Y-%m-%d"))
ds_speedy = make_ds_time_dim(ds_speedy_sim,6,startdate_speedy)

'''
ds_old_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.01_dp_beta_model_1_prior_0.0_overlap1_era5_6hrtimestep_tisr_full_years_0.00noise_speedy_states_climate_simtrial_12_31_1999_00.nc')
startdate_speedy = datetime(2010,1,1,0)
enddate_speedy = datetime(2018,12,31,23)
time_slice_speedy = slice(startdate_speedy.strftime("%Y-%m-%d"),enddate_speedy.strftime("%Y-%m-%d"))
ds_old_hybrid = make_ds_time_dim(ds_old_hybrid_sim,6,startdate_speedy)
'''
#qbo_plot_darpa_talk(ds_era5,ds_hybrid.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice),ds_old_hybrid.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice),ds_speedy.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice_speedy))


qbo_plot(ds_era5,ds_hybrid.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice),ds_speedy.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice_speedy))

'''

lat_slice = slice(55,65)
lon_slice = slice(0,365)
sigma = 0

startdate = datetime(1981,1,1,0)
enddate = datetime(2009,12,31,0)

ds_era5 = get_obs_era5_timeseries(startdate,enddate,6,lat_slice,lon_slice,sigma)
print(ds_era5.Timestep)

#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era9000_15_30_30_beta_res0.1_beta_model_1_prior_0.0_overlap1_era5_6hrtimestep_sst_true_climo_inputtrial_01_02_2000_00.nc')
#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_full_test_climate_all_tisr_longertrial_12_31_1999_00.nc')
ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_ml_only_ocean_beta_0.0001_sigma0.6_6diff_speedytrial_12_29_2006_00.nc')

#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era5000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_2_vertlap_2_full_test_climate_all_tisr_longertrial_12_28_2009_00.nc')
#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era4000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_4_vertlap_2_full_test_climate_all_tisr_longertrial_12_28_2009_00.nc')
startdate = datetime(2000,12,29,0)
enddate = datetime(2025,1,10,23)

ds_speedy_sim = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/climo_speedy/speedy_era_start_climo_sst_decade_sim12_31_1999_00.nc')
#ds_speedy_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.01_dp_beta_model_1_prior_0.0_overlap1_era5_6hrtimestep_tisr_full_years_0.00noise_speedy_states_climate_simtrial_12_31_1999_00.nc')
#ds_hybrid = ds_speedy_sim

time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"))

ds_hybrid = make_ds_time_dim(ds_hybrid,6,startdate)

startdate_speedy = datetime(2000,1,1,0)
enddate_speedy = datetime(2012,1,1,23)
time_slice_speedy = slice(startdate_speedy.strftime("%Y-%m-%d"),enddate_speedy.strftime("%Y-%m-%d"))
ds_speedy = make_ds_time_dim(ds_speedy_sim,6,startdate_speedy)

print(ds_hybrid.Timestep)
print(ds_era5)
print(ds_hybrid.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice))
print(ds_speedy.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice_speedy))

zonal_wind_mean_plot(ds_era5,ds_hybrid.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice),ds_speedy.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice_speedy))
'''
