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

        var = ['U-wind','logp']
        ds_era = ds_era['U-wind'].sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma)
        if start_year == currentdate.year:
           ds_merged = ds_era
        else:
           ds_merged = xr.merge([ds_merged,ds_era])

        currentdate = currentdate + timedelta(hours=ds_era.sizes['Timestep'])

    time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"),timestep)
    return  ds_merged.sel(Timestep=time_slice)

def get_6hr_precip_era5_timeseries(startdate,enddate,timestep,lat_slice,lon_slice):
    start_year = startdate.year
    end_year = enddate.year

    currentdate = startdate
    counter = 0
    while currentdate.year <= enddate.year:
        print(currentdate.year)
        try:
           ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{currentdate.year}/era_5_y{currentdate.year}_precip_regridded_mpi_fixed_var_gcc.nc')
        except:
           ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{currentdate.year}/era_5_y{currentdate.year}_regridded_mpi_fixed_var.nc')

        begin_year = datetime(currentdate.year,1,1,0)
        begin_year_str = begin_year.strftime("%Y-%m-%d")
        attrs = {"units": f"hours since {begin_year_str} "}
        ds_era = ds_era.assign_coords({"Timestep": ("Timestep", ds_era.Timestep.values, attrs)})
        ds_era = xr.decode_cf(ds_era)

        var = ['tp']
        ds_era = ds_era[var].sel(Lat=lat_slice,Lon=lon_slice)

        if start_year == currentdate.year:
           ds_merged = ds_era
        else:
           ds_merged = xr.merge([ds_merged,ds_era])

        currentdate = currentdate + timedelta(hours=ds_era.sizes['Timestep'])

    time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"),timestep)
    return  ds_merged.resample(Timestep = "6H").sum() #ds_merged.sel(Timestep=time_slice)

def make_ds_time_dim(ds,timestep,startdate):
    begin_year_str = startdate.strftime("%Y-%m-%d")

    attrs = {"units": f"hours since {begin_year_str} "}

    ds = ds.assign_coords({"Timestep": ("Timestep", ds.Timestep.values*timestep, attrs)})
    ds = xr.decode_cf(ds)
 
    return ds

def log_binning(data):
    #data unsorted
    copy = np.copy(data)

    copy = np.sort(copy)
    print(np.shape(copy))

    print(copy[-10:-1])

    
    percent_range = np.logspace(-4,2,100) #[80,90,95,99,99.9,99.95, 99.99,99.999, 99.9999] 
    power = percent_range
    percent_range = 100 - percent_range
    # log-scaled bins
    #bins = np.logspace(0, 3, 50)
    # Calculate histogram
    #hist = np.histogram(data, bins=bins)
    # normalize by bin width
    #hist_norm = hist[0]/widths

    # plot it!
    #plt.bar(bins[:-1], hist_norm, widths)
    #plt.xscale('log')
    #plt.yscale('log')
    return np.percentile(copy,percent_range), power #percent_range, power

#def mslp

def sigma_to_p(ds,variable):
    target_pressures = [25,95,200,350,500,680,850,950]
   
    test = lin_interp(ds[variable].values,ds['logp'].values,target_pressures)
    print(np.shape(test))
    print(np.shape(ds[variable].values))
    ds[variable][:] = test
    return ds

def uv_to_windspeedy(ds):
    return np.sqrt(np.power(ds['U-wind'].values,2) + np.power(ds['V-wind'].values, 2)).flatten()

def extreme_value_plot(ds_era,ds_hybrid,ds_speedy): 


    print(np.max(ds_hybrid['p6hr'].values))
    era_uwind = ds_era['tp'].to_numpy().flatten() * 1000.0 #uv_to_windspeedy(ds_era) #abs(ds_era['U-wind'].to_numpy().flatten())
    hybrid_uwind = ds_hybrid['p6hr'].to_numpy().flatten() * 25.4#uv_to_windspeedy(ds_hybrid) #abs(ds_hybrid['U-wind'].to_numpy().flatten())
    speedy_uwind = ds_speedy['rain'].to_numpy().flatten()/(3.6*4.0*6.0) #uv_to_windspeedy(ds_speedy) #abs(ds_speedy['U-wind'].to_numpy().flatten())

    print(np.shape(era_uwind))
    print(np.shape(hybrid_uwind))
    print(np.shape(speedy_uwind))

    print('max era',np.max(era_uwind))
    print('max hybrid',np.max(hybrid_uwind))
    print('max speedy',np.max(speedy_uwind))

    era_uwind_extreme, percent_range = log_binning(era_uwind) 
    hybrid_uwind_extreme, percent_range = log_binning(hybrid_uwind)
    speedy_uwind_extreme, percent_range = log_binning(speedy_uwind)
  
    percent_range = percent_range
    era_uwind_extreme = era_uwind_extreme
    hybrid_uwind_extreme = hybrid_uwind_extreme
    speedy_uwind_extreme = speedy_uwind_extreme

    fig1, ax1 = plt.subplots()

    ax1.semilogx(percent_range,era_uwind_extreme, linewidth=2,label='ERA5')
    ax1.semilogx(percent_range,hybrid_uwind_extreme,linewidth=2,label='Hybrid')
    ax1.semilogx(percent_range,speedy_uwind_extreme,linewidth=2,label="SPEEDY")

    ax1.set_xlim([10**-4,10])
    #ax1.set_xscale('log')
 
    x_labels_vals = [10**0,10**-1,10**-2,10**-3,10**-4]
    x_tick_labels = ['90%','99%','99.9%','99.99%','99.999%']
  
    ax1.tick_params(axis='y', which='major', labelsize=16)
    ax1.invert_xaxis()
    ax1.set_xticks(x_labels_vals)
    ax1.set_xticklabels(x_tick_labels,fontsize=16)  

    ax1.set_xlabel("Percentile",fontsize=20)

    ax1.set_ylabel("Total Precipitation (mm/6hr)",fontsize=20)
    

    #ax1.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    plt.legend(fontsize=20)
    plt.show()
    
def histo_precip():
    startdate = datetime(2000,1,1,0)
    enddate = datetime(2010,12,31,0)

    startdate_climo = datetime(2001,1,1,0)
    enddate_climo = datetime(2010,12,31,0)


    ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_ocean_model_false_log_preciptrial_12_31_1999_00.nc') #hybrid_prediction_era6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_ocean_model_false_log_precip_noise_add_raw_epsilon_0.0005trial_12_31_1999_00.nc')#

    timestep = 6

    ds_speedy = xr.open_dataset('/scratch/user/troyarcomano/letkf-hybrid-speedy/experiments/ERA_5_cycle_2_year_test/nature.nc')
    ds_speedy = speedy_total_precip(ds_speedy)
    ds_speedy = ds_speedy.sel(time=slice(startdate_climo.strftime("%Y-%m-%d"),enddate_climo.strftime("%Y-%m-%d")))
    print(ds_speedy)
    ds_speedy_annual_climo = ds_speedy.groupby('time.year').sum('time')['tp']
    ds_speedy_annual_climo = ds_speedy_annual_climo * 0.0393701
    print(np.shape(ds_speedy_annual_climo))

    ds_observed = get_era5_precip_timeseries(startdate_climo,enddate_climo,1)
    print(ds_observed)
    ds_observed_annual_climo = ds_observed.groupby('Timestep.year').sum('Timestep')['tp']
    ds_observed_annual_climo = ds_observed_annual_climo #* 39.3701
    print(ds_observed_annual_climo)

    ds_hybrid =  make_ds_time_dim(ds_hybrid,timestep,startdate)
    ds_hybrid = ds_hybrid.sel(Timestep=slice(startdate_climo.strftime("%Y-%m-%d"),enddate_climo.strftime("%Y-%m-%d")))
    ds_hybrid_annual_climo = ds_hybrid.groupby('Timestep.year').sum('Timestep')['p6hr']
    print(ds_hybrid_annual_climo)

    total_era5 = np.average(ds_observed_annual_climo.values, axis=0) * 39.37 #* 24 #ds_speedy_annual_climo[1,:,:]*365       #ds_observed_annual_climo[1,:,:]#ds_hybrid_annual_climo[1,:,:]
    total_speedy = np.average(ds_speedy_annual_climo.values, axis=0) * 365 * 0.0393701 * 2
    total_hybrid = np.average(ds_hybrid_annual_climo.values, axis=0)

    n_bins = 50

    bins = np.linspace(0,300,n_bins)

    sns.histplot(total_hybrid.flatten(), bins=bins,
                 multiple="layer",label='Hybrid',color='#377eb8',alpha=0.5)

    sns.histplot(total_speedy.flatten(), bins=bins,
                 multiple="layer",label='SPEEDY',color='#4daf4a',alpha=0.5)

    ist = sns.histplot(total_era5.flatten(), bins=bins,
                 multiple="layer",label='ERA5',color='#e41a1c',alpha=0.5)

    hist.set_ylabel("Cumulative Probability", fontsize = 14)
    hist.set_xlabel("Annual Precipitation (Inches)", fontsize = 14)

    plt.legend()

    plt.xlim([0,100])
    #plt.ylim([0,1])

    plt.show()
     
lat_slice = slice(-90,90)
lon_slice = slice(0,365)
sigma = 7

startdate = datetime(2000,1,1,0)
enddate = datetime(2011,1,1,0)

#ds_era5 = get_obs_era5_timeseries(startdate,enddate,6,lat_slice,lon_slice,sigma)
ds_era5 = get_6hr_precip_era5_timeseries(startdate,enddate,6,lat_slice,lon_slice)
#ds_era5 = sigma_to_p(ds_era5,'U-wind')
#ds_era5 = sigma_to_p(ds_era5,'V-wind')

ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_ml_only_ocean_beta_0.0001_sigma0.6_6diff_speedytrial_12_29_2006_00.nc') #hybrid_prediction_era6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_ocean_model_false_log_precip_noise_add_raw_epsilon_0.0005trial_12_31_1999_00.nc')

#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_full_test_climate_all_tisr_longertrial_12_31_1999_00.nc')
#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era5000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_2_vertlap_2_full_test_climate_all_tisr_longertrial_12_28_2009_00.nc')
#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era4000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_4_vertlap_2_full_test_climate_all_tisr_longertrial_12_28_2009_00.nc')

time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"))

ds_hybrid = make_ds_time_dim(ds_hybrid,6,startdate)

print(ds_hybrid)
#ds_hybrid = sigma_to_p(ds_hybrid,'U-wind')
#ds_hybrid = sigma_to_p(ds_hybrid,'V-wind')


#ds_speedy_sim = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/speedy_era_start_climo_sst_decade_sim12_31_1999_00.nc')
ds_speedy_sim = xr.open_dataset('/scratch/user/troyarcomano/letkf-hybrid-speedy/experiments/ERA_5_cycle_2_year_test/nature.nc')
ds_speedy = ds_speedy_sim


startdate_speedy = datetime(2000,1,1,0)
enddate_speedy = datetime(2010,1,1,23)
time_slice_speedy = slice(startdate_speedy.strftime("%Y-%m-%d"),enddate_speedy.strftime("%Y-%m-%d"))
#ds_speedy = make_ds_time_dim(ds_speedy_sim,6,startdate_speedy)
#ds_speedy = sigma_to_p(ds_speedy,'U-wind')
#ds_speedy = sigma_to_p(ds_speedy,'V-wind')

extreme_value_plot(ds_era5.sel(Lat=lat_slice,Lon=lon_slice,Timestep=time_slice),ds_hybrid.sel(Lat=lat_slice,Lon=lon_slice,Timestep=time_slice),ds_speedy.sel(lat=lat_slice,lon=lon_slice,time=time_slice_speedy))
