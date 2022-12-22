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

def sigma_to_pressure(ds,var):

    target_pressures = [25,95,200,350,500,680,850,950] 
    ds[var][:] = lin_interp(ds[var].values,ds['logp'].values,target_pressures)
   
    return ds 

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
           ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{current_year}/era_5_y{current_year}_regridded_mpi_fixed_var_gcc.nc')

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
           ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{currentdate.year}/era_5_y{currentdate.year}_regridded_mpi_fixed_var_gcc.nc')

        begin_year = datetime(currentdate.year,1,1,0)
        begin_year_str = begin_year.strftime("%Y-%m-%d")
        attrs = {"units": f"hours since {begin_year_str} "}
        ds_era = ds_era.assign_coords({"Timestep": ("Timestep", ds_era.Timestep.values, attrs)})
        ds_era = xr.decode_cf(ds_era)

        ds_era = ds_era.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma)

        if start_year == currentdate.year:
           ds_merged = ds_era
        else:
           #ds_merged = xr.merge([ds_merged,ds_era])
           ds_merged = xr.concat([ds_merged,ds_era],dim="Timestep") 

        currentdate = currentdate + timedelta(hours=ds_era.sizes['Timestep'])

    time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"),timestep)
    return ds_merged.sel(Timestep=time_slice)

def zonal_wind_mean_plot(ds_era,ds_hybrid,ds_speedy,ds_era_temp,ds_hybrid_temp,ds_speedy_temp):
    ds_zmean_era = ds_era.mean(dim='Lon')
    ds_zmean_era = ds_zmean_era.mean(dim='Lat') 
 
    ds_zmean_era_ndjfm = ds_zmean_era.where(ds_zmean_era.Timestep.dt.month.isin([1,2,3,11,12]))#[5,6,7,8,9,10]))#1,2,3,11,12]))
    plt.plot(ds_zmean_era_ndjfm['U-wind'].values)
    plt.show()
    ssw_events = ds_zmean_era_ndjfm.where(ds_zmean_era_ndjfm['U-wind'] < 0.0)
    print(np.shape(ssw_events['U-wind'].values)[0])
    time = np.arange(0,np.shape(ssw_events['U-wind'].values)[0],1)
    time = time/1460
    #plt.plot(ssw_events.Timestep,ssw_events['U-wind'].values)
    #plt.show()
    print(ssw_events) 
 
    mean_era = ds_zmean_era.groupby(ds_era.Timestep.dt.dayofyear).mean("Timestep")['U-wind']
    std_era = ds_zmean_era.groupby(ds_era.Timestep.dt.dayofyear).std("Timestep")['U-wind']

    ds_zmean_hybrid = ds_hybrid.mean(dim='Lon')
    ds_zmean_hybrid = ds_zmean_hybrid.mean(dim='Lat')

    ds_zmean_hybrid_ndjfm = ds_zmean_hybrid.where(ds_zmean_hybrid.Timestep.dt.month.isin([1,2,3,11,12]))#[5,6,7,8,9,10]))
    ssw_events = ds_zmean_hybrid_ndjfm.where(ds_zmean_hybrid_ndjfm['U-wind'] < 0.0)
    print(np.shape(ssw_events['U-wind'].values)[0])
    time = np.arange(0,np.shape(ssw_events['U-wind'].values)[0],1)
    time = time/1460
    #plt.plot(ssw_events.Timestep,ssw_events['U-wind'].values)
    plt.show()
    print(ssw_events)

    mean_hybrid = ds_zmean_hybrid.groupby(ds_hybrid.Timestep.dt.dayofyear).mean("Timestep")['U-wind']
    std_hybrid = ds_zmean_hybrid.groupby(ds_hybrid.Timestep.dt.dayofyear).std("Timestep")['U-wind']

    ds_zmean_speedy = ds_speedy.mean(dim='Lon')
    ds_zmean_speedy = ds_zmean_speedy.mean(dim='Lat')

    mean_speedy = ds_zmean_speedy.groupby(ds_speedy.Timestep.dt.dayofyear).mean("Timestep")['U-wind']
    std_speedy = ds_zmean_speedy.groupby(ds_speedy.Timestep.dt.dayofyear).std("Timestep")['U-wind']

    timestep = mean_era.dayofyear.values
    timestep = timestep[0:365]

    mean_era = np.roll(mean_era[0:365],-182)
    mean_hybrid = np.roll(mean_hybrid[0:365],-182)
    mean_speedy = np.roll(mean_speedy[0:365],-182)
    
    std_era = np.roll(std_era[0:365],-182)
    std_hybrid = np.roll(std_hybrid[0:365],-182)
    std_speedy = np.roll(std_speedy[0:365],-182)

    ###Temperature### 

    ds_zmean_era_temp = ds_era_temp.mean(dim='Lon')
    ds_zmean_era_temp = ds_zmean_era_temp.mean(dim='Lat')

    ds_zmean_hybrid_temp = ds_hybrid_temp.mean(dim='Lon')
    ds_zmean_hybrid_temp = ds_zmean_hybrid_temp.mean(dim='Lat')

    ds_zmean_speedy_temp = ds_speedy_temp.mean(dim='Lon')
    ds_zmean_speedy_temp = ds_zmean_speedy_temp.mean(dim='Lat')


    mean_hybrid_temp = ds_zmean_hybrid_temp.groupby(ds_hybrid.Timestep.dt.dayofyear).mean("Timestep")['Temperature']
    std_hybrid_temp = ds_zmean_hybrid_temp.groupby(ds_hybrid.Timestep.dt.dayofyear).std("Timestep")['Temperature']
    
    mean_speedy_temp = ds_zmean_speedy_temp.groupby(ds_speedy.Timestep.dt.dayofyear).mean("Timestep")['Temperature']
    std_speedy_temp = ds_zmean_speedy_temp.groupby(ds_speedy.Timestep.dt.dayofyear).std("Timestep")['Temperature']

    mean_era_temp = ds_zmean_era_temp.groupby(ds_era.Timestep.dt.dayofyear).mean("Timestep")['Temperature']
    std_era_temp = ds_zmean_era_temp.groupby(ds_era.Timestep.dt.dayofyear).std("Timestep")['Temperature']

    mean_era_temp = np.roll(mean_era_temp[0:365],-182)
    mean_hybrid_temp = np.roll(mean_hybrid_temp[0:365],-182)
    mean_speedy_temp = np.roll(mean_speedy_temp[0:365],-182)

    std_era_temp = np.roll(std_era_temp[0:365],-180)
    std_hybrid_temp = np.roll(std_hybrid_temp[0:365],-182)
    std_speedy_temp = np.roll(std_speedy_temp[0:365],-182)

   
    ###Example ERA
    start_date_era = datetime(2012,7,1,0)
    end_date_era = datetime(2013,6,30,23)
    time_slice_era = slice(start_date_era,end_date_era)

    resample_1d_era = ds_zmean_era.resample(Timestep="1D").mean("Timestep")
    era_2013 = resample_1d_era.sel(Timestep=time_slice_era)

    resample_1d_era_temp = ds_zmean_era_temp.resample(Timestep="1D").mean("Timestep")
    era_2013_temp = resample_1d_era_temp.sel(Timestep=time_slice_era)

    ###Example Hybrid
    start_date_hybrid = datetime(2018,7,1,0)
    end_date_hybrid = datetime(2019,6,30,23)
    time_slice_hybrid = slice(start_date_hybrid,end_date_hybrid)

    resample_1d_hybrid = ds_zmean_hybrid.resample(Timestep="1D").mean("Timestep")
    hybrid_2024 = resample_1d_hybrid.sel(Timestep=time_slice_hybrid)

    resample_1d_hybrid_temp = ds_zmean_hybrid_temp.resample(Timestep="1D").mean("Timestep")
    hybrid_2024_temp = resample_1d_hybrid_temp.sel(Timestep=time_slice_hybrid)

    print('shape(era_2013)',np.shape(era_2013['Temperature'].values))
  
    x_labels_vals = [15,45,76,107,138,168,199,229,257,287,317,348]
    x_tick_labels = ['July','Aug','Sept','Oct','Nov','Dec','Jan','Feb','March','April','May','June']

    fig, axes = plt.subplots(2,3,figsize=(18,10),sharex=False,sharey='row',constrained_layout=True)

    ####
    #ERA Strato
    ####
    ax = axes[0,0]

    print(mean_era,std_era)
    ax.plot(timestep,mean_era)
    ax.plot(timestep,era_2013['U-wind'].values,color='r')
    ax.fill_between(timestep,mean_era+std_era*2,mean_era-std_era*2, facecolor='grey', alpha=0.4) 

    ax.axhline(y = 0.0, color = 'k', linestyle = '--',linewidth=2)
    ax.axvline(x = 182, color = 'k')
    ax.axvline(x = 245, color = 'k')

    ax.text(5, 15, 'westerly\nflow',fontsize=18)
    ax.text(5, -12, 'easterly flow',fontsize=18)
    ax.text(185, 50, 'SSW',fontsize=18,color='r')
   
    ax.set_title('ERA5',fontsize=20,fontweight='bold')
    ax.set_xlim([np.min(timestep),np.max(timestep)]) 
    ax.set_ylim([-20,60])
    #ax.set_ylim([-30,80])
    #ax.set_xlabel('Day of Year')
    ax.set_ylabel('25 hPa Zonal Mean\nWind (m/s)',fontsize=16)

    ax.set_xticks(x_labels_vals)
    ax.set_xticklabels(x_tick_labels, rotation=30)

  
    ####
    #Hybrid Strato
    ####
    ax = axes[0,1]

    ax.plot(timestep,mean_hybrid)
    ax.plot(timestep,hybrid_2024['U-wind'].values,color='r')
    ax.fill_between(timestep,mean_hybrid+std_hybrid*2,mean_hybrid-std_hybrid*2, facecolor='grey', alpha=0.4) 

    ax.axhline(y = 0.0, color = 'k', linestyle = '--',linewidth=2)
    ax.axvline(x = 205, color = 'k')
    ax.axvline(x = 255, color = 'k')

    ax.text(5, 15, 'westerly\nflow',fontsize=18)
    ax.text(5, -12, 'easterly flow',fontsize=18)
    ax.text(206, 50, 'SSW',fontsize=18,color='r')

    ax.set_title('Hybrid Model',fontsize=20,fontweight='bold')
    ax.set_xlim([np.min(timestep),np.max(timestep)])
    ax.set_ylim([-20,60])
    #ax.set_ylim([-30,80])    
    #ax.set_xlabel('Day of Year')
    #ax.set_ylabel('Hybrid 20 hPa Zonal Mean Wind (2000-2011)')
    ax.set_xticks(x_labels_vals)
    ax.set_xticklabels(x_tick_labels, rotation=30)

    ####
    #SPEEDY Strato
    ####
    ax = axes[0,2]

    plot = ax.plot(timestep,mean_speedy,label='Mean Wind')
    fill = ax.fill_between(timestep,mean_speedy+std_speedy*2,mean_speedy-std_speedy*2, facecolor='grey', alpha=0.4,label='2 Standard Deviations')
    ax.axhline(y = 0.0, color = 'k', linestyle = '--',linewidth=2) 
    plot2 = ax.plot(timestep,mean_speedy+1000,label='SSW Event',color='r')
    ax.text(5, 15, 'westerly\nflow',fontsize=18)
    ax.text(5, -12, 'easterly flow',fontsize=18)

    ax.legend()
    ax.set_title('SPEEDY',fontsize=20,fontweight='bold')
    ax.set_xlim([np.min(timestep),np.max(timestep)])
    ax.set_ylim([-20,60])
    #ax.set_ylim([-30,80])
    #ax.set_xlabel('Day of Year')
    #ax.set_ylabel('SPEEDY 20 hPa Zonal Mean Wind (2000-2011)')
    ax.set_xticks(x_labels_vals)
    ax.set_xticklabels(x_tick_labels, rotation=30)


 
    ####
    #ERA Strato
    ####

    
    ax = axes[1,0]

    ax.plot(timestep,mean_era_temp)
    ax.plot(timestep,era_2013_temp['Temperature'].values,color='r')
    ax.fill_between(timestep,mean_era_temp+std_era_temp*2,mean_era_temp-std_era_temp*2, facecolor='grey', alpha=0.4)
    ax.axvline(x = 182, color = 'k')
    ax.axvline(x = 245, color = 'k')

    ax.text(185, 50, 'SSW',fontsize=18,color='r')

    #ax.set_title('ERA5',fontsize=20,fontweight='bold')
    ax.set_xlim([np.min(timestep),np.max(timestep)])
    ax.set_ylim([185,235])
    #ax.set_ylim([-30,80])
    #ax.set_xlabel('Day of Year')
    ax.set_ylabel('25 hPa Zonal Mean\nTemperature (K)',fontsize=16)

    ax.set_xticks(x_labels_vals)
    ax.set_xticklabels(x_tick_labels, rotation=30)


    ####
    #Hybrid Strato
    ####
    ax = axes[1,1]

    ax.plot(timestep,mean_hybrid_temp)
    ax.plot(timestep,hybrid_2024_temp['Temperature'].values,color='r')
    ax.fill_between(timestep,mean_hybrid_temp+std_hybrid_temp*2,mean_hybrid_temp-std_hybrid_temp*2, facecolor='grey', alpha=0.4)
    ax.axvline(x = 205, color = 'k')
    ax.axvline(x = 255, color = 'k')
    ax.text(206, 50, 'SSW',fontsize=18,color='r')

    #ax.set_title('Hybrid Model',fontsize=20,fontweight='bold')
    ax.set_xlim([np.min(timestep),np.max(timestep)])
    ax.set_ylim([185,235])
    #ax.set_ylim([-30,80])
    #ax.set_xlabel('Day of Year')
    #ax.set_ylabel('Hybrid 20 hPa Zonal Mean Wind (2000-2011)')
    ax.set_xticks(x_labels_vals)
    ax.set_xticklabels(x_tick_labels, rotation=30)

 
    ####
    #SPEEDY Strato
    ####
    ax = axes[1,2]

    plot = ax.plot(timestep,mean_speedy_temp,label='Mean Temp')
    fill = ax.fill_between(timestep,mean_speedy_temp+std_speedy_temp*2,mean_speedy_temp - std_speedy_temp*2, facecolor='grey', alpha=0.4,label='2 Standard Deviations')
    lot2 = ax.plot(timestep,mean_speedy_temp+1000,label='SSW Event',color='r')

    ax.legend()
    #ax.set_title('SPEEDY',fontsize=20,fontweight='bold')
    ax.set_xlim([np.min(timestep),np.max(timestep)])
    ax.set_ylim([185,235])
    #ax.set_ylim([-30,80])
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

def ssw_stuff():
    lat_slice_uwind = slice(55,65) #-65,-55)
    lat_slice_temp = slice(60,90)
    lat_slice = slice(55,90)

    lon_slice = slice(0,365)
    sigma_all = slice(0,8)
    sigma = 0

    startdate = datetime(1981,1,1,0)
    #startdate = datetime(2010,1,1,0)
    enddate = datetime(2017,12,31,0)
    #enddate = datetime(2015,12,31,0)

    ds_era5 = get_obs_era5_timeseries(startdate,enddate,6,lat_slice,lon_slice,sigma_all)

    ds_era5 = sigma_to_pressure(ds_era5,'Temperature')
    ds_era5 = sigma_to_pressure(ds_era5,'U-wind')

    ds_era5 = ds_era5.sel(Sigma_Level=sigma)

    ds_era5_uwind = ds_era5.sel(Lat=lat_slice_uwind)
    ds_era5_temp = ds_era5.sel(Lat=lat_slice_temp)

    print(ds_era5.Timestep)

    ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_sigma0.5_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_multi_gaussian_noise_newest_version_32_processors_root_ssttrial_12_29_2006_00.nc')

    startdate = datetime(2000,1,1,0)
    enddate = datetime(2045,1,1,0)

    ds_speedy_sim = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/speedy_1981_2020_era5_ics_6hourly_output.nc')#'/scratch/user/troyarcomano/temp_storage/climo_speedy/speedy_era_start_climo_sst_decade_sim12_31_1999_00.nc')

    time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"))

    ds_hybrid = make_ds_time_dim(ds_hybrid,6,startdate)

    startdate_speedy = datetime(2000,1,1,0)
    enddate_speedy = datetime(2040,1,1,23)

    time_slice_speedy = slice(startdate_speedy.strftime("%Y-%m-%d"),enddate_speedy.strftime("%Y-%m-%d"))
    ds_speedy = make_ds_time_dim(ds_speedy_sim,6,startdate_speedy)

    ds_speedy = sigma_to_pressure(ds_speedy,'Temperature')
    ds_speedy = sigma_to_pressure(ds_speedy,'U-wind')

    ds_hybrid = ds_hybrid.sel(Timestep=time_slice)

    ds_hybrid = sigma_to_pressure(ds_hybrid,'Temperature')
    ds_hybrid = sigma_to_pressure(ds_hybrid,'U-wind')

    print('shape(ds_hybrid)',np.shape(ds_hybrid['U-wind'].values))

    ds_hybrid_uwind = ds_hybrid.sel(Lat=lat_slice_uwind,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice)
    ds_speedy_uwind = ds_speedy.sel(Lat=lat_slice_uwind,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice_speedy)

    ds_hybrid_temp = ds_hybrid.sel(Lat=lat_slice_temp,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice)
    ds_speedy_temp = ds_speedy.sel(Lat=lat_slice_temp,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice_speedy)

    zonal_wind_mean_plot(ds_era5_uwind,ds_hybrid_uwind,ds_speedy_uwind,ds_era5_temp,ds_hybrid_temp,ds_speedy_temp)


'''
lat_slice = slice(-45,45)
lon_slice = slice(0,365)
sigma = 0

startdate = datetime(2007,1,1,0)
enddate = datetime(2018,12,31,0)

ds_era5 = get_obs_era5_timeseries(startdate,enddate,6,lat_slice,lon_slice,sigma)

#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era9000_15_30_30_beta_res0.1_beta_model_1_prior_0.0_overlap1_era5_6hrtimestep_sst_true_climo_inputtrial_01_02_2000_00.nc')
#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era5000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_2_vertlap_2_full_test_climate_all_tisr_longertrial_12_28_2009_00.nc')
ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_sigma0.5_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_multi_gaussian_noise_newest_version_32_processors_root_ssttrial_12_29_2006_00.nc')
#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era4000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_4_vertlap_2_full_test_climate_all_tisr_longertrial_12_28_2009_00.nc')

time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"))

ds_hybrid = make_ds_time_dim(ds_hybrid,6,startdate)


ds_speedy_sim = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/climo_speedy/speedy_era_start_climo_sst_decade_sim12_31_1999_00.nc')

#startdate_speedy = datetime(2000,1,1,0)
#enddate_speedy = datetime(2012,1,1,23)
startdate_speedy = datetime(2007,1,1,0)
enddate_speedy = datetime(2018,12,31,23)

time_slice_speedy = slice(startdate_speedy.strftime("%Y-%m-%d"),enddate_speedy.strftime("%Y-%m-%d"))
ds_speedy = make_ds_time_dim(ds_speedy_sim,6,startdate_speedy)

#ds_old_hybrid_sim = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.01_dp_beta_model_1_prior_0.0_overlap1_era5_6hrtimestep_tisr_full_years_0.00noise_speedy_states_climate_simtrial_12_31_1999_00.nc')
#startdate_speedy = datetime(2010,1,1,0)
#enddate_speedy = datetime(2018,12,31,23)
#time_slice_speedy = slice(startdate_speedy.strftime("%Y-%m-%d"),enddate_speedy.strftime("%Y-%m-%d"))
#ds_old_hybrid = make_ds_time_dim(ds_old_hybrid_sim,6,startdate_speedy)

#qbo_plot_darpa_talk(ds_era5,ds_hybrid.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice),ds_old_hybrid.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice),ds_speedy.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice_speedy))


qbo_plot(ds_era5['U-wind'],ds_hybrid.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice),ds_speedy.sel(Lat=lat_slice,Lon=lon_slice,Sigma_Level=sigma,Timestep=time_slice_speedy))
'''
ssw_stuff()
