import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from netCDF4 import Dataset
import cartopy as cart
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.geoaxes import GeoAxes
import xarray as xr
import glob
from datetime import datetime, timedelta
from dateutil.relativedelta import *
from numba import jit
import calendar
from mpl_toolkits.axes_grid1 import AxesGrid
import seaborn as sns


def make_ds_time_dim(ds,timestep,startdate):
    begin_year_str = startdate.strftime("%Y-%m-%d")

    attrs = {"units": f"hours since {begin_year_str} "}

    ds = ds.assign_coords({"Timestep": ("Timestep", ds.Timestep.values*timestep, attrs)})
    ds = xr.decode_cf(ds)

    return ds

def climo_sst(start_year,end_year,lon_slice,region_slice):
    root_path = '/scratch/user/troyarcomano/ERA/'
    hours_in_year = 24*365
    xgrid = 96
    ygrid = 48
    average_sst = np.zeros((hours_in_year,ygrid,xgrid))

    for current_year in range(start_year,end_year + 1):
        ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{current_year}/era_5_y{current_year}_sst_regridded_mpi_fixed_var.nc')

        if current_year == start_year:
           shape = np.shape(ds_era['sst'].sel(Lon=lon_slice,Lat=region_slice))
           xgrid = shape[2]
           ygrid = shape[1]
           average_sst = np.zeros((hours_in_year,ygrid,xgrid))

        if calendar.isleap(current_year):
           time_slice = slice(0,240*6)
           average_sst[0:240*6,:,:] += ds_era['sst'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice)

           time_slice = slice(244*6,1464*6)
           average_sst[240*6:1460*6,:,:] += ds_era['sst'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice)
        else:
           time_slice = slice(0,hours_in_year)
           average_sst += ds_era['sst'].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice)

    return average_sst/(end_year-start_year+1)

def get_obs_sst_timeseries(startdate,enddate,timestep):
    start_year = startdate.year
    end_year = enddate.year

    currentdate = startdate
    counter = 0
    while currentdate.year <= enddate.year:
        ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{currentdate.year}/era_5_y{currentdate.year}_sst_regridded_mpi_fixed_var.nc')

        begin_year = datetime(currentdate.year,1,1,0)
        begin_year_str = begin_year.strftime("%Y-%m-%d")
        attrs = {"units": f"hours since {begin_year_str} "}
        ds_era = ds_era.assign_coords({"Timestep": ("Timestep", ds_era.Timestep.values, attrs)})
        ds_era = xr.decode_cf(ds_era)

        if start_year == currentdate.year:
           ds_merged = ds_era
        else:
           ds_merged = xr.merge([ds_merged,ds_era])

        currentdate = currentdate + timedelta(hours=ds_era.sizes['Timestep'])

    time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"),timestep)
    return ds_merged.sel(Timestep=time_slice)

def monthly_temp_movie(ds_era,ds_hybrid):
    monthly_era = ds_era.resample(Timestep="1MS").mean(dim="Timestep") 
    monthly_hybrid = ds_hybrid.resample(Timestep="1MS").mean(dim="Timestep") 
    monthly_hybrid = monthly_hybrid['SST']
   
    climatology = monthly_era.groupby("Timestep.month").mean("Timestep")

    monthly_era_anoms = monthly_era.groupby("Timestep.month") - climatology  
    monthly_hybrid_anoms = monthly_hybrid.groupby("Timestep.month") - climatology
  
    print(monthly_hybrid_anoms) 
    startdate = datetime(2007,1,1,0)

    for i in range(12*11):#12*10):
        #single_sst_map(monthly_era,monthly_hybrid,i)
        currentdate = startdate + relativedelta(months=i)
        date_string = currentdate.strftime("%b - %Y")
        single_sst_anom_map(monthly_era_anoms,monthly_hybrid_anoms,i,date_string)
        #single_sst_map(monthly_era,monthly_hybrid,i,date_string) 

def single_sst_anom_map(monthly_era,monthly_hybrid,timestep,date):

    #temp_colormap = sns.color_palette("Spectral", as_cmap=True)
    #temp_colormap = temp_colormap.reversed()
    temp_colormap = 'seismic'

    print(monthly_era)
    era_data = monthly_era['sst'].values
    era_data = era_data[timestep,:,:] #- 273.15

    print(monthly_hybrid)
    hybrid_data = monthly_hybrid['sst'].values
    hybrid_data = hybrid_data[timestep,:,:] #- 273.15

    lons = monthly_era.Lon.values
    lats = monthly_era.Lat.values

    clevs = np.arange(-5,5.25,0.25)

    cyclic_data, cyclic_lons = add_cyclic_point(hybrid_data, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
    projection = ccrs.PlateCarree()
    axes_class = (GeoAxes,dict(map_projection=projection))
    plt.rc('font', family='serif')

    fig = plt.figure(figsize=(10,6))

    plt.suptitle(date, fontsize=20)

    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(1, 2),
                    axes_pad=0.9,
                    share_all=True,
                    #cbar_location='right',
                    #cbar_mode='each',
                    cbar_pad=0.1,
                    cbar_size='3%',
                    cbar_location="bottom",
                    cbar_mode="single",
                    label_mode='')  # note the empty label_mode


    ax1 = axgr[0]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax1.coastlines()

    ax1.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax1.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
 
    ax1.xaxis.set_major_formatter(lon_formatter)
    ax1.yaxis.set_major_formatter(lat_formatter)
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')

    #ax1.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    #ax1.set_yticks(np.linspace(-80, 80, 10), crs=projection)
    #lon_formatter = LongitudeFormatter(zero_direction_label=True)
    #lat_formatter = LatitudeFormatter()
    #ax1.xaxis.set_major_formatter(lon_formatter)
    #ax1.yaxis.set_major_formatter(lat_formatter)

    plot = ax1.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels=clevs,cmap=temp_colormap,extend='both')

    cbar = axgr.cbar_axes[0].colorbar(plot, extend='both')
    #cbar.set_label_text(f'hPa',fontsize=16)
    cbar.set_label(f'C',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax1.add_feature(cart.feature.LAND, zorder=100, edgecolor='k')

    ax1.set_title(f"Hybrid Monthly Averaged\nSST Anomalies",fontsize=18,fontweight="bold")

    ###########

    cyclic_data, cyclic_lons = add_cyclic_point(era_data, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

    print(np.shape(cyclic_data))
    print(np.shape(lons2d))

    ax2 = axgr[1] #plt.subplot(222,projection=ccrs.PlateCarree())

    ax2.coastlines()

    ax2.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax2.xaxis.set_major_formatter(lon_formatter)
    ax2.yaxis.set_major_formatter(lat_formatter)
    gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')

    #ax2.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    #ax2.set_yticks(np.linspace(-80, 80, 10), crs=projection)
    #lon_formatter = LongitudeFormatter(zero_direction_label=True)
    #lat_formatter = LatitudeFormatter()
    #ax2.xaxis.set_major_formatter(lon_formatter)
    #ax2.yaxis.set_major_formatter(lat_formatter)

    plot_2 = ax2.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels=clevs,cmap=temp_colormap,extend='both')
    cbar2 = axgr.cbar_axes[1].colorbar(plot_2, extend='both')
    #cbar2.set_label_text(f'hPa',fontsize=16)
    cbar2.set_label(f'C',fontsize=16)
    cbar2.ax.tick_params(labelsize=16)

    ax2.set_title(f"ERA Monthly Averaged\nSST Anomalies",fontsize=18,fontweight="bold")
    ax2.add_feature(cart.feature.LAND, zorder=100, edgecolor='k')

    plt.tight_layout()
    #plt.show()
    plt.savefig(f'/home/troyarcomano/FortranReservoir/vert_loc_hybridspeedy_leakage/animation/anom_sst_map_t{timestep:03}')
    plt.close("all")


def single_sst_map(monthly_era,monthly_hybrid,timestep,date):

    temp_colormap = sns.color_palette("Spectral", as_cmap=True)
    temp_colormap = temp_colormap.reversed()
 
    print(monthly_era) 
    era_data = monthly_era['sst'].values
    era_data = era_data[timestep,:,:] - 273.15

    print(monthly_hybrid)
    hybrid_data = monthly_hybrid['SST'].values
    hybrid_data = hybrid_data[timestep,:,:] - 273.15

    lons = monthly_era.Lon.values
    lats = monthly_era.Lat.values

    clevs = np.arange(0,34,2)

    cyclic_data, cyclic_lons = add_cyclic_point(hybrid_data, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
    projection = ccrs.PlateCarree()
    axes_class = (GeoAxes,dict(map_projection=projection))
    plt.rc('font', family='serif')

    fig = plt.figure(figsize=(10,6))
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(1, 2),
                    axes_pad=0.9,
                    share_all=True,
                    #cbar_location='right',
                    #cbar_mode='each',
                    cbar_pad=0.1,
                    cbar_size='3%',
                    cbar_location="bottom",
                    cbar_mode="single",
                    label_mode='')  # note the empty label_mode


    plt.suptitle(date, fontsize=20)

    ax1 = axgr[0]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax1.coastlines()
 
    ax1.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax1.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax1.xaxis.set_major_formatter(lon_formatter)
    ax1.yaxis.set_major_formatter(lat_formatter)
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')

    #ax1.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    #ax1.set_yticks(np.linspace(-80, 80, 10), crs=projection)
    #lon_formatter = LongitudeFormatter(zero_direction_label=True)
    #lat_formatter = LatitudeFormatter()
    #ax1.xaxis.set_major_formatter(lon_formatter)
    #ax1.yaxis.set_major_formatter(lat_formatter)

    plot = ax1.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels=clevs,cmap=temp_colormap,extend='both')

    cbar = axgr.cbar_axes[0].colorbar(plot, extend='both')
    #cbar.set_label_text(f'hPa',fontsize=16)
    cbar.set_label(f'C',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax1.add_feature(cart.feature.LAND, zorder=100, edgecolor='k')

    ax1.set_title(f"Hybrid Monthly Averaged SST",fontsize=18,fontweight="bold")

    ###########

    cyclic_data, cyclic_lons = add_cyclic_point(era_data, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

    print(np.shape(cyclic_data))
    print(np.shape(lons2d))

    ax2 = axgr[1] #plt.subplot(222,projection=ccrs.PlateCarree())

    ax2.coastlines()

    ax2.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax2.xaxis.set_major_formatter(lon_formatter)
    ax2.yaxis.set_major_formatter(lat_formatter)
    gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')

    #ax2.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    #ax2.set_yticks(np.linspace(-80, 80, 10), crs=projection)
    #lon_formatter = LongitudeFormatter(zero_direction_label=True)
    #lat_formatter = LatitudeFormatter()
    #ax2.xaxis.set_major_formatter(lon_formatter)
    #ax2.yaxis.set_major_formatter(lat_formatter)

    plot_2 = ax2.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels=clevs,cmap=temp_colormap,extend='both')
    cbar2 = axgr.cbar_axes[1].colorbar(plot_2, extend='both')
    #cbar2.set_label_text(f'hPa',fontsize=16)
    cbar2.set_label(f'C',fontsize=16)
    cbar2.ax.tick_params(labelsize=16)

    ax2.set_title(f"ERA Monthly Averaged SST",fontsize=18,fontweight="bold")
    ax2.add_feature(cart.feature.LAND, zorder=100, edgecolor='k')

    plt.tight_layout()
    #plt.show()
    plt.savefig(f'/home/troyarcomano/FortranReservoir/vert_loc_hybridspeedy_leakage/animation/sst_map_t{timestep:03}')
 

#ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_slab_ocean_model_7d_0.1rho_10noise_beta0.001_20yearstrial_01_09_2010_00.nc')

ds_hybrid = xr.open_dataset('//scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_ml_only_ocean_beta_0.0001_sigma0.6_6diff_speedytrial_12_29_2006_00.nc')

startdate = datetime(2007,1,1,0)
enddate = datetime(2018,1,1,0)
timestep = 6
ds_observed = get_obs_sst_timeseries(startdate,enddate,timestep)

ds_hybrid =  make_ds_time_dim(ds_hybrid,6,startdate)
monthly_temp_movie(ds_observed,ds_hybrid)
