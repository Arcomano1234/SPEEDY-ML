import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
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
from numba import jit
import calendar
from mpl_toolkits.axes_grid1 import AxesGrid
import fiona
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.prepared import prep
from scipy.ndimage.filters import uniform_filter1d
from scipy.signal import find_peaks, welch
import seaborn as sns

import pycwt as wavelet
from pycwt.helpers import find

#geoms = fiona.open(shpreader.natural_earth(resolution='10m', category='physical', name='land'))
#land_geom = sgeom.MultiPolygon([sgeom.shape(geom['geometry']) for geom in geoms])
#land = prep(land_geom)

#def is_land(lon,lat#):
#    return land.contains(sgeom.Point(lon, lat)) 

def get_e3sm_diff_colormap_greyer():
    e3sm_colors_grayer = [
        '#1b1c40',
        '#2a41a1',
        '#237cbc',
        '#75a9be',
        '#cad3d9',
        '#f1eceb',
        '#e5ccc4',
        '#d08b73',
        '#b9482d',
        '#860e27',
        '#3d0712'
    ]

    cmap = mpl.colors.ListedColormap(e3sm_colors_grayer)
    return cmap

def get_e3sm_diff_colormap():
    file_cm = 'e3sm.rgb'
    rgb_arr = np.loadtxt(file_cm)
    rgb_arr = rgb_arr / 255.0

    colormap = 'test_cmap'
    cmap = LinearSegmentedColormap.from_list(name=colormap, colors=rgb_arr)
    return cmap

@jit()
def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

@jit()
def rms(true,prediction):
    return np.sqrt(np.nanmean((prediction-true)**2))

@jit()
def anomly_cc(forecast,observed,climo):
    top = np.mean((forecast-climo)*(observed-climo))
    bottom = np.sqrt(np.mean((observed-climo)**2)*np.mean((observed-climo)**2))
    ACC = top/bottom
    return ACC

def climo_sst(start_year,end_year,lon_slice,region_slice):
    root_path = '/scratch/user/troyarcomano/ERA/'
    hours_in_year = 24*365
    xgrid = 96
    ygrid = 48
    average_sst = np.zeros((hours_in_year,ygrid,xgrid)) 

    for current_year in range(start_year,end_year + 1):
        ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{current_year}/era_5_y{current_year}_sst_regridded_fixed_var_gcc.nc')

        if current_year == start_year:
           shape = np.shape(ds_era['sst'].sel(lon=lon_slice,lat=region_slice))
           xgrid = shape[2]
           ygrid = shape[1]
           average_sst = np.zeros((hours_in_year,ygrid,xgrid)) 

        if calendar.isleap(current_year):
           startdate = datetime(current_year,1,1,0)
           enddate = datetime(current_year,2,29,0)
           time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"))
           average_sst[0:240*6,:,:] += ds_era['sst'].sel(time=time_slice,lon=lon_slice,lat=region_slice).values

           startdate = datetime(current_year,3,2,0)
           enddate = datetime(current_year,12,31,23)
           time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"))
           average_sst[240*6:1460*6,:,:] += ds_era['sst'].sel(time=time_slice,lon=lon_slice,lat=region_slice).values
        else:
           startdate = datetime(current_year,1,1,0)
           enddate = datetime(current_year,12,31,23)
           time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"))
           average_sst += ds_era['sst'].sel(time=time_slice,lon=lon_slice,lat=region_slice).values

    return average_sst/(end_year-start_year+1)

def climo_atmo_var(start_year,end_year,lon_slice,region_slice,varname):
    root_path = '/scratch/user/troyarcomano/ERA/'
    hours_in_year = 24*365
    xgrid = 96
    ygrid = 48
    average_var = np.zeros((hours_in_year,ygrid,xgrid))

    for current_year in range(start_year,end_year + 1):
        ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{current_year}/era_5_y{current_year}_regridded_mpi_fixed_var_gcc.nc')

        if current_year == start_year:
           shape = np.shape(ds_era[varname].sel(Lon=lon_slice,Lat=region_slice))
           xgrid = shape[2]
           ygrid = shape[1]
           average_var = np.zeros((hours_in_year,ygrid,xgrid))

        if calendar.isleap(current_year):
           time_slice = slice(0,240*6)
           average_var[0:240*6,:,:] += ds_era[varname].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice)

           time_slice = slice(244*6,1464*6)
           average_var[240*6:1460*6,:,:] += ds_era[varname].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice)
        else:
           time_slice = slice(0,hours_in_year)
           average_var += ds_era[varname].sel(Timestep=time_slice,Lon=lon_slice,Lat=region_slice)

    return average_var/(end_year-start_year+1)

def nino_index(sst_data,region):
    if region == "1+2":
       lat_slice = slice(-10,0)
       lon_slice = slice(360-90,360-80)
    elif region == "3":
       lat_slice = slice(-5,5)
       lon_slice = slice(360-150,360-90)
    elif region == "3.4" :
       lat_slice = slice(-5,5)
       lon_slice = slice(360-170,360-120)
    elif region == "4" :
       lat_slice = slice(-5,5)
       lon_slice = slice(360-200,360-150)
  
    start_year = 1981
    end_year = 2010

    try: 
      climate_data = climo_sst(start_year,end_year,lon_slice,lat_slice)
    
      sst_data = sst_data.sel(Lon=lon_slice,Lat=lat_slice)
     
      sst_anom = np.zeros((sst_data.sizes['Timestep'],sst_data.sizes['Lat'],sst_data.sizes['Lon'] ))
    except: 
      climate_data = climo_sst(start_year,end_year,lon_slice,lat_slice)

      sst_data = sst_data.sel(lon=lon_slice,lat=lat_slice)

      sst_anom = np.zeros((sst_data.sizes['time'],sst_data.sizes['lat'],sst_data.sizes['lon']))

    size = np.shape(climate_data)[0]
    print(np.shape(climate_data))
    print(size)
    try:
       for i in range(sst_data.sizes['Timestep']-1): 
           sst_anom[i,:,:] = sst_data[i,:,:] - climate_data[(i*6)%size,:,:]
    except:
       for i in range(sst_data.sizes['time']-1):
           sst_anom[i,:,:] = sst_data[i,:,:] - climate_data[(i*6)%size,:,:]
 
    return np.mean(sst_anom,axis=(1,2))

def nino_index_monthly(sst_data,region):
    if region == "1+2":
       lat_slice = slice(-10,0)
       lon_slice = slice(360-90,360-80)
    elif region == "3":
       lat_slice = slice(-5,5)
       lon_slice = slice(360-150,360-90)
    elif region == "3.4" :
       lat_slice = slice(-5,5)
       lon_slice = slice(360-170,360-120)
    elif region == "4" :
       lat_slice = slice(-5,5)
       lon_slice = slice(360-200,360-150)

    start_year = 1981
    end_year = 2010
   
    try:
       sst_data_copy = sst_data.sel(Lon=lon_slice,Lat=lat_slice)

       sst_data_copy = sst_data_copy.resample(Timestep="1MS").mean(dim="Timestep")
       gb = sst_data_copy.groupby('Timestep.month')
       tos_nino34_anom = gb - gb.mean(dim='Timestep')
       index_nino34 = tos_nino34_anom.mean(dim=['Lat', 'Lon'])
    except:
       sst_data_copy = sst_data.sel(lon=lon_slice,lat=lat_slice)
       sst_data_copy = sst_data_copy.resample(time="1MS").mean(dim="time")
       gb = sst_data_copy.groupby('time.month')
       tos_nino34_anom = gb - gb.mean(dim='time')
       index_nino34 = tos_nino34_anom.mean(dim=['lat', 'lon'])
    return index_nino34 

def nino_index_speedy(sst_data,startdate,enddate,region):
    if region == "1+2":
       lat_slice = slice(-10,0)
       lon_slice = slice(360-90,360-80)
    elif region == "3":
       lat_slice = slice(-5,5)
       lon_slice = slice(360-150,360-90)
    elif region == "3.4" :
       lat_slice = slice(-5,5)
       lon_slice = slice(360-170,360-120)
    elif region == "4" :
       lat_slice = slice(-5,5)
       lon_slice = slice(360-200,360-150)

    climate_start_year = 1981
    climate_end_year = 2000

    climate_data = climo_sst(climate_start_year,climate_end_year,lon_slice,lat_slice)

    sst_data = sst_data.sel(time=slice(startdate,enddate),lon=lon_slice,lat=lat_slice)
    print(sst_data)

    sst_anom = np.zeros((sst_data.sizes['time'],sst_data.sizes['lat'],sst_data.sizes['lon'] ))

    size = np.shape(climate_data)[0]
    for i in range(sst_data.sizes['time']-1):
        sst_anom[i,:,:] = sst_data[i,:,:].values - climate_data[(i*730)%size,:,:]

    return np.mean(sst_anom,axis=(1,2))

def get_obs_atmo_timeseries(startdate,enddate,timestep):
    start_year = startdate.year
    end_year = enddate.year

    currentdate = startdate
    counter = 0
    while currentdate.year <= enddate.year:
        ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{currentdate.year}/era_5_y{currentdate.year}_regridded_mpi_fixed_var_gcc.nc')

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

def get_obs_atmo_timeseries_var(startdate,enddate,timestep,var):
    start_year = startdate.year
    end_year = enddate.year

    currentdate = startdate
    counter = 0
    while currentdate.year <= enddate.year:
        print(currentdate.year)
        ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{currentdate.year}/era_5_y{currentdate.year}_regridded_mpi_fixed_var_gcc.nc')

        begin_year = datetime(currentdate.year,1,1,0)
        begin_year_str = begin_year.strftime("%Y-%m-%d")
        attrs = {"units": f"hours since {begin_year_str} "}
        ds_era = ds_era.assign_coords({"Timestep": ("Timestep", ds_era.Timestep.values, attrs)})
        ds_era = xr.decode_cf(ds_era)

        ds_era = ds_era[var]

        if start_year == currentdate.year:
           ds_merged = ds_era
        else:
           ds_merged = xr.merge([ds_merged,ds_era])

        currentdate = currentdate + timedelta(hours=ds_era.sizes['Timestep'])

    time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"),timestep)
    return ds_merged.sel(Timestep=time_slice)

def southern_oscillation_index(ds_hybrid,ds_era,startdate,enddate):
    darwin_lat, darwin_long = -12.4637, 130.8444
    tahiti_lat, tahiti_long = -17.6509, 210.426

    
    lon_slice = slice(darwin_long - 10, tahiti_long + 10)
    lat_slice = slice(tahiti_lat - 5, darwin_lat + 5)

    climate_start = datetime(1981,1,1,0)
    climate_end = datetime(2006,12,31,23)

    climo_time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"))
   
    data_tahiti = ds_era.sel(Lat = tahiti_lat, Lon = tahiti_long, method = 'nearest')['logp']
    data_darwin = ds_era.sel(Lat = darwin_lat, Lon = darwin_long, method = 'nearest')['logp']
   
    #data_tahiti = data_tahiti.resample(Timestep="1MS").mean("Timestep")
    #data_darwin = data_darwin.resample(Timestep="1MS").std("Timestep")

    data_tahiti = np.exp(data_tahiti)*1000.0
    data_darwin = np.exp(data_darwin)*1000.0

    diff_data = data_tahiti - data_darwin
   
    diff_climo = diff_data.sel(Timestep=climo_time_slice)

    climatology_mean = diff_climo.groupby("Timestep.month").mean("Timestep")
    climatology_std  = diff_climo.groupby("Timestep.month").std("Timestep")

    print('climatology_std',climatology_std)
    print('climatology_mean',climatology_mean)

    print('ds_hybrid',ds_hybrid)
    hybrid_data_tahiti = ds_hybrid.sel(Lat = tahiti_lat, Lon = tahiti_long, method = 'nearest')['logp']
    hybrid_data_darwin = ds_hybrid.sel(Lat = darwin_lat, Lon = darwin_long, method = 'nearest')['logp']

    #hybrid_data_tahiti = hybrid_data_tahiti.resample(Timestep="1MS")
    #hybrid_data_darwin = hybrid_data_darwin.resample(Timestep="1MS")

    hybrid_data_tahiti = np.exp(hybrid_data_tahiti)*1000.0
    hybrid_data_darwin = np.exp(hybrid_data_darwin)*1000.0

    hybrid_diff_data = hybrid_data_tahiti - hybrid_data_darwin

    hybrid_soi = xr.apply_ufunc(
        lambda x, m, s: (x - m) / s,
        hybrid_diff_data.groupby("Timestep.month"),
        climatology_mean,
        climatology_std,
    )

    print('hybrid_soi',hybrid_soi)
    hybrid_soi = hybrid_soi * 10.0
    #hybrid_soi = 10 * (hybrid_diff_data.groupby("Timestep.month") - climatology_mean)/climatology_std
  
    print('np.shape(hybrid_soi)',np.shape(hybrid_soi)) 

    era_soi = xr.apply_ufunc(
        lambda x, m, s: (x - m) / s,
        diff_climo.groupby("Timestep.month"),
        climatology_mean,
        climatology_std,
    )
    era_soi = era_soi * 10.0
    return hybrid_soi

def get_obs_sst_timeseries(startdate,enddate,timestep):
    start_year = startdate.year
    end_year = enddate.year

    currentdate = startdate
    counter = 0
    while currentdate.year <= enddate.year:
        print('year',currentdate.year)
        ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{currentdate.year}/era_5_y{currentdate.year}_sst_regridded_fixed_var_gcc.nc') 

        #begin_year = datetime(currentdate.year,1,1,0)
        #begin_year_str = begin_year.strftime("%Y-%m-%d")
        #attrs = {"units": f"hours since {begin_year_str} "}
        #ds_era = ds_era.assign_coords({"Timestep": ("Timestep", ds_era.time.values, attrs)})
        #ds_era = xr.decode_cf(ds_era)

        if start_year == currentdate.year:
           ds_merged = ds_era
        else:    
           ds_merged = xr.merge([ds_merged,ds_era]) 

        currentdate = currentdate + timedelta(hours=ds_era.sizes['time'])

    time_slice = slice(startdate.strftime("%Y-%m-%d"),enddate.strftime("%Y-%m-%d"),timestep)
    return ds_merged.sel(time=time_slice)

def sst_annual_variability(ds_model,ds_era,ds_speedy):
    lat_slice = slice(-90,90)
    lon_slice = slice(0,360)

    startdate_era = datetime(1981,1,1,0)
    enddate_era = datetime(2020,12,31,23)

    startdate_hybrid = datetime(2007,1,1,0)
    enddate_hybrid = datetime(2047,12,31,23)

    startdate_speedy = datetime(1980,1,1,0)
    enddate_speedy = datetime(2008,12,31,23)

    ds_model = ds_model.sel(Lon=lon_slice,Lat=lat_slice,Timestep=slice(startdate_hybrid.strftime("%Y-%m-%d"),enddate_hybrid.strftime("%Y-%m-%d")))['SST']
    ds_era = ds_era.sel(lon=lon_slice,lat=lat_slice,time=slice(startdate_era.strftime("%Y-%m-%d"),enddate_era.strftime("%Y-%m-%d")))['sst']
    ds_speedy = ds_speedy.sel(lon=lon_slice,lat=lat_slice,time=slice(startdate_speedy.strftime("%Y-%m-%d"),enddate_speedy.strftime("%Y-%m-%d")))['sst']
  
    ds_model = ds_model.groupby("Timestep.month").std("Timestep")
    ds_era = ds_era.groupby("time.month").std("time")
    ds_speedy = ds_speedy.groupby("time.month").std("time")

    lons = ds_model.Lon.values 
    lats = ds_model.Lat.values
    print('lons',lons)

    mean_hybrid = np.mean(ds_model.values,axis=(0))
    mean_climo = np.mean(ds_era.values,axis=(0))
    mean_speedy = np.mean(ds_speedy.values,axis=(0))

    projection = ccrs.PlateCarree(central_longitude=-179)
    axes_class = (GeoAxes,dict(map_projection=projection))
    plt.rc('font', family='serif')
    plt.rcParams['figure.constrained_layout.use'] = True


    cmap_e3sm_diff_colormap_greyer = get_e3sm_diff_colormap_greyer()

    fig = plt.figure(figsize=(6,10))
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(3, 1),
                    axes_pad=0.7,
                    cbar_location='right',
                    cbar_mode='each',
                    cbar_pad=0.2,
                    cbar_size='3%',
                    label_mode='')  # note the empty label_mode


    cyclic_data, cyclic_lons = add_cyclic_point(mean_hybrid, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

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

    levels = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6]
    temp_colormap = sns.color_palette("Spectral", as_cmap=True)
    temp_colormap = temp_colormap.reversed()
    cmap = temp_colormap
    plot = ax1.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[0].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax1.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax1.set_title(f"Coupled Model SST Standard Deviation of Monthly Means",fontsize=18,fontweight="bold")
    ax1.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    cyclic_data, cyclic_lons = add_cyclic_point(mean_climo, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
    ax2 = axgr[1]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax2.coastlines()

    ax2.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax2.xaxis.set_major_formatter(lon_formatter)
    ax2.yaxis.set_major_formatter(lat_formatter)
    gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                        linewidth=2, color='gray', alpha=0.5, linestyle='--')

    levels = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6]
    plot = ax2.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[1].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax2.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax2.set_title(f"ERA5 SST Standard Deviation",fontsize=18,fontweight="bold")
    ax2.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    #### Bias
    cyclic_data, cyclic_lons = add_cyclic_point(mean_hybrid - mean_climo, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
    ax3 = axgr[2]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax3.coastlines()
   
    ax3.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax3.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax3.xaxis.set_major_formatter(lon_formatter)
    ax3.yaxis.set_major_formatter(lat_formatter)
    gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                       linewidth=2, color='gray', alpha=0.5, linestyle='--')

    levels = [-0.6,-0.4,-0.2,-0.05,0.05,0.2,0.4,0.6]

    plot = ax3.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap_e3sm_diff_colormap_greyer,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[2].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax3.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax3.set_title(f"Standard Deviation Bias (Model - ERA5)",fontsize=18,fontweight="bold")
    ax3.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    plt.show()

def sst_climatology_error(ds_model,ds_era,ds_speedy):  
    lat_slice = slice(-90,90)
    lon_slice = slice(0,360) 

    startdate_era = datetime(1981,1,1,0)
    enddate_era = datetime(2020,12,31,23) 
    
    startdate_hybrid = datetime(2007,1,1,0)
    enddate_hybrid = datetime(2047,12,31,23)

    startdate_speedy = datetime(1980,1,1,0)
    enddate_speedy = datetime(2008,12,31,23)
   
    mean_hybrid = ds_model.sel(Lon=lon_slice,Lat=lat_slice,Timestep=slice(startdate_hybrid.strftime("%Y-%m-%d"),enddate_hybrid.strftime("%Y-%m-%d"))).mean("Timestep")['SST']
    mean_climo = ds_era.sel(lon=lon_slice,lat=lat_slice,time=slice(startdate_era.strftime("%Y-%m-%d"),enddate_era.strftime("%Y-%m-%d"))).mean("time")['sst']
    mean_speedy = ds_speedy.sel(lon=lon_slice,lat=lat_slice,time=slice(startdate_speedy.strftime("%Y-%m-%d"),enddate_speedy.strftime("%Y-%m-%d"))).mean("time")['sst']

    lons = mean_hybrid.Lon.values
    lats = mean_hybrid.Lat.values

    mean_hybrid = mean_hybrid.where(mean_hybrid > 272.0)

    mean_speedy = mean_speedy.where(mean_speedy > 272.0)

    mean_climo = mean_climo.where(mean_climo > 272.0)

    mean_hybrid = mean_hybrid.values
    mean_climo = mean_climo.values
    mean_speedy = mean_speedy.values

    print('mean_climo')
    print(mean_climo) 
    print(mean_hybrid)
    print(mean_speedy)
    #mean_hybrid = mean_hybrid.where(mean_hybrid > 273.0)
    #mean_speedy = mean_speedy.where(mean_speedy > 273.0)


    print('lons',lons) 
    diff_max = 5

    cmap = get_e3sm_diff_colormap()
 
    cmap_e3sm_diff_colormap_greyer = get_e3sm_diff_colormap_greyer()

    cyclic_data, cyclic_lons = add_cyclic_point(mean_hybrid - mean_climo, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
    projection = ccrs.PlateCarree(central_longitude=-179)
    axes_class = (GeoAxes,dict(map_projection=projection))
    plt.rc('font', family='serif')

    fig = plt.figure()
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(2, 1),
                    axes_pad=0.5,
                    cbar_location='right',
                    cbar_mode='each',
                    cbar_pad=0.2,
                    cbar_size='3%',
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

    levels = [-5,-3,-2,-1,-0.1,0.1,1,2,3,5]
    plot = ax1.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[0].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax1.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax1.set_title(f"Hybrid Sea Surface Temp. 40 yr. Avg. Bias",fontsize=18,fontweight="bold")
    ax1.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    print('Hybrid Sea Surface Temperature Bias',np.nanmin(cyclic_data),np.nanmax(cyclic_data),np.nanmean(cyclic_data))
    rms_hybrid = rms(mean_climo,mean_hybrid)
    print('Hybrid SST RMS',np.nanmin(rms_hybrid),np.nanmax(rms_hybrid),np.nanmean(rms_hybrid))
    ###########

    diff_max = 5
 
    cyclic_data, cyclic_lons = add_cyclic_point(mean_speedy - mean_climo, coord=lons)
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

    levels = [-5,-3,-2,-1,-0.1,0.1,1,2,3,5]
    plot2 = ax2.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels=levels,cmap=cmap,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[1].colorbar(plot2, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax2.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax2.set_title(f"SPEEDY Slab Ocean SST 40 yr. Avg. Bias",fontsize=18,fontweight="bold")
    ax2.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    print('SPEEDY Sea Surface Temperature Bias',np.nanmin(cyclic_data),np.nanmax(cyclic_data),np.nanmean(cyclic_data))
    rms_hybrid = rms(mean_climo,mean_speedy)
    print('SPEEDY SST RMS',np.nanmin(rms_hybrid),np.nanmax(rms_hybrid),np.nanmean(rms_hybrid))

    plt.show()

    fig = plt.figure(figsize=(6,10))
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(3, 1),
                    axes_pad=0.7,
                    cbar_location='right',
                    cbar_mode='each',
                    cbar_pad=0.2,
                    cbar_size='3%',
                    label_mode='')  # note the empty label_mode


    cyclic_data, cyclic_lons = add_cyclic_point(mean_hybrid - 273.15, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

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

    levels = np.arange(0,33,4)
    temp_colormap = sns.color_palette("Spectral", as_cmap=True)
    temp_colormap = temp_colormap.reversed() 
    cmap = temp_colormap
    plot = ax1.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[0].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax1.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax1.set_title(f"Coupled Model Mean SST",fontsize=18,fontweight="bold")
    ax1.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
  
 
    cyclic_data, cyclic_lons = add_cyclic_point(mean_climo - 273.15, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats) 
    ax2 = axgr[1]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax2.coastlines()

    ax2.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax2.xaxis.set_major_formatter(lon_formatter)
    ax2.yaxis.set_major_formatter(lat_formatter)
    gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                       linewidth=2, color='gray', alpha=0.5, linestyle='--')

    #plot = ax1.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=-1*diff_max,vmax=diff_max,cmap='seismic')
    levels = np.arange(0,33,4)
    plot = ax2.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[1].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax2.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax2.set_title(f"ERA5 Mean SST",fontsize=18,fontweight="bold")
    ax2.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
 
    #### Bias 
    cyclic_data, cyclic_lons = add_cyclic_point(mean_hybrid - mean_climo, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats) 
    ax3 = axgr[2]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax3.coastlines()

    ax3.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax3.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax3.xaxis.set_major_formatter(lon_formatter)
    ax3.yaxis.set_major_formatter(lat_formatter)
    gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                       linewidth=2, color='gray', alpha=0.5, linestyle='--')

    #plot = ax1.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=-1*diff_max,vmax=diff_max,cmap='seismic')
    levels = [-5,-3,-2,-1,-0.1,0.1,1,2,3,5]
    plot = ax3.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap_e3sm_diff_colormap_greyer,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[2].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax3.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax3.set_title(f"Sea Surface Temp. Bias (Model - ERA5)",fontsize=18,fontweight="bold")
    ax3.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    plt.show()

    ###########

def sst_annual_climo_and_var_grl2022_paper(ds_model,ds_era,ds_speedy):
    lat_slice = slice(-90,90)
    lon_slice = slice(0,360)

    startdate_era = datetime(1981,1,1,0)
    enddate_era = datetime(2020,12,31,23)

    startdate_hybrid = datetime(2007,1,1,0)
    enddate_hybrid = datetime(2047,12,31,23)

    startdate_speedy = datetime(1980,1,1,0)
    enddate_speedy = datetime(2008,12,31,23)

    mean_hybrid = ds_model.sel(Lon=lon_slice,Lat=lat_slice,Timestep=slice(startdate_hybrid.strftime("%Y-%m-%d"),enddate_hybrid.strftime("%Y-%m-%d"))).mean("Timestep")['SST']
    mean_climo = ds_era.sel(lon=lon_slice,lat=lat_slice,time=slice(startdate_era.strftime("%Y-%m-%d"),enddate_era.strftime("%Y-%m-%d"))).mean("time")['sst']
    mean_speedy = ds_speedy.sel(lon=lon_slice,lat=lat_slice,time=slice(startdate_speedy.strftime("%Y-%m-%d"),enddate_speedy.strftime("%Y-%m-%d"))).mean("time")['sst']

    lons = mean_hybrid.Lon.values
    lats = mean_hybrid.Lat.values

    mean_hybrid = mean_hybrid.where(mean_hybrid > 272.0)

    mean_speedy = mean_speedy.where(mean_speedy > 272.0)

    mean_climo = mean_climo.where(mean_climo > 272.0)

    mean_hybrid = mean_hybrid.values
    mean_climo = mean_climo.values
    mean_speedy = mean_speedy.values

    ds_model = ds_model.sel(Lon=lon_slice,Lat=lat_slice,Timestep=slice(startdate_hybrid.strftime("%Y-%m-%d"),enddate_hybrid.strftime("%Y-%m-%d")))['SST']
    ds_era = ds_era.sel(lon=lon_slice,lat=lat_slice,time=slice(startdate_era.strftime("%Y-%m-%d"),enddate_era.strftime("%Y-%m-%d")))['sst']
    ds_speedy = ds_speedy.sel(lon=lon_slice,lat=lat_slice,time=slice(startdate_speedy.strftime("%Y-%m-%d"),enddate_speedy.strftime("%Y-%m-%d")))['sst']

    ds_model_std = ds_model.groupby("Timestep.month").std("Timestep")
    ds_era_std = ds_era.groupby("time.month").std("time")
    ds_speedy_std = ds_speedy.groupby("time.month").std("time")

    mean_hybrid_std = np.mean(ds_model_std.values,axis=(0))
    mean_climo_std = np.mean(ds_era_std.values,axis=(0))
    mean_speedy_std = np.mean(ds_speedy_std.values,axis=(0))

    print('mean_climo')
    print(mean_climo)
    print(mean_hybrid)
    print(mean_speedy)
    #mean_hybrid = mean_hybrid.where(mean_hybrid > 273.0)
    #mean_speedy = mean_speedy.where(mean_speedy > 273.0)


    print('lons',lons)
    diff_max = 5

    cmap = get_e3sm_diff_colormap()

    cmap_e3sm_diff_colormap_greyer = get_e3sm_diff_colormap_greyer()

    cyclic_data, cyclic_lons = add_cyclic_point(mean_hybrid - mean_climo, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
    projection = ccrs.PlateCarree(central_longitude=-179)
    axes_class = (GeoAxes,dict(map_projection=projection))
    plt.rc('font', family='serif')
    plt.rcParams['figure.constrained_layout.use'] = True

    fig = plt.figure(figsize=(10,18))
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(3, 2),
                    axes_pad=(1.55,0.9),
                    cbar_location='right',
                    cbar_mode='each',
                    cbar_pad=0.2,
                    cbar_size='3%',
                    label_mode='')  # note the empty label_mode


    cyclic_data, cyclic_lons = add_cyclic_point(mean_hybrid - 273.15, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

    ax1 = axgr[2]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax1.coastlines()

    ax1.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax1.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax1.xaxis.set_major_formatter(lon_formatter)
    ax1.yaxis.set_major_formatter(lat_formatter)
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                       linewidth=2, color='gray', alpha=0.5, linestyle='--')

    levels = np.arange(0,33,4)
    temp_colormap = sns.color_palette("Spectral", as_cmap=True)
    temp_colormap = temp_colormap.reversed()
    cmap = temp_colormap
    plot = ax1.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[0].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax1.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax1.set_title(f"Hybrid Model Mean SST",fontsize=18,fontweight="bold")
    ax1.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())


    cyclic_data, cyclic_lons = add_cyclic_point(mean_climo - 273.15, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
    ax2 = axgr[0]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax2.coastlines()

    ax2.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax2.xaxis.set_major_formatter(lon_formatter)
    ax2.yaxis.set_major_formatter(lat_formatter)
    gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                       linewidth=2, color='gray', alpha=0.5, linestyle='--')

    #plot = ax1.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=-1*diff_max,vmax=diff_max,cmap='seismic')
    levels = np.arange(0,33,4)
    plot = ax2.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[2].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax2.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax2.set_title(f"ERA5 Mean SST",fontsize=18,fontweight="bold")
    ax2.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    #### Bias
    cyclic_data, cyclic_lons = add_cyclic_point(mean_hybrid - mean_climo, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
    ax3 = axgr[4]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax3.coastlines()

    ax3.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax3.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax3.xaxis.set_major_formatter(lon_formatter)
    ax3.yaxis.set_major_formatter(lat_formatter)
    gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                       linewidth=2, color='gray', alpha=0.5, linestyle='--')

    #plot = ax1.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=-1*diff_max,vmax=diff_max,cmap='seismic')
    levels = [-5,-3,-2,-1,-0.1,0.1,1,2,3,5]
    plot = ax3.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap_e3sm_diff_colormap_greyer,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[4].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax3.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax3.set_title(f"Sea Surface Temp. Bias\n(Model - ERA5)",fontsize=18,fontweight="bold")
    ax3.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    ###Monthly mean std###

    cyclic_data, cyclic_lons = add_cyclic_point(mean_hybrid_std, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

    ax4 = axgr[3]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax4.coastlines()

    ax4.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax4.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax4.xaxis.set_major_formatter(lon_formatter)
    ax4.yaxis.set_major_formatter(lat_formatter)
    gl = ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                       linewidth=2, color='gray', alpha=0.5, linestyle='--')

    levels = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6]
    temp_colormap = sns.color_palette("Spectral", as_cmap=True)
    temp_colormap = temp_colormap.reversed()
    cmap = temp_colormap
    plot = ax4.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap,extend='max')#'seismic',extend='both')

    cbar = axgr.cbar_axes[1].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax4.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax4.set_title(f"Hybrid Model SST Standard Deviation",fontsize=18,fontweight="bold")
    ax4.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    cyclic_data, cyclic_lons = add_cyclic_point(mean_climo_std, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
    ax5 = axgr[1]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax5.coastlines()

    ax5.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax5.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax5.xaxis.set_major_formatter(lon_formatter)
    ax5.yaxis.set_major_formatter(lat_formatter)
    gl = ax5.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                        linewidth=2, color='gray', alpha=0.5, linestyle='--')

    levels = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6]
    plot = ax5.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap,extend='max')#'seismic',extend='both')

    cbar = axgr.cbar_axes[3].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax5.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax5.set_title(f"ERA5 SST Standard Deviation",fontsize=18,fontweight="bold")
    ax5.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    #### Bias
    cyclic_data, cyclic_lons = add_cyclic_point(mean_hybrid_std - mean_climo_std, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
    ax6 = axgr[5]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax6.coastlines()

    ax6.set_xticks([-180,-120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax6.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax6.xaxis.set_major_formatter(lon_formatter)
    ax6.yaxis.set_major_formatter(lat_formatter)
    gl = ax6.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                       linewidth=2, color='gray', alpha=0.5, linestyle='--')

    levels = [-0.6,-0.4,-0.2,-0.05,0.05,0.2,0.4,0.6]

    plot = ax6.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap_e3sm_diff_colormap_greyer,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[5].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax6.add_feature(cart.feature.LAND, zorder=100, edgecolor='k',facecolor='#808080')

    ax6.set_title(f"Standard Deviation Difference\n(Model - ERA5)",fontsize=18,fontweight="bold")
    ax6.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    plt.show()

def make_ds_time_dim(ds,timestep,startdate):
    begin_year_str = startdate.strftime("%Y-%m-%d")

    attrs = {"units": f"hours since {begin_year_str} "}

    ds = ds.assign_coords({"Timestep": ("Timestep", ds.Timestep.values*timestep, attrs)})
    ds = xr.decode_cf(ds)

    return ds

def count_enso_peaks(array,distance,height):
    print(np.shape(array))
    peaks, _ = find_peaks(array, height=height,distance=distance)

    plt.plot(array)
    plt.plot(peaks,array[peaks], "x")
    plt.show()

def power_spectra_enso(data):
    time_step = 1/1460#4.0#1/14540
    power_spectra_era = np.zeros((np.shape(data)))

    hamming = np.hamming(np.shape(data)[0])

    fft_data = np.fft.fft(data*hamming)
    power_spectra_data = np.abs(fft_data)**2.0

    freq = np.fft.fftfreq(len(data),time_step)

    idx = np.argsort(freq)
    idx = idx[int(len(idx)/2)::]

    return power_spectra_data, freq, idx
    
def sst_trend(ds):
    result = ds.groupby("time.month") - ds.groupby("time.month").mean("time")  
    temp = result['sst'][:].values


    average_temp = np.nanmean(temp,axis=(1,2))
    print(np.shape(average_temp))

    time = np.arange(1,np.shape(average_temp)[0]+1)
    time = time/1460
    plt.plot(time,uniform_filter1d(average_temp,size=1460))

    ticks = np.arange(1,31)
    plt.xticks(ticks)

    plt.show()

def acf(x, length=20):
    return np.array([1]+[np.corrcoef(x[:-i], x[i:])[0,1]  \
        for i in range(1, length)])

def autocorr_plot(ds_hybrid,ds_era):
    hybrid_nino3_4 = nino_index_monthly(ds_hybrid,"3.4")
    
    observed_nino3_4 = nino_index_monthly(ds_observed,"3.4")
    print(hybrid_nino3_4)
    print(observed_nino3_4)

     
    hybrid_nino = uniform_filter1d(hybrid_nino3_4['SST'].values, size=3)
    era5_nino = uniform_filter1d(observed_nino3_4['sst'].values, size=3)

    hybrid_acf = acf(hybrid_nino,length=49)
    era_acf = acf(era5_nino,length=49)

    plt.rc('font', family='serif')

    x = np.arange(0,49)
    ticks = np.arange(0,49,6)
    plt.plot(x,hybrid_acf,label='Coupled Model')
    plt.plot(x,era_acf,label='ERA5') 
    plt.hlines(0.0,x[0],x[-1],linewidth=0.5,color='tab:gray',ls='--')
   
    plt.title('Nino3.4 Autocorrelation',fontsize=18) 
    plt.ylim([-1,1])
    plt.xlim([0,48])
    plt.yticks(fontsize=16)
    plt.xticks(ticks,fontsize=16)

    plt.xlabel('lag (months)',fontsize=16)
    plt.ylabel('autocorrelation',fontsize=16)
    plt.legend(fontsize=16)
    plt.show()
    

   
def oni_soi_timeseries(ds_hybrid,ds_observed,ds_era,ds_speedy):
    hybrid_soi = southern_oscillation_index(ds_hybrid,ds_era,startdate,enddate)

    hybrid_nino3_4 = nino_index(ds_hybrid['SST'],"3.4")
    observed_nino3_4 = nino_index(ds_observed['sst'],"3.4")
    speedy_nino = nino_index_speedy(ds_speedy['sst'],startdate,enddate,"3.4")
    #print(hybrid_nino3_4)
    #print(np.shape(hybrid_nino3_4))

    #print(np.shape(speedy_nino))
    time = np.arange(0,np.shape(hybrid_nino3_4)[0])
    time = time/1461
    time = time + startdate_hybrid.year

    time_soi = np.arange(0,np.shape(hybrid_soi)[0])
    time_soi = time_soi/1461
    time_soi = time_soi + startdate_hybrid.year

    speedy_time = np.arange(0,np.shape(speedy_nino)[0])
    speedy_time = speedy_time/12
    speedy_time = speedy_time + startdate_hybrid.year

    era_time = np.arange(0,np.shape(observed_nino3_4)[0])
    era_time = era_time/1461
    era_time = era_time + startdate_hybrid.year

    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()

    p1, = ax1.plot(time,uniform_filter1d(hybrid_nino3_4, size=360),color='k',ls='-',linewidth=2.0,label="Coupled Model ONI 3 Month Average")
    ax1.fill_between(time, 0.8, uniform_filter1d(hybrid_nino3_4, size=360), where = uniform_filter1d(hybrid_nino3_4, size=360) >= 0.8, alpha = 0.2, color = 'red')
    ax1.fill_between(time,uniform_filter1d(hybrid_nino3_4, size=360), -0.8, where = uniform_filter1d(hybrid_nino3_4, size=360) <= -0.8, alpha = 0.2, color = 'blue')
    #p2, = ax1.plot(era_time,uniform_filter1d(observed_nino3_4,size=360),color='#377eb8',ls='-',linewidth=2.0,label="ERA5 3 Month Average")
    #ax1.plot(speedy_time,uniform_filter1d(speedy_nino,size=3),color='#4daf4a',ls='-',linewidth=2.0,label="SPEEDY 3 Month Average")

    #ax1.hlines(-0.5,time[0],time[-1],color='b',ls='--')
    #ax1.hlines(0.5,time[0],time[-1],color='r',ls='--')

    #p2, = ax2.plot(time_soi,uniform_filter1d(hybrid_soi,size=120),color='#4daf4a',ls='-',linewidth=1.0,label="Hybrid SOI 1 Month Average")
    #p3, = ax2.plot(time_soi,uniform_filter1d(hybrid_soi,size=600),color='#4daf4a',ls='--',linewidth=2.0,label="Hybrid SOI 5 Month Average")
    p2, = ax2.plot(time_soi,uniform_filter1d(hybrid_soi,size=600),color='#4daf4a',ls='--',linewidth=2.0,label="Coupled SOI 5 Month Average")
    print(hybrid_soi)

    ticks = np.arange(startdate_hybrid.year,2047,5)#enddate_hybrid.year+1)#[0,1460,2920,4380,5840,7300]
    labels_ticks = ticks - startdate_hybrid.year

    ax1.set_xlim([np.min(ticks),np.max(ticks)])
    ax1.set_ylim([-3.2,3.2])
    ax2.set_ylim([-20,20]) 

    plt.title("Oceanic Nino Index / Southern Oscillation",fontsize=18,fontweight="bold")
    #plt.title(r"Oceanic Niño Index",fontsize=28,fontweight="bold")

    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels_ticks,fontsize=18)
    ax1.set_xlabel("Years into Simulation",fontsize=20)

    ax1.set_ylabel("ONI",fontsize=20)
    ax2.set_ylabel("SOI",fontsize=16)

    ax1.tick_params(axis='both', which='major', labelsize=18)
    ax1.tick_params(axis='both', which='minor', labelsize=18)

    ax2.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='minor', labelsize=14)

    #ax1.legend(handles=[p1, p2, p3],fontsize=16)
    ax1.legend(handles=[p1, p2],fontsize=16)

    plt.show()
    
def enso_combined_plots(ds_observed_sst,ds_era,ds_hybrid):
     
    print('startdate,enddate',startdate,enddate)
    hybrid_soi = southern_oscillation_index(ds_hybrid,ds_era,startdate,enddate)

    hybrid_nino3_4 = nino_index(ds_hybrid['SST'],"3.4")
    observed_nino3_4 = nino_index(ds_observed_sst['sst'],"3.4")
    speedy_nino = nino_index_speedy(ds_speedy['sst'],startdate,enddate,"3.4")

    ##Monthly or seasonal data
    hybrid_nino3_4_monthly = nino_index_monthly(ds_hybrid,"3.4")

    observed_nino3_4_monthly = nino_index_monthly(ds_observed_sst,"3.4")

    hybrid_nino_3monthly = uniform_filter1d(hybrid_nino3_4_monthly['SST'].values, size=3)
    era5_nino_3monthly = uniform_filter1d(observed_nino3_4_monthly['sst'].values, size=3)

    hybrid_acf = acf(hybrid_nino_3monthly,length=49)
    era_acf = acf(era5_nino_3monthly,length=49)

    file_name = '/home/troyarcomano/FortranReservoir/vert_loc_hybridspeedy_leakage/psl.noaa.gov/gcos_wgsp/Timeseries/Data/nino34.long.anom.data'
    A = np.genfromtxt(file_name,dtype=None,usecols=np.arange(1,13),skip_header=1,skip_footer=8)
    dat = A.flatten()

    t0 = 1870.
    dt = 1/12.0  # In years

    f_Had, Pxx_spec_Had = welch(observed_nino3_4_monthly['sst'].values,nperseg=128*2, scaling='spectrum')

    f_hybrid, Pxx_spec_hybrid = welch(hybrid_nino3_4_monthly['SST'].values,nperseg=128*2, scaling='spectrum')

    print('f_hybrid',f_hybrid,1/f_hybrid)

    HadISST_glbl_power, HadISST_fft_power, HadISST_period, HadISST_fftfreqs = get_wavelet_fft_power(dat,dt,t0)

    print('HadISST_fftfreqs',HadISST_fftfreqs,1/HadISST_fftfreqs)

    hybrid_glbl_power, hybrid_fft_power, hybrid_period, hybrid_fftfreqs = get_wavelet_fft_power(hybrid_nino3_4_monthly['SST'].values,dt,t0)


    hybrid_fft_power = Pxx_spec_hybrid
    HadISST_fft_power = Pxx_spec_Had
    HadISST_fftfreqs = f_Had*12
    hybrid_fftfreqs = f_hybrid*12

    time = np.arange(0,np.shape(hybrid_nino3_4)[0])
    time = time/1461
    time = time + startdate_hybrid.year

    time_soi = np.arange(0,np.shape(hybrid_soi)[0])
    time_soi = time_soi/1461
    time_soi = time_soi + startdate_hybrid.year

    speedy_time = np.arange(0,np.shape(speedy_nino)[0])
    speedy_time = speedy_time/12
    speedy_time = speedy_time + startdate_hybrid.year

    era_time = np.arange(0,np.shape(observed_nino3_4)[0])
    era_time = era_time/1461
    era_time = era_time + startdate_hybrid.year


    axd = plt.figure(constrained_layout=True).subplot_mosaic(
                     """
                     AA
                     BC
                     """
                     )
    print(axd)
    ax1 = axd['A']

    ax1b = ax1.twinx()

    p1, = ax1.plot(time,uniform_filter1d(hybrid_nino3_4, size=360),color='k',ls='-',linewidth=2.0,label="Hybrid Model ONI 3 Month Average")
    ax1.fill_between(time, 0.8, uniform_filter1d(hybrid_nino3_4, size=360), where = uniform_filter1d(hybrid_nino3_4, size=360) >= 0.8, alpha = 0.2, color = 'red')
    ax1.fill_between(time,uniform_filter1d(hybrid_nino3_4, size=360), -0.8, where = uniform_filter1d(hybrid_nino3_4, size=360) <= -0.8, alpha = 0.2, color = 'blue')
    #p2, = ax1.plot(era_time,uniform_filter1d(observed_nino3_4.sel(Timestep=slice(,size=360),color='#377eb8',ls='-',linewidth=2.0,label="ERA5 3 Month Average")
    #ax1.plot(speedy_time,uniform_filter1d(speedy_nino,size=3),color='#4daf4a',ls='-',linewidth=2.0,label="SPEEDY 3 Month Average")

    #ax1.hlines(-0.5,time[0],time[-1],color='b',ls='--')
    #ax1.hlines(0.5,time[0],time[-1],color='r',ls='--')

    #p2, = ax2.plot(time_soi,uniform_filter1d(hybrid_soi,size=120),color='#4daf4a',ls='-',linewidth=1.0,label="Hybrid SOI 1 Month Average")
    #p3, = ax2.plot(time_soi,uniform_filter1d(hybrid_soi,size=600),color='#4daf4a',ls='--',linewidth=2.0,label="Hybrid SOI 5 Month Average")
    p2, = ax1b.plot(time_soi,uniform_filter1d(hybrid_soi,size=600),color='#4daf4a',ls='--',linewidth=2.0,label="Hybrid Model SOI 5 Month Average")
    print(hybrid_soi)

    ticks = np.arange(startdate_hybrid.year,2057,5)#enddate_hybrid.year+1)#[0,1460,2920,4380,5840,7300]
    labels_ticks = ticks - startdate_hybrid.year

    ax1.set_xlim([np.min(ticks),np.max(ticks)])
    ax1.set_ylim([-3.2,3.2])
    ax1b.set_ylim([-20,20])

    plt.title("Oceanic Nino Index / Southern Oscillation",fontsize=18,fontweight="bold")
    #plt.title(r"Oceanic Niño Index",fontsize=28,fontweight="bold")

    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels_ticks,fontsize=18)
    ax1.set_xlabel("Years into Simulation",fontsize=20)

    ax1.set_ylabel("ONI",fontsize=20)
    ax1b.set_ylabel("SOI",fontsize=16)

    ax1.tick_params(axis='both', which='major', labelsize=18)
    ax1.tick_params(axis='both', which='minor', labelsize=18)

    ax1b.tick_params(axis='both', which='major', labelsize=14)
    ax1b.tick_params(axis='both', which='minor', labelsize=14)

    #ax1.legend(handles=[p1, p2, p3],fontsize=16)
    ax1.legend(handles=[p1, p2],fontsize=16)

    #####AAutocorrelation#####

    ax2 = axd['B']
    x = np.arange(0,49)
    ticks = np.arange(0,49,6)
    ax2.plot(x,era_acf,label='ERA5')
    ax2.plot(x,hybrid_acf,label='Hybrid Model')
    ax2.hlines(0.0,x[0],x[-1],linewidth=0.5,color='tab:gray',ls='--')

    ax2.set_title('Nino3.4 Autocorrelation',fontsize=18)
    ax2.set_ylim([-1,1])
    ax2.set_xlim([0,48])
    ax2.set_xticks(ticks)
    ax2.tick_params(axis='both', which='major', labelsize=16)

    ax2.set_xlabel('lag (months)',fontsize=16)
    ax2.set_ylabel('autocorrelation',fontsize=16)
    ax2.legend(fontsize=16)


    xticks = np.arange(0,11,1)

    ax3 = axd['C']
    ax3.plot(1/HadISST_fftfreqs,HadISST_fft_power,linewidth=1.5,label='ERA5')
    ax3.plot(1/hybrid_fftfreqs,hybrid_fft_power,linewidth=1.5,label='Hybrid Model')

    ax3.set_title('Nino 3.4 Power Spectrum',fontsize=20)
    ax3.set_ylabel(r'Power[$\degree C^{2}$]',fontsize=16)
    ax3.set_xlabel('Period (Years)',fontsize=16)
    #cx.set_yticks(np.arange(0,40,5))
    #cx.set_yticklabels(np.arange(0,40,5),fontsize=16)
    #cx.set_ylim([0, hybrid_glbl_power.max()])

    ax3.set_xticks(xticks)
    ax3.set_xticklabels(xticks,fontsize=16)
    ax3.set_xlim([0,xticks.max()])
    ax3.legend(fontsize=16)

    plt.show()

def get_wavelet_fft_power(dat,dt,t0):
    N = dat.size
    t = np.arange(0, N) * dt + t0
    #We write the following code to detrend and normalize the input data by its
    # standard deviation. Sometimes detrending is not necessary and simply
    # removing the mean value is good enough. However, if your dataset has a well
    # defined trend, such as the Mauna Loa CO\ :sub:`2` dataset available in the
    # above mentioned website, it is strongly advised to perform detrending.
    # Here, we fit a one-degree polynomial function and then subtract it from the
    # original data.
    p = np.polyfit(t - t0, dat, 1)
    dat_notrend = dat - np.polyval(p, t - t0)
    std = dat_notrend.std()  # Standard deviation
    var = std ** 2  # Variance
    dat_norm = dat_notrend / std  # Normalized dataset

    # The next step is to define some parameters of our wavelet analysis. We
    # select the mother wavelet, in this case the Morlet wavelet with
    # :math:`\omega_0=6`.
    mother = wavelet.Morlet(6)
    #s0 = 2 * dt  # Starting scale, in this case 2 * 0.25 years = 6 months
    s0 = 6 * dt
    dj = 1 / 12  # Twelve sub-octaves per octaves
    J = 7 / dj  # Seven powers of two with dj sub-octaves
    alpha, _, _ = wavelet.ar1(dat)  # Lag-1 autocorrelation for red noise

    # The following routines perform the wavelet transform and inverse wavelet
    # transform using the parameters defined above. Since we have normalized our
    # input time-series, we multiply the inverse transform by the standard
    # deviation.
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(dat_norm, dt, dj, s0, J,
                                                          mother)
    iwave = wavelet.icwt(wave, scales, dt, dj, mother) * std

    # We calculate the normalized wavelet and Fourier power spectra, as well as
    # the Fourier equivalent periods for each wavelet scale.
    power = (np.abs(wave)) ** 2
    fft_power = np.abs(fft) ** 2
    period = 1 / freqs

    # We could stop at this point and plot our results. However we are also
    # interested in the power spectra significance test. The power is significant
    # where the ratio ``power / sig95 > 1``.
    signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha,
                                             significance_level=0.95,
                                             wavelet=mother)
    sig95 = np.ones([1, N]) * signif[:, None]
    sig95 = power / sig95

    # Then, we calculate the global wavelet spectrum and determine its
    # significance level.
    glbl_power = power.mean(axis=1)
    dof = N - scales  # Correction for padding at edges
    glbl_signif, tmp = wavelet.significance(var, dt, scales, 1, alpha,
                                            significance_level=0.95, dof=dof,
                                            wavelet=mother)

    # We also calculate the scale average between 2 years and 8 years, and its
    # significance level.
    sel = find((period >= 2) & (period < 8))
    Cdelta = mother.cdelta
    scale_avg = (scales * np.ones((N, 1))).transpose()
    scale_avg = power / scale_avg  # As in Torrence and Compo (1998) equation 24
    scale_avg = var * dj * dt / Cdelta * scale_avg[sel, :].sum(axis=0)
    scale_avg_signif, tmp = wavelet.significance(var, dt, scales, 2, alpha,
                                                 significance_level=0.90,#5,
                                                 dof=[scales[sel[0]],
                                                      scales[sel[-1]]],
                                                 wavelet=mother)


    return var * glbl_power, var * fft_power, period, fftfreqs


ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_sigma0.5_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_multi_gaussian_noise_newest_version_32_processors_root_ssttrial_12_29_2006_00.nc')

#hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_ocean_model7d_0.0001beta_sigma0.6_dp_noise10trial_12_29_2006_00.nc')

#hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_ocean_model7d_0.0001beta_sigma0.6_dp_noise10trial_12_29_2006_00.nc')

#hybrid_prediction_era6000_20_20_20_sigma0.5_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_multi_gaussian_noise_trial_12_29_2006_00.nc')

#hybrid_prediction_era6000_20_20_20_sigma0.5_beta_res0.001_beta_model_0.01_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_multi_gaussian_noisetrial_01_03_2007_00.nc')

#hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_ml_only_ocean_beta_0.0001_sigma0.6_6diff_speedytrial_12_29_2006_00.nc')
#hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_tstep6_ml_only_ocean_input_speedy_testtrial_12_31_1999_00.nc') #hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_tstep6_ml_only_ocean_leak_averaged_atmotrial_12_31_1999_00.nc')#hybrid_prediction_era6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_ocean_model14d_0.1rho_10noise_beta0.001_20years_speedyinput_falsetrial_12_31_1999_00.nc')

#ds_speedy = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/speedy_slab_ocean_original_sim.nc')
ds_speedy = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/speedy_slab_with_qocn.nc')
print(ds_speedy)

startdate = datetime(1981,1,1,0)
enddate = datetime(2021,1,1,0)
#enddate = datetime(1985,1,1,0)

enddate_soi_era = datetime(2007,1,1,0)
#enddate_soi_era = datetime(1989,1,1,0)

startdate_hybrid = datetime(2007,1,1,0)
enddate_hybrid = datetime(2021,1,1,0)

enddate_hybrid_future = datetime(2072,1,1,0)

timestep = 6

ds_observed = get_obs_sst_timeseries(startdate,enddate,timestep)

ds_era = get_obs_atmo_timeseries_var(startdate,enddate_soi_era,timestep,'logp')
#print(ds_era)

ds_hybrid =  make_ds_time_dim(ds_hybrid,6,startdate_hybrid)

enso_combined_plots(ds_observed,ds_era,ds_hybrid)

#sst_annual_climo_and_var_grl2022_paper(ds_hybrid,ds_observed,ds_speedy)
#oni_soi_timeseries(ds_hybrid,ds_observed,ds_era,ds_speedy)

#sst_annual_variability(ds_hybrid,ds_observed,ds_speedy)

#sst_climatology_error(ds_hybrid,ds_observed,ds_speedy)


#autocorr_plot(ds_hybrid,ds_observed)

#count_enso_peaks(uniform_filter1d(hybrid_nino3_4, size=360),1000,0.8)


