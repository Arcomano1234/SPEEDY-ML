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
from scipy.signal import find_peaks
import iris as i
import iris.plot as iplt
import iris.quickplot as qplt

#geoms = fiona.open(shpreader.natural_earth(resolution='10m', category='physical', name='land'))
#land_geom = sgeom.MultiPolygon([sgeom.shape(geom['geometry']) for geom in geoms])
#land = prep(land_geom)

#def is_land(lon,lat#):
#    return land.contains(sgeom.Point(lon, lat)) 


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
    end_year = 1999

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
    '''
    plt.style.use('ggplot')

    data_file = '../output/exp_103/attm103.nc'

    cube = i.load(data_file, 'mean-sea-level pressure         [hPa]')[0]

    tahiti = cube.extract(i.Constraint(latitude=-16.7, longitude=210.0))
    darwin = cube.extract(i.Constraint(latitude=-12.989, longitude=131.25))

    diff = tahiti - darwin
    std = diff.collapsed('time', i.analysis.STD_DEV)
    soi = 10 * diff / std
    soi.rename('Southern Oscillation Index')

    qplt.plot(soi)
    plt.show()
    '''

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
    plt.plot(hybrid_soi)
    plt.show() 

    era_soi = xr.apply_ufunc(
        lambda x, m, s: (x - m) / s,
        diff_climo.groupby("Timestep.month"),
        climatology_mean,
        climatology_std,
    )
    era_soi = era_soi * 10.0
    plt.plot(era_soi)
    plt.show()
    return hybrid_soi

def get_obs_sst_timeseries(startdate,enddate,timestep):
    start_year = startdate.year
    end_year = enddate.year

    currentdate = startdate
    counter = 0
    while currentdate.year <= enddate.year:
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

def sst_climatology_error(ds_model,ds_era,ds_speedy):  
    lat_slice = slice(-90,90)
    lon_slice = slice(0,360) 

    startdate_era = datetime(2001,1,1,0)
    enddate_era = datetime(2010,12,31,23) 
    
    startdate_hybrid = datetime(2001,1,1,0)
    enddate_hybrid = datetime(2010,12,31,23)

    startdate_speedy = datetime(1997,1,1,0)
    enddate_speedy = datetime(2008,12,31,23)
   
    mean_hybrid = ds_model.sel(Lon=lon_slice,Lat=lat_slice,Timestep=slice(startdate_hybrid.strftime("%Y-%m-%d"),enddate_hybrid.strftime("%Y-%m-%d"))).mean("Timestep")['SST']
    mean_climo = ds_era.sel(lon=lon_slice,lat=lat_slice,time=slice(startdate_era.strftime("%Y-%m-%d"),enddate_era.strftime("%Y-%m-%d"))).mean("time")['sst']
    mean_speedy = ds_speedy.sel(lon=lon_slice,lat=lat_slice,time=slice(startdate_speedy.strftime("%Y-%m-%d"),enddate_speedy.strftime("%Y-%m-%d"))).mean("time")['sst']

    lons = mean_hybrid.Lon.values
    lats = mean_hybrid.Lat.values

    mean_hybrid = mean_hybrid.where(mean_hybrid > 273.0)

    mean_hybrid = mean_hybrid.values
    mean_climo = mean_climo.values
    mean_speedy = mean_speedy.values

    print('mean_climo')
    print(mean_climo) 
    print(mean_hybrid)
    print(mean_speedy)
    #mean_hybrid = mean_hybrid.where(mean_hybrid > 273.0)
    #mean_speedy = mean_speedy.where(mean_speedy > 273.0)

 
    diff_max = 5

    cmap = get_e3sm_diff_colormap()

    cyclic_data, cyclic_lons = add_cyclic_point(mean_hybrid - mean_climo, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
    projection = ccrs.PlateCarree()
    axes_class = (GeoAxes,dict(map_projection=projection))
    plt.rc('font', family='serif')

    fig = plt.figure()
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(2, 1),
                    axes_pad=0.9,
                    cbar_location='right',
                    cbar_mode='each',
                    cbar_pad=0.2,
                    cbar_size='3%',
                    label_mode='')  # note the empty label_mode


    ax1 = axgr[0]#plt.subplot(221,projection=ccrs.PlateCarree())
    ax1.coastlines()

    ax1.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    ax1.xaxis.set_major_formatter(lon_formatter)

    #plot = ax1.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=-1*diff_max,vmax=diff_max,cmap='seismic')
    levels = [-5,-3,-2,-1,-0.1,0.1,1,2,3,5]
    plot = ax1.contourf(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),levels = levels,cmap=cmap,extend='both')#'seismic',extend='both')

    cbar = axgr.cbar_axes[0].colorbar(plot, extend='both')
    cbar.set_ticks(levels)
    cbar.set_label(r'$\degree C$',fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax1.add_feature(cart.feature.LAND, zorder=100, edgecolor='k')

    ax1.set_title(f"Hybrid SST Bias",fontsize=18,fontweight="bold")
    ax1.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    print('Hybrid SST Bias',np.nanmin(cyclic_data),np.nanmax(cyclic_data),np.nanmean(cyclic_data))
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
    ax2.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    ax2.set_yticks(np.linspace(-90, 90, 5), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax2.xaxis.set_major_formatter(lon_formatter)
    ax2.yaxis.set_major_formatter(lat_formatter)

    plot_2 = ax2.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=-1*diff_max,vmax=diff_max,cmap='seismic')
    cbar2 = axgr.cbar_axes[1].colorbar(plot_2, extend='both')
    #cbar2.set_label_text(f'hPa',fontsize=16)
    cbar2.set_label(r'$\degree C$',fontsize=16)
    cbar2.ax.tick_params(labelsize=16)

    ax2.set_title(f"SPEEDY SST Bias",fontsize=18,fontweight="bold")
    ax2.add_feature(cart.feature.LAND, zorder=100, edgecolor='k')

    ax2.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
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
    
    
    
#def rmse_sst(truth,prediction,region):
    
           
ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_ml_only_ocean_beta_0.0001_sigma0.6_6diff_speedytrial_12_29_2006_00.nc')
#hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_tstep6_ml_only_ocean_input_speedy_testtrial_12_31_1999_00.nc') #hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_tstep6_ml_only_ocean_leak_averaged_atmotrial_12_31_1999_00.nc')#hybrid_prediction_era6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_ocean_model14d_0.1rho_10noise_beta0.001_20years_speedyinput_falsetrial_12_31_1999_00.nc')

#ds_speedy = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/speedy_slab_ocean_original_sim.nc')
ds_speedy = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/speedy_slab_with_qocn.nc')
print(ds_speedy)

startdate = datetime(1981,1,1,0)
enddate = datetime(2021,1,1,0)

enddate_soi_era = datetime(2006,1,1,0)

startdate_hybrid = datetime(2007,1,1,0)
enddate_hybrid = datetime(2021,1,1,0)

enddate_hybrid_future = datetime(2030,1,1,0)

timestep = 6
ds_observed = get_obs_sst_timeseries(startdate_hybrid,enddate_hybrid,timestep)

ds_era = get_obs_atmo_timeseries_var(startdate,enddate_soi_era,timestep,'logp')
print(ds_era)

ds_hybrid =  make_ds_time_dim(ds_hybrid,6,startdate_hybrid)

hybrid_soi = southern_oscillation_index(ds_hybrid,ds_era,startdate,enddate)
print(hybrid_soi)
sst_climatology_error(ds_hybrid,ds_observed,ds_speedy)

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

p1, = ax1.plot(time,uniform_filter1d(hybrid_nino3_4, size=360),color='k',ls='-',linewidth=2.0,label="Hybrid Model ONI 3 Month Average")
ax1.fill_between(time, 0.8, uniform_filter1d(hybrid_nino3_4, size=360), where = uniform_filter1d(hybrid_nino3_4, size=360) >= 0.8, alpha = 0.2, color = 'red')
ax1.fill_between(time,uniform_filter1d(hybrid_nino3_4, size=360), -0.8, where = uniform_filter1d(hybrid_nino3_4, size=360) <= -0.8, alpha = 0.2, color = 'blue')
#ax1.plot(era_time,uniform_filter1d(observed_nino3_4,size=360),color='#377eb8',ls='-',linewidth=2.0,label="ERA5 3 Month Average")
#ax1.plot(speedy_time,uniform_filter1d(speedy_nino,size=3),color='#4daf4a',ls='-',linewidth=2.0,label="SPEEDY 3 Month Average")

#ax1.hlines(-0.5,time[0],time[-1],color='b',ls='--')
#ax1.hlines(0.5,time[0],time[-1],color='r',ls='--')

p2, = ax2.plot(time_soi,uniform_filter1d(hybrid_soi,size=120),color='#4daf4a',ls='-',linewidth=1.0,label="Hybrid SOI 1 Month Average")
p3, = ax2.plot(time_soi,uniform_filter1d(hybrid_soi,size=600),color='#4daf4a',ls='--',linewidth=2.0,label="Hybrid SOI 5 Month Average")
print(hybrid_soi)

ticks = np.arange(startdate_hybrid.year,2035)#enddate_hybrid.year+1)#[0,1460,2920,4380,5840,7300]

ax1.set_xlim([np.min(ticks),np.max(ticks)])
ax2.set_ylim([-20,20])
plt.title("Oceanic Nino Index / Southern Oscillation")
ax1.set_xticks(ticks)
ax1.set_xticklabels(ticks)
plt.xlabel("Time")
ax1.set_ylabel("ONI")
ax2.set_ylabel("SOI")
ax1.legend(handles=[p1, p2, p3])
plt.show()

count_enso_peaks(uniform_filter1d(hybrid_nino3_4, size=360),1000,0.8)
