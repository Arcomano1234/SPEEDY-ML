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
import cmocean as ocean
from mpl_toolkits.axes_grid1 import AxesGrid
import seaborn as sns

class StretchOutNormalize(plt.Normalize):
    def __init__(self, vmin=None, vmax=None, low=None, up=None, clip=False):
        self.low = low
        self.up = up
        plt.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.low, self.up, self.vmax], [0, 0.5-1e-9, 0.5+1e-9, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def get_colorbars(fig):
    cbs = []
    for ax in fig.axes:
        cbs.extend(ax.findobj(lambda obj: hasattr(obj, "colorbar") and obj.colorbar))
    return [a.colorbar for a in cbs]

@jit()
def rms(true,prediction):
    return np.sqrt(np.nanmean((prediction-true)**2))

def latituded_weighted_rmse(true,prediction,lats):
    diff = prediction-true

    weights = np.cos(np.deg2rad(lats))

    weights2d = np.zeros(np.shape(diff))

    diff_squared = diff**2.0
    #weights = np.ones((10,96))

    weights2d = np.tile(weights,(96,1))
    weights2d = np.transpose(weights2d)

    masked = np.ma.MaskedArray(diff_squared, mask=np.isnan(diff_squared))
    weighted_average = np.ma.average(masked,weights=weights2d)

    return np.sqrt(weighted_average)

def make_ds_time_dim(ds,timestep,startdate):
    begin_year_str = startdate.strftime("%Y-%m-%d")

    attrs = {"units": f"hours since {begin_year_str} "}

    ds = ds.assign_coords({"Timestep": ("Timestep", ds.Timestep.values*timestep, attrs)})
    ds = xr.decode_cf(ds)

    return ds

def get_era5_precip_timeseries(startdate,enddate,timestep):
    start_year = startdate.year
    end_year = enddate.year

    currentdate = startdate
    counter = 0
    while currentdate.year <= enddate.year:
        ds_era = xr.open_dataset(f'/scratch/user/troyarcomano/ERA_5/{currentdate.year}/era_5_y{currentdate.year}_precip_regridded_mpi.nc')

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
    return ds_merged.sel(Timestep=time_slice) #['tp'] * 39.3701

def speedy_total_precip(ds): 
    ds = ds.assign(tp=(ds["precls"] + ds["precnv"]))
    return ds

def log_binning(data):
    #data unsorted
    copy = np.copy(data)

    copy = np.sort(copy)
    print(np.shape(copy))

    print(copy[-10:-1])


    percent_range = np.logspace(-3,2,100)#(-4,2,100) #[80,90,95,99,99.9,99.95, 99.99,99.999, 99.9999]
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
    print(percent_range)
    return np.nanpercentile(copy,percent_range), power #percent_range, power

def extreme_value_plot(ds_era,ds_hybrid,ds_speedy):

    print(np.max(ds_hybrid['p6hr'].values))
    ds_era = ds_era.where(ds_era['tp'] >= 0.0,drop=True)
    ds_speedy = ds_speedy.where(ds_speedy['rain'] >= 0.0,drop=True)

    era_uwind = ds_era['tp'].to_numpy().flatten() * 1000.0 #uv_to_windspeedy(ds_era) #abs(ds_era['U-wind'].to_numpy().flatten())
    hybrid_uwind = ds_hybrid['p6hr'].to_numpy().flatten() * 25.4#uv_to_windspeedy(ds_hybrid) #abs(ds_hybrid['U-wind'].to_numpy().flatten())
    speedy_uwind = ds_speedy['rain'].to_numpy().flatten()/(3.6*4.0*6.0) #uv_to_windspeedy(ds_speedy) #abs(ds_speedy['U-wind'].to_numpy().flatten())

    print(np.shape(era_uwind))
    print(np.shape(hybrid_uwind))
    print(np.shape(speedy_uwind))

    print('max era',np.max(era_uwind))
    print('max hybrid',np.max(hybrid_uwind))
    print('max speedy',np.max(speedy_uwind))

    print('min era',np.min(era_uwind))
    print('min hybrid',np.min(hybrid_uwind))
    print('min speedy',np.min(speedy_uwind))

    era_uwind_extreme, percent_range = log_binning(era_uwind)
    hybrid_uwind_extreme, percent_range = log_binning(hybrid_uwind)
    speedy_uwind_extreme, percent_range = log_binning(speedy_uwind)

    percent_range = percent_range
    era_uwind_extreme = era_uwind_extreme
    hybrid_uwind_extreme = hybrid_uwind_extreme
    speedy_uwind_extreme = speedy_uwind_extreme

    fig1, ax1 = plt.subplots()

    #ax1.semilogx(percent_range,era_uwind_extreme, linewidth=2,label='ERA5')
    #ax1.semilogx(percent_range,hybrid_uwind_extreme,linewidth=2,label='Hybrid')
    #ax1.semilogx(percent_range,speedy_uwind_extreme,linewidth=2,label="SPEEDY")
    ax1.loglog(percent_range,era_uwind_extreme, linewidth=2,label='ERA5')
    ax1.loglog(percent_range,hybrid_uwind_extreme,linewidth=2,label='Hybrid')
    ax1.loglog(percent_range,speedy_uwind_extreme,linewidth=2,label="SPEEDY")

    print(percent_range)
    ax1.set_xlim([10**-3,100])
    ax1.set_ylim([10**-2,10**2])#np.max(era_uwind_extreme)])
    #ax1.set_xscale('log')

    x_labels_vals = [10**1,10**0,10**-1,10**-2,10**-3]
    x_tick_labels = ['90%','99%','99.9%','99.99%','99.999%']

    ax1.tick_params(axis='y', which='major', labelsize=16)
    ax1.invert_xaxis()
    ax1.set_xticks(x_labels_vals)
    ax1.set_xticklabels(x_tick_labels,fontsize=16)

    ax1.set_xlabel("Percentile",fontsize=20)

    ax1.set_ylabel("Total Precipitation (mm/6 h)",fontsize=20)


    #ax1.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    plt.legend(fontsize=20)
    plt.show()

def histo_precip():
    startdate = datetime(2000,1,1,0)
    enddate = datetime(2010,12,31,0)

    startdate_climo = datetime(2001,1,1,0)
    enddate_climo = datetime(2010,12,31,0)


    ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_tstep6_ml_only_ocean_leak_averaged_atmotrial_12_31_1999_00.nc')#

    timestep = 6

    ds_speedy = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/speedy_slab_ocean_original_sim.nc')
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

    hist = sns.histplot(total_era5.flatten(), bins=bins,
                 multiple="layer",label='ERA5',color='#e41a1c',alpha=0.5)

    hist.set_ylabel("Cumulative Probability", fontsize = 14)
    hist.set_xlabel("Annual Precipitation (Inches)", fontsize = 14)

    plt.legend()

    plt.xlim([0,300])
    #plt.ylim([0,1])

    '''
    n_bins = 1000
    fig, ax = plt.subplots(figsize=(8, 4))
    n, bins, patches = ax.hist(total_speedy.flatten(), n_bins, density=True, histtype='step',
                               cumulative=True, label='SPEEDY')

    n, bins, patches = ax.hist(total_hybrid.flatten(), n_bins, density=True, histtype='step',
                               cumulative=True, label='Hybrid')
    n, bins, patches = ax.hist(total_era5.flatten(), n_bins, density=True, histtype='step',
                               cumulative=True, label='ERA5')

    plt.legend()

    plt.xlim([0,300])
    plt.ylim([0,1])
    '''
    plt.show()

  
def cdf_precip():
    startdate = datetime(2000,1,1,0)
    enddate = datetime(2010,12,31,0)

    startdate_climo = datetime(2001,1,1,0)
    enddate_climo = datetime(2010,12,31,0)


    ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_tstep6_ml_only_ocean_leak_averaged_atmotrial_12_31_1999_00.nc')#

    timestep = 6

    ds_speedy = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/speedy_slab_ocean_original_sim.nc')
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

    n_bins = 1000
    sns.histplot(total_speedy.flatten(), bins=n_bins,
                 log_scale=False, element="step", fill=False,
                 cumulative=True, stat="density", common_norm=False,label='SPEEDY',color='#4daf4a')
     
    sns.histplot(total_hybrid.flatten(), bins=n_bins,
                 log_scale=False, element="step", fill=False,
                 cumulative=True, stat="density", common_norm=False,label='Hybrid',color='#377eb8')

    hist = sns.histplot(total_era5.flatten(), bins=n_bins,
                 log_scale=False, element="step", fill=False,
                 cumulative=True, stat="density", common_norm=False,label='ERA5',color='#e41a1c')

    hist.set_ylabel("Cumulative Probability", fontsize = 14)
    hist.set_xlabel("Annual Precipitation (Inches)", fontsize = 14)

    plt.legend()

    plt.xlim([0,300])
    plt.ylim([0,1])
    '''
    n_bins = 1000
    fig, ax = plt.subplots(figsize=(8, 4))
    n, bins, patches = ax.hist(total_speedy.flatten(), n_bins, density=True, histtype='step',
                               cumulative=True, label='SPEEDY')

    n, bins, patches = ax.hist(total_hybrid.flatten(), n_bins, density=True, histtype='step',
                               cumulative=True, label='Hybrid')
    n, bins, patches = ax.hist(total_era5.flatten(), n_bins, density=True, histtype='step',
                               cumulative=True, label='ERA5')

    plt.legend()

    plt.xlim([0,300])
    plt.ylim([0,1])
    '''
    plt.show()

def total_precip_bias_plot():
    startdate = datetime(1981,1,1,0)
    enddate = datetime(2020,12,31,0)

    startdate_climo = datetime(1981,1,1,0)
    enddate_climo = datetime(2018,12,31,0)

    startdate_hybrid = datetime(2007,1,1,0)
    enddate_hybrid  = datetime(2047,12,31,0)
    #ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_sigma0.5_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precipclimo_2kbias_10_year_then_platue_speedy_bc_atmo_no_ice_2k_sst_mean_20std_increase_trial_12_29_2006_00.nc')
    ds_hybrid = xr.open_dataset('/scratch/user/troyarcomano/Predictions/Hybrid/hybrid_prediction_era6000_20_20_20_sigma0.5_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_multi_gaussian_noise_newest_version_32_processors_root_ssttrial_12_29_2006_00.nc') #hybrid_prediction_era6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_ocean_model_false_log_preciptrial_12_31_1999_00.nc')#

    timestep = 6

    ds_speedy = xr.open_dataset('/scratch/user/troyarcomano/temp_storage/speedy_slab_ocean_original_sim.nc')
    ds_speedy = speedy_total_precip(ds_speedy)

    ds_speedy = ds_speedy.sel(time=slice(startdate_climo.strftime("%Y-%m-%d"),enddate_climo.strftime("%Y-%m-%d")))
    ds_speedy_annual_climo = ds_speedy.groupby('time.year').mean('time')['tp']

    ds_observed = get_era5_precip_timeseries(startdate_climo,enddate_climo,1)
    ds_observed_annual_climo = ds_observed.groupby('Timestep.year').sum('Timestep')['tp']

    ds_hybrid =  make_ds_time_dim(ds_hybrid,timestep,startdate)
    ds_hybrid = ds_hybrid.sel(Timestep=slice(startdate_hybrid.strftime("%Y-%m-%d"),enddate_hybrid.strftime("%Y-%m-%d"))) #NOTE change back to start_climot end climato
    #ds_hybrid = ds_hybrid.sel(Timestep=slice(startdate_climo.strftime("%Y-%m-%d"),enddate_climo.strftime("%Y-%m-%d")))

    ds_hybrid_annual_climo = ds_hybrid.groupby('Timestep.year').sum('Timestep')['p6hr']
    print(ds_hybrid_annual_climo)

    total_era5 = np.average(ds_observed_annual_climo.values, axis=0) * (1000.0 / 365.25) #Units meters per year to mm per day
    total_speedy = np.average(ds_speedy_annual_climo.values, axis=0) #/ 12.0 #mm/day to  
    total_hybrid = np.average(ds_hybrid_annual_climo.values, axis=0) * (25.39998628 / 365.25) #Inches per year to mm/day

    speedy_bias = total_speedy - total_era5
    hybrid_bias = total_hybrid - total_era5

    weights = np.cos(np.deg2rad(ds_observed_annual_climo.Lat))

    weights2d = np.tile(weights,(96,1))
    weights2d = np.transpose(weights2d)

    print('ERA5',np.max(total_era5),np.min(total_era5),np.average(total_era5,weights=weights2d),np.std(total_era5),np.corrcoef(total_era5.flatten(),total_era5.flatten()))
    print('SPEEDY',np.max(total_speedy),np.min(total_speedy),np.average(total_speedy,weights=weights2d),np.std(total_speedy),np.corrcoef(total_era5.flatten(),total_speedy.flatten()))
    print('Hybrid',np.max(total_hybrid),np.min(total_hybrid),np.average(total_hybrid,weights=weights2d),np.std(total_hybrid),np.corrcoef(total_era5.flatten(),total_hybrid.flatten()))

    print('SPEEDY Bias',np.max(speedy_bias),np.min(speedy_bias),np.average(speedy_bias,weights=weights2d),np.std(speedy_bias))
    print('Hybrid Bias',np.max(hybrid_bias),np.min(hybrid_bias),np.average(hybrid_bias,weights=weights2d),np.std(hybrid_bias))

    print('SPEEDY abs error',np.max(abs(speedy_bias)),np.min(abs(speedy_bias)),np.average(abs(speedy_bias),weights=weights2d),np.std(abs(speedy_bias)))
    print('Hybrid abs error',np.max(abs(hybrid_bias)),np.min(abs(hybrid_bias)),np.average(abs(hybrid_bias),weights=weights2d),np.std(abs(hybrid_bias)))

    print('SPEEDY RMSE', latituded_weighted_rmse(total_era5,total_speedy,weights))
    print('Hybrid RMSE', latituded_weighted_rmse(total_era5,total_hybrid,weights))

    lons = ds_hybrid_annual_climo.Lon
    lats = ds_hybrid_annual_climo.Lat


    ####Precip extreme
    ds_era = ds_observed.resample(Timestep = "6H").sum()
    ds_era = ds_era.where(ds_era['tp'] >= 0.0,drop=True)

    era_uwind = ds_era['tp'].to_numpy().flatten() * 1000.0 #uv_to_windspeedy(ds_era) #abs(ds_era['U-wind'].to_numpy().flatten())
    hybrid_uwind = ds_hybrid['p6hr'].to_numpy().flatten() * 25.4#uv_to_windspeedy(ds_hybrid) #abs(ds_hybrid['U-wind'].to_numpy().flatten())

    print(np.shape(era_uwind))
    print(np.shape(hybrid_uwind))

    print('max era',np.max(era_uwind))
    print('max hybrid',np.max(hybrid_uwind))

    print('min era',np.min(era_uwind))
    print('min hybrid',np.min(hybrid_uwind))

    era_uwind_extreme, percent_range = log_binning(era_uwind)
    hybrid_uwind_extreme, percent_range = log_binning(hybrid_uwind)

    percent_range = percent_range
    era_uwind_extreme = era_uwind_extreme
    hybrid_uwind_extreme = hybrid_uwind_extreme

    nws_precip_colors = [
        "#04e9e7",  # 0.01 - 0.10 inches
        "#019ff4",  # 0.10 - 0.25 inches
        "#0300f4",  # 0.25 - 0.50 inches
        "#02fd02",  # 0.50 - 0.75 inches
        "#01c501",  # 0.75 - 1.00 inches
        "#008e00",  # 1.00 - 1.50 inches
        "#fdf802",  # 1.50 - 2.00 inches
        "#e5bc00",  # 2.00 - 2.50 inches
        "#fd9500",  # 2.50 - 3.00 inches
        "#fd0000",  # 3.00 - 4.00 inches
        "#d40000",  # 4.00 - 5.00 inches
        "#bc0000",  # 5.00 - 6.00 inches
        "#f800fd",  # 6.00 - 8.00 inches
        "#9854c6",  # 8.00 - 10.00 inches
        "#fdfdfd"   # 10.00+ 
    ]

    min_color = 0#0
    max_color = 16#320

    precip_colormap = mpl.colors.ListedColormap(nws_precip_colors)
 
    precip_colormap =  ocean.cm.rain

    precip_total_ticks = [0,2,4,6,8,10,12,14,16]#[0,50,100,150,200,250,300]
    precip_total_ticks_labels = list(map(str,precip_total_ticks)) #['0','50','100','150','200','250','300']

    bias_ticks = np.arange(-5,6,1)#[-50,-40,-30,-20,-10,0,10,20,30,40,50]
    bias_ticks_labels = list(map(str, bias_ticks))#['-50','-40','-30','-20','-10','0','10','20','30','40','50']
    print(bias_ticks_labels)

    cyclic_data, cyclic_lons = add_cyclic_point(total_speedy, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)
    projection = ccrs.PlateCarree(central_longitude=180)
    axes_class = (GeoAxes,dict(map_projection=projection))
    plt.rc('font', family='serif')

    fig = plt.figure(figsize=(16,13),constrained_layout=True)#figsize=(20,6))
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(3, 3),
                    axes_pad=0.9,
                    cbar_location='right', #'right',
                    cbar_mode='each',
                    cbar_pad=0.1,
                    cbar_size='3%',
                    label_mode='')  # note the empty label_mode


    ax1 = axgr[2]#plt.subplot(221,trained_layout="constrained"rojection=ccrs.PlateCarree())
    ax1.coastlines()

    ax1.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    ax1.set_yticks(np.linspace(-90, 90, 5), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax1.xaxis.set_major_formatter(lon_formatter)
    ax1.yaxis.set_major_formatter(lat_formatter)

    plot = ax1.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=min_color,vmax=max_color,cmap=precip_colormap)

    cbar = axgr.cbar_axes[0].colorbar(plot, extend='both')
    cbar.set_label(f'mm/d',fontsize=12)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_ticks(precip_total_ticks)
    cbar.set_ticklabels(precip_total_ticks_labels)

    ax1.set_title(f"SPEEDY\nMean Daily Precipitation Rates",fontsize=16,fontweight="bold")

    ###########
 
    cyclic_data, cyclic_lons = add_cyclic_point(total_hybrid, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

    ax2 = axgr[1] #plt.subplot(222,projection=ccrs.PlateCarree())
    ax2.coastlines()

    ax2.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    ax2.set_yticks(np.linspace(-90, 90, 5), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax2.xaxis.set_major_formatter(lon_formatter)
    ax2.yaxis.set_major_formatter(lat_formatter)

    plot_2 = ax2.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=min_color,vmax=max_color,cmap=precip_colormap)
    cbar2 = axgr.cbar_axes[1].colorbar(plot_2, extend='both')
    cbar2.set_label(f'mm/d',fontsize=12)
    cbar2.set_ticks(precip_total_ticks)
    cbar2.set_ticklabels(precip_total_ticks_labels)
    cbar2.ax.tick_params(labelsize=12)

    ax2.set_title(f"Hybrid Model\nMean Daily Precipitation Rates",fontsize=16,fontweight="bold")

    #####

    cyclic_data, cyclic_lons = add_cyclic_point(total_era5, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

    ax3 = axgr[0] #plt.subplot(223,projection=ccrs.PlateCarree())
    ax3.coastlines()
    
    ax3.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    ax3.set_yticks(np.linspace(-90, 90, 5), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax3.xaxis.set_major_formatter(lon_formatter)
    ax3.yaxis.set_major_formatter(lat_formatter)

    plot_3 = ax3.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=min_color,vmax=max_color,cmap=precip_colormap)

    cbar3 = axgr.cbar_axes[2].colorbar(plot_3,extend='both')
    cbar3.set_label(f'mm/d',fontsize=12)
    cbar3.set_ticks(precip_total_ticks)
    cbar3.set_ticklabels(precip_total_ticks_labels)
    cbar3.ax.tick_params(labelsize=12)

    ax3.set_title(f"ERA5\nMean Daily Precipitation Rates",fontsize=16,fontweight="bold")

    ###########BIAS
    min_color = -5#-200
    max_color = 5#200

    tarn = ocean.cm.tarn

    cyclic_data, cyclic_lons = add_cyclic_point(speedy_bias, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

    ax4 = axgr[5]#plt.subplot(224,projection=ccrs.PlateCarree())
    ax4.coastlines()

    ax4.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    ax4.set_yticks(np.linspace(-90, 90, 5), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax4.xaxis.set_major_formatter(lon_formatter)
    ax4.yaxis.set_major_formatter(lat_formatter)

    plot_4 = ax4.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=min_color,vmax=max_color,cmap=tarn)
    cbar4 = axgr.cbar_axes[3].colorbar(plot_4,extend='both')
    cbar4.set_label(f'mm/d',fontsize=12)
    cbar4.set_ticks(bias_ticks)
    cbar4.set_ticklabels(bias_ticks_labels)
    cbar4.ax.tick_params(labelsize=12)

    ax4.set_title(f"SPEEDY Bias",fontsize=16,fontweight="bold")


    ####
    cyclic_data, cyclic_lons = add_cyclic_point(hybrid_bias, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

    ax5 = axgr[4]#plt.subplot(224,projection=ccrs.PlateCarree())
    ax5.coastlines()

    ax5.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    ax5.set_yticks(np.linspace(-90, 90, 5), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax5.xaxis.set_major_formatter(lon_formatter)
    ax5.yaxis.set_major_formatter(lat_formatter)

    plot_5 = ax5.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=min_color,vmax=max_color,cmap=tarn)
    cbar5 = axgr.cbar_axes[4].colorbar(plot_5,extend='both')
    cbar5.set_label(f'mm/d',fontsize=12)
    cbar5.set_ticks(bias_ticks)
    cbar5.set_ticklabels(bias_ticks_labels)
    cbar5.ax.tick_params(labelsize=12)

    ax5.set_title(f"Hybrid Model Bias",fontsize=16,fontweight="bold")

    ####

    balance = ocean.cm.balance
    abs_bias_diff = abs(hybrid_bias) - abs(speedy_bias)
    print('Bias different (lower == better hybrid)',np.max(abs_bias_diff),np.min(abs_bias_diff),np.average(abs_bias_diff),np.std(abs_bias_diff))

    cyclic_data, cyclic_lons = add_cyclic_point(abs_bias_diff, coord=lons)
    lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

    ax6 = axgr[3]#plt.subplot(224,projection=ccrs.PlateCarree())
    ax6.coastlines()

    ax6.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    ax6.set_yticks(np.linspace(-90, 90, 5), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax6.xaxis.set_major_formatter(lon_formatter)
    ax6.yaxis.set_major_formatter(lat_formatter)

    plot_6 = ax6.pcolormesh(lons2d,lats2d,cyclic_data,transform=ccrs.PlateCarree(),vmin=min_color,vmax=max_color,cmap=balance)
    cbar6 = axgr.cbar_axes[5].colorbar(plot_6,extend='both')
    cbar6.set_label(f'mm/d',fontsize=12)
    cbar6.set_ticks(bias_ticks)
    cbar6.set_ticklabels(bias_ticks_labels)
    cbar6.ax.tick_params(labelsize=12)

    ax6.set_title(f"Abs Bias Diff (Hybrid - SPEEDY)",fontsize=16,fontweight="bold")


    ax7 = axgr[7]
    ax7.remove() 
    ax7 = fig.add_subplot(3, 3, (8))#plt.subplot(224,projection=ccrs.PlateCarree())

    ax7.loglog(percent_range,era_uwind_extreme, linewidth=2,label='ERA5')
    ax7.loglog(percent_range,hybrid_uwind_extreme,linewidth=2,label='Hybrid')

    print(percent_range)
    ax7.set_xlim([10**-3,100])
    ax7.set_ylim([10**-2,10**2])#np.max(era_uwind_extreme)])
    #ax1.set_xscale('log')

    x_labels_vals = [10**1,10**0,10**-1,10**-2,10**-3]
    x_tick_labels = ['90%','99%','99.9%','99.99%','99.999%']

    ax7.tick_params(axis='y', which='major', labelsize=12)
    ax7.invert_xaxis()
    ax7.set_xticks(x_labels_vals)
    ax7.set_xticklabels(x_tick_labels,fontsize=12)

    ax7.set_xlabel("Percentile",fontsize=12)

    ax7.set_ylabel("Total Precipitation\n(mm/6 h)",fontsize=12)
   
    ax7.legend(fontsize=12)

    ax8 = axgr[6]
    cb = axgr.cbar_axes[6].colorbar(plot_6)
    cb.remove()
    ax8.remove()

    ax9 = axgr[8]
    cb = axgr.cbar_axes[8].colorbar(plot_6)
    cb.remove()
    ax9.remove()

    cb_axes = get_colorbars(fig)
    print(cb_axes)

    #fig.delaxes(axgr[6])
    #axgr[6].set_visible(False)
   
    #fig.delaxes(axgr[8]) 
    #axgr[8].set_visible(False)
    #plt.subplots_adjust(top=0.982,bottom=0.018,left=0.043,right=0.979,hspace=0.2,wspace=0.2)
   
    plotdir = '/home/troyarcomano/FortranReservoir/vert_loc_hybridspeedy_leakage/plots/'

    #plt.savefig(f'{plotdir}annual_precip_totals_and_bias_sst_7day_biggercbar_range_mmday.pdf')
    #plt.close("all")
    #plt.subplot_tool()
    plt.subplots_adjust(top=1.0,bottom=0.064,left=0.053,right=0.953,hspace=0.896,wspace=0.039)
    plt.show()

total_precip_bias_plot()
#cdf_precip()
#histo_precip()
