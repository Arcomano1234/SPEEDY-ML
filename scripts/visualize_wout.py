import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import glob

def get_files(path,pattern):
    files = glob.glob(path+pattern)
    return files

def plot_wout(filename):
    nc = Dataset(filename)

    wout = nc['wout'][:]
    print(np.shape(wout)) 
    chunk_size = np.shape(wout)[1]
    print(chunk_size)
    

    for i in range(chunk_size):
        print(i,wout[i,i])

    plt.title(filename[-150:-1],fontsize=8)
    plt.pcolormesh(wout[0:chunk_size,0:chunk_size],vmin=-1,vmax=1,cmap='seismic')
    plt.colorbar()
    plt.show() 
    #return wout
    #plt.savefig(filename[0:-2]+'png')
    #plt.close("all")


#files = get_files('/scratch/user/troyarcomano/FortranReservoir/hybridspeedy/wout_tests/','worker_954_wout_6000_10_10_10_model_noise_0.0_prior_transbeta_res1.000000000000000E-003*.nc')

files = ['region_218_level_01wout_6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_ocean_model_false_log_precip_timestep_1hour_leak_0.17_sigma0.083.nc','region_218_level_01wout_6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_ocean_model_false_log_precip_timestep_1hour_leak_0.17.nc','region_218_level_01wout_6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_1_vertlap_0_ocean_model_false_log_precip_timestep_1hour_leak_0.17.nc']
wout = np.zeros((len(files),6116, 132))
i = 0 
for filename in files:
    #wout[i,:,:] = plot_wout(filename)
    plot_wout(filename)
    i += 1 

chunk_size = 132

diff = wout[0,0:chunk_size,0:chunk_size] - wout[1,0:chunk_size,0:chunk_size]

plt.title('Diff')
plt.pcolormesh(diff,vmin=-0.3,vmax=0.3,cmap='seismic')
plt.colorbar()
plt.show()
'''   
filename = 'worker_954_wout_6000_10_10_10_model_noise_0.0001_beta_res 1.00000000000000_beta_model 10.0000000000000_prior0.000000000000000E+000.nc'

nc = Dataset(filename)

wout = nc['wout'][:]

chunk_size = 132

for i in range(chunk_size):
    print(i,wout[i,i])
plt.pcolormesh(wout[0:chunk_size,0:chunk_size],vmin=-2,vmax=2,cmap='seismic')
plt.colorbar()
plt.show()

filename1 = 'res_prediction_era6000node_25_25_20_noise_beta_0.0001_degree6_cylcing12hr_sigma0.5_yes_regional_stand_yes_regional_radius0.3_0.7_all_heights_and_logp_standardized_persistencetrial_11_01_1986_00.nc'

filename2 = 'res_prediction_era6000node_25_25_20_noise_beta_0.0001_degree6_cylcing12hr_sigma0.5_yes_regional_stand_yes_regional_radius0.3_0.7_noraml_logp_standardizedtrial_11_01_1986_00.nc'

filename3 = 'era_truth6000node_25_25_20_noise_beta_0.0001_degree6_cylcing12hr_sigma0.5_yes_regional_stand_yes_regional_radius0.3_0.7_persistence_beta_10000_logp_standardizedtrial_11_01_1986_00.nc'

nc_persistence = Dataset(filename1)

nc_normal = Dataset(filename2)

nc_truth = Dataset(filename3)

temp_persistence = nc_persistence['Temperature'][:]

temp_normal = nc_normal['Temperature'][:]

temp_truth = nc_truth['Temperature'][:]

error_persistence = temp_persistence - temp_truth

error_normal = temp_normal - temp_truth

diff_in_error = abs(error_persistence) - abs(error_normal)

lats = nc_truth['Lat'][:]
lons = nc_truth['Lon'][:]

cyclic_data, cyclic_lons = add_cyclic_point(diff_in_error[1,2,:,:], coord=lons)
lons2d,lats2d = np.meshgrid(cyclic_lons,lats)

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()

plt.pcolormesh(lons2d,lats2d,cyclic_data,vmin=-1,vmax=1,cmap='seismic',transform=ccrs.PlateCarree())
plt.colorbar()
plt.title('Error Difference Persistence Wout - Control\n200 hPa V-wind')
plt.show()
'''
