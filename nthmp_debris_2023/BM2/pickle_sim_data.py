
from pylab import *
import os
import netCDF4
from scipy.interpolate import RegularGridInterpolator
import pickle

sim_data_dir = 'sim_data_BM2'

eta_fname = os.path.join(sim_data_dir, 'zeta.nc')
f = netCDF4.Dataset(eta_fname, 'r')
time_numerical = array(f.variables['time'])
x_numerical = array(f.variables['x'])
zeta = squeeze(array(f.variables['zeta']))  # squeeze to get rid of y-dimension index
depth = array(f.variables['depth'])

# numerical shift values - numerical domain was not same as physical domain
x_shift=-10.  # added 10-m to offshore boundary to create wave
t_shift=13.1  # cut off first 13 ish seconds since nothing happened

time=time_numerical+t_shift
x=x_numerical+x_shift


velo_fname = os.path.join(sim_data_dir, 'velo.nc')
f = netCDF4.Dataset(velo_fname, 'r')
u_vel = squeeze(array(f.variables['u_velo']))

blvs_fname = os.path.join(sim_data_dir, 'blvs.nc')
f = netCDF4.Dataset(blvs_fname, 'r')
wet_dry = squeeze(array(f.variables['boundary_boolean']))
nu_breaking = squeeze(array(f.variables['eddy_viscosity']))

nu_breaking = where(wet_dry==99, nan, nu_breaking)
zeta = where(wet_dry==99, nan, zeta)
u_vel = where(wet_dry==99, nan, u_vel)

zeta_fcn = RegularGridInterpolator((time,x), zeta, method='linear',
                bounds_error=False)
u_vel_fcn = RegularGridInterpolator((time,x), u_vel, method='linear',
                bounds_error=False)


picklefile = 'sim_data.pickle'
data = {'time':time, 'zeta_fcn':zeta_fcn, 'u_vel_fcn':u_vel_fcn}
with open(picklefile, 'wb') as f:
    pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

print('Created ',picklefile)
