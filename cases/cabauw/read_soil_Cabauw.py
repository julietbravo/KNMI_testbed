import matplotlib.pyplot as pl
import xarray as xr
import numpy as np
import datetime
from scipy import interpolate

def interpolate_soil(ds, z_obs, z_model, name, offset=0, plot=False):
    obs = np.zeros_like(z_obs, dtype=np.float)
    for k,z in enumerate(z_obs):
        obs[k] = ds['{0:}{1:02d}'.format(name, z)].values + offset

    print(z_obs, obs)

    # Remove nans..
    mask = np.invert(np.isnan(obs))
    z_obs = z_obs[mask]
    obs   = obs  [mask]

    print(z_obs, obs)

    # Interpolation function
    ip = interpolate.interp1d(z_obs, obs, kind='slinear', fill_value='extrapolate')

    # Interpolate
    model = ip(z_model)

    if (plot):
        pl.figure()
        pl.plot(obs, z_obs, '-x', label='obs')
        pl.plot(model, z_model, '-o', label='model')
        pl.legend()

    return model


def get_Tsoil_Cabauw(date, path, plot=False):
    # Read Cabauw soil data and select nearest time
    ds = xr.open_dataset('{0:}/cesar_soil_heat_lb1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path, date.year, date.month))
    ds = ds.sel(time=date, method='nearest')

    # Interpolate observations to model grid
    z_obs   = np.array([0,2,4,6,8,12,20,30,50])
    z_model = np.array([3.5, 17.5, 64, 195])
    name    = 'TS'

    T_model = interpolate_soil(ds, z_obs, z_model, name, offset=273.15, plot=plot)

    return T_model


def get_phisoil_Cabauw(date, path, plot=False):
    # Read Cabauw soil data and select nearest time
    ds = xr.open_dataset('{0:}/cesar_soil_water_lb1_t10_v1.1_{1:04d}{2:02d}.nc'.format(path, date.year, date.month))
    ds = ds.sel(time=date, method='nearest')

    # Interpolate observations to model grid
    z_obs   = np.array([5,19,33,40,56])
    z_model = np.array([3.5, 17.5, 64, 195])
    name    = 'TH'

    phi_model = interpolate_soil(ds, z_obs, z_model, name, offset=0, plot=plot)

    return phi_model

def array_to_string(arr):
    return str(arr).replace('[', '').replace(']', '')


if __name__ == '__main__':
    pl.close('all')

    #path = '/nobackup/users/stratum/Cabauw'
    path = '/scratch/ms/nl/nkbs/DOWA/Cabauw'
    date = datetime.datetime(year=2016, month=3, day=1, hour=0)

    Tsoil   = get_Tsoil_Cabauw  (date, path, True)
    phisoil = get_phisoil_Cabauw(date, path, True)
