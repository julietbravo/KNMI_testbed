import xarray as xr
import numpy as np
import datetime

def open_soil(date, lon, lat, path):
    file = 'soil_{0:04d}{1:02d}.nc'.format(date.year, date.month)
    ds   = xr.open_dataset('{}/{}'.format(path, file))
    return ds.sel(longitude=lon, latitude=lat, time=date, method='nearest')

def get_Tsoil(date, lon, lat, path):
    ds = open_soil(date, lon, lat, path)
    return np.array([ds.stl1.values, ds.stl2.values, ds.stl3.values, ds.stl4.values])

def get_qsoil(date, lon, lat, path):
    ds = open_soil(date, lon, lat, path)
    return np.array([ds.swvl1.values, ds.swvl2.values, ds.swvl3.values, ds.swvl4.values])

def array_to_string(arr):
    return str(arr).replace('[', '').replace(']', '')


if __name__ == '__main__':

    start = datetime.datetime(year=2017, month=2, day=1, hour=12)
    path  = '/nobackup/users/stratum/ERA5/soil/'

    Tsoil = get_Tsoil(start, 4.9, 51.97, path)
    qsoil = get_qsoil(start, 4.9, 51.97, path)
