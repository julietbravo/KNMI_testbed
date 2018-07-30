import matplotlib.pyplot as pl
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

from read_DDH_netcdf import *

def read_LEG_EPL(path):
    # Read and clean header
    f     = open(path, 'r', encoding="ISO-8859-1")
    info  = f.readline(). split(';')
    names = f.readline().split(';')
    f.close()

    for i in range(len(names)):
        names[i] = names[i].replace('LEG_', '')
        names[i] = names[i].replace('EPL_', '')

    return pd.read_csv(path, skiprows=3, header=1, sep=';', index_col=0, names=names, parse_dates=True, na_values=9999)

def mask_nan(array):
    return np.ma.masked_invalid(array)

def abs_speed(u, v):
    return (u**2 + v**2)**0.5

def components(speed, direction):
    u = -speed * np.sin(np.deg2rad(direction))
    v = -speed * np.cos(np.deg2rad(direction))
    return u, v


if __name__ == '__main__':

    if 'data_is_loaded' not in locals():
        data_is_loaded=True

        # Period to study
        start = datetime.datetime(year=2017, month=1, day=11, hour=0)
        end   = datetime.datetime(year=2017, month=1, day=13, hour=0)

        # Harmonie data from DDH files
        # ----------------------------
        path = '/Users/bart/meteo/data/Harmonie_DDH/'
        ham  = read_DDH_netcdf(start, end, path)

        # Lichteiland Goeree observations
        # -------------------------------
        leg = read_LEG_EPL('/Users/bart/meteo/data/offshore_wind/LEG/statistics_2017.csv')
        leg = leg.loc[start:end]

        # Lichteiland Goeree observations
        # -------------------------------
        epl = read_LEG_EPL('/Users/bart/meteo/data/offshore_wind/EPL/statistics_2017.csv')
        epl = epl.loc[start:end]

        # DALES
        # ------------------
        da1 = xr.open_dataset('profiles.001.nc')
        da1_time = np.array([start + datetime.timedelta(seconds=int(t)) for t in da1.time.values])

        da2 = xr.open_dataset('profiles.002.nc')
        da2_time = np.array([start + datetime.timedelta(seconds=int(t)) for t in da2.time.values])

    # Harmonie domain
    iloc = 2
    z = 291     # in {63,91,116,141,166,191,216,241,266,291}

    pl.close('all')
    pl.figure()

    epl_u, epl_v = mask_nan(components(epl['H{0:03d}_WsHor_avg'.format(z)], epl['H{0:03d}_Wd'.format(z)]))

    # Interpolation functions
    da1_u = interp1d(da1.zt.values, da1.u.values, axis=1)
    da1_v = interp1d(da1.zt.values, da1.v.values, axis=1)

    da2_u = interp1d(da2.zt.values, da2.u.values, axis=1)
    da2_v = interp1d(da2.zt.values, da2.v.values, axis=1)

    ham_u = interp1d(ham.z[0,iloc,:].values, ham.u.values[:,iloc,:], axis=1)
    ham_v = interp1d(ham.z[0,iloc,:].values, ham.v.values[:,iloc,:], axis=1)

    pl.subplot(221)
    pl.plot(epl.index, epl_u, '.', label='EPL')
    pl.plot(ham.time,  ham_u(z), label='Harmonie')
    pl.plot(da1_time,  da1_u(z), label='DALES')
    pl.plot(da2_time,  da2_u(z), label='DALES-nn')
    pl.legend()

    pl.subplot(222)
    pl.plot(epl.index, epl_v, '.', label='EPL')
    pl.plot(ham.time,  ham_v(z), label='Harmonie')
    pl.plot(da1_time,  da1_v(z), label='DALES')
    pl.plot(da2_time,  da2_v(z), label='DALES-nn')
    pl.legend()

    pl.subplot(223)
    pl.plot(ham.time,  ham['dtu_dyn'][:,iloc_epl,k]*3600., label='dyn')
    pl.plot(ham.time,  ham['dtu_phy'][:,iloc_epl,k]*3600., label='phy')
    pl.legend()

    pl.subplot(224)
    pl.plot(ham.time,  ham['dtv_dyn'][:,iloc_epl,k]*3600., label='dyn')
    pl.plot(ham.time,  ham['dtv_phy'][:,iloc_epl,k]*3600., label='phy')
    pl.legend()
