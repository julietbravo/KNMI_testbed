import matplotlib.pyplot as pl
import xarray as xr
import numpy as np
import datetime

# Custom script to read the FINO1 ASCII files as a Pandas dataframe
from read_FINO1 import read_FINO1

pl.close('all')

def abs(u,v):
    return (u**2 + v**2)**0.5

def components(speed, direction):
    u = -speed * np.sin(np.deg2rad(direction))
    v = -speed * np.cos(np.deg2rad(direction))
    return u, v

def esat(T):
    return 0.611e3 * np.exp(17.2694 * (T - 273.16) / (T - 35.86))

def qsat(T,p):
    return 0.622 * esat(T) / p

def interpz(z1, z2, zg, v1, v2):
    f2 = (zg - z1) / (z2 - z1)
    f1 = 1 - f2
    return f1*v1 + f2*v2


if __name__ == '__main__':

    start = datetime.datetime(year=2017, month=6, day=11, hour=0)
    end   = datetime.datetime(year=2017, month=6, day=13, hour=0)

    if 'f1' not in locals():
        # Read the FINO1 observations
        path = '/nobackup/users/stratum/FINO1_obs/'
        f1   = read_FINO1(path, start, end)

        # Fix units.....
        for z in [30, 40, 50, 70, 100]:
            f1['T_{}'.format(z) ] += 273.15      # Celcius -> Kelvin
        for z in [33, 40, 50, 70, 90]:
            f1['RH_{}'.format(z)] /= 100         # Procent -> fraction
        for z in [20, 90]:
            f1['p_{}'.format(z)] *= 100          # hPa -> Pa

        # Some (linear) interpolation to get obs at the same heights...
        f1['T_90'] = interpz(70, 100, 90, f1['T_70'].values, f1['T_100'].values)
        f1['T_33'] = interpz(30, 40,  33, f1['T_30'].values, f1['T_40' ].values)

        # Estimate pressure, and calculate specific humidity
        for z in [33, 40, 50, 70]:
            f1['p_{}'.format(z)] = f1['p_20'] * np.exp(-0.0342/f1['T_30']*(z-20))
        for z in [33, 40, 50, 70, 90]:
            f1['q_{}'.format(z)] = f1['RH_{}'.format(z)] * qsat(f1['T_{}'.format(z)], f1['p_{}'.format(z)])

        # Wind speed + direction to component/scratch/ms/nl/nkbs/DALES/KNMI_testbed/FINO1_tests
        for z in [33, 40, 50, 60, 70, 80, 90]:
            f1['u_{}'.format(z)], f1['v_{}'.format(z)] = components(f1['U_{}'.format(z)], f1['Udir_{}'.format(z)])


    if 'd1' not in locals():
        # Read the LES experiments
        d1   = xr.open_dataset('profiles.001.nc')
        d1_time = [start + datetime.timedelta(seconds=time) for time in d1.time.values]

        d1['Ta'] = d1['thl'] * (d1['presh'] / 1e5)**(287.05/1004.)

        d1t   = xr.open_dataset('tmser.001.nc')
        d1t_time = [start + datetime.timedelta(seconds=time) for time in d1t.time.values]



    pl.figure()
    ax=pl.subplot(111)
    pl.plot(f1.index, f1['Udir_90'])
    lim = ax.get_xlim()
    pl.fill_between(lim, [235, 235], [285, 285], color='g', alpha=0.3)
    pl.fill_between(lim, [0, 0],     [180, 180], color='r', alpha=0.3)
    pl.ylim(0,360)
    pl.xlim(lim)


    pl.figure()

    pl.subplot(221)
    pl.plot(f1.index, f1['u_90'], 'o', color='k', markersize=3)
    pl.plot(d1_time, d1['u'][:,4], 'k-')

    pl.subplot(222)
    pl.plot(f1.index, f1['v_90'], 'o', color='k', markersize=3)
    pl.plot(d1_time, d1['v'][:,4], 'k-')

    pl.subplot(223)
    pl.plot(f1.index, f1['T_90'], 'o', color='k', markersize=3)
    pl.plot(d1_time, d1['Ta'][:,4], 'k-')

    pl.subplot(224)
    pl.plot(f1.index, f1['q_90'], 'o', color='k', markersize=3)
    pl.plot(d1_time, d1['qt'][:,4], 'k-')






#pl.subplot(132)
#pl.plot(f1.index, f1['T_30']+273.15, 'o', color='k')
#pl.plot(d1_time, d1['Ta'][:,1], 'k-')
#
#pl.plot(f1.index, f1['T_100']+273.15, 'o', color='r')
#pl.plot(d1_time, d1['Ta'][:,4], 'r-')

#pl.subplot(133)
#pl.plot(f1.index, f1['Q_33'], 'o', color='k')
#pl.plot(d1_time, d1['qt'][:,1], 'k-')

#pl.plot(f1.index, f1['Q_90']+273.15, 'o', color='r')
#pl.plot(d1_time, d1['qt'][:,4], 'r-')
