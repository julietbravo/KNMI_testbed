import matplotlib.pyplot as pl
import xarray as xr
import pandas as pd
import numpy as np
import datetime

from read_DDH_netcdf import *
from read_FINO1 import *

pl.close('all')

path_hm = '/Users/bart/meteo/data/Harmonie_LES_forcing/'
path_cb = '/Users/bart/meteo/observations/Cabauw'
path_f1 = '/Users/bart/meteo/data/offshore_wind/FINO1/'


#path_hm = '/nobackup/users/stratum/DOWA/LES_forcing/'
#path_f1 = '/nobackup/users/stratum/FINO1_obs/'

class Statistics:
    def __init__(self, obs, model, scale=1):
        self.mean_diff = (model - obs).mean() * scale
        self.max_diff  = (model - obs).max() * scale
        self.rmse      = np.sqrt(((model - obs)**2).mean()) * scale

def esat(T):
    return 0.611e3 * np.exp(17.2694 * (T - 273.16) / (T - 35.86))

def qsat(T,p):
    return 0.622 * esat(T) / p

def components(speed, direction):
    u = -speed * np.sin(np.deg2rad(direction))
    v = -speed * np.cos(np.deg2rad(direction))
    return u, v

def interpz(data, z, zg):
    """
    Interpolate `data` to goal_height (`zg`) for
    time dependent vertical profiles (`z`)
    """

    z = z.values
    data = data.values

    # Find height index of nearest level below `zg`
    k0 = np.abs(z - zg).argmin(axis=1)
    z0 = z[:,k0][:,0]
    k0[z0 > zg] -= 1

    # Heights above and below interpolation height
    z1 = z[:,k0+1][:,0]
    z0 = z[:,k0]  [:,0]
    dz = z1 - z0

    # Interpolation factors
    f1 = (zg - z0) / dz
    f0 = 1-f1

    # Interpolate!
    return data[:,k0][:,0]*f0 + data[:,k0+1][:,0]*f1

def interpz_simple(z1, z2, zg, v1, v2):
    f2 = (zg - z1) / (z2 - z1)
    f1 = 1 - f2
    return f1*v1 + f2*v2

# -----------------
# Period
# -----------------
#for m in range(1,13):

start = datetime.datetime(year=2017, month=6, day=1, hour=3)
end   = datetime.datetime(year=2017, month=7, day=1, hour=0)

# -----------------
# Read Harmonie statistics
# -----------------
iloc = 7
variables = ['time', 'z', 'q']
if 'hm' not in locals():
    hm  = read_DDH_netcdf(start, end, path_hm, variables)

    # Interpolate to observation heights Cabauw, and store in DataFrame
    data = {'hm_q33' : interpz(hm.q[:,iloc,:], hm.z[:,iloc,:], 33. )*1e3,
            'hm_q40' : interpz(hm.q[:,iloc,:], hm.z[:,iloc,:], 40. )*1e3,
            'hm_q50' : interpz(hm.q[:,iloc,:], hm.z[:,iloc,:], 50. )*1e3,
            'hm_q70' : interpz(hm.q[:,iloc,:], hm.z[:,iloc,:], 70. )*1e3,
            'hm_q90' : interpz(hm.q[:,iloc,:], hm.z[:,iloc,:], 90. )*1e3}

    df_hm = pd.DataFrame(data, index=hm.time)
    df_hm.index = df_hm.index.round('min')
    df_hm = df_hm.loc[start:end]

# -----------------
# Read FINO1 statistics
# -----------------
if 'f1' not in locals():
    f1   = read_FINO1(path_f1, start, end)

    # Fix units.....
    for z in [30, 40, 50, 70, 100]:
        f1['T_{}'.format(z) ] += 273.15      # Celcius -> Kelvin
    for z in [33, 40, 50, 70, 90]:
        f1['RH_{}'.format(z)] /= 100         # Procent -> fraction
    for z in [20, 90]:
        f1['p_{}'.format(z)] *= 100          # hPa -> Pa

    # Some (linear) interpolation to get obs at the same heights...
    f1['T_90'] = interpz_simple(70, 100, 90, f1['T_70'].values, f1['T_100'].values)
    f1['T_33'] = interpz_simple(30, 40,  33, f1['T_30'].values, f1['T_40' ].values)

    # Estimate pressure, and calculate specific humidity
    for z in [33, 40, 50, 70]:
        f1['p_{}'.format(z)] = f1['p_20'] * np.exp(-0.0342/f1['T_30']*(z-20))
    for z in [33, 40, 50, 70, 90]:
        f1['q_{}'.format(z)] = f1['RH_{}'.format(z)] * qsat(f1['T_{}'.format(z)], f1['p_{}'.format(z)])

    # Wind speed + direction to component/scratch/ms/nl/nkbs/DALES/KNMI_testbed/FINO1_tests
    for z in [33, 40, 50, 60, 70, 80, 90]:
        f1['u_{}'.format(z)], f1['v_{}'.format(z)] = components(f1['U_{}'.format(z)], f1['Udir_{}'.format(z)])


# -----------------
# Merge them!
# -----------------
df = pd.concat([df_hm, f1], axis=1, join='inner')

for z in [33, 40, 50, 70, 90]:
    stat = Statistics(df['q_{}'.format(z)]*1000, df['hm_q{}'.format(z)])
    print('z={0:5.0f} m: mean_diff={1:6.2f}, rmse={2:6.2f}'.format(z, stat.mean_diff, stat.rmse))







#



#pl.figure()
#pl.subplot(211)
#pl.plot(df.index, df['q_hm'], label='HAM')
#pl.plot(df.index, df['q_cb'], label='Cabauw')
#pl.legend()
#
#
#pl.subplot(212)
#pl.plot(df.index, df['q_hm']-df['q_cb'])
