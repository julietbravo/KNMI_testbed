import matplotlib.pyplot as pl
import matplotlib.dates as mdates

import xarray as xr
import pandas as pd
import numpy as np
import datetime

# Custom script to read the FINO1 ASCII files as a Pandas dataframe
from read_FINO1 import read_FINO1
# Custom script to read Harmonie NetCDF LEs forcings
from read_DDH_netcdf import read_DDH_netcdf

pl.close('all')

def windspeed(u,v):
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

class Statistics:
    def __init__(self, obs, model, scale=1):
        self.mean_diff = (model - obs).mean() * scale
        self.max_diff  = (model - obs).max() * scale
        self.rmse      = np.sqrt(((model - obs)**2).mean()) * scale

# Data axis formatting...
def format_date():
    days     = mdates.DayLocator()
    days_fmt = mdates.DateFormatter('%d %b')

    ax = pl.gca()
    ax.xaxis.set_major_locator(days)
    ax.xaxis.set_major_formatter(days_fmt)


if __name__ == '__main__':

    start = datetime.datetime(year=2017, month=6, day=11, hour=0)
    end   = datetime.datetime(year=2017, month=6, day=18, hour=0)

    # Read the FINO1 observations
    if 'f1' not in locals():
        path = '/nobackup/users/stratum/FINO1_obs/'
        #path = '/Users/bart/meteo/data/offshore_wind/FINO1/'
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


    # Read the LES experiments
    if 'd1' not in locals():
        d1   = xr.open_dataset('profiles.001.nc')
        d1_time = [start + datetime.timedelta(seconds=int(time)) for time in d1.time.values]

        d1['Ta'] = d1['thl'] * (d1['presh'] / 1e5)**(287.05/1004.)

        d1t   = xr.open_dataset('tmser.001.nc')
        d1t_time = [start + datetime.timedelta(seconds=int(time)) for time in d1t.time.values]

        d2   = xr.open_dataset('profiles.002.nc')
        d2_time = [start + datetime.timedelta(seconds=int(time)) for time in d2.time.values]

        d2['Ta'] = d2['thl'] * (d2['presh'] / 1e5)**(287.05/1004.)

        d2t   = xr.open_dataset('tmser.002.nc')
        d2t_time = [start + datetime.timedelta(seconds=int(time)) for time in d2t.time.values]


    # Read the Harmonie LES forcings
    if 'h1' not in locals():
        path = '/nobackup/users/stratum/DOWA/LES_forcing/'
        #path = '/Users/bart/meteo/data/Harmonie_LES_forcing/'
        h1 = read_DDH_netcdf(start, end, path)
        iloc = 0


    # Put main model variables in FINO1 data frame
    if 'df' not in locals():
        data_les = {'LES_u_90': d1['u'][:,4], 'LES_v_90': d1['v'][:,4], 'LES_T_90': d1['Ta'][:,4], 'LES_q_90': d1['qt'][:,4]}
        df_les = pd.DataFrame(data_les, index=d1_time)

        data_ham = {'HAM_u_90': h1['u'][:,iloc,3], 'HAM_v_90': h1['v'][:,iloc,3], 'HAM_T_90': h1['T'][:,iloc,3], 'HAM_q_90': h1['q'][:,iloc,3]}
        df_ham = pd.DataFrame(data_ham, index=h1.time)
        df_ham.index = df_ham.index.round('min')

        # Merge them, and drop missing records (e.g. LES has 5min output while Harmonie has 10min...)
        df = pd.concat([f1, df_les, df_ham], axis=1, join='inner')


    # Statistics!
    if True:
        u_90_LES = Statistics(df['u_90'], df['LES_u_90'])
        u_90_HAM = Statistics(df['u_90'], df['HAM_u_90'])

        v_90_LES = Statistics(df['v_90'], df['LES_v_90'])
        v_90_HAM = Statistics(df['v_90'], df['HAM_v_90'])

        T_90_LES = Statistics(df['T_90'], df['LES_T_90'])
        T_90_HAM = Statistics(df['T_90'], df['HAM_T_90'])

        q_90_LES = Statistics(df['q_90'], df['LES_q_90'], scale=1000)
        q_90_HAM = Statistics(df['q_90'], df['HAM_q_90'], scale=1000)

        # Absolute wind
        df['LES_U_90'] = windspeed(df['LES_u_90'], df['LES_v_90'])
        df['HAM_U_90'] = windspeed(df['HAM_u_90'], df['HAM_v_90'])
        U_90_LES = Statistics(df['U_90'], df['LES_U_90'])
        U_90_HAM = Statistics(df['U_90'], df['HAM_U_90'])

        print('Absolute wind speeds:')
        print('LES: mean diff={0:.2f}, rmse={1:.2f}'.format(U_90_LES.mean_diff, U_90_LES.rmse))
        print('HAM: mean diff={0:.2f}, rmse={1:.2f}'.format(U_90_HAM.mean_diff, U_90_HAM.rmse))





    if False:
        pl.figure()
        ax=pl.subplot(111)
        pl.plot(f1.index, f1['Udir_90'])
        lim = ax.get_xlim()
        pl.fill_between(lim, [235, 235], [285, 285], color='g', alpha=0.3)
        pl.fill_between(lim, [0, 0],     [180, 180], color='r', alpha=0.3)
        pl.ylim(0,360)
        pl.xlim(lim)



    if True:
        pl.figure(figsize=(11,7))

        pl.subplot(221)
        pl.title('z = 90 m: ', loc='left')
        pl.plot(f1.index, f1['u_90'], '+', color='r', markersize=2, label='FINO1 obs')
        pl.plot(d1_time, d1['u'][:,4], 'k-',       label='LES: mean diff={0:.2f}, rmse={1:.2f}'.format(u_90_LES.mean_diff, u_90_LES.rmse))
        #pl.plot(d2_time, d2['u'][:,4], 'k--')
        pl.plot(h1.time, h1['u'][:,iloc,3], 'b--', label='HAM: mean diff={0:.2f}, rmse={1:.2f}'.format(u_90_HAM.mean_diff, u_90_HAM.rmse))
        format_date()
        pl.legend(fontsize='small')
        pl.ylabel('u (m s-1)')

        pl.subplot(222)
        pl.plot(f1.index, f1['v_90'], 'x', color='r', markersize=2, label='FINO1 obs')
        pl.plot(d1_time, d1['v'][:,4], 'k-',       label='LES: mean diff={0:.2f}, rmse={1:.2f}'.format(v_90_LES.mean_diff, v_90_LES.rmse))
        #pl.plot(d2_time, d2['v'][:,4], 'k--')
        pl.plot(h1.time, h1['v'][:,iloc,3], 'b--', label='HAM: mean diff={0:.2f}, rmse={1:.2f}'.format(v_90_HAM.mean_diff, v_90_HAM.rmse))
        format_date()
        pl.legend(fontsize='small')
        pl.ylabel('v (m s-1)')

        pl.subplot(223)
        pl.plot(f1.index, f1['T_90'], 'x', color='r', markersize=2, label='FINO1 obs')
        pl.plot(d1_time, d1['Ta'][:,4], 'k-',      label='LES: mean diff={0:.2f}, rmse={1:.2f}'.format(T_90_LES.mean_diff, T_90_LES.rmse))
        #pl.plot(d2_time, d2['Ta'][:,4], 'k--')
        pl.plot(h1.time, h1['T'][:,iloc,3], 'b--', label='HAM: mean diff={0:.2f}, rmse={1:.2f}'.format(T_90_HAM.mean_diff, T_90_HAM.rmse))
        format_date()
        pl.legend(fontsize='small')
        pl.ylabel('T (K)')

        pl.subplot(224)
        pl.plot(f1.index, f1['q_90']*1e3, 'x', color='r', markersize=2, label='FINO1 obs')
        pl.plot(d1_time, d1['qt'][:,4]*1e3, 'k-',      label='LES: mean diff={0:.2f}, rmse={1:.2f}'.format(q_90_LES.mean_diff, q_90_LES.rmse))
        #pl.plot(d2_time, d2['qt'][:,4]*1e3, 'k--')
        pl.plot(h1.time, h1['q'][:,iloc,3]*1e3, 'b--', label='HAM: mean diff={0:.2f}, rmse={1:.2f}'.format(q_90_HAM.mean_diff, q_90_HAM.rmse))
        format_date()
        pl.legend(fontsize='small')
        pl.ylabel('q (g kg-1)')

        pl.tight_layout()

