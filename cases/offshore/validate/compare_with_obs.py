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

# Data axis formatting...
def format_date_day():
    days     = mdates.DayLocator()
    days_fmt = mdates.DateFormatter('%d %b')

    ax = pl.gca()
    ax.xaxis.set_major_locator(days)
    ax.xaxis.set_major_formatter(days_fmt)

def format_date_hour():
    hours     = mdates.HourLocator()
    hours_fmt = mdates.DateFormatter('%H:%M')

    ax = pl.gca()
    ax.xaxis.set_major_locator(hours)
    ax.xaxis.set_major_formatter(hours_fmt)

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

def moving_av(x, N):
    return np.convolve(x, np.ones((N,))/N, mode='same')

class Statistics:
    def __init__(self, obs, model, label='', scale=1):
        self.mean_diff = (model - obs).mean() * scale
        self.max_diff  = (model - obs).max() * scale
        self.rmse      = np.sqrt(((model - obs)**2).mean()) * scale

        print('{0:10s}: Bias = {1:.3f}, RMSE = {2:.3f}'.format(label, self.mean_diff, self.rmse))

class Read_LES:
    def __init__(self, profiles, tmser, start):

        # Read profiles and convert time from seconds to date
        self.prof = xr.open_dataset(profiles)
        self.time = [start + datetime.timedelta(seconds=int(time)) for time in self.prof.time.values]

        # Some derived quantities
        self.prof['Ta'] = self.prof['thl'] * (self.prof['presh'] / 1e5)**(287.05/1004.)

        # Read time series
        self.tmser = xr.open_dataset(tmser)

        # Put some variables in Pandas data frame for statistical comparison
        data = {'LES_u_90': self.prof['u'][:,4], 'LES_v_90': self.prof['v'][:,4], 'LES_T_90': self.prof['Ta'][:,4], 'LES_q_90': self.prof['qt'][:,4]}
        self.df = pd.DataFrame(data, index=self.time)
        self.df['LES_U_90'] = windspeed(self.df['LES_u_90'], self.df['LES_v_90'])

    def calc_stat(self, obs, harm):
        # Merge all data frames
        df = pd.concat([obs, self.df, harm], axis=1, join='inner')

        self.U_90 = Statistics(df['Uc_90'], df['LES_U_90'], label='LES, U90')
        self.u_90 = Statistics(df['uc_90'], df['LES_u_90'], label='LES, u90')
        self.v_90 = Statistics(df['vc_90'], df['LES_v_90'], label='LES, v90')
        self.T_90 = Statistics(df['T_90'],  df['LES_T_90'], label='LES, v90')
        self.q_90 = Statistics(df['q_90'],  df['LES_q_90'], label='LES, q90', scale=1000)


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

    ## Read the Harmonie LES forcings
    if 'h1' not in locals():
        path = '/nobackup/users/stratum/DOWA/LES_forcing/'
        #path = '/Users/bart/meteo/data/Harmonie_LES_forcing/'
        h1   = read_DDH_netcdf(start, end, path)
        iloc = 0

        data_ham = {'HAM_u_90': h1['u'][:,iloc,3], 'HAM_v_90': h1['v'][:,iloc,3], 'HAM_T_90': h1['T'][:,iloc,3], 'HAM_q_90': h1['q'][:,iloc,3]}
        df_ham = pd.DataFrame(data_ham, index=h1.time)
        df_ham.index = df_ham.index.round('min')


    # Read the LES cases
    if 'e1' not in locals():
        e1 = Read_LES('profiles.001.nc',  'tmser.001.nc',  start)
        e2 = Read_LES('profiles.30m.nc',  'tmser.30m.nc',  start)
        e3 = Read_LES('profiles.60m.nc',  'tmser.60m.nc',  start)
        e4 = Read_LES('profiles.120m.nc', 'tmser.120m.nc', start)

        les_runs = [e1, e2, e3, e4]
        for r in les_runs:
            r.calc_stat(f1, df_ham)

    # Calculate Harmonie statistics
    df = pd.concat([f1, e1.df, df_ham], axis=1, join='inner')
    df['HAM_U_90'] = windspeed(df['HAM_u_90'], df['HAM_v_90'])
    U_90_HAM = Statistics(df['Uc_90'], df['HAM_U_90'], label='HAM, U90')
    u_90_HAM = Statistics(df['uc_90'], df['HAM_u_90'], label='HAM, u90')
    v_90_HAM = Statistics(df['vc_90'], df['HAM_v_90'], label='HAM, v90')
    T_90_HAM = Statistics(df['T_90'],  df['HAM_T_90'], label='HAM, T90')
    q_90_HAM = Statistics(df['q_90'],  df['HAM_q_90'], label='HAM, q90', scale=1000, )



    if False:
        pl.figure()
        ax=pl.subplot(111)
        pl.plot(f1.index, f1['Udir_90'])
        lim = ax.get_xlim()
        pl.fill_between(lim, [235, 235], [285, 285], color='g', alpha=0.3)
        pl.fill_between(lim, [0, 0],     [180, 180], color='r', alpha=0.3)
        pl.ylim(0,360)
        pl.xlim(lim)



    if False:
        pl.figure(figsize=(10,6))

        pl.subplot(221)
        pl.title('z = 90 m: ', loc='left')
        pl.plot(f1.index, f1['u_90'], '+', color='r', markersize=2, label='FINO1 obs')
        pl.plot(e1.time, e1.prof['u'][:,4], 'k-',  label='LES: mean diff={0:.2f} m/s, rmse={1:.2f} m/s'.format(e1.u_90.mean_diff, e1.u_90.rmse))
        pl.plot(h1.time, h1['u'][:,iloc,3], 'b--', label='HAM: mean diff={0:.2f} m/s, rmse={1:.2f} m/s'.format(u_90_HAM.mean_diff, u_90_HAM.rmse))
        format_date_day()
        pl.legend(fontsize='small', loc=3)
        pl.ylabel('u (m s-1)')
        pl.ylim(-15,20)

        pl.subplot(222)
        pl.plot(f1.index, f1['v_90'], '+', color='r', markersize=2, label='FINO1 obs')
        pl.plot(e1.time, e1.prof['v'][:,4], 'k-',  label='LES: mean diff={0:.2f} m/s, rmse={1:.2f} m/s'.format(e1.v_90.mean_diff, e1.v_90.rmse))
        pl.plot(h1.time, h1['v'][:,iloc,3], 'b--', label='HAM: mean diff={0:.2f} m/s, rmse={1:.2f} m/s'.format(v_90_HAM.mean_diff, v_90_HAM.rmse))
        format_date_day()
        pl.legend(fontsize='small')
        pl.ylabel('v (m s-1)')

        pl.subplot(223)
        pl.plot(f1.index, f1['T_90'], '+', color='r', markersize=2, label='FINO1 obs')
        pl.plot(e1.time, e1.prof['Ta'][:,4], 'k-',  label='LES: mean diff={0:.2f} K, rmse={1:.2f} K'.format(e1.T_90.mean_diff, e1.T_90.rmse))
        pl.plot(h1.time, h1['T'][:,iloc,3],  'b--', label='HAM: mean diff={0:.2f} K, rmse={1:.2f} K'.format(T_90_HAM.mean_diff, T_90_HAM.rmse))
        format_date_day()
        pl.legend(fontsize='small', loc=2)
        pl.ylabel('T (K)')

        pl.subplot(224)
        pl.plot(f1.index, f1['q_90']*1e3, '+', color='r', markersize=2, label='FINO1 obs')
        pl.plot(e1.time, e1.prof['qt'][:,4]*1e3, 'k-',  label='LES: mean diff={0:.2f} g/kg, rmse={1:.2f} g/kg'.format(e1.q_90.mean_diff, e1.q_90.rmse))
        pl.plot(h1.time, h1['q'][:,iloc,3]*1e3,  'b--', label='HAM: mean diff={0:.2f} g/kg, rmse={1:.2f} g/kg'.format(q_90_HAM.mean_diff, q_90_HAM.rmse))
        format_date_day()
        pl.legend(fontsize='small')
        pl.ylabel('q (g kg-1)')

        pl.tight_layout()
        pl.savefig('tser_full_LES10min.pdf')

    if False:

        pl.figure(figsize=(8,5))

        pl.subplot(111)
        pl.title('z = 90 m: ', loc='left')
        pl.plot(f1.index, f1['u_90'], 'o', color='r', markersize=4, label='FINO1 obs')
        pl.plot(e1.time, e1.prof['u'][:,4], 'k-', label='LES (10m): bias={0:.2f} m/s, rmse={1:.2f} m/s'.format(e1.u_90.mean_diff, e1.u_90.rmse))
        pl.plot(e2.time, e2.prof['u'][:,4], 'C3', dashes=[2,2], label='LES (30m): bias={0:.2f} m/s, rmse={1:.2f} m/s'.format(e2.u_90.mean_diff, e2.u_90.rmse))
        pl.plot(e3.time, e3.prof['u'][:,4], 'C2', dashes=[4,2], label='LES (60m): bias={0:.2f} m/s, rmse={1:.2f} m/s'.format(e3.u_90.mean_diff, e3.u_90.rmse))
        pl.plot(e4.time, e4.prof['u'][:,4], 'C0', dashes=[6,2], label='LES (120m): bias={0:.2f} m/s, rmse={1:.2f} m/s'.format(e4.u_90.mean_diff, e4.u_90.rmse))
        #pl.plot(h1.time, h1['u'][:,iloc,3], 'b:',  label='HAM: mean diff={0:.2f} m/s, rmse={1:.2f} m/s'.format(u_90_HAM.mean_diff, u_90_HAM.rmse))
        #format_date_day()
        format_date_hour()
        pl.xlim(datetime.datetime(2017,6,15,11), datetime.datetime(2017,6,15,19))
        pl.ylim(-5,15)
        pl.legend(fontsize='small')
        pl.ylabel('u (m s-1)')

        pl.tight_layout()
        pl.savefig('tser_u_LES_sensitivity_short.pdf')



    if True:
        dtu_dyn = np.abs(h1['dtu_dyn'][:,iloc,3])
        dtu_phy = np.abs(h1['dtu_phy'][:,iloc,3])
        dtu = dtu_dyn + dtu_phy

        pl.figure(figsize=(8,5))

        ax=pl.subplot(111)
        pl.title('z = 90 m: ', loc='left')
        pl.plot(f1.index, f1['u_90'], 'o', color='r', markersize=2, label='FINO1 obs')
        pl.plot(h1.time, h1['u'][:,iloc,3], 'k-',  label='HAM: mean diff={0:.2f} m/s, rmse={1:.2f} m/s'.format(u_90_HAM.mean_diff, u_90_HAM.rmse))
        format_date_day()
        #format_date_hour()
        #pl.xlim(datetime.datetime(2017,6,15,11), datetime.datetime(2017,6,15,19))
        #pl.ylim(-5,15)
        pl.legend(fontsize='small')
        pl.ylabel('u (m s-1)')

        ax=ax.twinx()
        pl.plot(h1.time, moving_av(dtu_dyn/dtu, 12), 'C2')
        pl.ylim(0,1)
        pl.ylabel('dtu_dyn / dtu (-)')
        ax.spines['right'].set_visible(True)

        pl.tight_layout()
