import matplotlib.pyplot as pl
import xarray as xr
import pandas as pd
import numpy as np
import datetime

pl.close('all')

def esat(T):
    return 0.611e3 * np.exp(17.2694 * (T - 273.16) / (T - 35.86))

def qsat(T,p):
    return 0.622 * esat(T) / p

def rmse(v1, v2):
    return np.sqrt(((v1 - v2)**2).mean())

def read_DOWA(i,j):
    # Path to the DOWA timeseries on Berts KNMI-pc
    path = '/net/mosuracmo/nobackup_2/users/ulft/DOWA/DOWA_40h12tg2_fERA5/timeseries/'

    # Read data with xarray
    ds = xr.open_dataset('{0:}/DOWA_40h12tg2_fERA5_NETHERLANDS.NL_ix{1:03d}_iy{2:03d}_2008010100-2018010100_v1.0.nc'.format(path, i, j))

    # Set invalid (runs is not yet completed...) values to nan


    # Calculate specific humidity
    ds['q'] = ds['hur'] * qsat(ds['ta'], ds['p'])

    return ds

def read_Cabauw():
    # Path to Cabauw data
    path = '/nobackup/users/stratum/Cabauw/'

    # Generate list of files to read
    files = []
    for year in range(2008, 2018):
        for month in range(1,13):
            files.append('{0:}/cesar_tower_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path, year, month))

    # Exclude `valid_dates`; xarray doesn't like multiple time dimensions....
    exclude = ['valid_dates']
    ds = xr.open_mfdataset(files, drop_variables=exclude)

    return ds



if __name__ == '__main__':
    # Read DOWA and Cabauw NetCDF data:
    if 'dowa_ds' not in locals():
        dowa_ds   = read_DOWA(98, 74)
        cabauw_ds = read_Cabauw()

    # Heights to compare
    heights = np.array([10, 80, 200])

    # Gather time series from DOWA and Cabauw data, and put in Pandas dataframe
    data_dw = {}
    data_cb = {}
    for z in heights:
        # Find height index:
        k_dw = int(np.abs(dowa_ds.height - z).argmin())
        k_cb = int(np.abs(cabauw_ds.z    - z).argmin())

        # Store moisture & temperature time series in dictionary
        data_dw['dw_q_{0:03d}'.format(z)] = dowa_ds  ['q'][:,k_dw,0,0]
        data_cb['cb_q_{0:03d}'.format(z)] = cabauw_ds['Q'][:,k_cb] * 1e-3

        data_dw['dw_ta_{0:03d}'.format(z)] = dowa_ds  ['ta'][:,k_dw,0,0]
        data_cb['cb_ta_{0:03d}'.format(z)] = cabauw_ds['TA'][:,k_cb]

        data_dw['dw_rh_{0:03d}'.format(z)] = dowa_ds  ['hur'][:,k_dw,0,0]
        data_cb['cb_rh_{0:03d}'.format(z)] = cabauw_ds['RH'][:,k_cb]

    # Create Pandas dataframes
    dowa_df   = pd.DataFrame(data_dw, index=dowa_ds  .time)
    cabauw_df = pd.DataFrame(data_cb, index=cabauw_ds.time)

    # Mask invalid (not-yet filled) part of DOWA data
    dowa_df = dowa_df.mask(dowa_df['dw_ta_010'] > 1e12)

    # Cabauw time is a bit inaccurate; round it...
    cabauw_df.index = cabauw_df.index.round('min')

    # Merge them; this drops times which don't exists in both time series (Harmonie = hourly, Cabauw = 10 min)
    df = pd.concat([dowa_df, cabauw_df], axis=1, join='inner')
    df.dropna(inplace=True)

    # Season selections
    #df = df.loc[(df.index.month >= 5)  & (df.index.month <= 8)]
    #df = df.loc[(df.index.month >= 11) | (df.index.month <= 2)]

    rmse_q = np.zeros((24, heights.size))
    rmse_T = np.zeros((24, heights.size))

    bias_q = np.zeros((24, heights.size))
    bias_T = np.zeros((24, heights.size))

    for t in range(24):
        for k in range(heights.size):
            df_hour = df.loc[df.index.hour == t]

            q1 = 'dw_q_{0:03d}'.format(heights[k])
            q2 = 'cb_q_{0:03d}'.format(heights[k])

            T1 = 'dw_ta_{0:03d}'.format(heights[k])
            T2 = 'cb_ta_{0:03d}'.format(heights[k])

            bias_q[t,k] = np.mean(df_hour[q1] - df_hour[q2])
            bias_T[t,k] = np.mean(df_hour[T1] - df_hour[T2])

            rmse_q[t,k] = rmse(df_hour[q1], df_hour[q2])
            rmse_T[t,k] = rmse(df_hour[T1], df_hour[T2])


    pl.figure(figsize=(8,6))
    pl.subplot(221)
    pl.plot(np.arange(24)+1, bias_q[:,0]*1000., label='10 m')
    pl.plot(np.arange(24)+1, bias_q[:,1]*1000., label='80 m')
    pl.plot(np.arange(24)+1, bias_q[:,2]*1000., label='200 m')
    pl.legend()
    pl.xlabel('time UTC (h)')
    pl.ylabel('<q model-obs> (g kg-1)')

    pl.subplot(222)
    pl.plot(np.arange(24)+1, bias_T[:,0])
    pl.plot(np.arange(24)+1, bias_T[:,1])
    pl.plot(np.arange(24)+1, bias_T[:,2])
    pl.xlabel('time UTC (h)')
    pl.ylabel('<T model-obs> (K)')

    pl.subplot(223)
    pl.plot(np.arange(24)+1, rmse_q[:,0]*1000.)
    pl.plot(np.arange(24)+1, rmse_q[:,1]*1000.)
    pl.plot(np.arange(24)+1, rmse_q[:,2]*1000.)
    pl.xlabel('time UTC (h)')
    pl.ylabel('RMSE q (g kg-1)')

    pl.subplot(224)
    pl.plot(np.arange(24)+1, rmse_T[:,0])
    pl.plot(np.arange(24)+1, rmse_T[:,1])
    pl.plot(np.arange(24)+1, rmse_T[:,2])
    pl.ylabel('RMSE T (K)')

    pl.tight_layout()






# -----------------
# --- DOWA DATA ---
# -----------------

# ---------------------------
# --- Cabauw observations ---
# ---------------------------
# Path to Cabauw observations




#from read_DDH_netcdf import *
#
#pl.close('all')
#
##path_hm = '/Users/bart/meteo/data/Harmonie_DDH/'
##path_cb = '/Users/bart/meteo/observations/Cabauw'
#
#path_hm = '/nobackup/users/stratum/DOWA/LES_forcing/'
#path_cb = '/nobackup/users/stratum/Cabauw/'
#
#class Statistics:
#    def __init__(self, obs, model, scale=1):
#        self.mean_diff = (model - obs).mean() * scale
#        self.max_diff  = (model - obs).max() * scale
#        self.rmse      = np.sqrt(((model - obs)**2).mean()) * scale
#
#def interpz(data, z, zg):
#    """
#    Interpolate `data` to goal_height (`zg`) for
#    data where the height levels change in time
#    """
#
#    # Needed for xarray
#    z = z.values
#    data = data.values
#
#    # Find height index of nearest level below `zg`
#    k0 = np.abs(z - zg).argmin(axis=1)
#    z0 = z[:,k0][:,0]
#    k0[z0 > zg] -= 1
#
#    # Heights above and below interpolation height
#    z1 = z[:,k0+1][:,0]
#    z0 = z[:,k0]  [:,0]
#    dz = z1 - z0
#
#    # Interpolation factors
#    f1 = (zg - z0) / dz
#    f0 = 1-f1
#
#    # Interpolate!
#    return f0*data[:,k0][:,0] + f1*data[:,k0+1][:,0]
#
#
## -----------------
## Period
## -----------------
#start = datetime.datetime(year=2017, month=1, day=1, hour=0)
#end   = datetime.datetime(year=2018, month=1, day=1, hour=0)
#
## -----------------
## Read Harmonie statistics
## -----------------
#iloc = 7    # 7 = single column Cabauw
#variables = ['time', 'z', 'q']
#if 'hm' not in locals():
#    hm  = read_DDH_netcdf(start, end, path_hm, variables)
#
## Interpolate to observation heights Cabauw, and store in DataFrame
#data = {'hm_q20' : interpz(hm.q[:,iloc,:], hm.z[:,iloc,:], 20. )*1e3,
#        'hm_q200': interpz(hm.q[:,iloc,:], hm.z[:,iloc,:], 200.)*1e3}
#
#df_hm = pd.DataFrame(data, index=hm.time)
#df_hm.index = df_hm.index.round('min')
#df_hm = df_hm.loc[start:end]
#
##    # -----------------
##    # Read Cabauw observations
##    # -----------------
##    files = []
##    for m in range(start.month, end.month+1):
##        files.append('{0}/cesar_tower_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path_cb, start.year, m))
##    #if 'cb' not in locals():
##    exclude = ['valid_dates']
##    cb = xr.open_mfdataset(files, concat_dim='time', drop_variables=exclude)
##
##    # Put in DataFrame
##    data = {'cb_q10' : cb['Q'][:,5],
##            'cb_q20' : cb['Q'][:,4],
##            'cb_q40' : cb['Q'][:,3],
##            'cb_q80' : cb['Q'][:,2],
##            'cb_q140': cb['Q'][:,1],
##            'cb_q200': cb['Q'][:,0]}
##
##    df_cb = pd.DataFrame(data, index=cb.time)
##    df_cb.index = df_cb.index.round('min')
##    df_cb = df_cb.loc[start:end]
##
##    # -----------------
##    # Merge them!
##    # -----------------
##    df = pd.concat([df_hm, df_cb], axis=1, join='inner')
##
##    for z in [10, 20, 40, 80, 140, 200]:
##        stat = Statistics(df['cb_q{}'.format(z)], df['hm_q{}'.format(z)])
##        print('z={0:5.0f} m: mean_diff={1:6.2f}, rmse={2:6.2f}'.format(z, stat.mean_diff, stat.rmse))
#
#
#
#
#
#
#
##
#
#
#
##pl.figure()
##pl.subplot(211)
##pl.plot(df.index, df['q_hm'], label='HAM')
##pl.plot(df.index, df['q_cb'], label='Cabauw')
##pl.legend()
##
##
##pl.subplot(212)
##pl.plot(df.index, df['q_hm']-df['q_cb'])
