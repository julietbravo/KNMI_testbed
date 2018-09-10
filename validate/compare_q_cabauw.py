import matplotlib.pyplot as pl
import xarray as xr
import pandas as pd
import numpy as np
import datetime

pl.close('all')

def calc_esat(T):
    return 0.611e3 * np.exp(17.2694 * (T - 273.16) / (T - 35.86))

def calc_qsat(T,p):
    return 0.622 * calc_esat(T) / p

def calc_rmse(v1, v2):
    return np.sqrt(((v1 - v2)**2).mean())



def read_DOWA(i,j):
    # Path to the DOWA timeseries on Berts KNMI-pc
    path = '/net/mosuracmo/nobackup_2/users/ulft/DOWA/DOWA_40h12tg2_fERA5/timeseries/'

    # Read data with xarray
    #ds = xr.open_dataset('{0:}/DOWA_40h12tg2_fERA5_NETHERLANDS.NL_ix{1:03d}_iy{2:03d}_2008010100-2018010100_v1.0.nc'.format(path, i, j))
    ds = xr.open_dataset('test.nc')

    # Calculate specific humidity
    ds['q'] = ds['hur'] * calc_qsat(ds['ta'], ds['p'])

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

    # Fix units...
    ds['RH'] /= 100.        # %    -> fraction
    df['Q']  /= 1000.       # g/kg -> kg/kg

    return ds



if __name__ == '__main__':
    # Read DOWA and Cabauw NetCDF data:
    if 'dw_ds' not in locals():
        dw_ds = read_DOWA(98, 74)
        cb_ds = read_Cabauw()

    # Heights to compare
    heights = np.array([10, 200])

    # Variable names in the DOWA and Cabauw files
    dw_vars = {'q': 'q', 'ta': 'ta', 'rh': 'hur', 'u': 'wspeed', 'udir': 'wdir'}
    cb_vars = {'q': 'Q', 'ta': 'TA', 'rh': 'RH',  'u': 'F',      'udir': 'D'}

    # Gather time series from DOWA and Cabauw data, and put in Pandas dataframe
    data_dw = {}
    data_cb = {}
    for z in heights:
        # Find height index:
        k_dw = int(np.abs(dw_ds.height - z).argmin())
        k_cb = int(np.abs(cb_ds.z      - z).argmin())

        for var in dw_vars.keys():
            data_dw['dw_{0:}_{1:03d}'.format(var, z)] = dw_ds[dw_vars[var]][:,k_dw,0,0]
            data_cb['cb_{0:}_{1:03d}'.format(var, z)] = cb_ds[cb_vars[var]][:,k_cb    ]

    # Create Pandas dataframes
    dw_df = pd.DataFrame(data_dw, index=dw_ds.time)
    cb_df = pd.DataFrame(data_cb, index=cb_ds.time)

    # Mask invalid (not-yet filled) part of DOWA data
    dw_df = dw_df.mask(dw_df['dw_ta_010'] > 1e12)

    # Cabauw time is a bit inaccurate; round it...
    cb_df.index = cb_df.index.round('min')

    # Merge them; this drops times which don't exists in both time series (Harmonie = hourly, Cabauw = 10 min)
    df = pd.concat([dw_df, cb_df], axis=1, join='inner')
    df.dropna(inplace=True)

    # Season selections
    #df = df.loc[(df.index.month >= 4)  & (df.index.month <= 9)]
    #df = df.loc[(df.index.month >= 11) | (df.index.month <= 2)]

    # Calculate statistics
    rmse = {}
    bias = {}
    for var in dw_vars.keys():
        rmse[var] = np.zeros((24, heights.size))
        bias[var] = np.zeros((24, heights.size))

    for t in range(24):
        df_hour = df.loc[df.index.hour == t]

        for var in dw_vars.keys():
            for k in range(heights.size):
                z = heights[k]

                v_dw = 'dw_{0:}_{1:03d}'.format(var, z)
                v_cb = 'cb_{0:}_{1:03d}'.format(var, z)

                bias[var][t,k] = np.mean  (df_hour[v_dw] - df_hour[v_cb])
                rmse[var][t,k] = calc_rmse(df_hour[v_dw],  df_hour[v_cb])



    pl.figure(figsize=(12,7))
    pl.subplot(241)
    pl.plot(np.arange(24), bias['q'][:,0]*1000., label='10 m')
    pl.plot(np.arange(24), bias['q'][:,1]*1000., label='200 m')
    pl.legend()
    pl.xlabel('time UTC (h)')
    pl.ylabel('<q model-obs> (g kg-1)')

    pl.subplot(242)
    pl.plot(np.arange(24), bias['ta'][:,0])
    pl.plot(np.arange(24), bias['ta'][:,1])
    pl.xlabel('time UTC (h)')
    pl.ylabel('<T model-obs> (K)')

    pl.subplot(243)
    pl.plot(np.arange(24), bias['rh'][:,0]*100)
    pl.plot(np.arange(24), bias['rh'][:,1]*100)
    pl.xlabel('time UTC (h)')
    pl.ylabel('<RH model-obs> (%)')

    pl.subplot(244)
    pl.plot(np.arange(24), bias['u'][:,0])
    pl.plot(np.arange(24), bias['u'][:,1])
    pl.xlabel('time UTC (h)')
    pl.ylabel('<u model-obs> (m/s)')

    pl.subplot(245)
    pl.plot(np.arange(24), rmse['q'][:,0]*1000., label='10 m')
    pl.plot(np.arange(24), rmse['q'][:,1]*1000., label='200 m')
    pl.legend()
    pl.xlabel('time UTC (h)')
    pl.ylabel('rmse q (g kg-1)')

    pl.subplot(246)
    pl.plot(np.arange(24), rmse['ta'][:,0])
    pl.plot(np.arange(24), rmse['ta'][:,1])
    pl.xlabel('time UTC (h)')
    pl.ylabel('rmse T (K)')

    pl.subplot(247)
    pl.plot(np.arange(24), rmse['rh'][:,0]*100)
    pl.plot(np.arange(24), rmse['rh'][:,1]*100)
    pl.xlabel('time UTC (h)')
    pl.ylabel('rmse RH (%)')

    pl.subplot(248)
    pl.plot(np.arange(24), rmse['u'][:,0])
    pl.plot(np.arange(24), rmse['u'][:,1])
    pl.xlabel('time UTC (h)')
    pl.ylabel('rmse u (m/s)')

    pl.tight_layout()

