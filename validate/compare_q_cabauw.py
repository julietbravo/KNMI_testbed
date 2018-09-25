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


def read_DOWA_2m(i,j):
    path = '/nobackup/users/stratum/DOWA/validation_2m/'

    ds_rh = xr.open_mfdataset('{}/hurs*'.format(path))
    ds_ta = xr.open_mfdataset('{}/tas*' .format(path))
    ds_ps = xr.open_mfdataset('{}/ps*'  .format(path))

    # Select Cabauw gridpoint
    x = ds_rh.x[i]
    y = ds_rh.y[j]

    ds_rh = ds_rh.sel(x=x, y=y)
    ds_ta = ds_ta.sel(x=x, y=y)
    ds_ps = ds_ps.sel(x=x, y=y)

    ds_rh = ds_rh.sortby('time')
    ds_ta = ds_ta.sortby('time')
    ds_ps = ds_ps.sortby('time')

    ds_q  = ds_rh.hurs * calc_qsat(ds_ta.tas, ds_ps.ps)

    data = {'hm_rh_002': ds_rh.hurs, 'hm_ta_002': ds_ta.tas, 'hm_q_002': ds_q.values}
    df = pd.DataFrame(data, index=ds_rh.time)

    return df


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
    ds['Q']  /= 1000.       # g/kg -> kg/kg

    return ds


def read_Cabauw_sfc():
    # Path to Cabauw data
    path = '/nobackup/users/stratum/Cabauw/'

    # Generate list of files to read
    files = []
    for year in range(2008, 2018):
        for month in range(1,13):
            files.append('{0:}/cesar_surface_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path, year, month))

    # Exclude `valid_dates`; xarray doesn't like multiple time dimensions....
    exclude = ['valid_dates']
    ds = xr.open_mfdataset(files, drop_variables=exclude)

    # Fix units...
    ds['RH002'] /= 100.        # %    -> fraction
    ds['Q002']  /= 1000.       # g/kg -> kg/kg

    return ds



if __name__ == '__main__':
    # Read DOWA and Cabauw NetCDF data:
    if 'hm_ds' not in locals():
        hm_ds    = read_DOWA(98, 74)
    if 'hm_df_2m' not in locals():
        hm_df_2m = read_DOWA_2m(97, 73)  # Indexing differs from hm_ds; Fortran vs Python...
    if 'cb_ds' not in locals():
        cb_ds     = read_Cabauw()
    if 'cb_ds_sfc' not in locals():
        cb_ds_sfc = read_Cabauw_sfc()

    # Heights to compare
    heights = np.array([10, 200])

    # Variable names in the DOWA and Cabauw files
    hm_vars = {'q': 'q', 'ta': 'ta', 'rh': 'hur', 'u': 'wspeed', 'udir': 'wdir'}
    cb_vars = {'q': 'Q', 'ta': 'TA', 'rh': 'RH',  'u': 'F',      'udir': 'D'}

    # Gather time series from DOWA and Cabauw data, and put in Pandas dataframe
    data_hm = {}
    data_cb = {}
    for z in heights:
        # Find height index:
        k_hm = int(np.abs(hm_ds.height - z).argmin())
        k_cb = int(np.abs(cb_ds.z      - z).argmin())

        for var in hm_vars.keys():
            data_hm['hm_{0:}_{1:03d}'.format(var, z)] = hm_ds[hm_vars[var]][:,k_hm,0,0]
            data_cb['cb_{0:}_{1:03d}'.format(var, z)] = cb_ds[cb_vars[var]][:,k_cb    ]

    # Create Pandas dataframes
    hm_df = pd.DataFrame(data_hm, index=hm_ds.time)
    cb_df = pd.DataFrame(data_cb, index=cb_ds.time)

    data_cb_sfc = {'cb_ta_002': cb_ds_sfc.TA002, 'cb_q_002': cb_ds_sfc.Q002, 'cb_rh_002': cb_ds_sfc.RH002}
    cb_df_sfc = pd.DataFrame(data_cb_sfc, index=cb_ds_sfc.time)
    cb_df_sfc.index = cb_df_sfc.index.round('min')

    # Mask invalid (not-yet filled) part of DOWA data
    hm_df = hm_df.mask(hm_df['hm_ta_010'] > 1e12)

    # Cabauw time is a bit inaccurate; round it...
    cb_df.index = cb_df.index.round('min')

    # Merge them; this drops times which don't exists in both time series (Harmonie = hourly, Cabauw = 10 min)
    df = pd.concat([hm_df, cb_df, hm_df_2m, cb_df_sfc], axis=1, join='inner')
    df.dropna(inplace=True)

    # Season selections
    #df = df.loc[(df.index.month >= 4)  & (df.index.month <= 9)]
    #df = df.loc[(df.index.month >= 11) | (df.index.month <= 2)]

    # Calculate statistics : atmosphere
    rmse = {}
    bias = {}
    for var in hm_vars.keys():
        rmse[var] = np.zeros((24, heights.size))
        bias[var] = np.zeros((24, heights.size))

    for t in range(24):
        df_hour = df.loc[df.index.hour == t]

        for var in hm_vars.keys():
            for k in range(heights.size):
                z = heights[k]

                v_hm = 'hm_{0:}_{1:03d}'.format(var, z)
                v_cb = 'cb_{0:}_{1:03d}'.format(var, z)

                bias[var][t,k] = np.mean  (df_hour[v_hm] - df_hour[v_cb])
                rmse[var][t,k] = calc_rmse(df_hour[v_hm],  df_hour[v_cb])

    # Calculate statistics : surface / 2m
    vars = ['ta002', 'q002', 'rh002']
    for var in vars:
        rmse[var] = np.zeros(24)
        bias[var] = np.zeros(24)

    for t in range(24):
        df_hour = df.loc[df.index.hour == t]

        bias['ta002'][t] = np.mean(df_hour['hm_ta_002'] - df_hour['cb_ta_002'])
        bias['rh002'][t] = np.mean(df_hour['hm_rh_002'] - df_hour['cb_rh_002'])
        bias['q002' ][t] = np.mean(df_hour['hm_q_002' ] - df_hour['cb_q_002' ])

        rmse['ta002'][t] = calc_rmse(df_hour['hm_ta_002'], df_hour['cb_ta_002'])
        rmse['rh002'][t] = calc_rmse(df_hour['hm_rh_002'], df_hour['cb_rh_002'])
        rmse['q002' ][t] = calc_rmse(df_hour['hm_q_002' ], df_hour['cb_q_002' ])



    pl.figure(figsize=(9,7))
    pl.subplot(221)
    pl.plot(np.arange(24), bias['q'][:,0]*1000., '--', label='10 m')
    pl.plot(np.arange(24), bias['q'][:,1]*1000., '-.', label='200 m')
    pl.plot(np.arange(24), bias['q002']*1000.,         label='2 m')
    pl.legend()
    pl.xlabel('time UTC (h)')
    pl.ylabel('<q model-obs> (g kg-1)')

    pl.subplot(222)
    pl.plot(np.arange(24), bias['ta'][:,0], '--')
    pl.plot(np.arange(24), bias['ta'][:,1], '-.')
    pl.plot(np.arange(24), bias['ta002'])
    pl.xlabel('time UTC (h)')
    pl.ylabel('<T model-obs> (K)')

    pl.subplot(223)
    pl.plot(np.arange(24), bias['rh'][:,0]*100, '--')
    pl.plot(np.arange(24), bias['rh'][:,1]*100, '-.')
    pl.plot(np.arange(24), bias['rh002']*100)
    pl.xlabel('time UTC (h)')
    pl.ylabel('<RH model-obs> (%)')

    pl.subplot(224)
    pl.plot(np.arange(24), bias['u'][:,0], '--')
    pl.plot(np.arange(24), bias['u'][:,1], '-.')
    pl.xlabel('time UTC (h)')
    pl.ylabel('<u model-obs> (m/s)')

    pl.tight_layout()




    pl.figure(figsize=(9,7))
    pl.subplot(221)
    pl.plot(np.arange(24), rmse['q'][:,0]*1000., '--')
    pl.plot(np.arange(24), rmse['q'][:,1]*1000., '-.')
    pl.plot(np.arange(24), rmse['q002']*1000.)
    pl.xlabel('time UTC (h)')
    pl.ylabel('rmse q (g kg-1)')

    pl.subplot(222)
    pl.plot(np.arange(24), rmse['ta'][:,0], '--')
    pl.plot(np.arange(24), rmse['ta'][:,1], '-.')
    pl.plot(np.arange(24), rmse['ta002'])
    pl.xlabel('time UTC (h)')
    pl.ylabel('rmse T (K)')

    pl.subplot(223)
    pl.plot(np.arange(24), rmse['rh'][:,0]*100, '--')
    pl.plot(np.arange(24), rmse['rh'][:,1]*100, '-.')
    pl.plot(np.arange(24), rmse['rh002']*100)
    pl.xlabel('time UTC (h)')
    pl.ylabel('rmse RH (%)')

    pl.subplot(224)
    pl.plot(np.arange(24), rmse['u'][:,0], '--')
    pl.plot(np.arange(24), rmse['u'][:,1], '-.')
    pl.xlabel('time UTC (h)')
    pl.ylabel('rmse u (m/s)')

    pl.tight_layout()

