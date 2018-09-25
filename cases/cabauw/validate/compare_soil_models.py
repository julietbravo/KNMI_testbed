import matplotlib.pyplot as pl
import pandas as pd
import xarray as xr
import numpy as np
import calendar

pl.close('all')

def interpz(z1, z2, zg, v1, v2):
    f2 = (zg - z1) / (z2 - z1)
    f1 = 1 - f2

    print(f1,f2)
    return f1*v1 + f2*v2

year = 2017

# Read ERA5 soil
pwd = '/nobackup/users/stratum/ERA5/soil/'
files = ['{0:}/soil_{1:04d}{2:02d}.nc'.format(pwd, year, m) for m in range(1,13)]
e5 = xr.open_mfdataset(files)
e5 = e5.sel(longitude=4.9, latitude=51.97, method='nearest')
e5 = e5.to_dataframe()
cols = ['e5_{}'.format(name) for name in e5.columns]
e5.columns = cols

# Read Harmonie soil
pwd = '/nobackup/users/stratum/DOWA/soil/'
files = ['{0:}/soil_{1:04d}{2:02d}.nc'.format(pwd, year, m) for m in range(1,12)]
hm = xr.open_mfdataset(files)
hm = hm.to_dataframe()
cols = ['hm_{}'.format(name) for name in hm.columns]
hm.columns = cols

# Read Cabauw observations
pwd = '/nobackup/users/stratum/Cabauw'
files1 = ['{0:}/cesar_soil_water_lb1_t10_v1.1_{1:04d}{2:02d}.nc'.format(pwd, year, m) for m in range(1,13)]
files2 = ['{0:}/cesar_soil_heat_lb1_t10_v1.0_{1:04d}{2:02d}.nc' .format(pwd, year, m) for m in range(1,13)]
files = files1 + files2
exclude = ['valid_dates', 'time_bnds', 'iso_dataset', 'product', 'station_details']
cb = xr.open_mfdataset(files, concat_dim='time', drop_variables=exclude)
cb = cb.to_dataframe()
cb.index = cb.index.round('min')
cols = ['cb_{}'.format(name) for name in cb.columns]
cb.columns = cols


# Merge!
df = pd.concat([e5, cb, hm], axis=1, join='inner')

# Interpolate Cabauw obs to ERA5 levels.....
df['cb_TS035i'] = interpz(2.0, 4.0, 3.5, df['cb_TS02'].values, df['cb_TS04'].values)

if False:
    pl.figure()
    pl.subplot(221)
    pl.plot(df.index, df['e5_stl1'], label='-3.5 cm')
    pl.plot(df.index, df['e5_stl2'], label='-17.5 cm')
    pl.plot(df.index, df['e5_stl3'], label='-64 cm')
    pl.plot(df.index, df['e5_stl4'], label='-195 cm')
    pl.legend()

    pl.subplot(222)
    pl.plot(df.index, df['cb_TS02'])
    pl.plot(df.index, df['cb_TS04'])
    pl.plot(df.index, df['cb_TS06'])
    pl.plot(df.index, df['cb_TS08'])
    pl.plot(df.index, df['cb_TS12'])
    pl.plot(df.index, df['cb_TS20'])
    pl.plot(df.index, df['cb_TS30'])
    pl.plot(df.index, df['cb_TS50'])
    pl.legend()

    pl.subplot(223)
    pl.plot(df.index, df['e5_swvl1'], label='-3.5 cm')
    pl.plot(df.index, df['e5_swvl2'], label='-17.5 cm')
    pl.plot(df.index, df['e5_swvl3'], label='-64 cm')
    pl.plot(df.index, df['e5_swvl4'], label='-195 cm')
    pl.legend()

    pl.subplot(224)
    #pl.plot(df.index, df['df_TH03'])
    pl.plot(df.index, df['cb_TH05'])
    #pl.plot(df.index, df['df_TH08'])
    pl.plot(df.index, df['cb_TH19'])
    #pl.plot(df.index, df['df_TH20'])
    pl.plot(df.index, df['cb_TH33'])
    pl.plot(df.index, df['cb_TH40'])
    pl.plot(df.index, df['cb_TH56'])
    pl.legend()

pl.figure()
pl.plot(cb.index, cb['cb_GWLEBN'])
pl.plot(cb.index, cb['cb_GWLEBS'])

pl.figure()

pl.subplot(221)
cc = pl.cm.hsv(np.linspace(0,1,12))
for m in range(12):
    df_m = df.loc[df.index.month == m+1]
    pl.scatter(df_m['cb_TH05'], df_m['e5_swvl1'], c=cc[m], s=5, label=calendar.month_name[m+1])
pl.legend(ncol=3)
pl.plot([0.15,0.60], [0.15,0.60], 'k--')
pl.xlabel('phi -5 cm Cabauw (-)')
pl.ylabel('phi -3.5 cm ERA5 (-)')

pl.subplot(222)
for m in range(12):
    df_m = df.loc[df.index.month == m+1]
    pl.scatter(df_m['cb_TH19'], df_m['e5_swvl2'], c=cc[m], s=5, label=str(m+1))
pl.plot([0.15,0.60], [0.15,0.60], 'k--')
pl.xlabel('phi -19 cm Cabauw (-)')
pl.ylabel('phi -17.5 cm ERA5 (-)')

pl.subplot(223)
cc = pl.cm.hsv(np.linspace(0,1,12))
for m in range(12):
    df_m = df.loc[df.index.month == m+1]
    pl.scatter(df_m['cb_TH05'], df_m['hm_wsa_L01'], c=cc[m], s=5, label=calendar.month_name[m+1])
pl.plot([0.15,0.60], [0.15,0.60], 'k--')
pl.xlabel('phi -5 cm Cabauw (-)')
pl.ylabel('phi -0.5 cm Harmonie (-)')

pl.subplot(224)
for m in range(12):
    df_m = df.loc[df.index.month == m+1]
    pl.scatter(df_m['cb_TH56'], df_m['hm_wsa_L02'], c=cc[m], s=5, label=str(m+1))
pl.plot([0.15,0.60], [0.15,0.60], 'k--')
pl.xlabel('phi -56 cm Cabauw (-)')
pl.ylabel('phi -50 cm Harmonie (-)')
