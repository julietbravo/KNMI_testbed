import matplotlib.pyplot as pl
import xarray as xr
import pandas as pd
import numpy as np
import datetime

from read_DDH_netcdf import *

pl.close('all')

#path_hm = '/Users/bart/meteo/data/Harmonie_DDH/'
#path_cb = '/Users/bart/meteo/observations/Cabauw'

path_hm = '/nobackup/users/stratum/DOWA/LES_forcing/'
path_cb = '/nobackup/users/stratum/Cabauw/'

class Statistics:
    def __init__(self, obs, model, scale=1):
        self.mean_diff = (model - obs).mean() * scale
        self.max_diff  = (model - obs).max() * scale
        self.rmse      = np.sqrt(((model - obs)**2).mean()) * scale

def interpz(data, z, zg):
    """
    Interpolate `data` to goal_height (`zg`) for
    data where the height levels change in time
    """

    # Needed for xarray
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
    return f0*data[:,k0][:,0] + f1*data[:,k0+1][:,0]


# -----------------
# Period
# -----------------
start = datetime.datetime(year=2017, month=1, day=1, hour=0)
end   = datetime.datetime(year=2018, month=1, day=1, hour=0)

# -----------------
# Read Harmonie statistics
# -----------------
iloc = 7    # 7 = single column Cabauw
variables = ['time', 'z', 'q']
if 'hm' not in locals():
    hm  = read_DDH_netcdf(start, end, path_hm, variables)

# Interpolate to observation heights Cabauw, and store in DataFrame
data = {'hm_q20' : interpz(hm.q[:,iloc,:], hm.z[:,iloc,:], 20. )*1e3,
        'hm_q200': interpz(hm.q[:,iloc,:], hm.z[:,iloc,:], 200.)*1e3}

df_hm = pd.DataFrame(data, index=hm.time)
df_hm.index = df_hm.index.round('min')
df_hm = df_hm.loc[start:end]

#    # -----------------
#    # Read Cabauw observations
#    # -----------------
#    files = []
#    for m in range(start.month, end.month+1):
#        files.append('{0}/cesar_tower_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path_cb, start.year, m))
#    #if 'cb' not in locals():
#    exclude = ['valid_dates']
#    cb = xr.open_mfdataset(files, concat_dim='time', drop_variables=exclude)
#
#    # Put in DataFrame
#    data = {'cb_q10' : cb['Q'][:,5],
#            'cb_q20' : cb['Q'][:,4],
#            'cb_q40' : cb['Q'][:,3],
#            'cb_q80' : cb['Q'][:,2],
#            'cb_q140': cb['Q'][:,1],
#            'cb_q200': cb['Q'][:,0]}
#
#    df_cb = pd.DataFrame(data, index=cb.time)
#    df_cb.index = df_cb.index.round('min')
#    df_cb = df_cb.loc[start:end]
#
#    # -----------------
#    # Merge them!
#    # -----------------
#    df = pd.concat([df_hm, df_cb], axis=1, join='inner')
#
#    for z in [10, 20, 40, 80, 140, 200]:
#        stat = Statistics(df['cb_q{}'.format(z)], df['hm_q{}'.format(z)])
#        print('z={0:5.0f} m: mean_diff={1:6.2f}, rmse={2:6.2f}'.format(z, stat.mean_diff, stat.rmse))







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
