import matplotlib.pyplot as pl
import xarray as xr
import pandas as pd
import numpy as np
import datetime

from read_DDH_netcdf import *

pl.close('all')

def components_to_direction(u, v):
    return np.rad2deg( np.arctan2(u, v) ) + 180

# -- Path to LES forcings --
#path = '/Users/bart/meteo/data/Harmonie_DDH/'
path = '/nobackup/users/stratum/DOWA/LES_forcing/'

# -- Period --
start = datetime.datetime(year=2017, month=5, day=1, hour=0)
end   = datetime.datetime(year=2017, month=9, day=1, hour=0)

# -- Read Harmonie statistics --
variables = ['u', 'v', 'z', 'time']
if 'ham' not in locals():
    ham  = read_DDH_netcdf(start, end, path, variables)


iloc = 0+12     # FINO1 10x10km
k    = 4        # 0 ~ 12m, 1 ~ 37m, 2 ~62m, 3~88m, 4~115m
n    = 12

pl.figure()
ax = pl.subplot(111)

if 'wdir' not in locals():
    wdir = components_to_direction(ham.u[:,iloc,k], ham.v[:,iloc,k])
mask = (wdir>=235)&(wdir<=285)



pl.plot(ham.time, wdir, 'k-', linewidth=0.5)
#pl.plot(ham.time[mask], wdir[mask], 'k-')

# Plot clear and disturbed sectors
lim = ax.get_xlim()
pl.fill_between(lim, [235, 235], [285, 285], color='g', alpha=0.3)
pl.fill_between(lim, [0, 0],     [180, 180], color='r', alpha=0.3)
pl.xlim(start, end)



