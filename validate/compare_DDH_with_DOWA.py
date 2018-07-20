import numpy as np
import matplotlib.pyplot as pl
import xarray as xr

from spatial_tools import *

pl.close('all')

pwd  = '/nobackup/users/stratum/DOWA/LES_forcing/2010/02/28/06'
pwd2 = '/nobackup/users/stratum/DOWA/LES_forcing/2010/02/28/09'

# Read the processed DDH files
ddh  = xr.open_dataset('{}/LES_forcing_2010022806.nc'.format(pwd))
ddh2 = xr.open_dataset('{}/LES_forcing_2010022809.nc'.format(pwd2))

# Read the normal DOWA NetCDF files
tsk = xr.open_dataset('{}/ts.his.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc'.format(pwd))
sst = xr.open_dataset('{}/sst.sfx.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc'.format(pwd))

# Coordinates DDH output domains
locations = [dict(name='Cabauw',   lat=51.97, lon=4.90),
             dict(name='IJmuiden', lat=52.85, lon=3.44)]
sizes = [-1, 10000, 30000]

pl.figure()

loc = locations[1]  # IJmuiden
j,i = find_nearest_latlon(sst.lat.values, sst.lon.values, loc['lat'], loc['lon'])

pl.subplot(121)
pl.title('IJmuiden', loc='left')
pl.plot(sst['time'], sst['sst'][:,j,i], '-x', label='sst, 2D field')
pl.plot(tsk['time'], tsk['ts' ][:,j,i], '-+', label='tsk, 2D field')
pl.plot(ddh['time'], ddh['Tsk'][:,1], label='DDH, single column')
pl.plot(ddh2['time'], ddh2['Tsk'][:,1], label='DDH, single column')
#pl.plot(ddh['time'], ddh['Tsk'][:,3], label='DDH, 10x10km')
#pl.plot(ddh['time'], ddh['Tsk'][:,5], label='DDH, 30x30km')
pl.ylim(sst['sst'][:,j,i].mean()-0.5, sst['sst'][:,j,i].mean()+0.5)
pl.legend()

loc = locations[0]  # Cabauw
j,i = find_nearest_latlon(sst.lat.values, sst.lon.values, loc['lat'], loc['lon'])

pl.subplot(122)
pl.title('Cabauw', loc='left')
pl.plot(sst['time'], sst['sst'][:,j,i], '-x', label='sst, 2D field')
pl.plot(tsk['time'], tsk['ts' ][:,j,i], '-+', label='tsk, 2D field')
pl.plot(ddh['time'], ddh['Tsk'][:,0], label='DDH, single column')
pl.plot(ddh2['time'], ddh2['Tsk'][:,0], label='DDH, single column')
#pl.plot(ddh['time'], ddh['Tsk'][:,2], label='DDH, 10x10km')
#pl.plot(ddh['time'], ddh['Tsk'][:,4], label='DDH, 30x30km')
pl.legend()

