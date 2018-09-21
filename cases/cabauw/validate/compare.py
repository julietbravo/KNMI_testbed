import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import xarray as xr
import pandas as pd
import numpy as np
import datetime
import sys
import os

def format_date_hour(interval):
    hours     = mdates.HourLocator(interval=interval)
    hours_fmt = mdates.DateFormatter('%H:%M')

    ax = pl.gca()
    ax.xaxis.set_major_locator(hours)
    ax.xaxis.set_major_formatter(hours_fmt)

from scipy import interpolate

# Add src directory to Python path, and import DALES specific tools
src_dir = os.path.abspath('{}/../../../src/'.format(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(src_dir)
from read_DDH_netcdf import *

pl.close('all')

# -- Period --
#start = datetime.datetime(year=2017, month=9, day=8, hour=0)
#end   = datetime.datetime(year=2017, month=9, day=10, hour=3)

start = datetime.datetime(year=2017, month=6, day=9, hour=0)
end   = datetime.datetime(year=2017, month=6, day=10, hour=3)

# -- Read Harmonie statistics --
iloc = 7+12
path  = '/nobackup/users/stratum/DOWA/LES_forcing'
#path = '/Users/bart/meteo/data/Harmonie_DDH/'
ham  = read_DDH_netcdf(start, end, path)

# -- Read DALES statistics --
daa  = xr.open_dataset('../profiles.001.nc')
das  = xr.open_dataset('../tmser.001.nc')

daa_time = [start + datetime.timedelta(seconds = int(t)) for t in daa.time.values]
das_time = [start + datetime.timedelta(seconds = int(t)) for t in das.time.values]

# -- Read Cabauw observations --
#path = '/Users/bart/meteo/observations/Cabauw'
path = '/nobackup/users/stratum/Cabauw'
files = ['{0}/cesar_surface_radiation_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path, start.year, start.month),
         '{0}/cesar_surface_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'    .format(path, start.year, start.month),
         '{0}/cesar_surface_flux_lc1_t10_v1.0_{1:04d}{2:02d}.nc'     .format(path, start.year, start.month),
         '{0}/cesar_soil_heat_lb1_t10_v1.0_{1:04d}{2:02d}.nc'        .format(path, start.year, start.month),
         '{0}/cesar_soil_water_lb1_t10_v1.1_{1:04d}{2:02d}.nc'       .format(path, start.year, start.month)]
cb = xr.open_mfdataset(files)

# -- Colors --
cd = 'k'
ch = 'r'
co = 'C2'


if True:
    # Surface radiation balance
    # -------------------------
    pl.figure(figsize=(12,8))
    pl.subplot(221)
    pl.plot(daa_time, -daa.lwd[:,0], color=cd, label='DALES')
    pl.plot(cb.time, cb.LWD[:], '+', ms=3, color=co, label='Cabauw')
    pl.legend()
    pl.xlim(start, end)
    pl.ylabel('LW_in (W m-2)')
    format_date_hour(3)

    pl.subplot(222)
    pl.plot(daa_time, daa.lwu[:,0], color=cd)
    pl.plot(cb.time, cb.LWU[:], '+', ms=3, color=co)
    pl.xlim(start, end)
    pl.ylabel('LW_out (W m-2)')
    format_date_hour(3)

    pl.subplot(223)
    pl.plot(daa_time, -daa.swd[:,0], color=cd)
    pl.plot(cb.time, cb.SWD[:], '+', ms=3, color=co)
    pl.xlim(start, end)
    pl.ylabel('SW_in (W m-2)')
    format_date_hour(3)

    pl.subplot(224)
    pl.plot(daa_time, daa.swu[:,0], color=cd)
    pl.plot(cb.time, cb.SWU[:], '+', ms=3, color=co)
    pl.xlim(start, end)
    pl.ylabel('SW_out (W m-2)')
    format_date_hour(3)

    pl.tight_layout()

if True:
    # Soil temperature & moisture
    # ----------------
    pl.figure(figsize=(12,8))
    pl.subplot(121)
    pl.plot(das_time, das.thlskin[:], '-',  color=cd, linewidth=2, label='DALES Ts')
    pl.plot(daa_time, daa.tsoil[:,0], '-',  color=cd, label='DALES {0:.1f} cm'.format(-daa.zts[0].values*100))
    pl.plot(daa_time, daa.tsoil[:,1], '--', color=cd, label='DALES {0:.1f} cm'.format(-daa.zts[1].values*100))
    pl.plot(daa_time, daa.tsoil[:,2], '-.', color=cd, label='DALES {0:.1f} cm'.format(-daa.zts[2].values*100))
    pl.plot(daa_time, daa.tsoil[:,3], ':',  color=cd, label='DALES {0:.1f} cm'.format(-daa.zts[3].values*100))

    pl.plot(cb.time, cb.TS00+273.15, '-',   linewidth=2,  color=ch, label='Cabauw Ts')
    pl.plot(cb.time, cb.TS04+273.15, '-',   color=ch, label='Cabauw 4 cm')
    pl.plot(cb.time, cb.TS20+273.15, '--',  color=ch, label='Cabauw 20 cm')
    pl.plot(cb.time, cb.TS50+273.15, '-.',  color=ch, label='Cabauw 50 cm')

    pl.legend()
    pl.xlim(start, end)
    format_date_hour(3)

    pl.subplot(122)
    pl.plot(daa_time, daa.phiw[:,0], '-',  color=cd, label='DALES {0:.1f} cm'.format(-daa.zts[0].values*100))
    pl.plot(daa_time, daa.phiw[:,1], '--', color=cd, label='DALES {0:.1f} cm'.format(-daa.zts[1].values*100))
    pl.plot(daa_time, daa.phiw[:,2], '-.', color=cd, label='DALES {0:.1f} cm'.format(-daa.zts[2].values*100))
    pl.plot(daa_time, daa.phiw[:,3], ':',  color=cd, label='DALES {0:.1f} cm'.format(-daa.zts[3].values*100))

    pl.plot(cb.time, cb.TH05, '-',   color=ch, label='Cabauw 5 cm')
    pl.plot(cb.time, cb.TH19, '--',  color=ch, label='Cabauw 19 cm')
    pl.plot(cb.time, cb.TH33, '-.',  color=ch, label='Cabauw 33 cm')
    pl.plot(cb.time, cb.TH40, ':',   color=ch, label='Cabauw 40 cm')

    pl.legend()
    pl.xlim(start, end)
    format_date_hour(3)


if False:
    # Soil temperature/moisture profiles
    cb1    = cb.sel(time=start, method='nearest')
    z_cb_q = np.array([5,19,33,40,56])
    z_cb_T = np.array([0,2,4,6,8,12,20,30,50])

    cb_q = np.zeros_like(z_cb_q, dtype=np.float)
    cb_T = np.zeros_like(z_cb_T, dtype=np.float)

    for k,z in enumerate(z_cb_q):
        cb_q[k] = cb1['TH{0:02d}'.format(z)].values

    for k,z in enumerate(z_cb_T):
        cb_T[k] = cb1['TS{0:02d}'.format(z)].values

    # yikes
    ip_T = interpolate.interp1d(z_cb_T, cb_T+273.15, kind='slinear', fill_value='extrapolate')
    ip_q = interpolate.interp1d(z_cb_q, cb_q       , kind='slinear', fill_value='extrapolate')

    T_dales = ip_T(-daa.zts*100)
    q_dales = ip_q(-daa.zts*100)

    pl.figure()
    pl.subplot(121)
    pl.plot(cb_T+273.15, -z_cb_T, '-x', label='Cabauw')
    pl.plot(daa.tsoil[0,:-1], daa.zts[:-1]*100, '-s', label='DALES')
    pl.plot(T_dales[:-1], daa.zts[:-1]*100, '-o', label='DALES new')
    pl.legend()
    pl.ylabel('z (cm)')
    pl.ylabel('T (K)')

    pl.subplot(122)
    pl.plot(cb_q, -z_cb_q, '-x', label='Cabauw')
    pl.plot(daa.phiw[0,:-1], daa.zts[:-1]*100, '-s', label='DALES')
    pl.plot(q_dales, daa.zts*100, '-o', label='DALES new')
    pl.legend()
    pl.ylabel('z (cm)')
    pl.ylabel('phi (-)')


if True:
    # Precipitation
    # --------------
    pl.figure(figsize=(12,8))
    pl.plot(cb.time, cb.RAIN*6, 'x', ms=3, color=co, label='Cabauw H')
    pl.plot(daa_time, daa.rainrate[:,0]/(daa.rhof[:,0]*2.45e6)*3600, '-', label='DALES')
    pl.xlim(start, end)
    pl.xlabel('date')
    pl.ylabel('rainrate (mm h-1)')
    format_date_hour(3)

if True:
    # Surface fluxes
    # --------------
    pl.figure(figsize=(12,8))
    pl.subplot(221)
    pl.plot(das_time, das.H, color=cd, label='DALES H')
    pl.plot(ham.time, -ham.H[:,iloc], color=ch, label='Harmonie H')
    pl.plot(cb.time, cb.H, 'x', ms=3, color=co, label='Cabauw H')
    pl.legend()
    pl.xlim(start, end)
    pl.ylabel('H (W m-2)')
    format_date_hour(3)

    pl.subplot(222)
    pl.plot(das_time, das.LE, color=cd)
    pl.plot(ham.time, -ham.LE[:,iloc], color=ch)
    pl.plot(cb.time, cb.LE, 'x', ms=3, color=co)
    pl.xlim(start, end)
    pl.ylabel('LE (W m-2)')
    format_date_hour(3)

    pl.subplot(223)
    pl.plot(das_time, das.G0, color=cd)
    pl.plot(cb.time, cb.G0, 'x', ms=3, color=co)
    pl.xlim(start, end)
    pl.ylabel('G (W m-2)')
    format_date_hour(3)

    pl.subplot(224)
    pl.plot(das_time, das.ustar, color=cd, label='DALES H')
    pl.plot(cb.time, cb.UST, 'x', ms=3, color=co)
    pl.xlim(start, end)
    pl.ylabel('u* (m s-1)')
    format_date_hour(3)

    pl.tight_layout()
