import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import xarray as xr
import pandas as pd
import numpy as np
import datetime
import sys
import os
from scipy import interpolate

# Add src directory to Python path, and import DALES specific tools
src_dir = os.path.abspath('{}/../../../src/'.format(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(src_dir)
from read_DDH_netcdf import *

pl.close('all')

def format_date_hour(interval):
    hours     = mdates.HourLocator(interval=interval)
    hours_fmt = mdates.DateFormatter('%H:%M')

    ax = pl.gca()
    ax.xaxis.set_major_locator(hours)
    ax.xaxis.set_major_formatter(hours_fmt)

def filter(arr, N):
    return np.convolve(arr, np.ones((N,))/N, mode='same')


# -- Period --
start = datetime.datetime(year=2017, month=4, day=2, hour=0)
end   = datetime.datetime(year=2017, month=4, day=6, hour=0)

# -- Read DALES statistics --
daa  = xr.open_dataset('../profiles.001.nc')
das  = xr.open_dataset('../tmser.001.nc')

daa_time = [start + datetime.timedelta(seconds = int(t)) for t in daa.time.values]
das_time = [start + datetime.timedelta(seconds = int(t)) for t in das.time.values]

if 'ham' not in locals():
    # -- Read Harmonie statistics --
    iloc = 7 #+12
    #path  = '/nobackup/users/stratum/DOWA/LES_forcing'
    path = '/Users/bart/meteo/data/Harmonie_LES_forcing/'
    ham  = read_DDH_netcdf(start, end, path)

if 'cb' not in locals():
    # -- Read Cabauw observations --
    path = '/Users/bart/meteo/observations/Cabauw'
    #path = '/nobackup/users/stratum/Cabauw'
    files = ['{0}/cesar_surface_radiation_lc1_t10_v1.0_{1:04d}{2:02d}.nc'    .format(path, start.year, start.month),
             '{0}/cesar_surface_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'        .format(path, start.year, start.month),
             '{0}/cesar_surface_flux_lc1_t10_v1.0_{1:04d}{2:02d}.nc'         .format(path, start.year, start.month),
             '{0}/cesar_soil_heat_lb1_t10_v1.0_{1:04d}{2:02d}.nc'            .format(path, start.year, start.month),
             '{0}/cesar_soil_water_lb1_t10_v1.1_{1:04d}{2:02d}.nc'           .format(path, start.year, start.month),
             '{0}/cesar_tower_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'          .format(path, start.year, start.month)]
    cb = xr.open_mfdataset(files)

    fnubi = '{0}/cesar_nubiscope_cloudcover_la1_t10_v1.0_{1:04d}{2:02d}.nc' .format(path, start.year, start.month)
    if os.path.exists(fnubi):
        cbc = xr.open_dataset(fnubi)

    nubi_type = {'LF':'light_fog', 'DF':'dense_fog', 'HP':'heavy_precipitation', 'LC':'low_transparent_clouds',
                 'TC':'transparent_clouds', 'OC':'overcast', 'BC':'broken_clouds', 'CI':'cirrus', 'CS':'clear_sky'}
    cbc_ct = cbc['obscuration_type'].to_dataframe()
    cbc_ct['obscuration_type'] = cbc_ct['obscuration_type'].str.decode("utf-8")
    for type in nubi_type.keys():
        cbc_ct[type] = cbc_ct['obscuration_type'] == type

# -- Colors et al. --
cd   = 'k'    # DALES color
ch   = 'C3'   # Harmonie color
co   = '#e41a1c'   # Obs color
mo   = 'o'    # Obs marker
ms   = 1.5     # Obs marker size
lw   = 1.5    # Linewidth
xint = 24     # Interval of x-markers
dash = [4,2]  # Format of dashed lines

if True:

    # Surface radiation balance
    # -------------------------
    pl.figure(figsize=(12,7))
    pl.suptitle('Cabauw: {} to {} UTC'.format(start, end), fontsize='medium')

    pl.subplot(241)

    pl.plot(daa_time, -daa.lwd[:,0], color=cd, linewidth=lw, label='DALES')
    pl.plot(cb.time, cb.LWD[:], mo, ms=ms, color=co, label='Cabauw')
    pl.legend()
    pl.xlim(start, end)
    pl.ylabel(r'LW$_\mathrm{in}$ (W m$^{-2}$)')
    format_date_hour(xint)

    pl.subplot(242)
    pl.plot(daa_time, daa.lwu[:,0], color=cd, linewidth=lw)
    pl.plot(cb.time, cb.LWU[:], mo, ms=ms, color=co)
    pl.xlim(start, end)
    pl.ylabel(r'LW$_\mathrm{out}$ (W m$^{-2}$)')
    format_date_hour(xint)

    pl.subplot(245)
    pl.plot(daa_time, -daa.swd[:,0], color=cd, linewidth=lw)
    pl.plot(cb.time, cb.SWD[:], mo, ms=ms, color=co)
    pl.xlim(start, end)
    pl.ylabel(r'SW$_\mathrm{in}$ (W m$^{-2}$)')
    pl.xlabel('time UTC')
    format_date_hour(xint)

    pl.subplot(246)
    pl.plot(daa_time, daa.swu[:,0], color=cd, linewidth=lw)
    pl.plot(cb.time, cb.SWU[:], mo, ms=ms, color=co)
    pl.xlim(start, end)
    pl.ylabel(r'SW$_\mathrm{out}$ (W m$^{-2}$)')
    pl.xlabel('time UTC')
    format_date_hour(xint)

    # Surface fluxes
    # --------------
    pl.subplot(243)
    pl.plot(das_time, das.H, color=cd, label='DALES H')
    #pl.plot(ham.time, -ham.H[:,iloc], color=ch, linewidth=lw, label='Harmonie H')
    pl.plot(cb.time, cb.H, mo, ms=ms, color=co, linewidth=lw, label='Cabauw H')
    pl.legend()
    pl.xlim(start, end)
    pl.ylabel(r'H (W m$^{-2}$)')
    pl.ylim(-50,200)
    format_date_hour(xint)

    pl.subplot(244)
    pl.plot(das_time, das.LE, color=cd)
    #pl.plot(ham.time, -ham.LE[:,iloc], color=ch, linewidth=lw)
    pl.plot(cb.time, cb.LE, mo, ms=ms, color=co, linewidth=lw)
    pl.xlim(start, end)
    pl.ylabel(r'LE (W m$^{-2}$)')
    pl.ylim(-50,300)
    format_date_hour(xint)

    pl.subplot(247)
    pl.plot(das_time, das.G0, color=cd)
    pl.plot(cb.time, cb.G0, mo, ms=ms, color=co, linewidth=lw)
    pl.xlim(start, end)
    pl.ylabel(r'G (W m$^{-2}$)')
    pl.xlabel('time UTC')
    format_date_hour(xint)

    pl.subplot(248)
    pl.plot(das_time, das.ustar, color=cd, linewidth=lw, label='DALES H')
    pl.plot(cb.time, cb.UST, mo, ms=ms, color=co, linewidth=lw)
    pl.xlim(start, end)
    pl.ylabel(r'u* (m s$^{-1}$)')
    pl.xlabel('time UTC')
    format_date_hour(xint)

    pl.tight_layout()
    fig = pl.gcf()
    fig.subplots_adjust(top=0.95)



if False:
    # Precipitation
    # --------------
    pl.figure(figsize=(12,7))
    pl.suptitle('Cabauw: {} to {} UTC'.format(start, end), fontsize='medium')

    pl.subplot(131)
    pl.plot(cb.time, cb.RAIN*6, mo, ms=ms, color=co, label='Cabauw')
    pl.plot(daa_time, daa.rainrate[:,0]/(daa.rhof[:,0]*2.45e6)*3600, color=cd, linewidth=lw, label='DALES')
    pl.xlim(start, end)
    pl.xlabel('time UTC')
    pl.ylabel(r'rainrate (mm h$^{-1}$)')
    format_date_hour(xint)

    pl.subplot(132)
    pl.plot(cbc.time, cbc.cldcover_total, mo, ms=ms, color=co)
    pl.plot(daa_time, np.max(daa.cfrac*100, axis=1), color=cd, linewidth=lw)
    pl.xlim(start, end)
    pl.xlabel('time UTC')
    pl.ylabel('cloud fraction (%)')
    format_date_hour(xint)

    pl.subplot(133)
    for type in nubi_type.keys():
        pl.plot(cbc_ct.index, filter(cbc_ct[type], 6), label=nubi_type[type])
    pl.xlim(start, end)
    pl.xlabel('time UTC')
    pl.ylabel('cloud type')
    pl.legend()
    format_date_hour(xint)



if True:
    # Soil temperature & moisture
    # ----------------
    exns = (daa.presh[:,0]/1e5)**(287/1004.)     # Not really surface; first model level... Ps is not in DALES output.....

    pl.figure(figsize=(12,6))

    pl.subplot(131)
    pl.plot(daa_time, exns*daa.thl[:,0  ], color='C1',  linewidth=lw, label='DALES Ta(10m)')
    pl.plot(das_time, exns*das.thlskin[:], color='k',   linewidth=lw, label='DALES Ts')

    pl.plot(cb.time, cb.TA[:,-2],                     color='C1',  linewidth=lw, label='Cabauw Ta(10m)', dashes=dash)
    pl.plot(cb.time, (cb.LWU / (0.98*5.67e-8))**0.25, color='k',   linewidth=lw, label='Cabauw Ts',  dashes=dash)

    pl.xlabel('time UTC')
    pl.ylabel('T (K)')
    pl.legend()
    pl.xlim(start, end)
    format_date_hour(xint)

    pl.subplot(132)
    pl.plot(daa_time, daa.tsoil[:,0],      color='C1',   linewidth=lw, label='DALES {0:.1f} cm'.format(-daa.zts[0].values*100))
    pl.plot(daa_time, daa.tsoil[:,1],      color='C2',   linewidth=lw, label='DALES {0:.1f} cm'.format(-daa.zts[1].values*100))
    pl.plot(daa_time, daa.tsoil[:,2],      color='C3',   linewidth=lw, label='DALES {0:.1f} cm'.format(-daa.zts[2].values*100))

    pl.plot(cb.time, cb.TS04+273.15, color='C1', linewidth=lw, label='Cabauw 4 cm',  dashes=dash)
    pl.plot(cb.time, cb.TS20+273.15, color='C2', linewidth=lw, label='Cabauw 20 cm', dashes=dash)
    pl.plot(cb.time, cb.TS50+273.15, color='C3', linewidth=lw, label='Cabauw 50 cm', dashes=dash)

    pl.xlabel('time UTC')
    pl.ylabel('T (K)')
    pl.legend()
    pl.xlim(start, end)
    format_date_hour(xint)

    pl.subplot(133)
    pl.plot(daa_time, daa.phiw[:,0], color='k',  linewidth=lw, label='DALES {0:.1f} cm'.format(-daa.zts[0].values*100))
    pl.plot(daa_time, daa.phiw[:,1], color='C1', linewidth=lw, label='DALES {0:.1f} cm'.format(-daa.zts[1].values*100))
    pl.plot(daa_time, daa.phiw[:,2], color='C2', linewidth=lw, label='DALES {0:.1f} cm'.format(-daa.zts[2].values*100))
    #pl.plot(daa_time, daa.phiw[:,3], color='C3', linewidth=lw, label='DALES {0:.1f} cm'.format(-daa.zts[3].values*100))

    pl.plot(cb.time, cb.TH05, color='k',  linewidth=lw, label='Cabauw 5 cm', dashes=dash)
    pl.plot(cb.time, cb.TH19, color='C1', linewidth=lw, label='Cabauw 19 cm', dashes=dash)
    #pl.plot(cb.time, cb.TH33, color='C2',linewidth=lw, label='Cabauw 33 cm', dashes=dash)
    pl.plot(cb.time, cb.TH40, color='C2', linewidth=lw, label='Cabauw 40 cm', dashes=dash)

    pl.xlabel('time UTC')
    pl.ylabel(r'$\phi_w$ (m$^3$ m$^{-3}$)')
    pl.legend()
    pl.xlim(start, end)
    format_date_hour(xint)
    pl.tight_layout()


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



#if True:
