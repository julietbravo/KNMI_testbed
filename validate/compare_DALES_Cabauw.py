import matplotlib.pyplot as pl
import xarray as xr
import pandas as pd
import numpy as np
import datetime

from read_DDH_netcdf import *

pl.close('all')

# -- Period --
start = datetime.datetime(year=2017, month=3, day=27, hour=12)
end   = datetime.datetime(year=2017, month=3, day=28, hour=21)

# -- Read Harmonie statistics --
iloc = 7+12
path = '/Users/bart/meteo/data/Harmonie_DDH/'
ham  = read_DDH_netcdf(start, end, path)

# -- Read DALES statistics --
daa  = xr.open_dataset('../cases/cabauw/profiles.001.nc')
das  = xr.open_dataset('../cases/cabauw/tmser.001.nc')

daa_time = [start + datetime.timedelta(seconds = int(t)) for t in daa.time.values]
das_time = [start + datetime.timedelta(seconds = int(t)) for t in das.time.values]

# -- Read Cabauw observations --
path = '/Users/bart/meteo/observations/Cabauw'
cbr  = xr.open_dataset('{0}/cesar_surface_radiation_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path, start.year, start.month))
cbs  = xr.open_dataset('{0}/cesar_surface_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path, start.year, start.month))
cbf  = xr.open_dataset('{0}/cesar_surface_flux_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path, start.year, start.month))
cbso = xr.open_dataset('{0}/cesar_soil_heat_lb1_t10_v1.0_{1:04d}{2:02d}.nc'.format(path, start.year, start.month))


# -- Colors --
cd = 'k'
ch = 'r'
co = '0.2'

# Surface/soil
# -------------
#pl.figure()
#pl.plot(das.

if True:
    # Surface radiation balance
    # -------------------------
    pl.figure()
    pl.subplot(221)
    pl.plot(daa_time, -daa.lwd[:,0], color=cd, label='DALES')
    pl.plot(cbr.time, cbr.LWD[:], '+', color=co, label='Cabauw')
    pl.legend()
    pl.xlim(start, end)
    pl.ylabel('LW_in (W m-2)')

    pl.subplot(222)
    pl.plot(daa_time, daa.lwu[:,0], color=cd)
    pl.plot(cbr.time, cbr.LWU[:], '+', color=co)
    pl.xlim(start, end)
    pl.ylabel('LW_out (W m-2)')

    pl.subplot(223)
    pl.plot(daa_time, -daa.swd[:,0], color=cd)
    pl.plot(cbr.time, cbr.SWD[:], '+', color=co)
    pl.xlim(start, end)
    pl.ylabel('SW_in (W m-2)')

    pl.subplot(224)
    pl.plot(daa_time, daa.swu[:,0], color=cd)
    pl.plot(cbr.time, cbr.SWU[:], '+', color=co)
    pl.xlim(start, end)
    pl.ylabel('SW_out (W m-2)')

    pl.tight_layout()

if True:
    # Soil temperature
    # ----------------
    pl.figure()
    pl.plot(das_time, das.thlskin[:], '-',  color=cd, linewidth=2, label='DALES 0 cm')
    pl.plot(daa_time, daa.tsoil[:,0], '-',  color=cd, label='DALES {0:.1f} cm'.format(-daa.zts[0].values*100))
    pl.plot(daa_time, daa.tsoil[:,1], '--', color=cd, label='DALES {0:.1f} cm'.format(-daa.zts[1].values*100))
    pl.plot(daa_time, daa.tsoil[:,2], '-.', color=cd, label='DALES {0:.1f} cm'.format(-daa.zts[2].values*100))
    pl.plot(daa_time, daa.tsoil[:,3], ':',  color=cd, label='DALES {0:.1f} cm'.format(-daa.zts[3].values*100))

    pl.plot(cbso.time, cbso.TS00+273.15, '-', linewidth=2,  color=ch, label='Cabauw 0 cm')
    pl.plot(cbso.time, cbso.TS04+273.15, '-',  color=ch, label='Cabauw 4 cm')
    pl.plot(cbso.time, cbso.TS20+273.15, '--',  color=ch, label='Cabauw 20 cm')
    pl.plot(cbso.time, cbso.TS50+273.15, '-.',  color=ch, label='Cabauw 50 cm')

    pl.legend()
    pl.xlim(start, end)

    cbso_t0 = cbso.sel(time=start)
    z_cb = np.array([0,2,4,8,12,20,50])
    T_cb = np.ma.zeros(z_cb.size)
    for k,z in enumerate(z_cb):
        try:
            T_cb[k] = cbso_t0['TS{0:02d}'.format(z)].values
        except:
            T_cb[k] = np.ma.masked

    pl.figure()
    pl.plot(T_cb+273.15, -z_cb, '-x')
    pl.plot(daa.tsoil[0,:], daa.zts*100, '-o')


if True:
    # Surface fluxes
    # --------------
    pl.figure()
    pl.subplot(221)
    pl.plot(das_time, das.H, color=cd, label='DALES H')
    pl.plot(ham.time, -ham.H[:,iloc], color=ch, label='Harmonie H')
    pl.plot(cbf.time, cbf.H, 'x', color=co, label='Cabauw H')
    pl.legend()
    pl.xlim(start, end)
    pl.ylabel('H (W m-2)')

    pl.subplot(222)
    pl.plot(das_time, das.LE, color=cd)
    pl.plot(ham.time, -ham.LE[:,iloc], color=ch)
    pl.plot(cbf.time, cbf.LE, 'x', color=co)
    pl.xlim(start, end)
    pl.ylabel('LE (W m-2)')

    pl.subplot(223)
    pl.plot(das_time, das.G0, color=cd)
    pl.plot(cbf.time, cbf.G0, 'x', color=co)
    pl.xlim(start, end)
    pl.ylabel('G (W m-2)')

    pl.subplot(224)
    pl.plot(das_time, das.ustar, color=cd, label='DALES H')
    pl.plot(cbf.time, cbf.UST, 'x', color=co)
    pl.xlim(start, end)
    pl.ylabel('u* (m s-1)')

    pl.tight_layout()
