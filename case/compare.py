#
# Compare DALES with Harmonie, ERA5, Cabauw observations
#

import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import matplotlib as mpl
import numpy as np
import xarray as xr

import datetime
import glob

from scipy import interpolate

pl.close('all')

# LS2D scripts, to calculate ERA5 tendencies
import sys; sys.path.append('/usr/people/stratum/meteo/models/LS2D/src')
from read_ERA5 import Read_ERA

pl.close('all')

def format_h_since(hours, start):
    return datetime.timedelta(hours=float(hours)) + start

def wind_to_components(speed, direction):
    u = -speed * np.sin(np.deg2rad(direction))
    v = -speed * np.cos(np.deg2rad(direction))
    return u,v

def interp(data, height, goal_height):
    # First index below goal height:
    klow = np.abs(height-goal_height).argmin()

    if (height[klow] > goal_height):
        klow -= 1
    khigh = klow+1

    if (klow >= 0):
        # Interpolate
        fhigh = (goal_height - height[klow]) / (height[khigh] - height[klow])
        flow  = 1-fhigh
        return flow*data[:,klow] + fhigh*data[:,khigh]
    else:
        # Extrapolate
        klow  += 1
        khigh += 1
        grad = (data[:,khigh] - data[:,klow]) / (height[khigh] - height[klow])
        return data[:,klow] - grad * (height[klow]-goal_height)

def format_ax(ax=None, interval=2):
    if ax is None:
        ax = pl.gca()

    hours = mdates.HourLocator(interval=interval)
    hours_fmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_locator(hours)
    ax.xaxis.set_major_formatter(hours_fmt)

def integrate_tend(u, v, dtu, dtv):
    ui = np.zeros_like(dtu)
    vi = np.zeros_like(dtu)

    ui[0,:] = u[:]
    vi[0,:] = v[:]

    for t in range(dtu.shape[0]-1):
        ui[t+1,:] = ui[t,:] + dtu[t,:] * 60.
        vi[t+1,:] = vi[t,:] + dtv[t,:] * 60.

    return ui, vi


if __name__ == '__main__':

    z = 200.     # Height to compare (m)

    # Fixed colors per model/obs type
    ch = 'C1'   # Harmonie
    cd = 'k'    # Dales
    ce = 'C2'   # ERA5
    cc = 'C3'   # Cabauw

    mpl.rcParams['lines.linewidth'] = 1.3

    # Read Cabauw observations
    # ========================
    pwd  = '/nobackup/users/stratum/Cabauw'
    cb   = xr.open_dataset('{}/cesar_tower_meteo_lc1_t10_v1.0_201002.nc'.format(pwd))
    k_cb = np.abs(cb.z.values-z).argmin()

    if (np.abs(cb.z[k_cb] - z) > 0.1):
        print('Height is not an observation height at Cabauw!!')

    # ERA5 data & tendencies
    # ======================
    settings = {
        'central_lat' : 51.971,
        'central_lon' : 4.927,
        'area_size'   : 2,
        'case_name'   : 'cabauw',
        'ERA5_path'   : '/nobackup/users/stratum/ERA5/LS2D/',
        'start_date'  : datetime.datetime(year=2010, month=2, day=28, hour=0),
        'end_date'    : datetime.datetime(year=2010, month=2, day=28, hour=23)
        }

    e5 = Read_ERA(settings)
    e5.calculate_forcings(n_av=0)
    t0 = datetime.datetime(1900, 1, 1)
    e5.datetime = [format_h_since(h, t0) for h in e5.time]

    # Interpolate wind to 10m
    e5_u = interp(e5.u_mean, e5.z_mean[0,:], z)
    e5_v = interp(e5.v_mean, e5.z_mean[0,:], z)

    # Interpolate dynamic tendency
    tmp = e5.dtu_advec+e5.dtu_coriolis
    e5_dtu_dyn = interp(tmp, e5.z_mean[0,:], z)
    tmp = e5.dtv_advec+e5.dtv_coriolis
    e5_dtv_dyn = interp(tmp, e5.z_mean[0,:], z)

    # Harmonie data & tendencies
    # ==========================
    pwd   = '/nobackup/users/stratum/DOWA/LES_forcing/'
    files = glob.glob('{}LES_forcings_20100228*'.format(pwd))
    files.sort()
    hm    = xr.open_mfdataset(files)

    # Interpolate wind to 10m
    hm_u = interp(hm.u.values, hm.zg.values[0,:], z)
    hm_v = interp(hm.v.values, hm.zg.values[0,:], z)

    # Interpolate tendencies
    hm_dtu_dyn = interp(hm.dtu_dyn.values, hm.zg.values[0,:], z)
    hm_dtu_phy = interp(hm.dtu_phy.values, hm.zg.values[0,:], z)

    hm_dtv_dyn = interp(hm.dtv_dyn.values, hm.zg.values[0,:], z)
    hm_dtv_phy = interp(hm.dtv_phy.values, hm.zg.values[0,:], z)

    # DALES runs
    # ==========
    da   = xr.open_dataset('profiles.001.nc')
    t0   = datetime.datetime(2010, 2, 28, 6)
    time = [format_h_since(s/3600., t0) for s in da.time]

    da_u = interp(da.u.values, da.zt.values, z)
    da_v = interp(da.v.values, da.zt.values, z)

    if (True):
        # ----------------
        # Overview & tendencies Harmonie / ERA5
        # ----------------
        xlim = [datetime.datetime(2010,2,28,6), datetime.datetime(2010,2,28,18)]

        pl.figure(figsize=(8,6))

        gs = gridspec.GridSpec(2, 2, height_ratios=[1,2])

        pl.subplot(gs[0,0])
        pl.title('2010-02-28, Cabauw, {}m'.format(z), loc='left')
        pl.plot(cb.time, wind_to_components(cb.F, cb.D)[0][:,k_cb], '+', label='Cabauw')
        pl.plot(e5.datetime, e5_u, color=ce, label='ERA5', dashes=[2,2])
        pl.plot(hm.time,     hm_u, color=ch, label='Harmonie')
        format_ax()
        pl.legend()
        pl.xlim(xlim)
        pl.ylabel('$u$ (m s$^{-1}$)')

        pl.subplot(gs[0,1])
        pl.plot(cb.time, wind_to_components(cb.F, cb.D)[1][:,k_cb], '+', label='Cabauw')
        pl.plot(e5.datetime, e5_v, color=ce, label='ERA5', dashes=[2,2])
        pl.plot(hm.time,     hm_v, color=ch, label='Harmonie')
        format_ax()
        pl.legend()
        pl.xlim(xlim)
        pl.ylabel('$v$ (m s$^{-1}$)')

        pl.subplot(gs[1,0])
        pl.plot(hm.time,     hm_dtu_dyn*3600, color=ch, label='Harmonie dynamics')
        pl.plot(hm.time,     hm_dtu_phy*3600, color=ch, label='Harmonie physics', dashes=[1,1])
        pl.plot(e5.datetime, e5_dtu_dyn*3600, color=ce, label='ERA5 dynamics')
        format_ax()
        pl.legend()
        pl.xlim(xlim)
        pl.ylabel('$\partial_t u$ (m s$^{-1} h^{-1}$)')

        pl.subplot(gs[1,1])
        pl.plot(hm.time,     hm_dtv_dyn*3600, color=ch, label='Harmonie dynamics')
        pl.plot(hm.time,     hm_dtv_phy*3600, color=ch, label='Harmonie physics', dashes=[1,1])
        pl.plot(e5.datetime, e5_dtv_dyn*3600, color=ce, label='ERA5 dynamics')
        format_ax()
        pl.legend()
        pl.xlim(xlim)
        pl.ylabel('$\partial_t v$ (m s$^{-2}$)')


    if (True):
        # ----------------
        # Compare DALES with others
        # ----------------

        pl.figure(figsize=(8,4))
        pl.subplot(121)
        pl.title('2010-02-28, Cabauw', loc='left')
        pl.plot(cb.time,     wind_to_components(cb.F, cb.D)[0][:,k_cb], '+', label='Cabauw')
        pl.plot(e5.datetime, e5_u, label='ERA5', dashes=[2,2])
        pl.plot(hm.time,     hm_u, label='Harmonie')
        pl.plot(time,        da_u, label='DALES', color='k')
        format_ax()
        pl.legend()
        pl.xlim(xlim)
        pl.ylabel('$u_{{{0:.0f} m}}$ (m s$^{{-1}}$)'.format(z))
        pl.xlabel('time (UTC)')

        pl.subplot(122)
        pl.plot(cb.time,     wind_to_components(cb.F, cb.D)[1][:,k_cb], '+', label='Cabauw')
        pl.plot(e5.datetime, e5_v, label='ERA5', dashes=[2,2])
        pl.plot(hm.time,     hm_v, label='Harmonie')
        pl.plot(time,        da_v, label='DALES', color='k')
        format_ax()
        pl.xlim(xlim)
        pl.ylabel('$v_{{{0:.0f} m}}$ (m s$^{{-1}}$)'.format(z))
        pl.xlabel('time (UTC)')

        pl.tight_layout()
        pl.savefig('figures/uv_{0:.0f}m.pdf'.format(z))
        pl.savefig('figures/uv_{0:.0f}m.png'.format(z))


    if (False):

        # Integrate dynamic tendencies
        #udi, vdi = integrate_tend(hm.u.values[0], hm.v.values[0], hm.dtu_dyn.values, hm.dtv_dyn.values)
        #upi, vpi = integrate_tend(hm.u.values[0], hm.v.values[0], hm.dtu_phy.values, hm.dtv_phy.values)
        #ui, vi   = integrate_tend(hm.u.values[0], hm.v.values[0], hm.dtu_phy.values+hm.dtu_dyn.values, hm.dtv_phy.values+hm.dtv_dyn.values)

        #hm_u10_di = interp(udi, hm.zg.values[0,:], 10.)
        #hm_v10_di = interp(vdi, hm.zg.values[0,:], 10.)

        #hm_u10_pi = interp(upi, hm.zg.values[0,:], 10.)
        #hm_v10_pi = interp(vpi, hm.zg.values[0,:], 10.)

        #hm_u10_i = interp(ui, hm.zg.values[0,:], 10.)
        #hm_v10_i = interp(vi, hm.zg.values[0,:], 10.)

        # ----------------
        # Budget of u,v tendencies
        # ----------------
        xlim = [datetime.datetime(2010,2,28,5), datetime.datetime(2010,2,28,18)]

        pl.figure(figsize=(8,4))
        pl.subplot(121)
        pl.title('2010-02-28, Cabauw', loc='left')
        pl.plot(hm.time,     hm_u10, label='Harmonie')
        pl.plot(hm.time,     hm_u10_di, label='$\int \partial_t u$ dyn')
        pl.plot(hm.time,     hm_u10_pi, label='$\int \partial_t u$ phy')
        pl.plot(hm.time,     hm_u10_i, label='$\int \partial_t u$')
        format_ax()
        pl.legend()
        pl.xlim(xlim)
        pl.ylabel('$u_{10m}$ (m s$^{-1}$)')
        pl.xlabel('time (UTC)')

        pl.subplot(122)
        pl.plot(hm.time,     hm_v10, label='Harmonie')
        pl.plot(hm.time,     hm_v10_di, label='Harmonie dyn')
        pl.plot(hm.time,     hm_v10_pi, label='Harmonie phy')
        pl.plot(hm.time,     hm_v10_i, label='Harmonie tot')
        format_ax()
        pl.xlim(xlim)
        pl.ylabel('$v_{10m}$ (m s$^{-1}$)')
        pl.xlabel('time (UTC)')

        pl.tight_layout()


