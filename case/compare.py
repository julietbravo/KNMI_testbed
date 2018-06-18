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

#-----------------------------------
#import matplotlib as mpl
#import matplotlib.font_manager as fm
#font_dirs = ['/Users/bart/Library/Fonts/']
#font_files = fm.findSystemFonts(fontpaths=font_dirs)
#font_list = fm.createFontList(font_files)
#fm.fontManager.ttflist.extend(font_list)
#
#prop = fm.FontProperties(fname='/Users/bart/Library/Fonts/RO-SerifWeb-Regular.ttf')
#prop = fm.FontProperties(fname='/Users/bart/Library/Fonts/rijksoverheidsanstext-regular.ttf')
#mpl.rcParams['font.family'] = prop.get_name()
#-----------------------------------

pl.close('all')

# LS2D scripts, to calculate ERA5 tendencies
import sys; sys.path.append('/usr/people/stratum/meteo/models/LS2D/src')
import sys; sys.path.append('/Users/bart/meteo/models/LS2D/src')
from read_ERA5 import Read_ERA

pl.close('all')

def format_h_since(hours, start):
    return datetime.timedelta(hours=float(hours)) + start

def wind_to_components(speed, direction):
    u = -speed * np.sin(np.deg2rad(direction))
    v = -speed * np.cos(np.deg2rad(direction))
    return u,v

def components_to_direction(u, v):
    return np.rad2deg(np.arctan2(u, v)) + 180

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

    z = 10.     # Height to compare (m)

    # Fixed colors per model/obs type
    ch = 'C3'   # Harmonie
    cd = 'k'    # Dales
    ce = 'C0'   # ERA5
    cc = 'C7'   # Cabauw

    mpl.rcParams['lines.linewidth'] = 1.3


    if 'cb' not in locals():

        # Read Cabauw observations
        # ========================
        #pwd  = '/nobackup/users/stratum/Cabauw'
        pwd  = '/Users/bart/meteo/data/Cabauw'
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
            #'ERA5_path'   : '/nobackup/users/stratum/ERA5/LS2D/',
            'ERA5_path'   : '/Users/bart/meteo/data/ERA5/LS2D/',
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
        #pwd   = '/nobackup/users/stratum/DOWA/LES_forcing/'
        pwd   = '/Users/bart/meteo/data/Harmonie_DDH/'
        files = glob.glob('{}LES_forcings_20100228*'.format(pwd))
        files.sort()
        hm    = xr.open_mfdataset(files)

        # Interpolate wind to 10m
        hm_u = interp(hm.u.values, hm.zg.values[0,:], z)
        hm_v = interp(hm.v.values, hm.zg.values[0,:], z)

        hm_T = interp(hm['T'].values, hm.zg.values[0,:], z)
        hm_q = interp(hm.q.values, hm.zg.values[0,:], z)

        # Interpolate tendencies
        hm_dtu_dyn = interp(hm.dtu_dyn.values, hm.zg.values[0,:], z)
        hm_dtu_phy = interp(hm.dtu_phy.values, hm.zg.values[0,:], z)

        hm_dtv_dyn = interp(hm.dtv_dyn.values, hm.zg.values[0,:], z)
        hm_dtv_phy = interp(hm.dtv_phy.values, hm.zg.values[0,:], z)

        hm_dtT_dyn = interp(hm.dtT_dyn.values, hm.zg.values[0,:], z)
        hm_dtT_phy = interp(hm.dtT_phy.values, hm.zg.values[0,:], z)

        hm_dtq_dyn = interp(hm.dtq_dyn.values, hm.zg.values[0,:], z)
        hm_dtq_phy = interp(hm.dtq_phy.values, hm.zg.values[0,:], z)

        # DALES runs
        # ==========
        da   = xr.open_dataset('profiles.force10min.nc')
        t0   = datetime.datetime(2010, 2, 28, 6)
        time = [format_h_since(s/3600., t0) for s in da.time]

        da_u = interp(da.u.values, da.zt.values, z)
        da_v = interp(da.v.values, da.zt.values, z)

        da_p = interp(da.presh, da.zt.values, z)
        da_exn = (da_p/1e5)**(287/1004)

        da_T = interp(da.thl.values, da.zt.values, z) * da_exn

        da_q = interp(da.qt.values, da.zt.values, z)


    if (False):
        # ------------------------
        # Demonstration tendencies
        # ------------------------

        k = 0           # Level
        height = 10.    # Height label in plot
        ev = 2          # Plot every `ev`'th time step

        xlim = [datetime.datetime(2010,2,28,14,45), datetime.datetime(2010,2,28,16,15)]

        # u, v components
        # ------------------------
        pl.figure(figsize=(8,5))
        gs = gridspec.GridSpec(2, 2, height_ratios=[1,2])

        pl.subplot(gs[0,0])
        pl.plot(hm.time[::ev], hm.u[::ev,k], color='C3', label='$u$')
        pl.plot(hm.time[::ev], hm.v[::ev,k], color='C0', label='$v$')
        pl.legend()
        pl.xlim(xlim)
        format_ax(interval=1)
        pl.ylabel('Speed ($\mathrm{m s^{-1}}$)')

        pl.subplot(gs[0,1])
        pl.plot(hm.time[::ev], components_to_direction(hm.u[::ev,k], hm.v[::ev,k]))
        pl.xlim(xlim)
        format_ax(interval=1)
        pl.ylim(170,300)
        pl.ylabel('Direction (deg) ')

        pl.subplot(gs[1,0])
        pl.plot(hm.time[::ev], hm.dtu_dyn[::ev,k]*3600, color='C3', label='dynamics')
        pl.plot(hm.time[::ev], hm.dtu_phy[::ev,k]*3600, color='C0', label='physics')
        pl.plot(hm.time[::ev], hm.dtu_tot[::ev,k]*3600, color='k', label='total', dashes=[2,1])
        pl.plot(e5.datetime,   (e5.dtu_advec+e5.dtu_coriolis)[:,0]*3600, '-x', label='ERA5', color='C9')
        pl.legend()
        pl.xlim(xlim)
        format_ax(interval=1)
        pl.ylabel('$\partial_t u$ ($\mathrm{m s^{-1} h^{-1}}$)')
        pl.xlabel('time (UTC)')

        pl.subplot(gs[1,1])
        pl.plot(hm.time[::ev], hm.dtv_dyn[::ev,k]*3600, color='C3')
        pl.plot(hm.time[::ev], hm.dtv_phy[::ev,k]*3600, color='C0')
        pl.plot(hm.time[::ev], hm.dtv_tot[::ev,k]*3600, color='k', dashes=[2,1])
        pl.plot(e5.datetime,   (e5.dtv_advec+e5.dtv_coriolis)[:,0]*3600, '-x', label='ERA5', color='C9')
        pl.xlim(xlim)
        format_ax(interval=1)
        pl.ylabel('$\partial_t v$ ($\mathrm{m s^{-1} h^{-1}}$)')
        pl.xlabel('time (UTC)')

        pl.tight_layout()


        # Temp & moisture
        # ------------------------
        pl.figure(figsize=(8,5))
        gs = gridspec.GridSpec(2, 2, height_ratios=[1,2])

        pl.subplot(gs[0,0])
        pl.plot(hm.time[::ev], hm['T'][::ev,k], color='C3')
        pl.xlim(xlim)
        format_ax(interval=1)
        pl.ylabel('$T$ ($\mathrm{K}$)')

        pl.subplot(gs[0,1])
        pl.plot(hm.time[::ev], hm['q'][::ev,k]*1000, color='C3')
        pl.xlim(xlim)
        format_ax(interval=1)
        pl.ylabel('$q_t$ ($\mathrm{kg kg^{-1}}$)')

        pl.subplot(gs[1,0])
        pl.plot(hm.time[::ev], hm.dtT_dyn[::ev,k]*3600, color='C3', label='dynamics')
        pl.plot(hm.time[::ev], hm.dtT_phy[::ev,k]*3600, color='C0', label='physics')
        pl.plot(hm.time[::ev], hm.dtT_tot[::ev,k]*3600, color='k', label='total', dashes=[2,1])
#        pl.plot(e5.datetime,   (e5.dtu_advec+e5.dtu_coriolis)[:,0]*3600, '-x', label='ERA5', color='C9')
        pl.legend()
        pl.xlim(xlim)
        format_ax(interval=1)
        pl.ylabel('$\partial_t T$ ($\mathrm{K h^{-1}}$)')
        pl.xlabel('time (UTC)')

        pl.subplot(gs[1,1])
        pl.plot(hm.time[::ev], hm.dtq_dyn[::ev,k]*3600*1000, color='C3')
        pl.plot(hm.time[::ev], hm.dtq_phy[::ev,k]*3600*1000, color='C0')
        pl.plot(hm.time[::ev], hm.dtq_tot[::ev,k]*3600*1000, color='k', dashes=[2,1])
#        pl.plot(e5.datetime,   (e5.dtv_advec+e5.dtv_coriolis)[:,0]*3600, '-x', label='ERA5', color='C9')
        pl.xlim(xlim)
        format_ax(interval=1)
        pl.ylabel('$\partial_t q_v$ ($\mathrm{kg kg^{-1} h^{-1}}$)')
        pl.xlabel('time (UTC)')

        pl.tight_layout()


    if (True):
        # ----------------
        # Compare DALES with others
        # ----------------

        xlim = [datetime.datetime(2010,2,28,6), datetime.datetime(2010,2,28,18)]

        pl.figure(figsize=(8,3))
        pl.subplot(121)
        pl.title('2010-02-28, Cabauw', loc='left')
        pl.plot(cb.time,     wind_to_components(cb.F, cb.D)[0][:,k_cb], '.', color='k', label='Cabauw', markersize=4)
        pl.plot(e5.datetime, e5_u, label='ERA5', dashes=[2,1], color=ce)
        pl.plot(hm.time,     hm_u, label='Harmonie', color=ch)
        pl.plot(time,        da_u, label='DALES', color=cd)
        format_ax()
        pl.legend(numpoints=3)
        pl.xlim(xlim)
        pl.ylabel('$u_{{{0:.0f} \mathrm{{m}}}}$ (m s$^{{-1}}$)'.format(z))
        pl.xlabel('time (UTC)')

        pl.subplot(122)
        pl.plot(cb.time,     wind_to_components(cb.F, cb.D)[1][:,k_cb], '.', color='k', label='Cabauw', markersize=4)
        pl.plot(e5.datetime, e5_v, label='ERA5', dashes=[2,1], color=ce)
        pl.plot(hm.time,     hm_v, color=ch)
        pl.plot(time,        da_v, color=cd)
        format_ax()
        pl.xlim(xlim)
        pl.ylabel('$v_{{{0:.0f} \mathrm{{m}}}}$ (m s$^{{-1}}$)'.format(z))
        pl.xlabel('time (UTC)')

        pl.tight_layout()
        pl.savefig('figures/uv_{0:.0f}m.pdf'.format(z))
        pl.savefig('figures/uv_{0:.0f}m.png'.format(z))



    if (False):
        # ----------------
        # Overview & tendencies Harmonie / ERA5
        # ----------------

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
