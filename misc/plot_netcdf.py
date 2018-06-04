import numpy as np
import xarray as xr
import datetime
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
from scipy import interpolate

# LS2D scripts, to calculate ERA5 tendencies
import sys; sys.path.append('/usr/people/stratum/meteo/models/LS2D/src')
from read_ERA5 import Read_ERA

pl.close('all')

def format_h_since(hours):
    return datetime.timedelta(hours=float(hours)) + datetime.datetime(1900, 1, 1)

def wind_to_components(speed, direction):
    u = -speed * np.sin(np.deg2rad(direction))
    v = -speed * np.cos(np.deg2rad(direction))
    return u,v

# ----------------------------
# Harmonie data and tendencies
# ----------------------------
f = xr.open_mfdataset(['20100228_06.nc','20100228_09.nc','20100228_12.nc'])


if(True):

    xlim  = [datetime.datetime(2010,2,28,9), datetime.datetime(2010,2,28,19)]
    xlimf = [datetime.datetime(2010,2,28,15), datetime.datetime(2010,2,28,16)]

    def format_ax(ax=None, interval=2):
        if ax is None:
            ax = pl.gca()

        hours = mdates.HourLocator(interval=interval)
        hours_fmt = mdates.DateFormatter('%H:%M')
        ax.xaxis.set_major_locator(hours)
        ax.xaxis.set_major_formatter(hours_fmt)

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


    # --------------------------------
    # Temperature wiggles
    # --------------------------------
    if (True):
        pl.figure()
        pl.subplot(221)
        pl.plot(f.time, f['T'][:,0], label='level 65')
        pl.plot(f.time, f['T'][:,10], label='level 55')
        pl.legend()
        pl.ylabel('T (K)')
        format_ax()

        pl.subplot(222)
        pl.plot(f.time, f['qv'][:,0], label='level 65')
        pl.plot(f.time, f['qv'][:,10], label='level 55')
        pl.legend()
        pl.ylabel('qv (kg kg-1)')
        format_ax()

        pl.subplot(223)
        pl.plot(f.time, f['u'][:,0], label='level 65')
        pl.plot(f.time, f['u'][:,10], label='level 55')
        pl.legend()
        pl.ylabel('u (m s-1')
        format_ax()

        pl.subplot(224)
        pl.plot(f.time, f['v'][:,0], label='level 65')
        pl.plot(f.time, f['v'][:,10], label='level 55')
        pl.legend()
        pl.ylabel('v (m s-1')
        format_ax()

        pl.tight_layout()


        pl.figure()
        pl.subplot(111)
        pl.plot(f.time, f['T'][:,0], label='level 65')
        pl.plot(f.time, f['T'][:,10], label='level 55')
        pl.legend()
        pl.ylabel('T (K)')
        format_ax(interval=1)

        pl.figure()
        pl.subplot(111)
        pl.plot(f.time, f['dtT_phy'][:,0]*3600., label='physics')
        pl.plot(f.time, f['dtT_dyn'][:,0]*3600., label='dynamics')
        pl.plot(f.time, f['dtT_tot'][:,0]*3600., label='total')
        pl.legend()
        pl.ylabel('dT/dt (K h-1)')
        format_ax(interval=1)


    # --------------------------------
    # u,v: Harmonie vs Cabauw vs ERA5
    # --------------------------------
    if (False):

        # ----
        # ERA5
        # ----
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
        e5.datetime = [format_h_since(h) for h in e5.time]

        # -------------------
        # Cabauw observations
        # -------------------
        pwd  = '/nobackup/users/stratum/Cabauw'
        cb_t = xr.open_dataset('{}/cesar_tower_meteo_lc1_t10_v1.0_201002.nc'.format(pwd))

        # Calculate wind components from speed/direction Cabauw
        cb_u, cb_v = wind_to_components(cb_t.F, cb_t.D)


        k = -2
        z = float(cb_t.z[k])
        print(z)

        # Interpolate Harmonie & ERA5 to Cabauw obs height
        u_h = interp(f['u'].values, f['zg'][0,:].values, z)
        u_e = interp(e5.u_mean, e5.z_mean[0,:], z)

        v_h = interp(f['v'].values, f['zg'][0,:].values, z)
        v_e = interp(e5.v_mean, e5.z_mean[0,:], z)

        ut_h = interp(f['dtu_dyn'].values, f['zg'][0,:].values, z)
        ut_e = interp(e5.dtu_advec+e5.dtu_coriolis, e5.z_mean[0,:], z)

        vt_h = interp(f['dtv_dyn'].values, f['zg'][0,:].values, z)
        vt_e = interp(e5.dtv_advec+e5.dtv_coriolis, e5.z_mean[0,:], z)

        ut_e_2 = interp(e5.dtu_advec_2+e5.dtu_coriolis_2, e5.z_mean[0,:], z)
        vt_e_2 = interp(e5.dtv_advec_2+e5.dtv_coriolis_2, e5.z_mean[0,:], z)


        pl.figure(figsize=(8,6))
        gs = gridspec.GridSpec(2, 2, height_ratios=[1,1])

        ax=pl.subplot(gs[0,0])
        pl.title('z={0:.0f} m'.format(z), loc='left')
        pl.plot(cb_t.time, cb_u[:,k], color='k',  label='Cabauw', marker='x', markersize=5, linewidth=0)
        pl.plot(f.time, u_h,           color='C1', label='Harmonie ',   linewidth=1.5)
        pl.plot(e5.datetime, u_e,      color='C2', label='ERA5',       marker='s', markersize=3)
        pl.ylabel('$u$ (m s$^{-1}$)')
        pl.legend()
        pl.xlim(xlim)
        pl.hlines(0, pl.xlim()[0], pl.xlim()[1], linestyles=':')
        format_ax()

        pl.subplot(gs[0,1])
        pl.plot(cb_t.time, cb_v[:,k], color='k', marker='x', markersize=5, linewidth=0)
        pl.plot(f.time, v_h,           color='C1', linewidth=1.5)
        pl.plot(e5.datetime, v_e,      color='C2', marker='s', markersize=3)
        pl.ylabel('$v$ (m s$^{-1}$)')
        pl.xlim(xlim)
        format_ax()

        pl.subplot(gs[1,0])
        pl.plot(f.time, ut_h*3600., color='C1', linewidth=1.5, label='dynamics')
        pl.plot(e5.datetime, ut_e*3600., color='C2', marker='s', markersize=3)
        pl.plot(e5.datetime, ut_e_2*3600., color='C3', marker='s', markersize=3)
        pl.ylabel('$\partial_t u$ (m s$^{-1}$ h$^{-1}$)')
        pl.xlabel('time (UTC)')
        pl.xlim(xlim)
        pl.hlines(0, pl.xlim()[0], pl.xlim()[1], linestyles=':')
        format_ax()

        pl.subplot(gs[1,1])
        pl.plot(f.time, vt_h*3600., color='C1', linewidth=1.5)
        pl.plot(e5.datetime, vt_e*3600., color='C2', marker='s', markersize=3)
        pl.plot(e5.datetime, vt_e_2*3600., color='C3', marker='s', markersize=3)
        pl.ylabel('$\partial_t v$ (m s$^{-1}$ h$^{-1}$)')
        pl.xlabel('time (UTC)')
        pl.xlim(xlim)
        pl.hlines(0, pl.xlim()[0], pl.xlim()[1], linestyles=':')
        format_ax()

        pl.tight_layout()

        # -------------------------------
        # Budget Harmonie tendencies
        # -------------------------------

        ut_dyn = interp(f['dtu_dyn'].values, f['zg'][0,:].values, z)
        ut_phy = interp(f['dtu_phy'].values, f['zg'][0,:].values, z)
        ut_tot = interp(f['dtu_tot'].values, f['zg'][0,:].values, z)

        vt_dyn = interp(f['dtv_dyn'].values, f['zg'][0,:].values, z)
        vt_phy = interp(f['dtv_phy'].values, f['zg'][0,:].values, z)
        vt_tot = interp(f['dtv_tot'].values, f['zg'][0,:].values, z)

        hours = mdates.MinuteLocator(interval=30)
        hours_fmt = mdates.DateFormatter('%H:%M')

        pl.figure(figsize=(3,6))

        ax=pl.subplot(211)
        pl.plot(f.time, ut_dyn*3600., color='C1', label='dyn', linewidth=1.5)
        pl.plot(f.time, ut_phy*3600., color='C2', label='phy', linewidth=1.5)
        pl.plot(f.time, ut_tot*3600., color='C3', label='tot', linewidth=1.5)
        pl.legend()
        pl.xlim(xlimf)
        pl.hlines(0, pl.xlim()[0], pl.xlim()[1], linestyles=':')
        pl.ylabel('$\partial_t u$ (m s$^{-1}$ h$^{-1}$)')
        pl.xlabel('time (UTC)')
        ax.xaxis.set_major_locator(hours)
        ax.xaxis.set_major_formatter(hours_fmt)


        ax=pl.subplot(212)
        pl.plot(f.time, vt_dyn*3600., color='C1', label='dyn', linewidth=1.5)
        pl.plot(f.time, vt_phy*3600., color='C2', label='phy', linewidth=1.5)
        pl.plot(f.time, vt_tot*3600., color='C3', label='tot', linewidth=1.5)
        pl.legend()
        pl.xlim(xlimf)
        pl.hlines(0, pl.xlim()[0], pl.xlim()[1], linestyles=':')
        pl.ylabel('$\partial_t v$ (m s$^{-1}$ h$^{-1}$)')
        pl.xlabel('time (UTC)')
        ax.xaxis.set_major_locator(hours)
        ax.xaxis.set_major_formatter(hours_fmt)

        pl.tight_layout()



