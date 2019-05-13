import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

import xarray as xr
import pandas as pd
import numpy as np

from datetime import datetime, timedelta
from scipy import interpolate

threehours = timedelta(hours=3)
oneday = timedelta(hours=24)

class Read_LES:
    def __init__(self, nc_path, start_date):
        print('Reading LES for {}'.format(start_date))

        date_str  = '{0:04d}{1:02d}{2:02d}'.format(start_date.year, start_date.month, start_date.day)

        # Domain mean statistics
        self.fp   = xr.open_dataset('{}/profiles_{}.nc'.format(nc_path, date_str))
        self.ft   = xr.open_dataset('{}/tmser_{}.nc'.format(nc_path, date_str))

        # Column sampled and averaged statistics
        self.fc   = xr.open_dataset('{}/column.i00097j00097_{}.nc'.format(nc_path, date_str))
        self.fc['Qnet'] = (('time'), -(self.fc['swd'][:,0] + self.fc['swu'][:,0] + self.fc['lwd'][:,0] + self.fc['lwu'][:,0]))

        # Calculate LWP
        dz  = self.fp['zm'][1:].values - self.fp['zm'][:-1].values
        dz  = np.append(dz, dz[-1])
        dims = self.fp['rhof'].shape
        self.fc['lwp'] = np.sum(self.fc['ql'] * np.broadcast_to(self.fp['rhof'], dims) * np.broadcast_to(dz, dims), axis=1)

        self.time = [start_date+timedelta(seconds=int(self.fp.time[t])) for t in range(self.fp.time.size)]


def read_all(start, end, LES_path):
    runs = []
    date  = start
    while date < end:
        l = Read_LES(LES_path, date)
        l.fp['rainrate'] /= (l.fp['rhof']*2.45e6)
        runs.append(l)
        date += timedelta(hours=24)
    return runs


def lim_and_line(vmin, vmax):
    pl.xlim(vmin, vmax)
    pl.ylim(vmin, vmax)
    pl.plot([vmin,vmax], [vmin,vmax], 'k:', linewidth=1)
    pl.plot([vmin,vmax], [0,0], 'k:', linewidth=1)
    pl.plot([0,0], [vmin,vmax], 'k:', linewidth=1)


def lim_and_line2(v1, v2):
    vmin = np.min((v1.min(), v2.min()))
    vmax = np.max((v1.max(), v2.max()))

    vmin = np.floor(vmin)
    vmax = np.ceil(vmax)

    ax=pl.gca()
    pl.plot([vmin,vmax], [vmin,vmax], 'k:', linewidth=1)
    pl.plot([vmin,vmax], [0,0], 'k:', linewidth=1)
    pl.plot([0,0], [vmin,vmax], 'k:', linewidth=1)
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)


def format_ax():
    ax = pl.gca()
    h24 = mdates.HourLocator(interval=48)
    fmt = mdates.DateFormatter('%d-%m')
    ax.xaxis.set_major_locator(h24)
    ax.xaxis.set_major_formatter(fmt)


def interp_z(array, heights, goal):
    """
    Interpolate `array` at heights `heights` to goal height `goal`
    """
    f = interpolate.interp1d(heights, array, fill_value='extrapolate')
    return f(goal)


def interp_zt(array, heights, goal):
    """
    Interpolate `array` at heights `heights` to goal height `goal`,
    when the input heights are time dependent (e.g. Harmonie, ERA5, ..)
    """
    nt  = array.shape[0]
    out = np.empty((nt,goal.size), dtype=np.float)
    for i in range(nt):
        if i%100 == 0:
            print(i,nt)
        out[i,:] = np.interp(goal, heights[i,:], array[i,:])
    return out


def calc_swd(lon, lat, hour, doy):
    """
    Calculate theoretical shortwave incoming radiation,
    from (time) index of a Pandas dataframe
    """
    lon    = -lon
    sda    = 0.409 * np.cos(2. * np.pi * (doy - 173.) / 365.)
    sinlea = np.sin(2. * np.pi * lat / 360.) * np.sin(sda) - \
             np.cos(2. * np.pi * lat / 360.) * np.cos(sda) * \
             np.cos(2. * np.pi * (hour*3600.) / 86400. - 2. * np.pi * lon / 360.)
    sinlea = np.maximum(sinlea, 1e-9)
    Tr     = (0.6 + 0.2 * sinlea)
    swin   = 1368. * Tr * sinlea

    return swin


def absval(a, b):
    return (a**2 + b**2)**0.5

def calc_rmse(a, b):
    return np.sqrt(np.mean((a-b)**2))

if __name__ == '__main__':
    pl.close('all')

    # Period to read/plot/..
    start = datetime(year=2016, month=8, day=4,  hour=0)
    end   = datetime(year=2016, month=8, day=18, hour=0)

    # Local file paths
    # ---- Macbook ----
    """
    LES_path  = '/Users/bart/meteo/data/KNMI_testbed/cabauw_20160804_20160818_lambda'
    CB_path   = '/Users/bart/meteo/observations/Cabauw'
    HM_path   = '/Users/bart/meteo/data/Harmonie_LES_forcing'
    E5_path   = '/Users/bart/meteo/data//LS2D/cabauw/ERA5'
    """

    # ---- KNMI Desktop ----
    LES_path  = '/nobackup/users/stratum/KNMI_testbed/cases/cabauw_aug2018'
    CB_path   = '/nobackup/users/stratum/Cabauw'
    HM_path   = '/nobackup/users/stratum/DOWA/LES_forcing'
    E5_path   = '/nobackup/users/stratum/ERA5/LS2D/cabauw/ERA5'


    # Read the LES data
    # -----------------
    if 'runs' not in locals():
        runs = read_all(start, end, LES_path)


    # Read HARMONIE data
    # ------------------
    if 'hm' not in locals():
        print('Reading HARMONIE')

        iloc = 7 #+12     # 7=2.5x2.5, 19=10x10km Cabauw
        files = []
        t = start
        while t <= end:
            files.append('{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/LES_forcing_{1:04d}{2:02d}{3:02d}{4:02d}.nc'\
                    .format(HM_path, t.year, t.month, t.day, t.hour))
            t += threehours
        hm = xr.open_mfdataset(files)
        hm = hm.loc[{'domain': iloc}]


    # Read ERA5 data
    # --------------
    if 'e5' not in locals():
        print('Reading ERA5')

        files = []
        t = start
        while t< end:
            files.append('{0:}/{1:04d}/{2:02d}/{3:02d}/model_an.nc'\
                    .format(E5_path, t.year, t.month, t.day))
            t += oneday
        e5 = xr.open_mfdataset(files)
        e5 = e5.sel(longitude=4.91, latitude=51.97, method='nearest')


    # Read Cabauw observations
    # ------------------------
    if 'cb_sf' not in locals():
        print('Reading Cabauw')

        def open_and_sel(name, start, end):
            f = xr.open_mfdataset(name, drop_variables=['valid_dates'])
            f = f.sel(time=slice(start, end))
            return f

        cb_sm  = open_and_sel('{0}/cesar_surface_mete*{1:04d}{2:02d}.nc'.format(CB_path, start.year, start.month), start, end)
        cb_sf  = open_and_sel('{0}/cesar_surface_flux*{1:04d}{2:02d}.nc'.format(CB_path, start.year, start.month), start, end)
        cb_sr  = open_and_sel('{0}/cesar_surface_radi*{1:04d}{2:02d}.nc'.format(CB_path, start.year, start.month), start, end)
        cb_ns  = open_and_sel('{0}/cesar_nubiscope*{1:04d}{2:02d}.nc'   .format(CB_path, start.year, start.month), start, end)
        cb_lwc = open_and_sel('{0}/{1:04d}{2:02d}*_cabauw_lwc-scaled-adiabatic.nc'.format(CB_path, start.year, start.month), start, end)
        cb_tm  = open_and_sel('{0}/cesar_tower_meteo*{1:04d}{2:02d}.nc'.format(CB_path, start.year, start.month), start, end)

        # Cabauw soil heat flux calculated as residual term
        cb_G_res = cb_sf.QN - (cb_sf.H + cb_sf.LE)

        # Indices of heights
        cb_k = {}
        for i in range(cb_tm.z.size):
            cb_k[int(cb_tm.z[i])] = i


    # Sync times in a Pandas dataframe
    # ------------------
    if 'dfs' not in locals():

        # Read selected LES variables in Pandas DataFrame's
        dfs = []
        for r in runs:
            exner = (r.fp['presh'] / 1e5)**(287.05/1004.)
            T     = (r.fc['thl'] + (2.45e6 / 1004.) * r.fc['ql']) * exner
            r.fc['T'] = T

            data = { 'LE_LES':   r.fc['LE'],
                     'H_LES':    r.fc['H'],
                     'G_LES':    r.fc['G'],
                     'Qn_LES':   r.fc['Qnet'],
                     'swd_LES':  r.fc['swd'][:,0],
                     'swu_LES':  r.fc['swu'][:,0],
                     'lwd_LES':  r.fc['lwd'][:,0],
                     'lwu_LES':  r.fc['lwu'][:,0],
                     #'cc_LES':   r.ft['cfrac'],     # Was ist loss mit cfrac?
                     'rr_LES':   r.fc['rainrate'][:,0]*r.fp.rhobh[:,0],
                     'lwp_LES':  r.fc['lwp'],
                     'U010_LES': absval(interp_z(r.fc['u'], r.fc['zt'], 10),  interp_z(r.fc['v'], r.fc['zt'], 10)),
                     'U200_LES': absval(interp_z(r.fc['u'], r.fc['zt'], 200), interp_z(r.fc['v'], r.fc['zt'], 200)),
                     'T010_LES': interp_z(r.fc['T'],  r.fc['zt'], 10),
                     'T200_LES': interp_z(r.fc['T'],  r.fc['zt'], 200),
                     'q010_LES': interp_z(r.fc['qt'], r.fc['zt'], 10),
                     'q200_LES': interp_z(r.fc['qt'], r.fc['zt'], 200)
                   }

            dfs.append( pd.DataFrame(data, index=r.time) )
        df_LES = pd.concat(dfs)

        # Put Cabauw observations in DataFrame
        data = { 'LE_CB':   cb_sf['LE'],
                 'LE2_CB':  cb_sf['LE2'],
                 'H_CB':    cb_sf['H'],
                 'G_CB':    cb_sf['G0'],
                 'Qn_CB':   cb_sf['QN'],
                 'swd_CB':  cb_sr['SWD'],
                 'swu_CB':  cb_sr['SWU'],
                 'lwd_CB':  cb_sr['LWD'],
                 'lwu_CB':  cb_sr['LWU'],
                 'rr_CB':   cb_sm['RAIN'],
                 'U010_CB': cb_tm['F'] [:,cb_k[10] ],
                 'U200_CB': cb_tm['F'] [:,cb_k[200]],
                 'T010_CB': cb_tm['TA'][:,cb_k[10] ],
                 'T200_CB': cb_tm['TA'][:,cb_k[200]],
                 'q010_CB': cb_tm['Q'] [:,cb_k[10] ],
                 'q200_CB': cb_tm['Q'] [:,cb_k[200]]
               }

        data2 = {'cc_CB':  cb_ns.cldcover_total/100.}

        df_CB = pd.DataFrame(data, index=cb_sf.time)
        df_CB.index = df_CB.index.round('1min')


        df_CB2 = pd.DataFrame(data2, index=cb_ns.time)
        df_CB2.index = df_CB2.index.round('1min')

        df_CB = pd.concat([df_CB, df_CB2], axis=1)
        df_CB.dropna(inplace=True)

        # Merge DataFrame's
        df = pd.concat([df_LES, df_CB], axis=1)
        df.dropna(inplace=True)

        # Theoretical shortwave incoming radiation
        hour = df.index.hour + df.index.minute/60.
        doy  = df.index.dayofyear
        df['swd_theory'] = calc_swd(4.9, 51.97, hour, doy)
        df['is_night'] = df['swd_theory'] < 0.1


    # Plot settings
    c1 = '#4d4d4d'   # Green
    c2 = '#4daf4a'    # Blue

    c_cb  = '#4daf4a'   # Green
    c_cb2 = '#377eb8'   # Blue
    c_da  = '#4d4d4d'   # Gray
    c_da2 = '#b2182b'   # DarkRed


    def scatter_stat(obs, model, label, xlabel, ylabel, night_mask, ax=None):
        # Statistics
        rmse_all = calc_rmse(obs, model)
        diff_all = (model-obs).mean()

        rmse_night = calc_rmse(obs[night_mask], model[night_mask])
        diff_night = (model[night_mask]-obs[night_mask]).mean()

        rmse_day = calc_rmse(obs[~night_mask], model[~night_mask])
        diff_day = (model[~night_mask]-obs[~night_mask]).mean()

        print('---------------------------')
        print(label)
        print('ALL:   RMSE = {0:6.1f}, DIFF = {1:6.1f}'.format(rmse_all, diff_all))
        print('DAY:   RMSE = {0:6.1f}, DIFF = {1:6.1f}'.format(rmse_day, diff_day))
        print('NIGHT: RMSE = {0:6.1f}, DIFF = {1:6.1f}'.format(rmse_night, diff_night))


        pl.scatter(obs, model, s=1, color=c2)
        lim_and_line2(obs, model)
        pl.xlabel(xlabel)
        pl.ylabel(ylabel)



    if False:
        # --------------
        # LES vs HARMONIE vs ERA5
        # --------------

        # ---- LES ----
        dfs = []
        for r in runs:
            data = { 'U010_LES': absval(interp_z(r.fp['u'], r.fp['zt'], 10),  interp_z(r.fp['v'], r.fp['zt'], 10)),
                     'U100_LES': absval(interp_z(r.fp['u'], r.fp['zt'], 100), interp_z(r.fp['v'], r.fp['zt'], 100)),
                     'U200_LES': absval(interp_z(r.fp['u'], r.fp['zt'], 200), interp_z(r.fp['v'], r.fp['zt'], 200)),
                     'T010_LES': interp_z(r.fp['T'],  r.fp['zt'], 10),
                     'T100_LES': interp_z(r.fp['T'],  r.fp['zt'], 100),
                     'T200_LES': interp_z(r.fp['T'],  r.fp['zt'], 200),
                     'q010_LES': interp_z(r.fp['qt'], r.fp['zt'], 10),
                     'q100_LES': interp_z(r.fp['qt'], r.fp['zt'], 100),
                     'q200_LES': interp_z(r.fp['qt'], r.fp['zt'], 200) }
            dfs.append( pd.DataFrame(data, index=r.time) )
        df_LES = pd.concat(dfs)

        # ---- HARMONIE ----
        if 'z' not in locals():
            z = np.array((10,100,200))
            u = interp_zt(hm['u'], hm['z'], z)
            v = interp_zt(hm['v'], hm['z'], z)
            T = interp_zt(hm['T'], hm['z'], z)
            q = interp_zt(hm['q'], hm['z'], z)

        data = { 'U010_HAM': absval(u[:,0], v[:,0]),
                 'U100_HAM': absval(u[:,1], v[:,1]),
                 'U200_HAM': absval(u[:,2], v[:,2]),
                 'T010_HAM': T[:,0],
                 'T100_HAM': T[:,1],
                 'T200_HAM': T[:,2],
                 'q010_HAM': q[:,0],
                 'q100_HAM': q[:,1],
                 'q200_HAM': q[:,2] }

        df_HAM = pd.DataFrame(data, index=hm.time)
        df_HAM.index = df_HAM.index.round('1min')

        # ---- ERA5 ----
        if 'z' not in locals():
            z = np.array((10,100,200))
            u = interp_zt(hm['u'], hm['z'], z)
            v = interp_zt(hm['v'], hm['z'], z)
            T = interp_zt(hm['T'], hm['z'], z)
            q = interp_zt(hm['q'], hm['z'], z)

        data = { 'U010_HAM': absval(u[:,0], v[:,0]),
                 'U100_HAM': absval(u[:,1], v[:,1]),
                 'U200_HAM': absval(u[:,2], v[:,2]),
                 'T010_HAM': T[:,0],
                 'T100_HAM': T[:,1],
                 'T200_HAM': T[:,2],
                 'q010_HAM': q[:,0],
                 'q100_HAM': q[:,1],
                 'q200_HAM': q[:,2] }

        df_HAM = pd.DataFrame(data, index=hm.time)
        df_HAM.index = df_HAM.index.round('1min')




    if False:
        # --------------
        # Atmospheric variables (wind, temp, moisture)
        # --------------

        # ---- Wind ----
        pl.figure(figsize=(10,5))
        gs = gridspec.GridSpec(2, 2, width_ratios=[3.8,1])

        ax=pl.subplot(gs[0,0])
        pl.plot(cb_tm.time, cb_tm['F'][:,cb_k[10]], 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, absval(interp_z(r.fc['u'], r.fc['zt'], 10),  interp_z(r.fc['v'], r.fc['zt'], 10)), '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'U$_\mathrm{10m}$ (m s$^{-1}$')
        format_ax()

        ax=pl.subplot(gs[1,0])
        pl.plot(cb_tm.time, cb_tm['F'][:,cb_k[200]], 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, absval(interp_z(r.fc['u'], r.fc['zt'], 200),  interp_z(r.fc['v'], r.fc['zt'], 200)), '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'U$_\mathrm{200m}$ (m s$^{-1}$')
        format_ax()

        pl.subplot(gs[0,1])
        scatter_stat(df['U010_CB'], df['U010_LES'], 'U_10m (m/s)', r'OBS (m s$^{-1}$)', r'LES (m s$^{-1}$)', df['is_night'])

        pl.subplot(gs[1,1])
        scatter_stat(df['U200_CB'], df['U200_LES'], 'U_200m (m/s)', r'OBS (m s$^{-1}$)', r'LES (m s$^{-1}$)', df['is_night'])

        pl.tight_layout()
        pl.savefig('wind_tser_scatter.pdf')


        # ---- Temperature ----
        pl.figure(figsize=(10,5))
        gs = gridspec.GridSpec(2, 2, width_ratios=[3.8,1])

        ax=pl.subplot(gs[0,0])
        pl.plot(cb_tm.time, cb_tm['TA'][:,cb_k[10]], 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, interp_z(r.fc['T'], r.fc['zt'], 10), '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'T$_\mathrm{10m}$ (K)')
        format_ax()

        ax=pl.subplot(gs[1,0])
        pl.plot(cb_tm.time, cb_tm['TA'][:,cb_k[200]], 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, interp_z(r.fc['T'], r.fc['zt'], 200), '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'T$_\mathrm{200m}$ (K)')
        format_ax()

        pl.subplot(gs[0,1])
        scatter_stat(df['T010_CB'], df['T010_LES'], 'T_10m (K)', r'OBS (K)', r'LES (K)', df['is_night'])

        pl.subplot(gs[1,1])
        scatter_stat(df['T200_CB'], df['T200_LES'], 'T_10m (K)', r'OBS (K)', r'LES (K)', df['is_night'])

        pl.tight_layout()
        pl.savefig('temperature_tser_scatter.pdf')


        # ---- Specific humidity ----
        pl.figure(figsize=(10,5))
        gs = gridspec.GridSpec(2, 2, width_ratios=[3.8,1])

        ax=pl.subplot(gs[0,0])
        pl.plot(cb_tm.time, cb_tm['Q'][:,cb_k[10]], 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, interp_z(r.fc['qt']*1e3, r.fc['zt'], 10), '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'q$_\mathrm{10m}$ (g kg$^{-1}$)')
        format_ax()

        ax=pl.subplot(gs[1,0])
        pl.plot(cb_tm.time, cb_tm['Q'][:,cb_k[200]], 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, interp_z(r.fc['qt']*1e3, r.fc['zt'], 200), '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'q$_\mathrm{200m}$ (g kg$^{-1}$)')
        format_ax()

        pl.subplot(gs[0,1])
        scatter_stat(df['q010_CB'], df['q010_LES']*1000, 'q_10m (g/kg)', r'OBS (g kg$^{-1}$)', r'LES (g kg$^{-1}$)', df['is_night'])

        pl.subplot(gs[1,1])
        scatter_stat(df['q200_CB'], df['q200_LES']*1000, 'q_200m (g/kg)', r'OBS (g kg$^{-1}$)', r'LES (g kg$^{-1}$)', df['is_night'])

        pl.tight_layout()
        pl.savefig('spechum_tser_scatter.pdf')


    if False:
        # --------------
        # Surface fluxes
        # --------------
        pl.figure(figsize=(10,8))
        gs = gridspec.GridSpec(4, 2, width_ratios=[3.8,1])

        # Time series
        # Time series
        ax=pl.subplot(gs[0,0])
        pl.plot(cb_sf.time.values, cb_sf.H, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.fc.H, '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'H (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[1,0], sharex=ax)
        pl.plot(cb_sr.time.values,  cb_sf.LE, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.fc.LE, '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'LE (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[2,0], sharex=ax)
        pl.plot(cb_sr.time.values, cb_sf.G0, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.fc.G, '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'G (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[3,0], sharex=ax)
        pl.plot(cb_sr.time.values, cb_sf.QN, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.fc.Qnet, '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'Q$_\mathrm{net}$ (W m$^{-2}$')
        format_ax()

        # Scatter plots
        pl.subplot(gs[0,1])
        scatter_stat(df['H_CB'], df['H_LES'], 'H (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', df['is_night'])

        pl.subplot(gs[1,1])
        scatter_stat(df['LE_CB'], df['LE_LES'], 'LE (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', df['is_night'])

        pl.subplot(gs[2,1])
        scatter_stat(df['G_CB'], df['G_LES'], 'G (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', df['is_night'])

        pl.subplot(gs[3,1])
        scatter_stat(df['Qn_CB'], df['Qn_LES'], 'Qnet (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', df['is_night'])

        pl.tight_layout()
        pl.savefig('surface_flux_tser_scatter.pdf')


    if False:
        # --------------
        # Surface radiation
        # --------------

        pl.figure(figsize=(10,8))

        gs = gridspec.GridSpec(4, 2, width_ratios=[3.8,1])

        ax=pl.subplot(gs[0,0])
        pl.plot(cb_sr.time.values, -cb_sr.SWD, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.fc.swd[:,0], '-', color=c1)   # Column
        pl.xlim(start, end)
        pl.ylabel(r'$SW_\mathrm{down}$ (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[1,0], sharex=ax)
        pl.plot(cb_sr.time.values,  cb_sr.SWU, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.fc.swu[:,0], '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'$SW_\mathrm{up}$ (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[2,0], sharex=ax)
        pl.plot(cb_sr.time.values, -cb_sr.LWD, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.fc.lwd[:,0], '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'$LW_\mathrm{down}$ (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[3,0], sharex=ax)
        pl.plot(cb_sr.time.values, cb_sr.LWU, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.fc.lwu[:,0], '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'$LW_\mathrm{up}$ (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[0,1])
        scatter_stat(-df['swd_CB'], df['swd_LES'], 'SWdown (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', df['is_night'])

        pl.subplot(gs[1,1])
        scatter_stat(df['swu_CB'], df['swu_LES'], 'SWup (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', df['is_night'])

        pl.subplot(gs[2,1])
        scatter_stat(-df['lwd_CB'], df['lwd_LES'], 'LWdown (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', df['is_night'])

        pl.subplot(gs[3,1])
        scatter_stat(df['lwu_CB'], df['lwu_LES'], 'LWup (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', df['is_night'])

        pl.tight_layout()
        pl.savefig('radiation_tser_scatter.pdf')


    if False:
        # --------------
        # Surface meteo
        # --------------

        pl.figure(figsize=(10,4.8))
        gs = gridspec.GridSpec(2, 2, width_ratios=[3.8,1])

        ax=pl.subplot(gs[0,0])
        for i,r in enumerate(runs):
            pl.plot(r.time, r.ft.cfrac, '-', color=c1)
        pl.plot(cb_ns.time.values, cb_ns.cldcover_total/100., 'o', mfc=c2, mec=c2, ms=2)
        pl.plot([start, end], [0,0], 'k:')
        pl.xlim(start, end)
        pl.ylim(0,1)
        pl.ylabel(r'$cc$ (-)')

        pl.subplot(gs[1,0], sharex=ax)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.fc.rainrate[:,0]*r.fp.rhobh[:,0]*3600, '-', color=c1)
        pl.plot(cb_sm.time.values, cb_sm.RAIN*6, 'o', mfc=c2, mec=c2, ms=2)
        pl.plot([start, end], [0,0], 'k:')
        pl.xlim(start, end)
        pl.ylim(0,1)
        pl.ylabel(r'$rr$ (mm h$^{-1}$)')

        ax=pl.subplot(gs[0,1])
        pl.scatter(df['cc_CB'], df['cc_LES'], s=1, color=c2)
        lim_and_line2(df['cc_CB'], df['cc_LES'])
        pl.xlabel(r'OBS (-)')
        pl.ylabel(r'LES (-)')

        ax=pl.subplot(gs[1,1])
        pl.scatter(df['rr_CB']*6, df['rr_LES']*3600, s=1, color=c2)
        lim_and_line(0,2)
        pl.xlabel(r'OBS (mm h$^{-1}$)')
        pl.ylabel(r'LES (mm h$^{-1}$)')

        pl.tight_layout()
        pl.savefig('clouds_rain_tser_scatter.pdf')


    if False:
        # --------------
        # Check SEB closure observations
        # --------------

        pl.figure()
        pl.scatter(cb_sf['QN'], cb_sf['H']+cb_sf['LE']+cb_sf['G0'], s=1, color='r')
        lim_and_line(-100,600)
        pl.plot([-50,0], [-100,0], 'k:')
        pl.xlabel(r'Q$_\mathrm{net}$ (W m$^{-2}$)')
        pl.ylabel(r'LE+H+G (W m$^{-2}$)')


