import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

import xarray as xr
import pandas as pd
import numpy as np

from datetime import datetime, timedelta
from scipy import interpolate, stats

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

# Enable LaTeX plotting
from matplotlib import rc
rc('text', usetex=True)
rc('font', size=13)
rc('legend', fontsize=11)
rc('text.latex', preamble=r'\usepackage{sansmathfonts}')

# Custom module (from LS2D)
from IFS_tools import IFS_tools

threehours = timedelta(hours=3)
oneday = timedelta(hours=24)

class Read_LES:
    def __init__(self, nc_path, column_i, column_j, start_date, start_hour):
        print('Reading LES for {}'.format(start_date))

        date_str  = '{0:04d}{1:02d}{2:02d}'.format(start_date.year, start_date.month, start_date.day)

        # Domain mean statistics
        self.fp   = xr.open_dataset('{}/profiles_{}.nc'.format(nc_path, date_str))
        self.ft   = xr.open_dataset('{}/tmser_{}.nc'.format(nc_path, date_str))

        exner = (self.fp['presh'] / 1e5)**(287.05/1004.)

        # Column sampled and averaged statistics
        self.fc   = xr.open_dataset('{0}/column.i{1:05d}j{2:05d}_{3}.nc'.format(nc_path, column_i, column_j, date_str))
        self.fc['Qnet'] = (('time'), -(self.fc['swd'][:,0] + self.fc['swu'][:,0] + self.fc['lwd'][:,0] + self.fc['lwu'][:,0]))
        self.fc['T'] = (self.fc['thl'] + (2.45e6 / 1004.) * self.fc['ql'])*exner

        # Calculate LWP
        dz  = self.fp['zm'][1:].values - self.fp['zm'][:-1].values
        dz  = np.append(dz, dz[-1])
        dims = self.fp['rhof'].shape
        self.fc['lwp'] = np.sum(self.fc['ql'] * np.broadcast_to(self.fp['rhof'], dims) * np.broadcast_to(dz, dims), axis=1)

        self.time = [start_date+timedelta(hours=start_hour)+timedelta(seconds=int(self.fp.time[t])) for t in range(self.fp.time.size)]

        # Cloud fraction in the time series is set at fill_value if cc=0 (why?) -> change NaNs to zeros
        self.ft['cfrac'] = self.ft['cfrac'].fillna(0)


def read_all(start, end, start_hour, LES_path, column_i, column_j):
    runs = []
    date  = start
    while date < end:
        l = Read_LES(LES_path, column_i, column_j, date, start_hour)
        l.fp['rainrate'] /= (l.fp['rhof']*2.45e6)
        runs.append(l)
        date += timedelta(hours=24)
    return runs


def absval(a, b):
    return (a**2 + b**2)**0.5


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

    print('Coffee time :-)')
    nt  = array.shape[0]
    out = np.empty((nt,goal.size), dtype=np.float)
    for i in range(nt):
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


def lim_and_line(vmin, vmax):
    """
    Set the x and y axis limits to `vmin, vmax`,
    and add 1:1  and x=0 and y=0 lines
    """

    pl.xlim(vmin, vmax)
    pl.ylim(vmin, vmax)
    pl.plot([vmin,vmax], [vmin,vmax], 'k:', linewidth=1)
    pl.plot([vmin,vmax], [0,0], 'k:', linewidth=1)
    pl.plot([0,0], [vmin,vmax], 'k:', linewidth=1)


def lim_and_line2(v1, v2, round_lims=True):
    """
    Automaticcaly set the x and y axis limits,
    and add 1:1  and x=0 and y=0 lines
    """

    vmin = np.min((v1.min(), v2.min()))
    vmax = np.max((v1.max(), v2.max()))

    if round_lims:
        vmin = np.floor(vmin)
        vmax = np.ceil(vmax)

    ax=pl.gca()
    pl.plot([vmin,vmax], [vmin,vmax], 'k:', linewidth=1)
    pl.plot([vmin,vmax], [0,0], 'k:', linewidth=1)
    pl.plot([0,0], [vmin,vmax], 'k:', linewidth=1)

    # Aarghh
    for i in range(2):
        ticks = ax.get_yticks()
        ax.set_xticks(ticks)

        ax.set_xlim(vmin, vmax)
        ax.set_ylim(vmin, vmax)




def format_ax():
    """
    Format date/time x-axis
    """

    ax = pl.gca()
    h24 = mdates.HourLocator(interval=24)
    h48 = mdates.HourLocator(interval=48)
    fmt = mdates.DateFormatter('%d-%m')
    ax.xaxis.set_minor_locator(h24)
    ax.xaxis.set_major_locator(h48)
    ax.xaxis.set_major_formatter(fmt)
    ax.grid(which='both', axis='x')


class Stats:
    """
    Simple statistics class (RMSE, mean diff, ...)
    """
    def __init__(self, obs, model):

        self.mse   = np.mean((model-obs)**2)
        self.rmse  = np.sqrt(self.mse)
        self.diff  = (model-obs).mean()
        self.slope, self.intercept, self.rvalue, self.pvalue, self.stderr = stats.linregress(obs, model)


def scatter_stat(obs, model, label, xlabel, ylabel, rasterized=False):
    """
    Scatter plot of `model` vs `obs`
    """

    pl.scatter(obs, model, s=1, color=c2, rasterized=rasterized)
    lim_and_line2(obs, model)
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)



def pretty_align(label, names, values):
    """
    Align legend items in rows at `=` character
    """

    str_out = '{0:^12s}\n'.format(label)
    for name, value in zip(names, values):
        str_out += '{0:<4s} = {1:<+7.2f}\n'.format(name, value)
    str_out = str_out[:-2]  # Trim the last '\n'
    return str_out


def scatter_stat2(obs, model, label, xlabel, ylabel, night_mask, c_day, c_nig, xlim=None, ylim=None):
    """
    Complex version of `scatter_stat` which calculates / adds statistics,
    calculated over full run, and day and night periods
    """

    # Get current axes and figure objects
    fig = pl.gcf()
    ax  = pl.gca()

    # Statistics (full day, day and night)
    full  = Stats(obs, model)
    day   = Stats(obs[~night_mask], model[~night_mask])
    night = Stats(obs[ night_mask], model[ night_mask])

    # Scatter day and night in different colors
    p1=pl.scatter(obs[~night_mask], model[~night_mask], s=1, color=c_day)
    p2=pl.scatter(obs[ night_mask], model[ night_mask], s=1, color=c_nig)

    lim_and_line2(obs, model)
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)

    # Add statistics to legend
    label1 = pretty_align('Day',   ['RMSE', 'diff'], [day.rmse, day.diff])
    label2 = pretty_align('Night', ['RMSE', 'diff'], [night.rmse, night.diff])

    l1 = pl.legend([p1], [label1], loc=2, handlelength=0, handletextpad=0, prop={'family': 'monospace', 'size': 8})
    l2 = pl.legend([p2], [label2], loc=4, handlelength=0, handletextpad=0, prop={'family': 'monospace', 'size': 8})
    fig.add_artist(l1)

    # Text color of legend
    for text in l1.get_texts():
        text.set_color(c_day)
    for text in l2.get_texts():
        text.set_color(c_nig)

    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    return ax.get_xlim(), ax.get_ylim()




if __name__ == '__main__':
    pl.close('all')
    pl.ion()

    # Period to read/plot/..
    start = datetime(year=2016, month=8, day=4,  hour=0)
    end   = datetime(year=2016, month=8, day=18, hour=0)
    start_hour = 0

    # ---- Macbook ----
    LES_path1  = '/Users/bart/meteo/data/KNMI_testbed/cabauw_20160804_20160818_ref_lr'
    LES_path2  = '/Users/bart/meteo/data/KNMI_testbed/cabauw_20160804_20160818_LS2D_lr'

    CB_path   = '/Users/bart/meteo/observations/Cabauw'
    HM_path   = '/Users/bart/meteo/data/Harmonie_LES_forcing'
    E5_path   = '/Users/bart/meteo/data//LS2D/cabauw/ERA5'

    fig_path = '/Users/bart/meteo/KNMI_git/DOWA/reports/LES_downscaling/figs/'

    # ---- KNMI Desktop ----
    """
    LES_path  = '/nobackup/users/stratum/KNMI_testbed/cases/{}_{}'.format(base, case)
    CB_path   = '/nobackup/users/stratum/Cabauw'
    HM_path   = '/nobackup/users/stratum/DOWA/LES_forcing'
    E5_path   = '/nobackup/users/stratum/ERA5/LS2D/cabauw/ERA5'
    """

    # Plot settings
    c1 =  'k'      # Green
    c2 = '#3579c1'      # Blue

    c_cb  = '#4daf4a'   # Green
    c_cb2 = '#377eb8'   # Blue
    c_da  = '#4d4d4d'   # Gray
    c_da2 = '#b2182b'   # DarkRed

    c_day = '#b2182b'   # DarkRed
    c_nig = c2   # Green

    lw = 1.0

    # --------------------------------
    #
    # Read all the data
    #
    # --------------------------------

    #
    # Read the LES data
    #
    if 'runs1' not in locals():
        runs1 = read_all(start, end, start_hour, LES_path1, 97, 97)  # Harmonie -> LES
        #runs1 = read_all(start, end, start_hour, LES_path1, 161, 161)  # Harmonie -> LES
        runs2 = read_all(start, end, start_hour, LES_path2, 97, 97)  # ERA5 -> LES

    #
    # Read HARMONIE data
    #
    if 'hm' not in locals():
        print('Reading HARMONIE')

        iloc = 7 #+12     # 7=2.5x2.5, 19=10x10km Cabauw
        files = []
        t = start
        while t <= end:
            files.append('{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/LES_forcing_{1:04d}{2:02d}{3:02d}{4:02d}.nc'\
                    .format(HM_path, t.year, t.month, t.day, t.hour))
            t += threehours
        hm = xr.open_mfdataset(files, combine='by_coords')
        hm = hm.loc[{'domain': iloc}]

    #
    # Read ERA5 data
    #
    if 'e5' not in locals():
        print('Reading ERA5')

        files = []
        t = start
        while t < end:
            files.append('{0:}/{1:04d}/{2:02d}/{3:02d}/model_an.nc'\
                    .format(E5_path, t.year, t.month, t.day))
            t += oneday
        e5 = xr.open_mfdataset(files, combine='by_coords')
        e5 = e5.sel(longitude=4.91, latitude=51.97, method='nearest')

        # Calculate heights of model levels
        e5['zf'] = (('time','level'), np.zeros((e5.dims['time'], e5.dims['level'])))
        for t in range(e5.dims['time']):
            Tv = IFS_tools.calc_virtual_temp(e5['t'][t,:].values,    e5['q'][t,:].values,
                                              e5['clwc'][t,:].values, e5['ciwc'][t,:].values)
            ps = np.exp(e5['lnsp'][t,0].values)
            ph = IFS_tools.calc_half_level_pressure(ps)
            z  = IFS_tools.calc_full_level_Zg(ph, Tv)
            e5['zf'][t,:] = z[::-1]

    #
    # Read Cabauw observations
    #
    if 'cb_sf' not in locals():
        print('Reading Cabauw')

        def open_and_sel(name, start, end):
            f = xr.open_mfdataset(name, combine='by_coords', drop_variables=['valid_dates'])
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


    # --------------------------------
    #
    # Comparison LES - Cabauw
    #
    # --------------------------------
    if True:

        #
        # Sync observastion and LES times in Pandas dataframe
        #
        if 'df' not in locals():

            # Read selected LES variables in Pandas DataFrame's
            #
            # 1. HARMONIE forced LES
            #
            dfs = []
            for r in runs1:
                exner = (r.fp['presh'] / 1e5)**(287.05/1004.)
                T     = (r.fc['thl'] + (2.45e6 / 1004.) * r.fc['ql']) * exner
                r.fc['T'] = T

                data = { 'LE_LES1':   r.fc['LE'],
                         'H_LES1':    r.fc['H'],
                         'G_LES1':    r.fc['G'],
                         'Qn_LES1':   r.fc['Qnet'],
                         'swd_LES1':  r.fc['swd'][:,0],
                         'swu_LES1':  r.fc['swu'][:,0],
                         'lwd_LES1':  r.fc['lwd'][:,0],
                         'lwu_LES1':  r.fc['lwu'][:,0],
                         'cc_LES1':   r.ft['cfrac'],
                         'rr_LES1':   r.fc['rainrate'][:,0]*r.fp.rhobh[:,0],
                         'lwp_LES1':  r.fc['lwp'],
                         'U010_LES1': absval(interp_z(r.fc['u'], r.fc['zt'], 10),  interp_z(r.fc['v'], r.fc['zt'], 10)),
                         'U200_LES1': absval(interp_z(r.fc['u'], r.fc['zt'], 200), interp_z(r.fc['v'], r.fc['zt'], 200)),
                         'T010_LES1': interp_z(r.fc['T'],  r.fc['zt'], 10),
                         'T200_LES1': interp_z(r.fc['T'],  r.fc['zt'], 200),
                         'q010_LES1': interp_z(r.fc['qt'], r.fc['zt'], 10),
                         'q200_LES1': interp_z(r.fc['qt'], r.fc['zt'], 200)
                       }

                dfs.append( pd.DataFrame(data, index=r.time) )
            df_LES1 = pd.concat(dfs)

            #
            # 2. EAR5 forced LES
            #
            dfs = []
            for r in runs2:
                exner = (r.fp['presh'] / 1e5)**(287.05/1004.)
                T     = (r.fc['thl'] + (2.45e6 / 1004.) * r.fc['ql']) * exner
                r.fc['T'] = T

                t0=1    # Remove first output time LES -> bug RRTMG

                data = { 'LE_LES2':   r.fc['LE'][t0:],
                         'H_LES2':    r.fc['H'][t0:],
                         'G_LES2':    r.fc['G'][t0:],
                         'Qn_LES2':   r.fc['Qnet'][t0:],
                         'swd_LES2':  r.fc['swd'][t0:,0],
                         'swu_LES2':  r.fc['swu'][t0:,0],
                         'lwd_LES2':  r.fc['lwd'][t0:,0],
                         'lwu_LES2':  r.fc['lwu'][t0:,0],
                         'cc_LES2':   r.ft['cfrac'][t0:],
                         'rr_LES2':   r.fc['rainrate'][t0:,0]*r.fp.rhobh[t0:,0],
                         'lwp_LES2':  r.fc['lwp'][t0:],
                         'U010_LES2': absval(interp_z(r.fc['u'][t0:], r.fc['zt'], 10),  interp_z(r.fc['v'][t0:], r.fc['zt'], 10)),
                         'U200_LES2': absval(interp_z(r.fc['u'][t0:], r.fc['zt'], 200), interp_z(r.fc['v'][t0:], r.fc['zt'], 200)),
                         'T010_LES2': interp_z(r.fc['T'][t0:],  r.fc['zt'], 10),
                         'T200_LES2': interp_z(r.fc['T'][t0:],  r.fc['zt'], 200),
                         'q010_LES2': interp_z(r.fc['qt'][t0:], r.fc['zt'], 10),
                         'q200_LES2': interp_z(r.fc['qt'][t0:], r.fc['zt'], 200)
                       }

                dfs.append( pd.DataFrame(data, index=r.time[t0:]) )
            df_LES2 = pd.concat(dfs)

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

            df_CB = pd.DataFrame(data, index=cb_sf.time.values)
            df_CB.index = df_CB.index.round('1min')

            df_CB2 = pd.DataFrame(data2, index=cb_ns.time.values)
            df_CB2.index = df_CB2.index.round('1min')

            df_CB = pd.concat([df_CB, df_CB2], axis=1)
            df_CB.dropna(inplace=True)

            # Merge DataFrame's
            df = pd.concat([df_LES1, df_LES2, df_CB], axis=1)
            df.dropna(inplace=True)

            # Theoretical shortwave incoming radiation
            hour = df.index.hour + df.index.minute/60.
            doy  = df.index.dayofyear
            df['swd_theory'] = calc_swd(4.9, 51.97, hour, doy)
            df['is_night'] = df['swd_theory'] < 0.1


        def add_stat(label, obs, model, x, unit):
            s = Stats(obs, model)
            pl.title('{0}) RMSE={1:.1f} {3:}, ME={2:.1f} {3:}'.format(label, s.rmse, s.diff, unit), loc='left', fontsize=10)


        figsize2 = (8,5)
        figsize4 = (8,8)
        ratio = (2.5,1)
        raster = True
        ms = 1

        if True: 
            #
            # Wind 
            #

            pl.figure(figsize=figsize2)
            gs = gridspec.GridSpec(2, 2, width_ratios=ratio)

            ax=pl.subplot(gs[0,0])
            pl.plot(df.index, df['U010_CB'], 'o', mec=c2, mfc=c2, ms=ms, alpha=1, rasterized=raster)
            pl.plot(df.index, df['U010_LES1'], '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylim(0,17)
            add_stat('a', df['U010_CB'], df['U010_LES1'], df.index[50], r'm s$^{-1}$')
            pl.ylabel(r'$U_\mathrm{10m}$ (m s$^{-1}$)')
            format_ax()
            
            ax=pl.subplot(gs[1,0])
            pl.plot(df.index, df['U200_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['U200_LES1'], '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylim(0,17)
            add_stat('c', df['U200_CB'], df['U200_LES1'], df.index[50], r'm s$^{-1}$')
            pl.ylabel(r'$U_\mathrm{200m}$ (m s$^{-1}$)')
            format_ax()

            pl.subplot(gs[0,1])
            pl.title('b)', loc='left', fontsize=10)
            scatter_stat(df['U010_CB'], df['U010_LES1'], 'U_10m (m/s)', r'OBS (m s$^{-1}$)', r'LES (m s$^{-1}$)', rasterized=raster)

            pl.subplot(gs[1,1])
            pl.title('d)', loc='left', fontsize=10)
            scatter_stat(df['U200_CB'], df['U200_LES1'], 'U_200m (m/s)', r'OBS (m s$^{-1}$)', r'LES (m s$^{-1}$)', rasterized=raster)

            pl.tight_layout()
            pl.savefig('{}/wind_tser_scatter.pdf'.format(fig_path))


        if True:
            #
            # Temperature
            #

            pl.figure(figsize=figsize2)
            gs = gridspec.GridSpec(2, 2, width_ratios=ratio)

            ax=pl.subplot(gs[0,0])
            pl.plot(df.index, df['T010_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster, label='LES')
            pl.plot(df.index, df['T010_LES1'], '-', color=c1, linewidth=lw, label='OBS')
            pl.xlim(start, end)
            pl.ylabel(r'T$_\mathrm{10m}$ (K)')
            add_stat('a', df['T010_CB'], df['T010_LES1'], df.index[50], r'K')
            format_ax()
            pl.legend()

            ax=pl.subplot(gs[1,0])
            pl.plot(df.index, df['T200_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['T200_LES1'], '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylabel(r'T$_\mathrm{200m}$ (K)')
            add_stat('c', df['T200_CB'], df['T200_LES1'], df.index[50], r'K')
            format_ax()

            pl.subplot(gs[0,1])
            pl.title('b)', loc='left', fontsize=10)
            scatter_stat(df['T010_CB'], df['T010_LES1'], 'T_10m (K)', r'OBS (K)', r'LES (K)', rasterized=raster)

            pl.subplot(gs[1,1])
            pl.title('d)', loc='left', fontsize=10)
            scatter_stat(df['T200_CB'], df['T200_LES1'], 'T_10m (K)', r'OBS (K)', r'LES (K)', rasterized=raster)

            pl.tight_layout()
            pl.savefig('{}/temperature_tser_scatter.pdf'.format(fig_path))


        if True:
            #
            # Specific humidity
            #

            pl.figure(figsize=figsize2)
            gs = gridspec.GridSpec(2, 2, width_ratios=ratio)

            ax=pl.subplot(gs[0,0])
            pl.plot(df.index, df['q010_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['q010_LES1']*1000, '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylabel(r'q$_\mathrm{10m}$ (g kg$^{-1}$)')
            add_stat('a', df['q010_CB'], df['q010_LES1']*1000, df.index[50], r'g kg$^{-1}$')
            format_ax()

            ax=pl.subplot(gs[1,0])
            pl.plot(df.index, df['q200_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['q200_LES1']*1000, '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylabel(r'q$_\mathrm{200m}$ (g kg$^{-1}$)')
            add_stat('c', df['q200_CB'], df['q200_LES1']*1000, df.index[50], r'g kg$^{-1}$')
            format_ax()

            pl.subplot(gs[0,1])
            pl.title('b)', loc='left', fontsize=10)
            scatter_stat(df['q010_CB'], df['q010_LES1']*1000, 'q_10m (g/kg)', r'OBS (g kg$^{-1}$)', r'LES (g kg$^{-1}$)', rasterized=raster)

            pl.subplot(gs[1,1])
            pl.title('d)', loc='left', fontsize=10)
            scatter_stat(df['q200_CB'], df['q200_LES1']*1000, 'q_200m (g/kg)', r'OBS (g kg$^{-1}$)', r'LES (g kg$^{-1}$)', rasterized=raster)

            pl.tight_layout()
            pl.savefig('{}/spechum_tser_scatter.pdf'.format(fig_path))


        if True:
            #
            # Surface fluxes
            #

            pl.figure(figsize=figsize4)
            gs = gridspec.GridSpec(4, 2, width_ratios=ratio)

            ax=pl.subplot(gs[0,0])
            pl.plot(df.index, df['H_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['H_LES1'], '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylabel(r'H (W m$^{-2}$')
            add_stat('a', df['H_CB'], df['H_LES1'], df.index[50], r'W m$^{-2}$')
            format_ax()

            pl.subplot(gs[1,0], sharex=ax)
            pl.plot(df.index, df['LE_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['LE_LES1'], '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylabel(r'LE (W m$^{-2}$')
            add_stat('c', df['LE_CB'], df['LE_LES1'], df.index[50], r'W m$^{-2}$')
            format_ax()

            pl.subplot(gs[2,0], sharex=ax)
            pl.plot(df.index, df['G_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['G_LES1'], '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylabel(r'G (W m$^{-2}$')
            add_stat('e', df['G_CB'], df['G_LES1'], df.index[50], r'W m$^{-2}$')
            format_ax()

            pl.subplot(gs[3,0], sharex=ax)
            pl.plot(df.index, df['Qn_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['Qn_LES1'], '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylabel(r'Q$_\mathrm{net}$ (W m$^{-2}$')
            add_stat('g', df['Qn_CB'], df['Qn_LES1'], df.index[50], r'W m$^{-2}$')
            format_ax()

            # Scatter plots
            pl.subplot(gs[0,1])
            pl.title('b)', loc='left', fontsize=10)
            scatter_stat(df['H_CB'], df['H_LES1'], 'H (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', rasterized=raster)

            pl.subplot(gs[1,1])
            pl.title('d)', loc='left', fontsize=10)
            scatter_stat(df['LE_CB'], df['LE_LES1'], 'LE (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', rasterized=raster)

            pl.subplot(gs[2,1])
            pl.title('f)', loc='left', fontsize=10)
            scatter_stat(df['G_CB'], df['G_LES1'], 'G (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', rasterized=raster)

            pl.subplot(gs[3,1])
            pl.title('h)', loc='left', fontsize=10)
            scatter_stat(df['Qn_CB'], df['Qn_LES1'], 'Qnet (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', rasterized=raster)

            pl.tight_layout()
            pl.savefig('{}/surface_flux_tser_scatter.pdf'.format(fig_path))


        if True:
            #
            # Surface radiation
            #

            pl.figure(figsize=figsize4)

            gs = gridspec.GridSpec(4, 2, width_ratios=ratio)

            ax=pl.subplot(gs[0,0])
            pl.plot(df.index, -df['swd_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['swd_LES1'], '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylabel(r'$SW_\mathrm{down}$ (W m$^{-2}$)')
            pl.ylim(-1000,20)
            add_stat('a', -df['swd_CB'], df['swd_LES1'], df.index[50], r'W m$^{-2}$')
            format_ax()

            pl.subplot(gs[1,0], sharex=ax)
            pl.plot(df.index, df['swu_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['swu_LES1'], '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylabel(r'$SW_\mathrm{up}$ (W m$^{-2}$)')
            add_stat('c', df['swu_CB'], df['swu_LES1'], df.index[50], r'W m$^{-2}$')
            format_ax()

            pl.subplot(gs[2,0], sharex=ax)
            pl.plot(df.index, -df['lwd_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['lwd_LES1'], '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylabel(r'$LW_\mathrm{down}$ (W m$^{-2}$)')
            add_stat('e', -df['lwd_CB'], df['lwd_LES1'], df.index[50], r'W m$^{-2}$')
            format_ax()

            pl.subplot(gs[3,0], sharex=ax)
            pl.plot(df.index, df['lwu_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['lwu_LES1'], '-', color=c1, linewidth=lw)
            pl.xlim(start, end)
            pl.ylabel(r'$LW_\mathrm{up}$ (W m$^{-2}$)')
            add_stat('g', df['lwu_CB'], df['lwu_LES1'], df.index[50], r'W m$^{-2}$')
            format_ax()

            pl.subplot(gs[0,1])
            pl.title('b)', loc='left', fontsize=10)
            scatter_stat(-df['swd_CB'], df['swd_LES1'], 'SWdown (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', rasterized=raster)

            pl.subplot(gs[1,1])
            pl.title('d)', loc='left', fontsize=10)
            scatter_stat(df['swu_CB'], df['swu_LES1'], 'SWup (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', rasterized=raster)

            pl.subplot(gs[2,1])
            pl.title('f)', loc='left', fontsize=10)
            scatter_stat(-df['lwd_CB'], df['lwd_LES1'], 'LWdown (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', rasterized=raster)

            pl.subplot(gs[3,1])
            pl.title('h)', loc='left', fontsize=10)
            scatter_stat(df['lwu_CB'], df['lwu_LES1'], 'LWup (W/m2)', r'OBS (W m$^{-2}$)', r'LES (W m$^{-2}$)', rasterized=raster)

            pl.tight_layout()
            pl.savefig('{}/radiation_tser_scatter.pdf'.format(fig_path))



        if False:
            #
            # Surface meteo
            #

            pl.figure(figsize=(10,4.8))
            gs = gridspec.GridSpec(2, 2, width_ratios=[3.8,1])

            ax=pl.subplot(gs[0,0])
            pl.plot(df.index, df['cc_CB'], 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['cc_LES1'], '-', color=c1, linewidth=lw)
            pl.plot(cb_ns.time.values, cb_ns.cldcover_total/100., 'o', mfc=c2, mec=c2, ms=ms)
            pl.plot([start, end], [0,0], 'k:')
            pl.xlim(start, end)
            pl.ylim(-0.01,1.01)
            pl.ylabel(r'$cc$ (-)')

            pl.subplot(gs[1,0], sharex=ax)
            pl.plot(df.index, df['rr_CB']*6, 'o', mec=c2, mfc=c2, ms=ms, rasterized=raster)
            pl.plot(df.index, df['rr_LES1']*3600, '-', color=c1, linewidth=lw)
            pl.plot(cb_sm.time.values, cb_sm.RAIN*6, 'o', mfc=c2, mec=c2, ms=ms)
            pl.plot([start, end], [0,0], 'k:')
            pl.xlim(start, end)
            pl.ylim(0,2.5)
            pl.ylabel(r'$rr$ (mm h$^{-1}$)')

            ax=pl.subplot(gs[0,1])
            pl.scatter(df['cc_CB'], df['cc_LES1'], s=1, color=c2)
            lim_and_line2(df['cc_CB'], df['cc_LES1'])
            pl.xlabel(r'OBS (-)')
            pl.ylabel(r'LES (-)')

            ax=pl.subplot(gs[1,1])
            pl.scatter(df['rr_CB']*6, df['rr_LES1']*3600, s=1, color=c2)
            lim_and_line(0,2)
            pl.xlabel(r'OBS (mm h$^{-1}$)')
            pl.ylabel(r'LES (mm h$^{-1}$)')

            pl.tight_layout()
            pl.savefig('clouds_rain_tser_scatter.pdf')



















    if False:

        # --------------------------------
        #
        # Sync Cabauw and LES + HARMONIE + ERA5 data in Pandas data frame
        #
        # --------------------------------
        if 'df2' not in locals():

            # Heights to consider. Note: Cabauw data is not automatically
            # interpolated to other heights...
            heights = np.array((10,20,40,80,140,200))

            # ---- LES ----
            dfs = []
            for r in runs1:
                data = {}
                for z in heights:
                    data['U{0:03d}_LES1'.format(z)] = absval(interp_z(r.fc['u'], r.fc['zt'], z), interp_z(r.fc['v'], r.fc['zt'], z))
                    data['T{0:03d}_LES1'.format(z)] = interp_z(r.fc['T'],  r.fc['zt'], z)
                    data['q{0:03d}_LES1'.format(z)] = interp_z(r.fc['qt'], r.fc['zt'], z)*1000
                dfs.append( pd.DataFrame(data, index=r.time) )
            df_LES1 = pd.concat(dfs)

            dfs = []
            for r in runs2:
                data = {}
                for z in heights:
                    data['U{0:03d}_LES2'.format(z)] = absval(interp_z(r.fc['u'], r.fc['zt'], z), interp_z(r.fc['v'], r.fc['zt'], z))
                    data['T{0:03d}_LES2'.format(z)] = interp_z(r.fc['T'],  r.fc['zt'], z)
                    data['q{0:03d}_LES2'.format(z)] = interp_z(r.fc['qt'], r.fc['zt'], z)*1000
                dfs.append( pd.DataFrame(data, index=r.time) )
            df_LES2 = pd.concat(dfs)

            # ---- HARMONIE ----
            # Interpolate model levels to fixed heights
            if 'u_tmp1' not in locals():
                u_tmp1 = interp_zt(hm['u'], hm['z'], heights)
                v_tmp1 = interp_zt(hm['v'], hm['z'], heights)
                T_tmp1 = interp_zt(hm['T'], hm['z'], heights)
                q_tmp1 = interp_zt(hm['q'], hm['z'], heights)

            data = {}
            for k,z in enumerate(heights):
                data['U{0:03d}_HAM'.format(z)] = absval(u_tmp1[:,k], v_tmp1[:,k])
                data['T{0:03d}_HAM'.format(z)] = T_tmp1[:,k]
                data['q{0:03d}_HAM'.format(z)] = q_tmp1[:,k]*1000.

            df_HAM = pd.DataFrame(data, index=hm.time.values)
            df_HAM.index = df_HAM.index.round('1min')

            # ---- ERA5 ----
            # Interpolate model levels to fixed heights
            if 'u_tmp2' not in locals():
                u_tmp2 = interp_zt(e5['u'][:,::-1], e5['zf'][:,::-1], heights)
                v_tmp2 = interp_zt(e5['v'][:,::-1], e5['zf'][:,::-1], heights)
                T_tmp2 = interp_zt(e5['t'][:,::-1], e5['zf'][:,::-1], heights)
                q_tmp2 = interp_zt(e5['q'][:,::-1], e5['zf'][:,::-1], heights)

            data = {}
            for k,z in enumerate(heights):
                data['U{0:03d}_ERA'.format(z)] = absval(u_tmp2[:,k], v_tmp2[:,k])
                data['T{0:03d}_ERA'.format(z)] = T_tmp2[:,k]
                data['q{0:03d}_ERA'.format(z)] = q_tmp2[:,k]*1000.

            df_ERA = pd.DataFrame(data, index=e5.time.values)
            df_ERA.index = df_ERA.index.round('1min')

            # ---- Cabauw ----
            data = {}
            cb_k = {200:0, 140:1, 80:2, 40:3, 20:4, 10:5}

            for k,z in enumerate(heights):
                kk = cb_k[int(z)]

                data['U{0:03d}_CB'.format(z)] = cb_tm['F'] [:,kk]
                data['T{0:03d}_CB'.format(z)] = cb_tm['TA'][:,kk]
                data['q{0:03d}_CB'.format(z)] = cb_tm['Q'] [:,kk]

            df_CB = pd.DataFrame(data, index=cb_tm.time.values)
            df_CB.index = df_CB.index.round('1min')

            # ---- Merge all data frames, and drop missing rows to sync times ----
            df2 = pd.concat([df_LES1, df_LES2, df_HAM, df_ERA, df_CB], axis=1)
            df2.dropna(inplace=True)

            # Add day/night flag based on incoming shortwave radiation
            hour = df2.index.hour + df2.index.minute/60.
            doy  = df2.index.dayofyear
            df2['swd_theory'] = calc_swd(4.9, 51.97, hour, doy)
            df2['is_night'] = df2['swd_theory'] < 0.1


    
        if False:
            #
            # Print statistics
            #
            vars = ['q','T','U']
            units = ['g kg-1', 'K', 'm s-1']

            for var, unit in zip(vars, units):
                print('--------------------------')
                for z in heights:
                    CB    = df2['{0:}{1:03d}_CB'  .format(var, z)]
                    vLES1 = df2['{0:}{1:03d}_LES1'.format(var, z)]
                    vLES2 = df2['{0:}{1:03d}_LES2'.format(var, z)]
                    vHAM  = df2['{0:}{1:03d}_HAM' .format(var, z)]
                    vERA  = df2['{0:}{1:03d}_ERA' .format(var, z)]

                    # Statistics
                    cLES1 = Stats(CB, vLES1)
                    cLES2 = Stats(CB, vLES2)
                    cHAM  = Stats(CB, vHAM)
                    cERA  = Stats(CB, vERA)

                    print('{0}_{1:03d} | LES-HM: RMSE={2:6.2f}, diff={3:6.2f} | LES-E5: RMSE={4:6.2f}, diff={5:6.2f}| HAM: RMSE={6:6.2f}, diff={7:6.2f} | ERA: RMSE={8:6.2f}, diff={9:6.2f}'\
                            .format(var, z, cLES1.rmse, cLES1.diff, cLES2.rmse, cLES2.diff, cHAM.rmse, cHAM.diff, cERA.rmse, cERA.diff ))



        if False:
            #
            # Scatter plots with statistics
            #

            vars = ['q','T','U']
            units = [r'g kg$^\mathrm{-1}$', 'K', r'm s$^\mathrm{-1}$']

            for var,unit in zip(vars, units):

                pl.figure(figsize=(9,8)); sp=1

                for z in [10,80,200]:
                    for model, name in zip(['HAM','ERA','LES1', 'LES2'], ['HARMONIE','ERA5','LES-HM', 'LES-ERA5']):

                        if (sp-1)%3 == 0:
                            xlim = None
                            ylim = None

                        pl.subplot(3,4,sp); sp+=1
                        pl.title(r'${}_\mathrm{{{}m}}$'.format(var,z), loc='left')
                        xlim, ylim = scatter_stat2(df2['{0:}{1:03d}_CB'.format(var,z)], df2['{0:}{1:03d}_{2:}'.format(var,z,model)],
                                    '{0:}_{1:}m ({2:})'.format(var, z, unit), r'Cabauw ({})'.format(unit), r'{0:} ({1:})'.format(name,unit),
                                    df2['is_night'], c_day, c_nig, xlim, ylim)



                pl.tight_layout()


        if True:
            #
            # Line plots of statistics vs height
            #

            heights = np.array((10,20,40,80,140,200))
            is_night = df2['is_night']

            fig,ax = pl.subplots(nrows=3, ncols=2, figsize=(8,6))

            cm = pl.cm.Paired
            colors = [cm(5), cm(1), cm(5), cm(1)]
            dashes = ['--', '--', '-', '-']

            for i,var in enumerate(['U','T','q']):

                for model,name,color,lt in zip(['HAM','ERA','LES1','LES2'], ['HARMONIE','ERA5','LES-HARMONIE','LES-ERA5'], colors, dashes):

                    rmse_all = np.zeros_like(heights, dtype=np.float)
                    diff_all = np.zeros_like(heights, dtype=np.float)

                    rmse_day = np.zeros_like(heights, dtype=np.float)
                    diff_day = np.zeros_like(heights, dtype=np.float)

                    rmse_night = np.zeros_like(heights, dtype=np.float)
                    diff_night = np.zeros_like(heights, dtype=np.float)

                    for k,z in enumerate(heights):
                        v1 = '{0}{1:03d}_CB'.format(var,z)
                        v2 = '{0}{1:03d}_{2}'.format(var,z,model)

                        stat = Stats(df2[v1], df2[v2])
                        rmse_all[k] = stat.rmse
                        diff_all[k] = stat.diff

                        stat = Stats(df2[v1][~is_night], df2[v2][~is_night])
                        rmse_day[k] = stat.rmse
                        diff_day[k] = stat.diff

                        stat = Stats(df2[v1][is_night], df2[v2][is_night])
                        rmse_night[k] = stat.rmse
                        diff_night[k] = stat.diff

                    ax[i,0].plot(rmse_all,   heights, lt, color=color, linewidth=1.4, label='{}'.format(name))
                    ax[i,1].plot(diff_all,   heights, lt, color=color, linewidth=1.4, label='{}'.format(name))

                    #ax[i,0].plot(rmse_day,   heights, '-', color=color, dashes=[1,1], label='{}-day'.format(name))
                    #ax[i,1].plot(diff_day,   heights, '-', color=color, dashes=[1,1], label='{}-day'.format(name))
                    #                                
                    #ax[i,0].plot(rmse_night, heights, '-', color=color, dashes=[4,2], label='{}-night'.format(name))
                    #ax[i,1].plot(diff_night, heights, '-', color=color, dashes=[4,2], label='{}-night'.format(name))


            ax[0,0].set_xlabel(r'RMSE U (m s$^{-1}$)')
            ax[1,0].set_xlabel(r'RMSE T (K)')
            ax[2,0].set_xlabel(r'RMSE q (g kg$^{-1}$)')
            ax[0,1].set_xlabel(r'diff U (m s$^{-1}$)')
            ax[1,1].set_xlabel(r'diff T (K)')
            ax[2,1].set_xlabel(r'diff q (g kg$^{-1}$)')

            ax[0,0].set_ylabel('z (m)') 
            ax[1,0].set_ylabel('z (m)') 
            ax[2,0].set_ylabel('z (m)') 

            ax[0,1].vlines(0, ymin=0, ymax=200, colors='k', linestyles='dotted') 
            ax[1,1].vlines(0, ymin=0, ymax=200, colors='k', linestyles='dotted') 
            ax[2,1].vlines(0, ymin=0, ymax=200, colors='k', linestyles='dotted') 

            ax[0,0].set_xlim(0,1.7)
            ax[1,0].set_xlim(0,1)
            ax[2,0].set_xlim(0,1)

            ax[0,0].legend(ncol=1, fontsize=9)

            pl.tight_layout()

