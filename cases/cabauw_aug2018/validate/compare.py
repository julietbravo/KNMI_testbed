import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

import xarray as xr
import pandas as pd
import numpy as np

from datetime import datetime, timedelta

class Read_LES:
    def __init__(self, nc_path, start_date):
        print('Reading LES for {}'.format(start_date))

        date_str  = '{0:04d}{1:02d}{2:02d}'.format(start_date.year, start_date.month, start_date.day)
        self.fp   = xr.open_mfdataset('{}/profiles_{}.nc'.format(nc_path, date_str))
        self.ft   = xr.open_mfdataset('{}/tmser_{}.nc'.format(nc_path, date_str))
        self.time = [start_date+timedelta(seconds=int(self.fp.time[t])) for t in range(self.fp.time.size)]


def lim_and_line(vmin, vmax):
    pl.xlim(vmin, vmax)
    pl.ylim(vmin, vmax)
    pl.plot([vmin,vmax], [vmin,vmax], 'k:', linewidth=1)
    pl.plot([vmin,vmax], [0,0], 'k:', linewidth=1)
    pl.plot([0,0], [vmin,vmax], 'k:', linewidth=1)


def lim_and_line2(v1, v2):
    vmin = np.min((v1.min(), v2.min()))
    vmax = np.max((v1.max(), v2.max()))

    pl.xlim(vmin, vmax)
    pl.ylim(vmin, vmax)
    pl.plot([vmin,vmax], [vmin,vmax], 'k:', linewidth=1)
    pl.plot([vmin,vmax], [0,0], 'k:', linewidth=1)
    pl.plot([0,0], [vmin,vmax], 'k:', linewidth=1)


def format_ax():
    ax = pl.gca()
    h24 = mdates.HourLocator(interval=48)
    fmt = mdates.DateFormatter('%d-%m')
    ax.xaxis.set_major_locator(h24)
    ax.xaxis.set_major_formatter(fmt)


if __name__ == '__main__':
    pl.close('all')


    # Period to read/plot/..
    start = datetime(year=2016, month=8, day=4,  hour=0)
    end   = datetime(year=2016, month=8, day=19, hour=0)

    # Local file paths
    LES_path = '/Users/bart/meteo/data/KNMI_testbed/cabauw_20160804_20160818_soil_scaled_ccn_sgs'
    CB_path  = '/Users/bart/meteo/observations/Cabauw'
    HM_path  = '/Users/bart/meteo/data/Harmonie_LES_forcing/'


    # Read the LES data
    # -----------------
    if 'runs' not in locals():

        runs = []
        date  = start
        while date < end:
            l = Read_LES(LES_path, date)
            l.fp['rainrate'] /= (l.fp['rhof']*2.45e6)
            runs.append(l)
            date += timedelta(hours=24)


    # Read HARMONIE data
    # -----------------
    if 'hm' not in locals():

        iloc = 7+12     # 10x10km Cabauw
        files = []
        t = start
        while t <= end:
            files.append('{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/LES_forcing_{1:04d}{2:02d}{3:02d}{4:02d}.nc'\
                    .format(HM_path, t.year, t.month, t.day, t.hour))
            t += timedelta(hours=3)
        hm = xr.open_mfdataset(files)
        hm = hm.loc[{'domain': iloc}]


    # Read Cabauw observations
    # ------------------------
    if 'cb_sf' not in locals():

        cb_sm = xr.open_mfdataset('{0}/cesar_surface_mete*{1:04d}{2:02d}.nc'.format(CB_path, start.year, start.month), drop_variables=['valid_dates'])
        cb_sf = xr.open_mfdataset('{0}/cesar_surface_flux*{1:04d}{2:02d}.nc'.format(CB_path, start.year, start.month), drop_variables=['valid_dates'])
        cb_sr = xr.open_mfdataset('{0}/cesar_surface_radi*{1:04d}{2:02d}.nc'.format(CB_path, start.year, start.month), drop_variables=['valid_dates'])
        cb_ns = xr.open_mfdataset('{0}/cesar_nubiscope*{1:04d}{2:02d}.nc'   .format(CB_path, start.year, start.month), drop_variables=['valid_dates'])
        cb_lwc = xr.open_mfdataset('{0}/{1:04d}{2:02d}*_cabauw_lwc-scaled-adiabatic.nc'.format(CB_path, start.year, start.month))

        # Cabauw soil heat flux calculated as residual term
        cb_G_res = cb_sf.QN - (cb_sf.H + cb_sf.LE)


    # Sync times in a Pandas dataframe
    # ------------------
    if 'dfs' not in locals():

        # Read selected LES variables in Pandas DataFrame's
        dfs = []
        for r in runs:
            data = { 'LE_LES':  r.ft.LE,
                     'H_LES':   r.ft.H,
                     'G_LES':   r.ft.G0,
                     'Qn_LES':  r.ft.Qnet,
                     'swd_LES': r.fp.swd[:,0],
                     'swu_LES': r.fp.swu[:,0],
                     'lwd_LES': r.fp.lwd[:,0],
                     'lwu_LES': r.fp.lwu[:,0],
                     'cc_LES':  r.ft.cfrac,
                     'rr_LES':  r.fp.rainrate[:,0],
                     'lwp_LES': r.ft.lwp_bar}
            dfs.append( pd.DataFrame(data, index=r.time) )
        df_LES = pd.concat(dfs)

        # Put Cabauw observations in DataFrame
        data = { 'LE_CB':  cb_sf.LE,
                 'H_CB':   cb_sf.H,
                 'G_CB':   cb_sf.G0,
                 'G2_CB':  cb_G_res,
                 'Qn_CB':  cb_sf.QN,
                 'swd_CB': cb_sr.SWD,
                 'swu_CB': cb_sr.SWU,
                 'lwd_CB': cb_sr.LWD,
                 'lwu_CB': cb_sr.LWU,
                 'rr_CB':  cb_sm.RAIN}

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


    # Plot settings
    c1 = '#4d4d4d'   # Green
    c2 = '#4daf4a'    # Blue

    c_cb  = '#4daf4a'   # Green
    c_cb2 = '#377eb8'   # Blue
    c_da  = '#4d4d4d'   # Gray
    c_da2 = '#b2182b'   # DarkRed


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

        pl.figure()
        pl.subplot(311)
        pl.plot(np.arange(cb_sf['QN'].size), cb_sf['QN'] - (cb_sf['H']+cb_sf['LE']+cb_sf['G0']))

        pl.subplot(312)
        pl.plot(np.arange(cb_sf['QN'].size), cb_sf['QN'] - (cb_sf['H']+cb_sf['LE']))

        pl.subplot(313)
        pl.plot(np.arange(cb_sf['QN'].size), cb_sf['IG0'])


    if False:
        # --------------
        # Cloud statistics
        # --------------
        def rmean(x, N=20):
            return np.convolve(x, np.ones((N,))/N, mode='same')


        pl.figure()
        pl.plot(df_LES.index, df_LES.lwp_LES, color=c1)
        pl.plot(cb_lwc.time, rmean(cb_lwc.lwp), color=c2)
        pl.xlim(start, end)


        pl.figure()
        ax=pl.subplot(311)
        for i,r in enumerate(runs):
            pl.pcolormesh(r.time, r.fp.zt, r.fp.cfrac.T, vmin=0, vmax=1, cmap=pl.cm.tab20c_r)
        pl.colorbar()

        pl.subplot(312, sharex=ax)
        for i,r in enumerate(runs):
            pl.pcolormesh(r.time, r.fp.zt, r.fp.ql.T*1000, vmin=0, vmax=0.25, cmap=pl.cm.tab20c_r)
        pl.colorbar()

        pl.subplot(313, sharex=ax)
        for i,r in enumerate(runs):
            pl.pcolormesh(r.time, r.fp.zt, np.log(r.fp.sv002.T+1e-12), cmap=pl.cm.tab20c_r)
        pl.colorbar()


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
            pl.plot(r.time, r.fp.rainrate[:,0]*3600, '-', color=c1)
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


        pl.figure()
        pl.plot(df.index, np.cumsum(df['rr_LES'])*600, color=c1, label='LES')
        pl.plot(df.index, np.cumsum(df['rr_CB']), color=c2, label='Cabauw')
        pl.ylabel('cumulative sum(rr) (mm)')
        pl.legend()


    if True:
        # --------------
        # Surface fluxes
        # --------------


        pl.figure(figsize=(10,8))

        gs = gridspec.GridSpec(4, 2, width_ratios=[3.8,1])

        # Time series
        ax=pl.subplot(gs[0,0])
        pl.plot(cb_sf.time.values, cb_sf.H, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.ft.H, '-', color=c1)
        #pl.plot(hm.time, -hm.H)
        pl.xlim(start, end)
        pl.ylabel(r'H (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[1,0], sharex=ax)
        pl.plot(cb_sr.time.values,  cb_sf.LE, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.ft.LE, '-', color=c1)
        #pl.plot(hm.time, -hm.LE)
        pl.xlim(start, end)
        pl.ylabel(r'LE (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[2,0], sharex=ax)
        pl.plot(cb_sr.time.values, cb_sf.G0, 'o', mfc=c2, mec=c2, ms=2)
        pl.plot(cb_sr.time.values, cb_G_res, 'o', mfc='r', mec='r', ms=1)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.ft.G0, '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'G (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[3,0], sharex=ax)
        pl.plot(cb_sr.time.values, cb_sf.QN, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.ft.Qnet, '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'Q$_\mathrm{net}$ (W m$^{-2}$')
        format_ax()

        # Scatter plots
        pl.subplot(gs[0,1])
        pl.scatter(df['H_CB'], df['H_LES'], s=1, color=c2)
        lim_and_line2(df['H_CB'], df['H_LES'])
        pl.xlabel(r'OBS (W m$^{-2}$)')
        pl.ylabel(r'LES (W m$^{-2}$)')

        pl.subplot(gs[1,1])
        pl.scatter(df['LE_CB'], df['LE_LES'], s=1, color=c2)
        lim_and_line2(df['LE_CB'], df['LE_LES'])
        pl.xlabel(r'OBS (W m$^{-2}$)')
        pl.ylabel(r'LES (W m$^{-2}$)')

        pl.subplot(gs[2,1])
        pl.scatter(df['G_CB'], df['G_LES'], s=1, color=c2)
        pl.scatter(df['G2_CB'], df['G_LES'], s=1, color='r')
        lim_and_line2(df['G_CB'], df['G_LES'])
        pl.xlabel(r'OBS (W m$^{-2}$)')
        pl.ylabel(r'LES (W m$^{-2}$)')

        pl.subplot(gs[3,1])
        pl.scatter(df['Qn_CB'], df['Qn_LES'], s=1, color=c2)
        lim_and_line2(df['Qn_CB'], df['Qn_LES'])
        pl.xlabel(r'OBS (W m$^{-2}$)')
        pl.ylabel(r'LES (W m$^{-2}$)')

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
            pl.plot(r.time, r.fp.swd[:,0], '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'$SW_\mathrm{down}$ (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[1,0], sharex=ax)
        pl.plot(cb_sr.time.values,  cb_sr.SWU, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.fp.swu[:,0], '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'$SW_\mathrm{up}$ (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[2,0], sharex=ax)
        pl.plot(cb_sr.time.values, -cb_sr.LWD, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.fp.lwd[:,0], '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'$LW_\mathrm{down}$ (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[3,0], sharex=ax)
        pl.plot(cb_sr.time.values, cb_sr.LWU, 'o', mfc=c2, mec=c2, ms=2)
        for i,r in enumerate(runs):
            pl.plot(r.time, r.fp.lwu[:,0], '-', color=c1)
        pl.xlim(start, end)
        pl.ylabel(r'$LW_\mathrm{up}$ (W m$^{-2}$')
        format_ax()

        pl.subplot(gs[0,1])
        pl.scatter(-df['swd_CB'], df['swd_LES'], s=1, color=c2)
        lim_and_line(-800,0)
        pl.xlabel(r'OBS (W m$^{-2}$)')
        pl.ylabel(r'LES (W m$^{-2}$)')

        pl.subplot(gs[1,1])
        pl.scatter(df['swu_CB'], df['swu_LES'], s=1, color=c2)
        lim_and_line(0,200)
        pl.xlabel(r'OBS (W m$^{-2}$)')
        pl.ylabel(r'LES (W m$^{-2}$)')

        pl.subplot(gs[2,1])
        pl.scatter(-df['lwd_CB'], df['lwd_LES'], s=1, color=c2)
        lim_and_line(-420,-280)
        pl.xlabel(r'OBS (W m$^{-2}$)')
        pl.ylabel(r'LES (W m$^{-2}$)')

        pl.subplot(gs[3,1])
        pl.scatter(df['lwu_CB'], df['lwu_LES'], s=1, color=c2)
        lim_and_line(350,460)
        pl.xlabel(r'OBS (W m$^{-2}$)')
        pl.ylabel(r'LES (W m$^{-2}$)')

        pl.tight_layout()
        pl.savefig('radiation_tser_scatter.pdf')


    if False:
        # --------------
        # Closure SEB
        # --------------

        pl.figure(figsize=(10,8))

        gs = gridspec.GridSpec(4, 2, width_ratios=[3.8,1])



