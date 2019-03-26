import matplotlib.pyplot as pl
import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime, timedelta



class Read_LES:
    def __init__(self, nc_path, start_date):
        print('Reading LES for {}'.format(start_date))

        date_str = '{0:04d}{1:02d}{2:02d}'.format(start_date.year, start_date.month, start_date.day)
        self.fp = xr.open_mfdataset('{}/profiles_{}.nc'.format(nc_path, date_str))
        self.ft = xr.open_mfdataset('{}/tmser_{}.nc'.format(nc_path, date_str))
        
        self.time = [start_date+timedelta(seconds=int(self.fp.time[t])) for t in range(self.fp.time.size)]

class Read_Cabauw:
    def __init__(self, path, date):
        print('Reading Cabauw {}'.format(date))

def lim_and_line(vmin, vmax):
    pl.xlim(vmin, vmax)
    pl.ylim(vmin, vmax)
    pl.plot([vmin,vmax], [vmin,vmax], 'k:')
    pl.plot([vmin,vmax], [0,0], 'k:')
    pl.plot([0,0], [vmin,vmax], 'k:')



if __name__ == '__main__':
    pl.close('all')


    # Period to read/plot/..
    start = datetime(year=2016, month=8, day=4,  hour=0)
    end   = datetime(year=2016, month=8, day=19, hour=0)

    # Local file paths
    LES_path  = '/Users/bart/meteo/data/KNMI_testbed/cabauw_20160804_20160818'
    CB_path   = '/Users/bart/meteo/observations/Cabauw'

    # Read the LES data
    # -----------------
    if 'runs' not in locals():

        runs = []
        date  = start
        while date < end:
            runs.append( Read_LES(LES_path, date) )
            date += timedelta(hours=24)

    # Read Cabauw observations
    # ------------------------
    if 'cb_sf' not in locals():

        cb_sm = xr.open_mfdataset('{0}/cesar_surface_mete*{1:04d}{2:02d}.nc'.format(CB_path, start.year, start.month), drop_variables=['valid_dates'])
        cb_sf = xr.open_mfdataset('{0}/cesar_surface_flux*{1:04d}{2:02d}.nc'.format(CB_path, start.year, start.month), drop_variables=['valid_dates'])
        cb_sr = xr.open_mfdataset('{0}/cesar_surface_radi*{1:04d}{2:02d}.nc'.format(CB_path, start.year, start.month), drop_variables=['valid_dates'])
        cb_ns = xr.open_mfdataset('{0}/cesar_nubiscope*{1:04d}{2:02d}.nc'.format(CB_path, start.year, start.month), drop_variables=['valid_dates'])


    # Plot settings
    c_cb  = '#4daf4a'
    c_da  = '#4d4d4d'
    c_da2 = '#b2182b'

    if True:
        # --------------
        # Surface meteo
        # --------------

        pl.figure(figsize=(10,7))
        pl.subplot(211)

        for i,r in enumerate(runs):
            label = 'LES' if i==0 else ''
            pl.plot(r.time, r.ft.cfrac, '-', color=c_da, label=label)

        pl.plot(cb_ns.time.values, cb_ns.cldcover_total/100., 'o', label='Cabauw', mfc=c_cb, mec=c_cb, ms=2)

        pl.plot([start, end], [0,0], 'k:')
        pl.xlim(start, end)
        pl.ylim(0,1)
        pl.legend()
        pl.ylabel(r'$cc$ (-)')

        pl.subplot(212)

        for i,r in enumerate(runs):
            label = 'LES' if i==0 else ''
            pl.plot(r.time, r.fp.rainrate[:,0]/(r.fp.rhof[:,0]*2.45e6)*600, '-', color=c_da, label=label)

        pl.plot(cb_sm.time.values, cb_sm.RAIN, 'o', label='Cabauw', mfc=c_cb, mec=c_cb, ms=2)

        pl.plot([start, end], [0,0], 'k:')
        pl.xlim(start, end)
        pl.ylim(0,1)
        pl.legend()
        pl.ylabel(r'$rr$ (-)')

        pl.tight_layout()


    if True:
        # --------------
        # Surface fluxes
        # --------------

        pl.figure(figsize=(10,7))
        pl.subplot(211)

        for i,r in enumerate(runs):
            label = 'LES' if i==0 else ''
            pl.plot(r.time, r.ft.H, '-', color=c_da, label=label)

        pl.plot(cb_sf.time.values, cb_sf.H, 'o', label='Cabauw', mfc=c_cb, mec=c_cb, ms=2)

        pl.plot([start, end], [0,0], 'k:')
        pl.xlim(start, end)
        pl.ylim(-75,200)
        pl.legend()
        pl.ylabel(r'$H$ (W m$^{-2}$')

        pl.subplot(212)

        for i,r in enumerate(runs):
            label = 'LES' if i==0 else ''
            pl.plot(r.time, r.ft.LE, '-', color=c_da, label=label)

        pl.plot(cb_sf.time.values, cb_sf.LE, 'o', label='Cabauw', mfc=c_cb, mec=c_cb, ms=2)

        pl.plot([start, end], [0,0], 'k:')
        pl.xlim(start, end)
        pl.legend()
        pl.ylabel(r'$LE$ (W m$^{-2}$')

        pl.tight_layout()


    if True:
        # --------------
        # Surface radiation
        # --------------

        pl.figure(figsize=(10,7))

        pl.subplot(211)

        for i,r in enumerate(runs):
            pl.plot(r.time, r.fp.swd[:,0], '-', color=c_da)
            pl.plot(r.time, r.fp.swu[:,0], '--', color=c_da)

        pl.plot(cb_sr.time.values, -cb_sr.SWD, 'o', mfc=c_cb, mec=c_cb, ms=2)
        pl.plot(cb_sr.time.values, cb_sr.SWU, 'o', mfc=c_cb, mec=c_cb, ms=2)

        pl.plot([start, end], [0,0], 'k:')
        pl.xlim(start, end)
        pl.ylabel(r'$SW$ (W m$^{-2}$')

        pl.subplot(212)

        for i,r in enumerate(runs):
            pl.plot(r.time, r.fp.lwd[:,0], '-', color=c_da)
            pl.plot(r.time, r.fp.lwu[:,0], '--', color=c_da)

        pl.plot(cb_sr.time.values, -cb_sr.LWD, 'o', mfc=c_cb, mec=c_cb, ms=2)
        pl.plot(cb_sr.time.values, cb_sr.LWU, 'o', mfc=c_cb, mec=c_cb, ms=2)

        pl.plot([start, end], [0,0], 'k:')
        pl.xlim(start, end)
        pl.ylabel(r'$LW$ (W m$^{-2}$')

        pl.tight_layout()




    if True:
        # ------------------
        # Compare LES vs OBS
        # ------------------

        # Read selected LES variables in Pandas DataFrame's
        dfs = []
        for r in runs:
            data = { 'LE_LES':  r.ft.LE, 
                     'H_LES':   r.ft.H,
                     'swd_LES': r.fp.swd[:,0],
                     'lwd_LES': r.fp.lwd[:,0],
                     'lwu_LES': r.fp.lwu[:,0],
                     'cc_LES':  r.ft.cfrac}
            dfs.append( pd.DataFrame(data, index=r.time) )
        df_LES = pd.concat(dfs)

        # Put Cabauw observations in DataFrame
        data = { 'LE_CB':  cb_sf.LE,
                 'H_CB':   cb_sf.H,
                 'swd_CB': cb_sr.SWD,
                 'lwd_CB': cb_sr.LWD,
                 'lwu_CB': cb_sr.LWU }

        data2 = {'cc_CB': cb_ns.cldcover_total/100.}

        df_CB = pd.DataFrame(data, index=cb_sf.time)
        df_CB.index = df_CB.index.round('1min')

        df_CB2 = pd.DataFrame(data2, index=cb_ns.time)
        df_CB2.index = df_CB2.index.round('1min')

        df_CB = pd.concat([df_CB, df_CB2], axis=1)
        df_CB.dropna(inplace=True) 

        # Merge DataFrame's
        df = pd.concat([df_LES, df_CB], axis=1)
        df.dropna(inplace=True) 



        pl.figure(figsize=(8,4))
        pl.subplot(121, aspect='equal')
        pl.scatter(df['H_CB'], df['H_LES'], s=1)
        lim_and_line(-75,200)
        pl.xlabel(r'H$_\mathrm{OBS}$ (W m$^{-2}$)')
        pl.ylabel(r'H$_\mathrm{LES}$ (W m$^{-2}$)')

        pl.subplot(122, aspect='equal')
        pl.scatter(df['LE_CB'], df['LE_LES'], s=1)
        lim_and_line(-50,500)
        pl.xlabel(r'LE$_\mathrm{OBS}$ (W m$^{-2}$)')
        pl.ylabel(r'LE$_\mathrm{LES}$ (W m$^{-2}$)')

        pl.tight_layout()


        pl.figure(figsize=(11.5,4))
        pl.subplot(131, aspect='equal')
        pl.scatter(-df['swd_CB'], df['swd_LES'], s=1)
        lim_and_line(-800,0)
        pl.xlabel(r'SWD$_\mathrm{OBS}$ (W m$^{-2}$)')
        pl.ylabel(r'SWD$_\mathrm{LES}$ (W m$^{-2}$)')

        pl.subplot(132, aspect='equal')
        pl.scatter(-df['lwd_CB'], df['lwd_LES'], s=1)
        lim_and_line(-400,-270)
        pl.xlabel(r'LWD$_\mathrm{OBS}$ (W m$^{-2}$)')
        pl.ylabel(r'LWD$_\mathrm{LES}$ (W m$^{-2}$)')

        pl.subplot(133, aspect='equal')
        pl.scatter(df['lwu_CB'], df['lwu_LES'], s=1)
        lim_and_line(330,475)
        pl.xlabel(r'LWU$_\mathrm{OBS}$ (W m$^{-2}$)')
        pl.ylabel(r'LWU$_\mathrm{LES}$ (W m$^{-2}$)')

        pl.tight_layout()


        pl.figure(figsize=(4,4))
        pl.subplot(111, aspect='equal')
        pl.scatter(df['cc_CB'], df['cc_LES'], s=1)
        lim_and_line(0,1)
        pl.xlabel(r'cc$_\mathrm{OBS}$ (-)')
        pl.ylabel(r'cc$_\mathrm{LES}$ (-)')
