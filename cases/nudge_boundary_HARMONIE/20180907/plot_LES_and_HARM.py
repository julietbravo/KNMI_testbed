import matplotlib.pyplot as pl
import xarray as xr
import numpy as np

import cartopy.crs as ccrs

# Custom code
from hybrid_sigma_grid import Sigma_grid
from helpers import interp_z, setup_map, plot_LES_domain, add_labels, minmax

# Enable LaTeX plotting
from matplotlib import rc
rc('text', usetex=True)
rc('font', size=13)
rc('legend', fontsize=11)
rc('text.latex', preamble=r'\usepackage{sansmathfonts}')

#proj_str = '+proj=lcc +lat_1=52.500000 +lat_2=52.500000 +lat_0=52.500000 +lon_0=.000000 +k_0=1.0 +x_0=649536.512574 +y_0=1032883.739533 +a=6371220.000000 +b=6371220.000000'


def f(fld):
    return gaussian_filter(fld, sigma=30)

def minmax(arr):
    return np.floor(arr.min()), np.ceil(arr.max())

def albedo(lwp):
    tau = 0.19 * lwp**(5./6.) * 300e6**(1./3.)
    return tau / (6.8 + tau)


if __name__ == '__main__':

    pl.close('all')
    pl.ion()

    #
    # Map projection
    #
    proj_str = '+proj=lcc +lat_1=52.500000 +lat_2=52.500000 +lat_0=52.500000 +lon_0=.000000 \
            +k_0=1.0 +x_0=649536.512574 +y_0=1032883.739533 +a=6371220.000000 +b=6371220.000000'

    proj = pyproj.Proj(proj_str)

    path = '/nobackup/users/stratum/KNMI_testbed/cases/nudge_boundary_HARMONIE_20180907/'
    #path = '/Users/bart/meteo/data/KNMI_testbed/nudge_boundary_HARMONIE_20180907/data/'

    #
    # Read HARMONIE data
    #
    if 'hm_u' not in locals():
        # Read HARMONIE 3D fields
        hm_u = xr.open_dataset('{}/HARMONIE/ua.Slev.his.MINI.DOWA_40h12tg2_fERA5_ptF.20180907.nc'.format(path))
        hm_v = xr.open_dataset('{}/HARMONIE/va.Slev.his.MINI.DOWA_40h12tg2_fERA5_ptF.20180907.nc'.format(path))
        hm_t = xr.open_dataset('{}/HARMONIE/ta.Slev.his.MINI.DOWA_40h12tg2_fERA5_ptF.20180907.nc'.format(path))
        hm_q = xr.open_dataset('{}/HARMONIE/hus.Slev.his.MINI.DOWA_40h12tg2_fERA5_ptF.20180907.nc'.format(path))
        hm_p = xr.open_dataset('{}/HARMONIE/ps.his.MINI.DOWA_40h12tg2_fERA5_ptF.20180907.nc'.format(path))

        # Calculate HARMONIE pressure/height levels
        grid = Sigma_grid('{}/H40_65lev.txt'.format(path))

        # Time in seconds since 00 UTC
        time = [float(hm_u['time'][i] - hm_u['time'][0])/1e9 for i in range(hm_u.dims['time'])]

    itot_hm = hm_u.dims['x']
    jtot_hm = hm_u.dims['y']
    ktot_hm = hm_u.dims['lev']

    #
    # Read DALES data
    #
    k_les = 42  # {1,4,15,26,42,66,89,111}
    les_u = xr.open_dataset('{0}/DALES/crossxy.{1:04d}.u.001.nc'  .format(path, k_les))
    les_v = xr.open_dataset('{0}/DALES/crossxy.{1:04d}.v.001.nc'  .format(path, k_les))
    les_q = xr.open_dataset('{0}/DALES/crossxy.{1:04d}.qt.001.nc' .format(path, k_les))
    les_t = xr.open_dataset('{0}/DALES/crossxy.{1:04d}.thl.001.nc'.format(path, k_les))
    k_les -= 1

    les_lwp = xr.open_dataset('{0}/DALES/crossxy.lwp.001.nc'.format(path, k_les))
    les_alb = albedo(les_lwp['lwpxy'].values)

    lon = np.load('{}/DALES/lon_LES.npy'.format(path))
    lat = np.load('{}/DALES/lat_LES.npy'.format(path))

    z_les = np.loadtxt('{}/DALES/prof.inp.001'.format(path), skiprows=2, usecols=0)

    #
    # Time indices in HARMONIE dan DALES data
    #
    t_les = 5
    t_hm  = t_les+9

    #
    # Calculate pressure/height levels HARMONIE
    #
    Tv    = grid.calc_virtual_temp(hm_t['ta'][t_hm,:,:,:].values, hm_q['hus'][t_hm,:,:,:].values)
    ph    = grid.calc_half_level_pressure(hm_p['ps'][t_hm,:,:].values)
    zh_hm = grid.calc_half_level_Zg(ph, Tv)[::-1,:,:]
    p_hm  = grid.calc_full_level_pressure(ph)[::-1,:,:]
    z_hm  = grid.calc_full_level_Zg(zh_hm)


    # Why did I need this?
    x0_les = 700000
    y0_les = 1200000

    print('z_LES={}'.format(z_les[k_les]))

    # Test test
    ztest = np.zeros((jtot_hm, itot_hm), dtype=np.float64)
    interp_z(ztest, z_hm[::-1,:,:], z_hm[::-1,:,:], z_les[k_les], itot_hm, jtot_hm, ktot_hm)
    print('z_HARMONIE=', ztest.min(), ztest.max())


    if True:
        #
        # Plot LES domain
        #

        fig = pl.figure(figsize=(5, 5.))
        fig.subplots_adjust(left=0.11, bottom=0.07, right=0.98, top=0.98, wspace=0.08)

        ax = setup_map(subplot=111, extent=[-5, 12, 50, 60])
        ax.outline_patch.set_visible(False)

        plot_LES_domain(lon, lat, 'r')
        add_labels()

        pl.savefig('large_les_domain.pdf')

    if True:
        #
        # Clouds
        #
        vmin = 0
        vmax = 1

        #
        # Plot!
        #
        fig = pl.figure(figsize=(6, 6))
        fig.subplots_adjust(left=0.03, bottom=0.05, right=0.98, top=0.9)

        ax = setup_map(subplot=111)
        ax.outline_patch.set_visible(False)

        pc = ax.pcolormesh(
                lon, lat, les_alb[t_les,:,:].T,
                vmin=vmin, vmax=vmax, cmap=pl.cm.Blues_r,
                rasterized=True, transform=ccrs.PlateCarree())

        plot_LES_domain(lon, lat, 'r')
        add_labels()

        cax = fig.add_axes([1-0.45, 0.96, 0.35, 0.015])
        cb = pl.colorbar(pc, orientation='horizontal', cax=cax)
        pl.figtext(1-0.46, 0.96, r'albedo (-)', ha='right', fontsize=10)

        pl.savefig('albedo.pdf'.format(z_les[k_les]))


    if True:
        #
        # Wind barbs
        #
        ua = np.zeros((jtot_hm, itot_hm))
        va = np.zeros((jtot_hm, itot_hm))

        interp_z(ua, hm_u['ua'][t_hm,::-1,:,:].values, z_hm[::-1,:,:], z_les[k_les], itot_hm, jtot_hm, ktot_hm)
        interp_z(va, hm_v['va'][t_hm,::-1,:,:].values, z_hm[::-1,:,:], z_les[k_les], itot_hm, jtot_hm, ktot_hm)

        vmin, vmax = minmax(calc_abs(ua, va))

        #
        # Plot!
        #
        fig = pl.figure(figsize=(8, 4.5))
        fig.subplots_adjust(left=0.07, bottom=0.05, right=0.98, top=0.88, wspace=0.08)

        ax = setup_map(subplot=121)
        pl.title(r'$z$ = {0:.1f} m'.format(z_les[k_les]), loc='left')
        ax.outline_patch.set_visible(False)

        ax.pcolormesh(
                hm_u['lon'], hm_u['lat'],
                calc_abs(ua, va),
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r,
                rasterized=True, transform=ccrs.PlateCarree())

        ev=8
        s = np.s_[ev//2::ev,ev//2::ev]
        pl.barbs(
                hm_u['lon'].values[s], hm_u['lat'].values[s], ua[s], va[s],
                length=5, pivot='middle', linewidth=0.8, transform=ccrs.PlateCarree())

        plot_LES_domain(lon, lat, 'r')
        add_labels()

        ax = setup_map(subplot=122)
        ax.outline_patch.set_visible(False)

        pc = ax.pcolormesh(
                lon, lat,
                calc_abs(les_u['uxy'][t_les,:,:].values, les_v['vxy'][t_les,:,:].values).T,
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r,
                rasterized=True, transform=ccrs.PlateCarree())

        ev=100
        s = np.s_[ev//2::ev,ev//2::ev]
        pl.barbs(
                lon[s], lat[s], les_u['uxy'][t_les,:].values[s].T, les_v['vxy'][t_les,:].values[s].T,
                length=5, pivot='middle', linewidth=0.8, transform=ccrs.PlateCarree())

        plot_LES_domain(lon, lat, 'r')
        add_labels(ylabels=False)

        cax = fig.add_axes([1-0.28, 0.96, 0.25, 0.02])
        cb = pl.colorbar(pc, orientation='horizontal', cax=cax)
        pl.figtext(1-0.29, 0.96, r'$U$ (m s$^{-1}$)', ha='right', fontsize=10)

        pl.savefig('U_{0:.0f}m.pdf'.format(z_les[k_les]))


    if True:
        #
        # Specific humidity
        #
        # Interpolate HARMONIE field to LES height
        hus = np.zeros((jtot_hm, itot_hm))

        interp_z(hus, hm_q['hus'][t_hm,::-1,:,:].values, z_hm[::-1,:,:], z_les[k_les], itot_hm, jtot_hm, ktot_hm)

        if k_les == 0:
            vmin = 5
            vmax = 10
        elif k_les == 41:
            vmin = 1
            vmax = 8
        else:
            vmin, vmax = minmax(hus*1e3)

        #
        # Plot!
        #
        fig = pl.figure(figsize=(6, 6))
        fig.subplots_adjust(left=0.03, bottom=0.05, right=0.98, top=0.9)

        ax = setup_map(subplot=111)
        pl.title(r'$z$ = {0:.1f} m'.format(z_les[k_les]), loc='left')
        ax.outline_patch.set_visible(False)

        ax.pcolormesh(
                hm_u['lon'], hm_u['lat'], hus*1e3,
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r,
                rasterized=True, transform=ccrs.PlateCarree())

        pc = ax.pcolormesh(
                lon, lat, les_q['qtxy'][t_les,:,:].values.T*1e3,
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r,
                rasterized=True, transform=ccrs.PlateCarree())

        plot_LES_domain(lon, lat, 'r')
        add_labels()

        cax = fig.add_axes([1-0.45, 0.96, 0.35, 0.015])
        cb = pl.colorbar(pc, orientation='horizontal', cax=cax)
        pl.figtext(1-0.46, 0.96, r'$q_\mathrm{t}$ (g kg$^{-1}$)', ha='right', fontsize=10)


        pl.savefig('qt_{0:.0f}m.pdf'.format(z_les[k_les]))






