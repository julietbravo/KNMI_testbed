import matplotlib.pyplot as pl
import xarray as xr
import numpy as np

import cartopy.crs as ccrs

# Custom scripts
import spectral_tools as st
import interpolate as ip
from hybrid_sigma_grid import Sigma_grid
from helpers import interp_z, setup_map, plot_LES_domain, calc_abs, minmax, Grid


if __name__ == '__main__':
    pl.close('all')
    pl.ion()


    path = '/nobackup/users/stratum/KNMI_testbed/cases/nudge_boundary_HARMONIE_20180907/'
    #path = '/Users/bart/meteo/data/KNMI_testbed/nudge_boundary_HARMONIE_20180907/data/'

    #
    # Settings LES domain
    #
    x0 = 700000
    y0 = 1200000
 
    # Domain size LES (m)
    xsize = 1680*200
    ysize = 1680*200
 
    # Number of grid points LES
    itot = 1680
    jtot = 1680
    ktot = 128

    # LES grid definition (simple version..)
    LES_grid = Grid(xsize, ysize, 1., itot, jtot, ktot)

    #
    # Read LES data:
    #
    if 'les_u' not in locals():
        k_les = 1    # in {1,4,15,26,42,66,89,111}
        les_u = xr.open_dataset('{0}/DALES/crossxy.{1:04d}.u.001.nc'  .format(path, k_les))
        les_v = xr.open_dataset('{0}/DALES/crossxy.{1:04d}.v.001.nc'  .format(path, k_les))
        les_q = xr.open_dataset('{0}/DALES/crossxy.{1:04d}.qt.001.nc' .format(path, k_les))
        les_t = xr.open_dataset('{0}/DALES/crossxy.{1:04d}.thl.001.nc'.format(path, k_les))
        k_les -= 1

        lon = np.load('{}/DALES/lon_LES.npy'.format(path))
        lat = np.load('{}/DALES/lat_LES.npy'.format(path))

        z_les = np.loadtxt('{}/DALES/prof.inp.001'.format(path), skiprows=2, usecols=0)

    #
    # Read HARMONIE data:
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

        # Grid dimensions
        itot_hm = hm_u.dims['x']
        jtot_hm = hm_u.dims['y']
        ktot_hm = hm_u.dims['lev']

    #
    # Time to analyse
    #
    t_les = 1
    t_hm  = t_les+9

    #
    # Calculate pressure/height levels HARMONIE
    #
    Tv    = grid.calc_virtual_temp(hm_t['ta'][t_hm,:,:,:].values, hm_q['hus'][t_hm,:,:,:].values)
    ph    = grid.calc_half_level_pressure(hm_p['ps'][t_hm,:,:].values)
    zh_hm = grid.calc_half_level_Zg(ph, Tv)[::-1,:,:]
    p_hm  = grid.calc_full_level_pressure(ph)[::-1,:,:]
    z_hm  = grid.calc_full_level_Zg(zh_hm)

    #
    # Interpolate HARMONIE to LES height
    #
    ua = np.zeros((jtot_hm, itot_hm))
    va = np.zeros((jtot_hm, itot_hm))
    qa = np.zeros((jtot_hm, itot_hm))

    interp_z(ua, hm_u['ua'] [t_hm,::-1,:,:].values, z_hm[::-1,:,:], z_les[k_les], itot_hm, jtot_hm, ktot_hm)
    interp_z(va, hm_v['va'] [t_hm,::-1,:,:].values, z_hm[::-1,:,:], z_les[k_les], itot_hm, jtot_hm, ktot_hm)
    interp_z(qa, hm_q['hus'][t_hm,::-1,:,:].values, z_hm[::-1,:,:], z_les[k_les], itot_hm, jtot_hm, ktot_hm)

    vmin, vmax = minmax(calc_abs(ua, va))

    #
    # Interpolate HARMONIE to LES grid
    #
    intp = ip.Grid_interpolator(
            hm_u['x'].values, hm_u['y'].values, None, 
            LES_grid.x, LES_grid.y, None, 
            LES_grid.xh, LES_grid.yh, None, 
            x0, y0)

    ua_at_les = intp.interpolate_2d(ua, 'x', 'y')
    va_at_les = intp.interpolate_2d(va, 'x', 'y')



    if True:
        # Plot!
        pl.figure(figsize=(10,4))

        ax=setup_map(subplot=131)
        plot_LES_domain(lon, lat, 'r')
 
        ax.pcolormesh(
                hm_u['lon'], hm_u['lat'], calc_abs(ua, va),
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r,
                rasterized=True, transform=ccrs.PlateCarree())


        ax=setup_map(subplot=132)
        plot_LES_domain(lon, lat, 'r')

        ax.pcolormesh(
                lon, lat,
                calc_abs(les_u['uxy'][t_les,:,:].values, les_v['vxy'][t_les,:,:].values).T,
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r,
                rasterized=True, transform=ccrs.PlateCarree())

        ax=setup_map(subplot=133)
        plot_LES_domain(lon, lat, 'r')

        dU = calc_abs(les_u['uxy'][t_les,:,:].values, les_v['vxy'][t_les,:,:].values).T - calc_abs(ua_at_les, va_at_les)

        pc=ax.pcolormesh(
                lon, lat, dU,
                vmin=-6, vmax=6, cmap=pl.cm.RdBu_r,
                rasterized=True, transform=ccrs.PlateCarree())
        pl.colorbar(pc)
