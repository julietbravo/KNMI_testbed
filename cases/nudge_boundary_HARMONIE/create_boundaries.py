#
# Script to create LES lateral boundaries from HARMONIE 3D fields
# Bart van Stratum (KNMI), Dec. 2018
#

import xarray as xr
import numpy as np

# Custom scripts
import hybrid_sigma_grid as hsg
import interpolate as ip

class Grid:
    """
    Simple version of LES grid, without ghost cells (not needed here),
    """

    def __init__(self, xsize, ysize, zsize, itot, jtot, ktot):

        # Store input/settings
        self.xsize = xsize
        self.ysize = ysize
        self.zsize = zsize

        self.itot = itot
        self.jtot = jtot
        self.ktot = ktot

        # Calculate grid
        self.dx = xsize / itot
        self.dy = ysize / jtot
        self.dz = zsize / ktot

        self.x = np.arange(0.5*self.dx, self.xsize, self.dx)
        self.y = np.arange(0.5*self.dy, self.ysize, self.dy)
        self.z = np.arange(0.5*self.dz, self.zsize, self.dz)

        self.xh = np.arange(0, self.xsize, self.dx)
        self.yh = np.arange(0, self.ysize, self.dy)
        self.zh = np.arange(0, self.zsize, self.dz)



if __name__ == '__main__':
    import matplotlib.pyplot as pl
    pl.close('all')

    # Settings
    x0 = 150000        # Lower left corner LES domain in Harmonie (m)
    y0 = 150000
    
    xsize = 50000    # Domain size LES (m)
    ysize = 50000
    zsize = 3000
    
    itot = 128    # Number of grid points LES (m)
    jtot = 128
    ktot = 128
    
    # LES grid
    grid = Grid(xsize, ysize, zsize, itot, jtot, ktot)
    
    # HARMONIE data
    u  = xr.open_dataset('/Users/bart/meteo/data/DOWA_fulldomain/2010/02/28/00/ua.Slev.his.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc')
    v  = xr.open_dataset('/Users/bart/meteo/data/DOWA_fulldomain/2010/02/28/00/va.Slev.his.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc')
    T  = xr.open_dataset('/Users/bart/meteo/data/DOWA_fulldomain/2010/02/28/00/ta.Slev.his.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc')
    q  = xr.open_dataset('/Users/bart/meteo/data/DOWA_fulldomain/2010/02/28/00/hus.Slev.his.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc')
    ps = xr.open_dataset('/Users/bart/meteo/data/DOWA_fulldomain/2010/02/28/00/ps.his.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc')

    # Select a sub-area around the LES domain, to speed up
    # calcuations done over the entire HARMONIE grid
    def sel_sub_area(ds, x0, y0, xsize, ysize, margin=5000):
        return ds.sel(x=slice(x0-margin, x0+xsize+margin), y=slice(y0-margin, y0+ysize+margin))

    u  = sel_sub_area(u,  x0, y0, xsize, ysize)
    v  = sel_sub_area(v,  x0, y0, xsize, ysize)
    T  = sel_sub_area(T,  x0, y0, xsize, ysize)
    q  = sel_sub_area(q,  x0, y0, xsize, ysize)
    ps = sel_sub_area(ps, x0, y0, xsize, ysize)
    
    # Time index
    t = 0
    
    # Calculate full and half pressure & height HARMONIE grid levels
    grid_sig = hsg.Sigma_grid('data/H40_65lev.txt')

    ph = grid_sig.calc_half_level_pressure(ps['ps'][t,:,:].values)
    zh = grid_sig.calc_half_level_Zg(ph, T['ta'][t,:,:,:].values)   # TO-DO: Tv instead of T!!
    
    p  = grid_sig.calc_full_level_pressure(ph)
    z  = grid_sig.calc_full_level_Zg(zh)

    # Create the interpolator for HARMONIE -> LES
    intp = ip.Grid_interpolator(u['x'].values, u['y'].values, z, grid.x, grid.y, grid.z, grid.xh, grid.yh, grid.zh, x0, y0)

    # Interpolate HARMONIE onto LES grid
    ul = intp.interpolate(u['ua' ][0,::-1,:,:].values, 'xh', 'y', 'z')
    vl = intp.interpolate(v['va' ][0,::-1,:,:].values, 'x', 'yh', 'z')
    Tl = intp.interpolate(T['ta' ][0,::-1,:,:].values, 'x', 'y',  'z')
    ql = intp.interpolate(q['hus'][0,::-1,:,:].values, 'x', 'y',  'z')
    



    if False:
        kLES = 0
        kHM = -1

        pl.figure()
        pl.subplot(321)
        pl.pcolormesh(u['x']/1000., u['y']/1000., u['ua'][t,kHM,:,:], vmin=ul[:,:,kLES].min(), vmax=ul[:,:,kLES].max(), cmap=pl.cm.magma)
        pl.colorbar()
        
        pl.subplot(322)
        pl.pcolormesh(grid.xh/1000., grid.y/1000, ul[:,:,kLES].T, vmin=ul[:,:,kLES].min(), vmax=ul[:,:,kLES].max(), cmap=pl.cm.magma)
        pl.colorbar()

        pl.subplot(323)
        pl.pcolormesh(u['x']/1000., u['y']/1000., v['va'][t,kHM,:,:], vmin=vl[:,:,kLES].min(), vmax=vl[:,:,kLES].max(), cmap=pl.cm.magma)
        pl.colorbar()
        
        pl.subplot(324)
        pl.pcolormesh(grid.x/1000., grid.yh/1000, vl[:,:,kLES].T, vmin=vl[:,:,kLES].min(), vmax=vl[:,:,kLES].max(), cmap=pl.cm.magma)
        pl.colorbar()

        pl.subplot(325)
        pl.pcolormesh(u['x']/1000., u['y']/1000., T['ta'][t,kHM,:,:], vmin=Tl[:,:,kLES].min(), vmax=Tl[:,:,kLES].max(), cmap=pl.cm.magma)
        pl.colorbar()
        
        pl.subplot(326)
        pl.pcolormesh(grid.x/1000., grid.y/1000, Tl[:,:,kLES].T, vmin=Tl[:,:,kLES].min(), vmax=Tl[:,:,kLES].max(), cmap=pl.cm.magma)
        pl.colorbar()
