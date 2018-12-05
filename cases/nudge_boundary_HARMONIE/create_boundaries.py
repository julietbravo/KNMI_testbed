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
    Simple version of LES grid, without ghost cells (not needed here)
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

#def write_LBCs(u, v, thl, qt



if __name__ == '__main__':
    import matplotlib.pyplot as pl
    pl.close('all')

    # Settings
    # Start and end time (index in HARMONIE files)
    t0 = 0
    t1 = 3

    # Lower left corner LES domain in Harmonie (m)
    x0 = 150000
    y0 = 150000

    # Domain size LES (m)
    xsize = 100000
    ysize = 100000
    zsize = 3000

    # Number of grid points LES
    itot = 100
    jtot = 100
    ktot = 128

    # Number of x,y MPI processes
    nprocx = 2
    nprocy = 2

    # Width of the lateral boundaries (# gridpoints)
    n_lbc = 10

    # DALES constants (modglobal.f90)
    cd = dict(p0=1.e5, Rd=287.04, cp=1004., Lv=2.53e6)

    # LES grid
    grid = Grid(xsize, ysize, zsize, itot, jtot, ktot)

    # Hybrid sigma grid tools
    grid_sig = hsg.Sigma_grid('data/H40_65lev.txt')

    # HARMONIE data
    u  = xr.open_dataset('/Users/bart/meteo/data/DOWA_fulldomain/2010/02/28/00/ua.Slev.his.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc')
    v  = xr.open_dataset('/Users/bart/meteo/data/DOWA_fulldomain/2010/02/28/00/va.Slev.his.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc')
    T  = xr.open_dataset('/Users/bart/meteo/data/DOWA_fulldomain/2010/02/28/00/ta.Slev.his.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc')
    q  = xr.open_dataset('/Users/bart/meteo/data/DOWA_fulldomain/2010/02/28/00/hus.Slev.his.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc')
    ql = xr.open_dataset('/Users/bart/meteo/data/DOWA_fulldomain/2010/02/28/00/clw.Slev.his.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc')
    ps = xr.open_dataset('/Users/bart/meteo/data/DOWA_fulldomain/2010/02/28/00/ps.his.NETHERLANDS.DOWA_40h12tg2_fERA5_ptA.20100228.nc')

    # Select a sub-area around the LES domain, to speed up
    # calculations done over the entire HARMONIE grid
    def sel_sub_area(ds, x0, y0, xsize, ysize, margin=2500):
        return ds.sel(x=slice(x0-margin, x0+xsize+margin), y=slice(y0-margin, y0+ysize+margin))

    u  = sel_sub_area(u,  x0, y0, xsize, ysize)['ua' ]
    v  = sel_sub_area(v,  x0, y0, xsize, ysize)['va' ]
    T  = sel_sub_area(T,  x0, y0, xsize, ysize)['ta' ]
    q  = sel_sub_area(q,  x0, y0, xsize, ysize)['hus']
    ql = sel_sub_area(ql, x0, y0, xsize, ysize)['clw']
    ps = sel_sub_area(ps, x0, y0, xsize, ysize)['ps' ]

    # Create hourly boundaries:
    for t in range(t0, t1):
        print('Processing t={0:>2d}/{1:<2d}'.format(t+1, (t1-t0)))

        # Load data from current time step (not really necessary, but otherwise
        # from here on some variable do have a time dimension, and some dont't....)
        # Also, nice time to drop the no longer necessary xarray stuff with `.values`
        u_t  = u [t,:,:,:].values
        v_t  = v [t,:,:,:].values
        T_t  = T [t,:,:,:].values
        q_t  = q [t,:,:,:].values
        ql_t = ql[t,:,:,:].values
        ps_t = ps[t,:,:  ].values

        # Calculate pressure and height on full and half HARMONIE grid levels
        ph = grid_sig.calc_half_level_pressure(ps_t)
        zh = grid_sig.calc_half_level_Zg(ph, T_t)   # TO-DO: Tv instead of T!!
        p  = grid_sig.calc_full_level_pressure(ph)
        z  = grid_sig.calc_full_level_Zg(zh)

        # Conversions HARMONIE quantities -> LES
        exner = (p/cd['p0'])**(cd['Rd']/cd['cp'])
        th_t  = T_t / exner
        thl_t = th_t - cd['Lv'] / (cd['cp'] * exner) * ql_t
        qt_t  = q_t + ql_t

        # Create the interpolator for HARMONIE -> LES
        intp = ip.Grid_interpolator(u['x'].values, u['y'].values, z, grid.x, grid.y, grid.z, grid.xh, grid.yh, grid.zh, x0, y0)

        # Interpolate HARMONIE onto LES grid
        # `::-1` reverses the vertical dimension (HARMONIE's data
        # is aranged from top-to-bottom, LES from bottom-to-top
        u_LES   = intp.interpolate(u_t  [::-1,:,:], 'xh', 'y', 'z')
        v_LES   = intp.interpolate(v_t  [::-1,:,:], 'x', 'yh', 'z')
        thl_LES = intp.interpolate(thl_t[::-1,:,:], 'x', 'y',  'z')
        ql_LES  = intp.interpolate(qt_t [::-1,:,:], 'x', 'y',  'z')


#        # Write lateral boundaries as binary files which can directly be read by DALES
#        # 1. Zonal boundaries
#        #block_x = int(itot / nprocx)
#        #block_y = int(jtot / nprocy)
#        #for j in range(nprocy):
#        #    # 1.1 Left (west) boundary
#        #    sy = np.s_[0:n_lbc, j*block_y:(j+1)*block_y, :]
#




    if False:
        kLES = 0
        kHM = -1

        pl.figure()
        pl.subplot(321)
        pl.pcolormesh((u['x']-x0)/1000., (u['y']-y0)/1000., u['ua'][t,kHM,:,:], vmin=u_LES[:,:,kLES].min(), vmax=u_LES[:,:,kLES].max(), cmap=pl.cm.magma)
        pl.colorbar()
        pl.xlim(0,xsize/1000)
        pl.ylim(0,ysize/1000)

        pl.subplot(322)
        pl.pcolormesh(grid.xh/1000., grid.y/1000, u_LES[:,:,kLES].T, vmin=u_LES[:,:,kLES].min(), vmax=u_LES[:,:,kLES].max(), cmap=pl.cm.magma)
        pl.colorbar()
        pl.xlim(0,xsize/1000)
        pl.ylim(0,ysize/1000)

        pl.subplot(323)
        pl.pcolormesh((u['x']-x0)/1000., (u['y']-y0)/1000., v['va'][t,kHM,:,:], vmin=v_LES[:,:,kLES].min(), vmax=v_LES[:,:,kLES].max(), cmap=pl.cm.magma)
        pl.colorbar()
        pl.xlim(0,xsize/1000)
        pl.ylim(0,ysize/1000)

        pl.subplot(324)
        pl.pcolormesh(grid.x/1000., grid.yh/1000, v_LES[:,:,kLES].T, vmin=v_LES[:,:,kLES].min(), vmax=v_LES[:,:,kLES].max(), cmap=pl.cm.magma)
        pl.colorbar()
        pl.xlim(0,xsize/1000)
        pl.ylim(0,ysize/1000)

        pl.subplot(325)
        pl.pcolormesh((u['x']-x0)/1000., (u['y']-y0)/1000., T['ta'][t,kHM,:,:], vmin=T_LES[:,:,kLES].min(), vmax=T_LES[:,:,kLES].max(), cmap=pl.cm.magma)
        pl.colorbar()
        pl.xlim(0,xsize/1000)
        pl.ylim(0,ysize/1000)

        pl.subplot(326)
        pl.pcolormesh(grid.x/1000., grid.y/1000, T_LES[:,:,kLES].T, vmin=T_LES[:,:,kLES].min(), vmax=T_LES[:,:,kLES].max(), cmap=pl.cm.magma)
        pl.colorbar()
        pl.xlim(0,xsize/1000)
        pl.ylim(0,ysize/1000)
