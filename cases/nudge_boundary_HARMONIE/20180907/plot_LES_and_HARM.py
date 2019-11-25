import matplotlib.pyplot as pl
import xarray as xr
import numpy as np
import pyproj
from numba import jit

# Custom code
from hybrid_sigma_grid import Sigma_grid

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

proj_str = '+proj=lcc +lat_1=52.500000 +lat_2=52.500000 +lat_0=52.500000 +lon_0=.000000 +k_0=1.0 +x_0=649536.512574 +y_0=1032883.739533 +a=6371220.000000 +b=6371220.000000'


def add_features(mask_land=False, ccolor='k'):
    ax = pl.gca() # Get current axis

    # Draw coast lines, resolution = {'10m','50m','110m'}
    ax.coastlines(resolution='10m', linewidth=0.8, color=ccolor)

    # Load country geometries and lakes (for the IJselmeer)
    # from www.naturalearthdata.com and add to axes
    if mask_land:
        land = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor=ccolor, facecolor='#DADEDF', hatch='//', zorder=99)
        ax.add_feature(land, edgecolor=ccolor, linewidth=0.8)
    
    countries = cfeature.NaturalEarthFeature(
            category='cultural', name='admin_0_boundary_lines_land', scale='10m', facecolor='none', zorder=100)
    ax.add_feature(countries, edgecolor=ccolor, linewidth=0.8)

    lakes = cfeature.NaturalEarthFeature(
            category='physical', name='lakes', scale='10m', facecolor='none', zorder=100)
    ax.add_feature(lakes, edgecolor=ccolor, linewidth=0.8)


def setup_map(mask_land=False, ccolor='k', extent=[-0.5, 7, 53.0, 57.5], subplot=111):

    # Map projection
    proj = ccrs.LambertConformal(central_longitude=4.9, central_latitude=51.967)

    # Create single axes object in Lambert projection
    #ax = pl.axes(projection=proj)
    ax = pl.subplot(subplot, projection=proj)

    # Add coast lines, country borders, ..
    add_features(mask_land, ccolor)

    # Set spatial extent of map, and add grid lines
    ax.set_extent(extent, ccrs.PlateCarree())

    # Remove minor ticks
    ax.minorticks_off()

    return ax


def calc_abs(a,b):
    return np.sqrt(a**2 + b**2)


def f(fld):
    #return fld
    return gaussian_filter(fld, sigma=30)


def plot_LES_domain(lon, lat, color='k'):

    bl = np.s_[0,0]
    br = np.s_[-1,0]
    tl = np.s_[0,-1]
    tr = np.s_[-1,-1]

    pl.plot([lon[bl], lon[br]], [lat[bl], lat[br]], '--', color=color, linewidth=1, transform=ccrs.PlateCarree())
    pl.plot([lon[bl], lon[tl]], [lat[bl], lat[tl]], '--', color=color, linewidth=1, transform=ccrs.PlateCarree())
    pl.plot([lon[tl], lon[tr]], [lat[tl], lat[tr]], '--', color=color, linewidth=1, transform=ccrs.PlateCarree())
    pl.plot([lon[tr], lon[br]], [lat[tr], lat[br]], '--', color=color, linewidth=1, transform=ccrs.PlateCarree())


def minmax(arr):
    return np.floor(arr.min()), np.ceil(arr.max())


#@jit(nopython=True, nogil=True)
def interp_z(fld_out, fld_in, z, zg, itot, jtot, ktot):
    """
    Interpolate 3D array `fld_in` on 3D heights `z` to 2D field fld_out at fixed height `zg`
    """

    for j in range(jtot):
        for i in range(itot):
            # Find grid level directly below `zg`
            for k in range(ktot):
                if z[k,j,i] > zg:
                    k0 = k-1
                    break

            # Extrapolate if requested height is below first grid level
            if k0 == -1:
                k0 = 0

            # Interpolation factor, and interpolation
            f0 = (z[k,j,i] - zg) / (z[k+1,j,i] - z[k,j,i])
            fld_out[j,i] = f0 * fld_in[k0,j,i] + (1.-f0) * fld_in[k0+1,j,i]
                




if __name__ == '__main__':

    pl.close('all')
    pl.ion()

    #
    # Map projection
    #
    proj_str = '+proj=lcc +lat_1=52.500000 +lat_2=52.500000 +lat_0=52.500000 +lon_0=.000000 \
            +k_0=1.0 +x_0=649536.512574 +y_0=1032883.739533 +a=6371220.000000 +b=6371220.000000'

    proj = pyproj.Proj(proj_str)


    #
    # Read HARMONIE data
    #
    if 'hm_u' not in locals():
        # Read HARMONIE 3D fields
        hm_u = xr.open_dataset('data/HARMONIE/ua.Slev.his.MINI.DOWA_40h12tg2_fERA5_ptF.20180907.nc')
        hm_v = xr.open_dataset('data/HARMONIE/va.Slev.his.MINI.DOWA_40h12tg2_fERA5_ptF.20180907.nc')
        hm_t = xr.open_dataset('data/HARMONIE/ta.Slev.his.MINI.DOWA_40h12tg2_fERA5_ptF.20180907.nc')
        hm_q = xr.open_dataset('data/HARMONIE/hus.Slev.his.MINI.DOWA_40h12tg2_fERA5_ptF.20180907.nc')
        hm_p = xr.open_dataset('data/HARMONIE/ps.his.MINI.DOWA_40h12tg2_fERA5_ptF.20180907.nc')

        # Calculate HARMONIE pressure/height levels
        grid = Sigma_grid('data/H40_65lev.txt')

        # Time in seconds since 00 UTC
        time = [float(hm_u['time'][i] - hm_u['time'][0])/1e9 for i in range(hm_u.dims['time'])]

    itot_hm = hm_u.dims['x']
    jtot_hm = hm_u.dims['y']
    ktot_hm = hm_u.dims['lev']



    #
    # Read DALES data
    #
    k_les = 1  # {1,4,15,26,42,66,89,111}
    f0001_u = xr.open_dataset('data/DALES/crossxy.{0:04d}.u.001.nc'.format(k_les))
    f0001_v = xr.open_dataset('data/DALES/crossxy.{0:04d}.v.001.nc'.format(k_les))
    f0001_q = xr.open_dataset('data/DALES/crossxy.{0:04d}.qt.001.nc'.format(k_les))
    f0001_t = xr.open_dataset('data/DALES/crossxy.{0:04d}.thl.001.nc'.format(k_les))
    k_les -= 1
        
    lon = np.load('data/DALES/lon_les.npy')
    lat = np.load('data/DALES/lat_les.npy')

    z_les = np.loadtxt('data/DALES/prof.inp.001', skiprows=2, usecols=0)

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


    print('z={}'.format(z_les[k_les]))

    # Tes test
    ztest = np.zeros((jtot_hm, itot_hm), dtype=np.float64)
    interp_z(ztest, z_hm[::-1,:,:], z_hm[::-1,:,:], 15., itot_hm, jtot_hm, ktot_hm)  


    if False:
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
        fig.subplots_adjust(left=0.03, bottom=0.05, right=0.98, top=0.88, wspace=0.08)

        ax = setup_map(subplot=121)
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
                length=5, pivot='middle', transform=ccrs.PlateCarree())

        plot_LES_domain(lon, lat, 'r')

        ax = setup_map(subplot=122)
        ax.outline_patch.set_visible(False)

        pc = ax.pcolormesh(
                lon, lat, 
                calc_abs(f0001_u['uxy'][t_les,:,:].values, f0001_v['vxy'][t_les,:,:].values).T, 
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r, 
                rasterized=True, transform=ccrs.PlateCarree()) 

        ev=100
        s = np.s_[ev//2::ev,ev//2::ev]
        pl.barbs(
                lon[s], lat[s], f0001_u['uxy'][t_les,:].values[s].T, f0001_v['vxy'][t_les,:].values[s].T,
                length=5, pivot='middle', transform=ccrs.PlateCarree())

        plot_LES_domain(lon, lat, 'r')

        cax = fig.add_axes([1-0.28, 0.96, 0.25, 0.02]) 
        cb = pl.colorbar(pc, orientation='horizontal', cax=cax) 
        pl.figtext(1-0.29, 0.96, r'$U$ (m s$^{-1}$)', ha='right', fontsize=10)

        pl.savefig('U_{0:.0f}m.pdf'.format(z_les[k_les]))

    
    if False:
        #
        # Specific humidity
        #
        # Interpolate HARMONIE field to LES height
        hus = np.zeros((jtot_hm, itot_hm))

        interp_z(hus, hm_q['hus'][t_hm,::-1,:,:].values, z_hm[::-1,:,:], z_les[k_les], itot_hm, jtot_hm, ktot_hm)  

        if k_les == 0:
            vmin = 5
            vmax=10
        else:
            vmin, vmax = minmax(hus*1e3)

        #
        # Plot!
        #
        fig = pl.figure(figsize=(6, 6))
        fig.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.9)

        ax = setup_map(subplot=111)
        ax.outline_patch.set_visible(False)

        ax.pcolormesh(
                hm_u['lon'], hm_u['lat'], hus*1e3, 
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r, 
                rasterized=True, transform=ccrs.PlateCarree()) 

        pc = ax.pcolormesh(
                lon, lat, f0001_q['qtxy'][t_les,:,:].values.T*1e3,
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r, 
                rasterized=True, transform=ccrs.PlateCarree()) 

        plot_LES_domain(lon, lat, 'k')

        cax = fig.add_axes([1-0.35, 0.96, 0.25, 0.015]) 
        cb = pl.colorbar(pc, orientation='horizontal', cax=cax) 
        pl.figtext(1-0.36, 0.96, r'$q_\mathrm{t}$ (g kg$^{-1}$)', ha='right', fontsize=10)

        pl.savefig('qt_{0:.0f}m.pdf'.format(z_les[k_les]))

    
    if False:
        #
        # Potential temperature
        #
        vmin = 284
        vmax = 289

        pl.figure(figsize=(8,8))
        ax = setup_map(subplot=111)

        exn   = (1e5 / hm_p['ps'][t_hm,:,:].values)**(287.05/1004.)
        hm_th = hm_t['ta'][t_hm,k_hm,:,:].values * exn
        
        ax.pcolormesh(
                hm_u['lon'], hm_u['lat'], hm_th, 
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r, transform=ccrs.PlateCarree()) 

        ax.pcolormesh(
                lon, lat, f0001_t['thlxy'][t,:,:].values.T,
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r, transform=ccrs.PlateCarree()) 

        plot_LES_domain(lon, lat)


    if False:
        #
        # Wind speed
        #
        ua = np.zeros((jtot_hm, itot_hm))
        va = np.zeros((jtot_hm, itot_hm))

        interp_z(ua, hm_u['ua'][t_hm,::-1,:,:].values, z_hm[::-1,:,:], z_les[k_les], itot_hm, jtot_hm, ktot_hm)  
        interp_z(va, hm_v['va'][t_hm,::-1,:,:].values, z_hm[::-1,:,:], z_les[k_les], itot_hm, jtot_hm, ktot_hm)  

        vmin, vmax = minmax(calc_abs(ua, va))

        #
        # Plot!
        #
        pl.figure(figsize=(8,8))
        ax = setup_map(subplot=111)

        ax.pcolormesh(
                hm_u['lon'], hm_u['lat'], 
                calc_abs(ua, va), 
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r, transform=ccrs.PlateCarree()) 

        ax.pcolormesh(
                lon, lat, 
                calc_abs(f0001_u['uxy'][t_les,:,:].values, f0001_v['vxy'][t_les,:,:].values).T, 
                vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r, transform=ccrs.PlateCarree()) 



    if False:
        pl.figure()
        #ax = setup_map(subplot=121)
        #ax.pcolormesh(
        #        hm_u['lon'], hm_u['lat'], 
        #        calc_abs(hm_u['ua'][t_hm,k_hm,:,:].values, hm_v['va'][t_hm,k_hm,:,:].values), 
        #        vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r, transform=ccrs.PlateCarree()) 

        ## Overlay LES domain
        #plot_LES_domain(lon, lat)

        #ax = setup_map(subplot=111)

        #pl.streamplot(
        #        hm_u['lon'].values, hm_u['lat'].values,
        #        hm_u['ua'][t_hm,k_hm,:,:].values, hm_v['va'][t_hm,k_hm,:,:].values,
        #        transform=ccrs.PlateCarree())
        

        #ax.pcolormesh(
        #        hm_u['lon'], hm_u['lat'], 
        #        calc_abs(hm_u['ua'][t_hm,k_hm,:,:].values, hm_v['va'][t_hm,k_hm,:,:].values), 
        #        vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r, transform=ccrs.PlateCarree()) 

        #ax.pcolormesh(
        #        lon, lat, 
        #        calc_abs(f0001_u['uxy'][t,:,:].values, f0001_v['vxy'][t,:,:].values).T, 
        #        vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r, transform=ccrs.PlateCarree()) 

        #plot_LES_domain(lon, lat)


    
    



    if False:
        pl.figure()
        pl.subplot(121)
        pl.streamplot(hm_u['x'].values, hm_u['y'].values, hm_u['ua'][t_hm,k_hm,:,:].values, hm_v['va'][t_hm,k_hm,:,:].values, density=1)
        pl.xlim(500000, 1300000)
        pl.ylim(1000000, 1800000)

        pl.subplot(122)
        pl.streamplot(f0001_u['xt'].values+x0_les, f0001_u['yt'].values+y0_les, f0001_u['uxy'][t_les,:,:].values, f0001_v['vxy'][t_les,:,:].values, density=1)
        pl.xlim(500000, 1300000)
        pl.ylim(1000000, 1800000)


