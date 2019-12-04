import matplotlib.pyplot as pl
import xarray as xr
import numpy as np
import pyproj
from numba import jit

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

# Custom scripts:
import tick_labels as tl


def minmax(arr):
    return np.floor(arr.min()), np.ceil(arr.max())

def calc_abs(a,b):
    return np.sqrt(a**2 + b**2)


def add_features(mask_land=False, ccolor='k'):
    ax = pl.gca() # Get current axis

    # Draw coast lines, resolution = {'10m','50m','110m'}
    ax.coastlines(resolution='10m', linewidth=0.8, color=ccolor)

    # Load country geometries and lakes (for the IJselmeer)
    # from www.naturalearthdata.com and add to axes
    #if mask_land:
    #    land = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor=ccolor, facecolor='#DADEDF', hatch='//', zorder=99)
    #    ax.add_feature(land, edgecolor=ccolor, linewidth=0.8)

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
    ax = pl.subplot(subplot, projection=proj)

    # Add coast lines, country borders, ..
    add_features(mask_land, ccolor)

    # Set spatial extent of map, and add grid lines
    ax.set_extent(extent, ccrs.PlateCarree())

    # Remove minor ticks
    ax.minorticks_off()

    return ax


def add_labels(gridlines=True, xlabels=True, ylabels=True):
    fig = pl.gcf()
    ax  = pl.gca()

    # *must* call draw in order to get the axis boundary used to add ticks:
    fig.canvas.draw()

    # Define gridline locations and draw the lines using cartopy's built-in gridliner:
    xticks = np.arange(-10, 20, 2)
    yticks = np.arange(45, 65, 2)
    if gridlines:
        ax.gridlines(xlocs=xticks, ylocs=yticks, linestyle=':', color='0.3')

    # Label the end-points of the gridlines using the custom tick makers:
    if xlabels:
        ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
        tl.lambert_xticks(ax, list(xticks))

    if ylabels:
        ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
        tl.lambert_yticks(ax, list(yticks))


def plot_LES_domain(lon, lat, color='k'):
    bl = np.s_[0,0]
    br = np.s_[-1,0]
    tl = np.s_[0,-1]
    tr = np.s_[-1,-1]

    pl.plot([lon[bl], lon[br]], [lat[bl], lat[br]], '--', color=color, linewidth=1, transform=ccrs.PlateCarree())
    pl.plot([lon[bl], lon[tl]], [lat[bl], lat[tl]], '--', color=color, linewidth=1, transform=ccrs.PlateCarree())
    pl.plot([lon[tl], lon[tr]], [lat[tl], lat[tr]], '--', color=color, linewidth=1, transform=ccrs.PlateCarree())
    pl.plot([lon[tr], lon[br]], [lat[tr], lat[br]], '--', color=color, linewidth=1, transform=ccrs.PlateCarree())


@jit(nopython=True, nogil=True)
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
            k1 = k0+1

            # Interpolation factors
            f0 = (z[k1,j,i] - zg) / (z[k1,j,i] - z[k0,j,i])
            f1 = 1.-f0

            # Interpolate!
            fld_out[j,i] = f0*fld_in[k0,j,i] + f1*fld_in[k1,j,i]


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
