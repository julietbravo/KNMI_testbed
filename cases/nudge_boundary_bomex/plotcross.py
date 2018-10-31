import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4 as nc4
import numpy as np

pl.close('all')

def roomyspine(spines=['left','bottom'], spacing=5):
    ax = pl.gca()
    for loc, spine in ax.spines.items():
        spine.set_position(('outward',spacing))

def add_colorbar(label=''):
    ax = pl.gca()
    fig = pl.gcf()
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size="1%", pad=0.1, pack_start=False)
    fig.add_axes(cax)
    cbar = pl.colorbar(cax=cax)
    cbar.ax.tick_params(labelsize=9) 
    cbar.ax.set_ylabel(label, rotation=90)


f = nc4.Dataset('crossxy.0002.001_full.nc')
time = f.variables['time'][:]
x = f.variables['xt'][:] / 1000.
y = f.variables['yt'][:] / 1000.
thl = f.variables['thlxy'][:]
w   = f.variables['wxy'][:]
u   = f.variables['uxy'][:]

pl.ioff()

for t in range(0,f.dimensions['time'].size,1):
    print(t,time.size)

    pl.figure(figsize=(12.8,7.2))
    ax=pl.subplot(311)
    pl.title(r'$t=$ {0:.1f} h'.format(time[t]/3600.), loc='left')
    pl.pcolormesh(x, y, thl[t], cmap=pl.cm.RdBu_r, vmin=298.65, vmax=298.9)
    pl.ylabel(r'$y$ (km)')
    roomyspine()
    add_colorbar(r'$\theta_l$ (K)')
    ax.set_xticklabels([])

    ax=pl.subplot(312)
    pl.pcolormesh(x, y, u[t], cmap=pl.cm.RdBu_r, vmin=4.5, vmax=9)
    pl.ylabel(r'$y$ (km)')
    roomyspine()
    add_colorbar(r'$u$ (m s$^{-1}$)')
    ax.set_xticklabels([])

    ax=pl.subplot(313)
    pl.pcolormesh(x, y, w[t], cmap=pl.cm.RdBu_r, vmin=-0.6, vmax=0.6)
    pl.xlabel(r'$x$ (km)')
    pl.ylabel(r'$y$ (km)')
    roomyspine()
    add_colorbar(r'$w$ (m s$^{-1}$)')

    pl.tight_layout()
    pl.savefig('figures/fig{0:04d}.png'.format(t))
    pl.close('all')
