import matplotlib.pyplot as pl
import xarray as xr
import numpy as np
from itertools import cycle

pl.close('all')

#path = '/nobackup/users/stratum/DOWA/nudge_boundary/nudge_boundary_bomex'
path = 'results'

class Linecycler:
    def __init__(self):
        #self.colors = ['#2A3132','#FF420E','#89DA59','#336B87','#FFBB00']
        self.colors = ['#FF420E','#89DA59','#336B87']
        self.dashes = [(1,0), (4,2), (2,2)]

        self.ic = 0
        self.id = 0

    def next(self):
        color = self.colors[self.ic]
        dash  = self.dashes[self.id]

        self.ic += 1

        if self.ic > len(self.colors)-1:
            self.ic = 0
            self.id += 1

        return color, dash

class Read_span:
    def __init__(self, path, file, t_start, t_end, x_start, x_end):
        print('Reading {}'.format(file))
        self.f = xr.open_dataset('{}/{}'.format(path, file))

        # Period to average
        self.t0 = int(np.abs(self.f.time - t_start).argmin())
        self.t1 = int(np.abs(self.f.time - t_end  ).argmin())

        # Area to average
        self.i0 = int(np.abs(self.f.xt - x_start).argmin())
        self.i1 = int(np.abs(self.f.xt - x_end  ).argmin())

        # Time averaged
        self.um    = self.f['uxz'  ][self.t0:self.t1+1,:,:].mean(axis=0)
        self.vm    = self.f['vxz'  ][self.t0:self.t1+1,:,:].mean(axis=0)
        self.wm    = self.f['wxz'  ][self.t0:self.t1+1,:,:].mean(axis=0)
        self.thlm  = self.f['thlxz'][self.t0:self.t1+1,:,:].mean(axis=0)
        self.qtm   = self.f['qtxz' ][self.t0:self.t1+1,:,:].mean(axis=0)
        self.qlm   = self.f['qlxz' ][self.t0:self.t1+1,:,:].mean(axis=0)

        self.u2m   = self.f['u2rxz'  ][self.t0:self.t1+1].mean(axis=0)
        self.v2m   = self.f['v2rxz'  ][self.t0:self.t1+1].mean(axis=0)
        self.w2m   = self.f['w2rxz'  ][self.t0:self.t1+1].mean(axis=0)
        self.thl2m = self.f['thl2rxz'][self.t0:self.t1+1].mean(axis=0)
        self.qt2m  = self.f['qt2rxz' ][self.t0:self.t1+1].mean(axis=0)
        self.ql2m  = self.f['ql2rxz' ][self.t0:self.t1+1].mean(axis=0)

        # Space (x) averaged
        self.ub    = self.um  [:,self.i0:self.i1+1].mean(axis=1)
        self.vb    = self.vm  [:,self.i0:self.i1+1].mean(axis=1)
        self.wb    = self.wm  [:,self.i0:self.i1+1].mean(axis=1)
        self.thlb  = self.thlm[:,self.i0:self.i1+1].mean(axis=1)
        self.qtb   = self.qtm [:,self.i0:self.i1+1].mean(axis=1)
        self.qlb   = self.qlm [:,self.i0:self.i1+1].mean(axis=1)

        self.u2b   = self.u2m  [:,self.i0:self.i1+1].mean(axis=1)
        self.v2b   = self.v2m  [:,self.i0:self.i1+1].mean(axis=1)
        self.w2b   = self.w2m  [:,self.i0:self.i1+1].mean(axis=1)
        self.thl2b = self.thl2m[:,self.i0:self.i1+1].mean(axis=1)
        self.qt2b  = self.qt2m [:,self.i0:self.i1+1].mean(axis=1)
        self.ql2b  = self.ql2m [:,self.i0:self.i1+1].mean(axis=1)

        self.x = self.f['xt']
        self.z = self.f['zt']


class Read_span_mean:
    def __init__(self, path, file, x_start, x_end, label):
        self.label = label

        print('Reading {}'.format(file))
        self.f = xr.open_dataset('{}/{}'.format(path, file))

        # Area to average
        self.i0 = int(np.abs(self.f.xt - x_start).argmin())
        self.i1 = int(np.abs(self.f.xt - x_end  ).argmin())

        # Time averaged
        self.um    = self.f['uxz'  ][:,:]
        self.vm    = self.f['vxz'  ][:,:]
        self.wm    = self.f['wxz'  ][:,:]
        self.thlm  = self.f['thlxz'][:,:]
        self.qtm   = self.f['qtxz' ][:,:]
        self.qlm   = self.f['qlxz' ][:,:]

        self.u2m   = self.f['u2rxz'  ][:,:]
        self.v2m   = self.f['v2rxz'  ][:,:]
        self.w2m   = self.f['w2rxz'  ][:,:]
        self.thl2m = self.f['thl2rxz'][:,:]
        self.qt2m  = self.f['qt2rxz' ][:,:]
        self.ql2m  = self.f['ql2rxz' ][:,:]

        # Space (x) averaged
        self.ub    = self.um  [:,self.i0:self.i1+1].mean(axis=1)
        self.vb    = self.vm  [:,self.i0:self.i1+1].mean(axis=1)
        self.wb    = self.wm  [:,self.i0:self.i1+1].mean(axis=1)
        self.thlb  = self.thlm[:,self.i0:self.i1+1].mean(axis=1)
        self.qtb   = self.qtm [:,self.i0:self.i1+1].mean(axis=1)
        self.qlb   = self.qlm [:,self.i0:self.i1+1].mean(axis=1)

        self.u2b   = self.u2m  [:,self.i0:self.i1+1].mean(axis=1)
        self.v2b   = self.v2m  [:,self.i0:self.i1+1].mean(axis=1)
        self.w2b   = self.w2m  [:,self.i0:self.i1+1].mean(axis=1)
        self.thl2b = self.thl2m[:,self.i0:self.i1+1].mean(axis=1)
        self.qt2b  = self.qt2m [:,self.i0:self.i1+1].mean(axis=1)
        self.ql2b  = self.ql2m [:,self.i0:self.i1+1].mean(axis=1)

        self.x = self.f['xt']
        self.z = self.f['zt']

if 'ref_sp' not in locals():
    #ref_sp  = Read_span(path, 'crossxzspan.001.nc', 6*3600, 12*3600, 0,     6400)
    #nud_sp  = Read_span(path, 'crossxzspan.002.nc', 6*3600, 12*3600, 45000, 55000)
    #nud_spl = Read_span(path, 'crossxzspan.003.nc', 6*3600, 12*3600, 55000, 75000)

    # Online time averaged + spanwise averaged cross-sections
    ref_sp   = Read_span_mean(path, 'crossxzspan.mean.001.nc', 0,     6400 , r'Reference')
    nud_sp0  = Read_span_mean(path, 'crossxzspan.mean.002.nc', 45000, 55000, r'Laminar')
    #nud_sp1  = Read_span_mean(path, 'crossxzspan.mean.004.nc', 45000, 55000, r'1x1 pert. ($\pm$0.25K)')
    #nud_sp2  = Read_span_mean(path, 'crossxzspan.mean.005.nc', 45000, 55000, r'5x5 pert. ($\pm$0.25k)')
    #nud_sp3  = Read_span_mean(path, 'crossxzspan.mean.006.nc', 45000, 55000, r'10x10 pert. ($\pm$0.25k)')
    nud_sp4  = Read_span_mean(path, 'crossxzspan.mean.007.nc', 45000, 55000, r'1x1 pert. ($\pm$0.05 K)')
    nud_sp5  = Read_span_mean(path, 'crossxzspan.mean.008.nc', 45000, 55000, r'5x5 pert. ($\pm$0.05 K)')
    nud_sp6  = Read_span_mean(path, 'crossxzspan.mean.009.nc', 45000, 55000, r'10x10 pert. ($\pm$0.05 K)')
    nud_sp7  = Read_span_mean(path, 'crossxzspan.mean.010.nc', 45000, 55000, r'1x1 pert. ($\pm \sigma$ K)')
    nud_sp8  = Read_span_mean(path, 'crossxzspan.mean.011.nc', 45000, 55000, r'5x5 pert. ($\pm \sigma$ K)')
    nud_sp9  = Read_span_mean(path, 'crossxzspan.mean.012.nc', 45000, 55000, r'10x10 pert. ($\pm \sigma$ K)')

    #sps = [nud_sp0, nud_sp4, nud_sp5, nud_sp6]      # small amplitude perturbations
    ##sps = [nud_sp0, nud_sp1, nud_sp2, nud_sp3]      # large amplitude perturbations
    #sps = [nud_sp0, nud_sp7, nud_sp8, nud_sp9]      # variance scaling perturbations
    sps = [nud_sp0, nud_sp4, nud_sp5, nud_sp6, nud_sp7, nud_sp8, nud_sp9]      # all

    ref_pr   = xr.open_dataset('{}/profiles.001.nc'.format(path))
    nud_pr0  = xr.open_dataset('{}/profiles.002.nc'.format(path))
    nud_pr1  = xr.open_dataset('{}/profiles.004.nc'.format(path))
    nud_pr2  = xr.open_dataset('{}/profiles.005.nc'.format(path))
    nud_pr3  = xr.open_dataset('{}/profiles.006.nc'.format(path))
    nud_pr4  = xr.open_dataset('{}/profiles.007.nc'.format(path))

    prs = [nud_pr0, nud_pr1, nud_pr2, nud_pr3, nud_pr4]



if (True):
    # ------------------------------
    # Change variance in streamwise direction
    # ------------------------------

    class Variable:
        def __init__(self, name, label, scale=1):
            self.name  = name
            self.label = label
            self.scale = scale

    u2  = Variable('u2',   r'$\sigma^2_{u}$ (m$^2$ s$^{-2}$)')
    v2  = Variable('v2',   r'$\sigma^2_{v}$ (m$^2$ s$^{-2}$)')
    w2  = Variable('w2',   r'$\sigma^2_{w}$ (m$^2$ s$^{-2}$)')
    th2 = Variable('thl2', r'$\sigma^2_{\theta}$ (K$^2$)')
    qt2 = Variable('qt2',  r'$\sigma^2_{q_t}$ (g$^2$ kg$^{-2}$)', 1e6)
    ql2 = Variable('ql2',  r'$\sigma^2_{q_l}$ (g$^2$ kg$^{-2}$)', 1e6)

    #variables = [u2, v2, w2, th2, qt2, ql2]
    variables = [u2, w2, th2]

    heights = np.array([100, 350, 850, 1400])
    indexes = np.zeros_like(heights)
    for i,z in enumerate(heights):
        indexes[i] = np.abs(nud_sp0.z - (z)).argmin()

    lw = 1.1

    x0_ref = 0
    x1_ref = 60

    for v in variables:
        pl.figure(figsize=(10,6))
        for k in range(heights.size):
            pl.subplot(2,2,k+1)
            pl.title('z={}m'.format(heights[k]), loc='left', fontsize='small')
            lines = Linecycler()

            for e,exp in enumerate(sps):
                if e==0:
                    color, dash = 'k', (1,0)
                else:
                    color, dash = lines.next()
                pl.plot(exp.x/1000, getattr(exp, '{}m'.format(v.name))[indexes[k],:]*v.scale, color=color, linewidth=lw, dashes=dash, label=exp.label)
            pl.plot(x1_ref, getattr(ref_sp, '{}b'.format(v.name))[indexes[k]]*v.scale, 'k+', label='Reference')

            pl.ylabel(v.label)
            pl.xlabel(r'$x$ (km)')
            if k == 0:
                pl.legend(fontsize='small', ncol=2, frameon=True, framealpha=0.8, edgecolor='w')
        pl.tight_layout()


if (False):
    # ------------------------------
    # Mixed layer scaling Sorbjan (1989)
    # ------------------------------

    # Average reference case
    ref_mean = ref_pr.sel(time=slice(6*3600, 12*3600)).mean(dim='time')
    z = ref_mean['zt'].values

    thetastar = -0.028571   # Fixed for BOMEX
    zi = 600

    var = (thetastar**2. * (2*(z/zi)**(-2/3.) * (1-(z/zi))**(4/3.) + 0.94*(z/zi)**(4/3.) * (1-(z/zi))**(-2/3.)))**0.5
    


    pl.figure()
    pl.plot(ref_mean['thl2r']**0.5, z, label='LES')
    pl.plot(var, z, label='LES')
    pl.legend()





















if False:
    # ------------------------------
    # Mean quantities from spanwise averaged cross sections
    # ------------------------------

    pl.figure(figsize=(10,6))
    pl.subplot(231)
    pl.plot(ref_sp.ub,  ref_sp.z,  label='Reference')
    pl.plot(nud_sp.ub,  nud_sp.z,  label='Small')
    pl.plot(nud_spl.ub, nud_spl.z, label='Large')
    pl.legend()
    pl.xlabel(r'$u$ (m s$^{-1}$)')
    pl.ylabel(r'$z$ (m)')

    pl.subplot(232)
    pl.plot(ref_sp.vb,  ref_sp.z)
    pl.plot(nud_sp.vb,  nud_sp.z)
    pl.plot(nud_spl.vb, nud_spl.z)
    pl.xlabel(r'$v$ (m s$^{-1}$)')

    pl.subplot(233)
    pl.plot(ref_sp.wb,  ref_sp.z)
    pl.plot(nud_sp.wb,  nud_sp.z)
    pl.plot(nud_spl.wb, nud_spl.z)
    pl.xlabel(r'$w$ (m s$^{-1}$)')

    pl.subplot(234)
    pl.plot(ref_sp.thlb,  ref_sp.z)
    pl.plot(nud_sp.thlb,  nud_sp.z)
    pl.plot(nud_spl.thlb, nud_spl.z)
    pl.xlabel(r'$\theta_l$ (K)')
    pl.ylabel(r'$z$ (m)')

    pl.subplot(235)
    pl.plot(ref_sp.qtb*1e3,  ref_sp.z)
    pl.plot(nud_sp.qtb*1e3,  nud_sp.z)
    pl.plot(nud_spl.qtb*1e3, nud_spl.z)
    pl.xlabel(r'$q_t$ (g kg$^{-1}$)')

    pl.subplot(236)
    pl.plot(ref_sp.qlb*1e3,  ref_sp.z)
    pl.plot(nud_sp.qlb*1e3,  nud_sp.z)
    pl.plot(nud_spl.qlb*1e3, nud_spl.z)
    pl.xlabel(r'$q_l$ (g kg$^{-1}$)')

    pl.tight_layout()


if False:
    # ------------------------------
    # Change variance in streamwise direction
    # ------------------------------
    heights = np.array([100, 350, 850, 1400])
    indexes = np.zeros_like(heights)
    for i,z in enumerate(heights):
        indexes[i] = np.abs(nud_sp.z - (z)).argmin()

    colors = pl.cm.Set1.colors

    x0_ref = 40
    x1_ref = 85

    pl.figure(figsize=(10,6))
    ax=pl.subplot(231)

    for i,k in enumerate(indexes):
        #pl.plot(nud_sp.x/1000,  nud_sp.u2m[k,:],  color=colors[i], label=r'$z$ = {0:4.0f} m'.format(heights[i]))
        pl.plot(nud_spl.x/1000, nud_spl.u2m[k,:], color=colors[i], label=r'$z$ = {0:4.0f} m'.format(heights[i]))
    for i,k in enumerate(indexes):
        label2 = 'reference' if i == 0 else ''
        pl.plot(x1_ref, ref_sp.u2b[k], '+', color=colors[i], markeredgewidth=1.5, label=label2)
    pl.legend()
    pl.ylim(-0.015, 0.4)
    pl.ylabel(r'$\sigma^2_{u}$ (m$^2$ s$^{-2}$)')

    pl.subplot(232)
    for i,k in enumerate(indexes):
        #pl.plot(nud_sp.x/1000,  nud_sp.v2m[k,:], color=colors[i])
        pl.plot(nud_spl.x/1000, nud_spl.v2m[k,:], color=colors[i])
        pl.plot(x1_ref, ref_sp.v2b[k], '+', color=colors[i], markeredgewidth=1.5)
    pl.ylim(-0.015, 0.3)
    pl.ylabel(r'$\sigma^2_{v}$ (m$^2$ s$^{-2}$)')

    pl.subplot(233)
    for i,k in enumerate(indexes):
        #pl.plot(nud_sp.x/1000,  nud_sp.w2m[k,:], color=colors[i])
        pl.plot(nud_spl.x/1000, nud_spl.w2m[k,:], color=colors[i])
        pl.plot(x1_ref, ref_sp.w2b[k], '+', color=colors[i], markeredgewidth=1.5)
    pl.ylim(-0.015, 0.4)
    pl.ylabel(r'$\sigma^2_{w}$ (m$^2$ s$^{-2}$)')

    pl.subplot(234)
    for i,k in enumerate(indexes):
        #pl.plot(nud_sp.x/1000,  nud_sp.thl2m[k,:], color=colors[i])
        pl.plot(nud_spl.x/1000, nud_spl.thl2m[k,:], color=colors[i])
        pl.plot(x1_ref, ref_sp.thl2b[k], '+', color=colors[i], markeredgewidth=1.5)
    pl.ylim(-0.01, 0.15)
    pl.ylabel(r'$\sigma^2_{\theta}$ (K$^2$)')
    pl.xlabel(r'$x$ (km)')

    pl.subplot(235)
    for i,k in enumerate(indexes):
        #pl.plot(nud_sp.x/1000,  nud_sp.qt2m[k,:]*1e6, color=colors[i])
        pl.plot(nud_spl.x/1000, nud_spl.qt2m[k,:]*1e6, color=colors[i])
        pl.plot(x1_ref, ref_sp.qt2b[k]*1e6, '+', color=colors[i], markeredgewidth=1.5)
    pl.ylim(-0.015, 0.5)
    pl.ylabel(r'$\sigma^2_{qt}$ (g$^2$ kg$^{-2}$)')
    pl.xlabel(r'$x$ (km)')

    pl.subplot(236)
    for i,k in enumerate(indexes):
        #pl.plot(nud_sp.x/1000, nud_sp.ql2m[k,:]*1e6, color=colors[i])
        pl.plot(nud_spl.x/1000, nud_spl.ql2m[k,:]*1e6, color=colors[i])
        pl.plot(x1_ref, ref_sp.ql2b[k]*1e6, '+', color=colors[i], markeredgewidth=1.5)
    pl.ylim(-0.0002, 0.02)
    pl.ylabel(r'$\sigma^2_{ql}$ (g$^2$ kg$^{-2}$)')
    pl.xlabel(r'$x$ (km)')

    pl.tight_layout()


if False:
    # ------------------------------
    # Vertical profiles variance
    # ------------------------------
    x_locs  = np.array([15, 20, 25, 30, 35, 40, 50, 60, 70, 80])
    indexes = np.zeros_like(x_locs, dtype=np.int)
    for i,x in enumerate(x_locs):
        indexes[i] = np.abs(nud_spl.x - (x*1000)).argmin()

    colors  = pl.cm.tab20.colors   #(np.linspace(0,1,indexes.size))
    lines   = [[1,0], [4,2]]
    lw = 1.5

    pl.figure(figsize=(10,5))
    pl.subplot(131)
    pl.plot(ref_sp.u2b, ref_sp.z, 'k-', label='Reference', linewidth=2, zorder=999)
    for i,index in enumerate(indexes):
        pl.plot(nud_spl.u2m[:,index], nud_spl.z, color=colors[i], linewidth=lw, dashes=lines[i%len(lines)], label='x ={0:4.1f} km'.format(x_locs[i]))
    pl.legend()
    pl.xlabel(r'$\sigma^2_{u}$ (m$^2$ s$^{-2}$)')
    pl.ylabel(r'$z$ (m)')

    pl.subplot(132)
    pl.plot(ref_sp.v2b, ref_sp.z, 'k-', linewidth=2, zorder=999)
    for i,index in enumerate(indexes):
        pl.plot(nud_spl.v2m[:,index], nud_spl.z, color=colors[i], linewidth=lw, dashes=lines[i%len(lines)])
    pl.xlabel(r'$\sigma^2_{v}$ (m$^2$ s$^{-2}$)')

    pl.subplot(133)
    pl.plot(ref_sp.w2b, ref_sp.z, 'k-', linewidth=2, zorder=999)
    for i,index in enumerate(indexes):
        pl.plot(nud_spl.w2m[:,index], nud_spl.z, color=colors[i], linewidth=lw, dashes=lines[i%len(lines)])
    pl.xlabel(r'$\sigma^2_{w}$ (m$^2$ s$^{-2}$)')

    pl.tight_layout()


    pl.figure(figsize=(10,5))
    pl.subplot(131)
    pl.plot(ref_sp.thl2b, ref_sp.z, 'k-', label='Reference', linewidth=2, zorder=999)
    for i,index in enumerate(indexes):
        pl.plot(nud_spl.thl2m[:,index], nud_spl.z, color=colors[i], linewidth=lw, dashes=lines[i%len(lines)], label='x ={0:4.1f} km'.format(x_locs[i]))
    pl.legend()
    pl.xlabel(r'$\sigma^2_{\theta}$ (K$^2$)')
    pl.ylabel(r'$z$ (m)')

    pl.subplot(132)
    pl.plot(ref_sp.qt2b*1e6, ref_sp.z, 'k-', linewidth=2, zorder=999)
    for i,index in enumerate(indexes):
        pl.plot(nud_spl.qt2m[:,index]*1e6, nud_spl.z, color=colors[i], linewidth=lw, dashes=lines[i%len(lines)])
    pl.xlabel(r'$\sigma^2_{qt}$ (g$^2$ kg$^{-2}$)')

    pl.subplot(133)
    pl.plot(ref_sp.ql2b*1e6, ref_sp.z, 'k-', linewidth=2, zorder=999)
    for i,index in enumerate(indexes):
        pl.plot(nud_spl.ql2m[:,index]*1e6, nud_spl.z, color=colors[i], linewidth=lw, dashes=lines[i%len(lines)])
    pl.xlabel(r'$\sigma^2_{ql}$ (g$^2$ kg$^{-2}$)')
    pl.xlim(0,0.008)

    pl.tight_layout()

if False:
    # ------------------------------
    # Time and spanwise averaged cross-sections
    # ------------------------------
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    def add_colorbar(label=''):
        ax = pl.gca()
        fig = pl.gcf()
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal(size="1%", pad=0.1, pack_start=False)
        fig.add_axes(cax)
        cbar = pl.colorbar(cax=cax)
        cbar.ax.tick_params(labelsize=9) 
        cbar.ax.set_ylabel(label, rotation=90)


    cmap = pl.cm.terrain_r

    pl.figure(figsize=(9,7))
    pl.subplot(511)
    pl.pcolormesh(nud_spl.x/1000, nud_spl.z/1000, nud_spl.qlm*1e3, cmap=cmap)
    pl.ylabel(r'$z$ (km)')
    add_colorbar(r'$q_l$ (g kg$^{-1}$)')

    pl.subplot(512)
    pl.pcolormesh(nud_spl.x/1000, nud_spl.z/1000, nud_spl.u2m, cmap=cmap)
    pl.ylabel(r'$z$ (km)')
    add_colorbar(r'$\sigma^2_u$ (g$^2$ kg$^{-2}$)')

    pl.subplot(513)
    pl.pcolormesh(nud_spl.x/1000, nud_spl.z/1000, nud_spl.v2m, cmap=cmap)
    pl.ylabel(r'$z$ (km)')
    add_colorbar(r'$\sigma^2_v$ (g$^2$ kg$^{-2}$)')

    pl.subplot(514)
    pl.pcolormesh(nud_spl.x/1000, nud_spl.z/1000, nud_spl.w2m, cmap=cmap)
    pl.ylabel(r'$z$ (km)')
    add_colorbar(r'$\sigma^2_w$ (g$^2$ kg$^{-2}$)')

    pl.subplot(515)
    pl.pcolormesh(nud_spl.x/1000, nud_spl.z/1000, nud_spl.thl2m, cmap=cmap)
    pl.ylabel(r'$z$ (km)')
    pl.xlabel(r'$x$ (km)')
    add_colorbar(r'$\sigma^2_\theta$ (K$^2$)')

    pl.tight_layout()
