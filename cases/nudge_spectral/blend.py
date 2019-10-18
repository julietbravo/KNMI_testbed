import numpy as np
import sys

import spectral_tools as st
import DALES_tools as dt

def read_binary(base):
    print('Reading {}'.format(base))

    u = np.zeros((itot, jtot, ktot))
    v = np.zeros((itot, jtot, ktot))
    T = np.zeros((itot, jtot, ktot))
    q = np.zeros((itot, jtot, ktot))
    
    for i in range(nprocx):
        for j in range(nprocy):
            fld = np.fromfile('{0}_x{1:03d}y{2:03d}.001'.format(base,i,j), count=n).reshape(
                    nfld, ktot, blocky, blockx).T

            # Indices of current sub-domain in full 3D fields
            slice = np.s_[i*blockx:(i+1)*blockx, j*blocky:(j+1)*blocky, :]
    
            # Read blocks of data
            u[slice] = fld[:,:,:,0]
            v[slice] = fld[:,:,:,1]
            T[slice] = fld[:,:,:,2]
            q[slice] = fld[:,:,:,3]

    return u,v,T,q


def write_binary(base, u, v, thl, qt):
    print('Writing {}'.format(base))

    for i in range(nprocx):
        for j in range(nprocy):

            name = '{0}_x{1:03d}y{2:03d}.001'.format(base,i,j)
            f = open(name, 'wb+')

            # Indices of current sub-domain in full 3D fields
            slice = np.s_[i*blockx:(i+1)*blockx, j*blocky:(j+1)*blocky, :]

            # Write data to binary file
            np.transpose(u  [slice]).tofile(f)
            np.transpose(v  [slice]).tofile(f)
            np.transpose(thl[slice]).tofile(f)
            np.transpose(qt [slice]).tofile(f)

            f.close()


if __name__ == '__main__':

    # Working directory (in/output of files)
    #path = '/nobackup/users/stratum/KNMI_testbed/cases/nudge_spectral/'
    path = '.'

    # Read namelist
    nl = dt.Read_namelist('namoptions.001')
    
    # Number of fields in binary files (u,v,thl,qt)
    nfld = 4        # = u, v, thl, qt
    
    # Grid dimensions
    itot = nl['domain']['itot']
    jtot = nl['domain']['jtot']
    ktot = nl['domain']['kmax']
    
    # MPI decomposition
    nprocx = nl['run']['nprocx']
    nprocy = nl['run']['nprocy']
    
    # Domain dimensions
    xsize = nl['domain']['xsize']
    ysize = nl['domain']['ysize']
    
    # Grid spacing
    dx = xsize/itot
    dy = ysize/jtot
    
    # Spatial grid (just for plotting)
    x = np.arange(dx/2, xsize, dx)
    y = np.arange(dy/2, ysize, dy)
    
    # Grid size per MPI task
    blockx = itot//nprocx
    blocky = jtot//nprocy
    
    n = blockx*blocky*ktot*nfld
    
    # Get requested time from command line
    time   = float(sys.argv[1])
    minute = int((time%3600)/60)
    hour   = int(time//3600)
    
    # Read 3D LES fields
    u_les, v_les, thl_les, qt_les = read_binary('{0}/fld{1:03d}h{2:02d}m'.format(path, hour, minute))
    
    # Read interpolated HARMONIE fields
    u_hrm, v_hrm, thl_hrm, qt_hrm = read_binary('{0}/lbc{1:03d}h{2:02d}m'.format(path, hour, minute))
    
    # Spectral blending of the fields (smallest wavenumbers HARMONIE + largest LES)
    u_blend   = np.zeros_like(u_les)
    v_blend   = np.zeros_like(u_les)
    thl_blend = np.zeros_like(u_les)
    qt_blend  = np.zeros_like(u_les)
    
    for k in range(ktot):
        print('Blending k={}'.format(k))
        u_blend  [:,:,k] = st.spectral_blend_2d(u_hrm  [:,:,k], u_les  [:,:,k], 5)
        v_blend  [:,:,k] = st.spectral_blend_2d(v_hrm  [:,:,k], v_les  [:,:,k], 5)
        thl_blend[:,:,k] = st.spectral_blend_2d(thl_hrm[:,:,k], thl_les[:,:,k], 5)
        qt_blend [:,:,k] = st.spectral_blend_2d(qt_hrm [:,:,k], qt_les [:,:,k], 5)
    
    # Calculate increments
    u_diff   = u_blend   - u_les
    v_diff   = v_blend   - v_les
    thl_diff = thl_blend - thl_les
    qt_diff  = qt_blend  - qt_les

    # Write increment to binary file
    write_binary('{0}/inc{1:03d}h{2:02d}m'.format(path, hour, minute), u_diff, v_diff, thl_diff, qt_diff)



    if True:
        #
        # Visualisation of blending / increment / ...
        #

        import matplotlib.pyplot as pl
        pl.close('all')
        pl.ion()


        def plot(label, hrm, les, blend, increment, k):

            vmin = min(min(hrm[:,:,k].min(), les[:,:,k].min()), blend[:,:,k].min())
            vmax = max(max(hrm[:,:,k].max(), les[:,:,k].max()), blend[:,:,k].max())

            pl.figure()
    
            ax=pl.subplot(221, aspect='equal')
            pl.title('{} HARMONIE'.format(label), loc='left')
            pl.pcolormesh(x/1e3, y/1e3, hrm[:,:,k].T, vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r)
            pl.colorbar()
            pl.xlabel('x (km)')
            pl.xlabel('y (km)')
    
            pl.subplot(222, aspect='equal', sharex=ax, sharey=ax)
            pl.title('{} LES'.format(label), loc='left')
            pl.pcolormesh(x/1e3, y/1e3, les[:,:,k].T, vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r)
            pl.colorbar()
            pl.xlabel('x (km)')
            pl.xlabel('y (km)')
    
            pl.subplot(223, aspect='equal', sharex=ax, sharey=ax)
            pl.title('{} blended'.format(label), loc='left')
            pl.pcolormesh(x/1e3, y/1e3, blend[:,:,k].T, vmin=vmin, vmax=vmax, cmap=pl.cm.RdBu_r)
            pl.colorbar()
            pl.xlabel('x (km)')
            pl.xlabel('y (km)')
    
            pl.subplot(224, aspect='equal', sharex=ax, sharey=ax)
            pl.title('{} increment'.format(label), loc='left')
            pl.pcolormesh(x/1e3, y/1e3, increment[:,:,k].T, cmap=pl.cm.RdBu_r)
            pl.colorbar()
            pl.xlabel('x (km)')
            pl.xlabel('y (km)')


        plot('thl', thl_hrm, thl_les, thl_blend, thl_diff, k=0)
        plot('qt',  qt_hrm,  qt_les,  qt_blend,  qt_diff,  k=0)
        plot('u',   u_hrm,   u_les,   u_blend,   u_diff,   k=0)
        plot('v',   v_hrm,   v_les,   v_blend,   v_diff,   k=0)

