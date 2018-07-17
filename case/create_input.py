import numpy as np
import matplotlib.pyplot as pl
import xarray as xr
from collections import OrderedDict as odict

import datetime
import sys
import os
import glob

# Add src directory to Python path
sys.path.append(os.path.abspath('{}/../src/'.format(os.path.dirname(os.path.abspath(__file__)))))

# DALES specific tools to read/write in- and output,
# plus some help functions in tools.py
from DALES_tools import *
from tools import *

if __name__ == '__main__':
    pl.close('all')

    # Location (domain) in NetCDF file
    iloc = 0+12

    # Start and endtime of experiments:
    start = datetime.datetime(year=2016, month=12, day=1, hour=6)
    end   = datetime.datetime(year=2016, month=12, day=1, hour=18)

    print('Creating LES forcings for {} to {}'.format(start, end))

    # Path of DDH data. Data structure below is expected to be in format "path/yyyy/mm/dd/hh/"
    #path  = '/nobackup/users/stratum/DOWA/LES_forcing'
    path  = '/Users/bart/meteo/data/Harmonie_DDH/'

    # Get list of NetCDF files which need to be processed
    nc_files = get_file_list(path, start, end)

    # Read Harmonie DDH NetCDF files
    nc_data = xr.open_mfdataset(nc_files)

    # Name of domain & information to add to ascii files
    domain    = nc_data.name[0][iloc].values
    docstring = 'Harmonie LES forcings, {}: {} to {}'.format(domain, start, end)
    print(docstring)

    # Get start and end indices in `nc_data`
    t0, t1 = get_start_end_indices(start, end, nc_data.time.values)

    # Read DALES namelist
    nl = Read_namelist('namoptions.001')

    # Create vertical grid
    grid = Grid_stretched(nl['domain']['kmax'], dz0=20, nloc1=80, nbuf1=20, dz1=130)
    #grid.plot()

    # Initial vertical profiles
    # -------------------------
    # Interpolate from Harmonie to LES grid
    p   = interpz( nc_data['z'][t0, iloc, :], grid.z, nc_data['p' ][t0, iloc, :] )
    T   = interpz( nc_data['z'][t0, iloc, :], grid.z, nc_data['T' ][t0, iloc, :] )
    qt  = interpz( nc_data['z'][t0, iloc, :], grid.z, nc_data['q' ][t0, iloc, :] )
    ql  = interpz( nc_data['z'][t0, iloc, :], grid.z, nc_data['ql'][t0, iloc, :] )
    u   = interpz( nc_data['z'][t0, iloc, :], grid.z, nc_data['u' ][t0, iloc, :] )
    v   = interpz( nc_data['z'][t0, iloc, :], grid.z, nc_data['v' ][t0, iloc, :] )
    tke = np.ones(grid.kmax) * 0.1

    # Conversions from Harmonie -> DALES
    exner  = (p / constants['p0'])**(constants['rd']/constants['cp'])
    theta  = T / exner
    thetal = theta - constants['lv'] / (constants['cp'] * exner) * ql

    # Write to prof.inp.001
    output = odict([('z (m)', grid.z), ('thl (K)', thetal), ('qt (kg kg-1)', qt), \
                    ('u (m s-1)', u), ('v (m s-1)', v), ('tke (m2 s-2)', tke)])
    write_profiles('prof.inp.001', output, docstring)

    # Initial scalar profiles (for microphysics) are zero
    zero = np.zeros(grid.kmax)
    output = odict([('z (m)', grid.z), ('qr (kg kg-1)', zero), ('nr (kg kg-1)', zero)])
    write_profiles('scalar.inp.001', output, docstring)

    # Surface and atmospheric forcings
    # --------------------------------
    time_sec    = (nc_data.time[t0:t1  ] - nc_data.time[t0]).values / 1e9

    # Forcings are accumulated, so "valid" time is half a dt earlier.
    # Read one extra time step and modify time for input. The forcings for
    # t==0 are corrected later
    time_sec_ls = (nc_data.time[t0:t1+1] - nc_data.time[t0]).values / 1e9
    dt = time_sec[1] - time_sec[0]
    time_sec_ls[1:] -= 0.5*dt

    # Surface variables
    ps     = nc_data['p_s'][t0:t1, iloc].values
    Ts     = nc_data['T_s'][t0:t1, iloc].values
    qs     = nc_data['q_s'][t0:t1, iloc].values
    zero_s = np.zeros_like(Ts)

    # Atmospheric forcings
    dtT    = interpz_time( nc_data['z'][t0:t1+1, iloc, :], grid.z, nc_data['dtT_dyn' ][t0:t1+1, iloc, :] )
    dtu    = interpz_time( nc_data['z'][t0:t1+1, iloc, :], grid.z, nc_data['dtu_dyn' ][t0:t1+1, iloc, :] )
    dtv    = interpz_time( nc_data['z'][t0:t1+1, iloc, :], grid.z, nc_data['dtv_dyn' ][t0:t1+1, iloc, :] )
    dtqv   = interpz_time( nc_data['z'][t0:t1+1, iloc, :], grid.z, nc_data['dtqv_dyn'][t0:t1+1, iloc, :] )
    dtql   = interpz_time( nc_data['z'][t0:t1+1, iloc, :], grid.z, nc_data['dtql_dyn'][t0:t1+1, iloc, :] )
    zero_a = np.zeros_like(dtT)

    # Conversions atmospheric forcings
    p_tmp  = interpz_time(nc_data['z'][t0:t1+1, iloc, :], grid.z, nc_data['p'][t0:t1+1, iloc, :])
    exner  = (p_tmp / constants['p0'])**(constants['rd']/constants['cp'])
    dtth   = dtT / exner
    dtthl  = dtth - constants['lv'] / (constants['cp'] * exner) * dtql

    # Interpolate the forcings for t==0
    dtthl[0, :] = 0.5*(dtthl[0, :] + dtthl[1, :])
    dtu  [0, :] = 0.5*(dtu  [0, :] + dtu  [1, :])
    dtv  [0, :] = 0.5*(dtv  [0, :] + dtv  [1, :])
    dtqv [0, :] = 0.5*(dtqv [0, :] + dtqv [1, :])
    dtql [0, :] = 0.5*(dtql [0, :] + dtql [1, :])

    # Write to ls_flux.inp
    output_sfc = odict([('time', time_sec), ('wthl_s', zero_s), ('wqt_s', zero_s), \
                        ('p_s', ps), ('T_s', Ts), ('qt_s', qs)])
    output_ls  = odict([('time', time_sec_ls), ('z', grid.z), ('ug', zero_a), ('vg', zero_a), \
                        ('dqtdt', dtqv), ('dthldt', dtthl), ('dudt', dtu), ('dvdt', dtv)])
    write_forcings('ls_flux.inp.001', output_sfc, output_ls, docstring)

    # Dummy forcings for the microphysics scalars
    write_dummy_forcings('ls_fluxsv.inp.001', 2, grid.z)

    # Also create non-time dependent file (lscale.inp), required by DALES (why?)
    zero = np.zeros_like(grid.z)
    output_ls  = odict([('height', grid.z), ('ug', zero), ('vg', zero), ('wfls', zero), \
                        ('dqtdxls', zero), ('dqtdyls', zero), ('dqtdtls', zero), ('dthldt', zero)])
    write_profiles('lscale.inp.001', output_ls, docstring)
