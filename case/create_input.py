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

# DALES specific tools to read/write in- and output:
from DALES_tools import *

def interp_z(z_Harmonie, z_LES, variable):
    return np.interp(z_LES, z_Harmonie, variable)

def interp_z_time(z_Harmonie, z_LES, variable):
    data = np.zeros((variable.shape[0], z_LES.size))
    for t in range(variable.shape[0]):
        data[t,:] = interp_z(z_Harmonie[t,:], z_LES, variable[t,:])
    return data

def add_zeros(array):
    return np.insert(array, 0, np.zeros(array.shape[1]), axis=0)

def get_files(path, starttime, endtime):

    #if starttime.minute != 0 or starttime.second != 0:
    #    sys.exit('Can only create forcings starting at a complete hour!')

    ## First file to read; if hour%3==0, we also need the previous cycle
    #if starttime.hour % 3 == 0 and starttime.minute == 0:
    #    start = starttime - datetime.timedelta(hours=3)
    #else:
    #    start = starttime

    # Number of cycles to convert
    n_cycles = int((endtime-starttime).total_seconds() / 3600. / 3.) + 1

    # Create list of cycles to convert
    files = []
    for i in range(n_cycles):
        date = start + i * datetime.timedelta(hours=3)
        in_file = '{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/LES_forcing_{1:04d}{2:02d}{3:02d}{4:02d}.nc'.format(path, date.year, date.month, date.day, date.hour)
        files.append(in_file)

    return files


if __name__ == '__main__':
    pl.close('all')

    # Start and endtime of experiments. Currently has to
    # match the 3-hourly Harmonie cycles. TO-DO: fix this :)
    start = datetime.datetime(year=2010, month=2, day=28, hour=9)
    end   = datetime.datetime(year=2010, month=2, day=28, hour=18)

    # Location in NetCDF file
    iloc = 0

    # Path of DDH data. Data structure below is expected to be in format "path/yyyy/mm/dd/hh/"
    path  = '/nobackup/users/stratum/DOWA/LES_forcing'

    # Get list of NetCDF files which need to be processed
    nc_files = get_files(path, start, end)

    # Read DALES namelist
    nl = Read_namelist('namoptions.001')

    # Create vertical grid
    # ====================
    dz0  = 20.   # Grid spacing near surface
    dz1  = 130   # Grid spacing near domain top
    kloc = 80    # Transition level from dz0->dz1
    nloc = 20    # Depth of transition
    grid = Grid_stretched(nl['domain']['kmax'], dz0, kloc, nloc, dz1)
    #grid.plot()

    # Read Harmonie initial conditions & forcings
    input = xr.open_mfdataset(nc_files)

    # Get nearest start/end time in NetCDF files
    # AAAARGGGGHHHH why are there 100 different time types in Python....?!
    # Select times manually for now.
    t0 = 17
    t1 = -1

    print(input.time[t0].values, input.time[t1].values)

    # Create initial profiles
    # =======================
    # Interpolate from Harmonie to LES grid
    p   = interp_z(input['z'][t0,iloc,:], grid.z, input['p' ][t0,iloc,:])
    T   = interp_z(input['z'][t0,iloc,:], grid.z, input['T' ][t0,iloc,:])
    qt  = interp_z(input['z'][t0,iloc,:], grid.z, input['q' ][t0,iloc,:])
    ql  = interp_z(input['z'][t0,iloc,:], grid.z, input['ql'][t0,iloc,:])
    u   = interp_z(input['z'][t0,iloc,:], grid.z, input['u' ][t0,iloc,:])
    v   = interp_z(input['z'][t0,iloc,:], grid.z, input['v' ][t0,iloc,:])
    tke = np.ones(grid.kmax) * 0.1

    # Conversions from Harmonie -> DALES
    exner  = (p / constants['p0'])**(constants['rd']/constants['cp'])
    theta  = T / exner
    thetal = theta - constants['lv'] / (constants['cp'] * exner) * ql

    # Write to prof.inp.001
    output = odict([('z (m)', grid.z), ('thl (K)', thetal), ('qt (kg kg-1)', qt), \
                    ('u (m s-1)', u), ('v (m s-1)', v), ('tke (m2 s-2)', tke)])
    write_profiles('prof.inp.001', output, 'DOWA Harmonie testbed')

    # Initial scalar profiles (for microphysics) are zero
    zero = np.zeros(grid.kmax)
    output = odict([('z (m)', grid.z), ('qr (kg kg-1)', zero), ('nr (kg kg-1)', zero)])
    write_profiles('scalar.inp.001', output, 'DOWA Harmonie testbed')

    # Surface and atmospheric forcings
    # ================================
    time_sec = (input.time[t0:t1] - input.time[t0]).values / 1e9

    # Surface variables
    ps   = input['p_s'][t0:t1,iloc].values
    Ts   = input['T_s'][t0:t1,iloc].values
    qs   = input['q_s'][t0:t1,iloc].values
    dums = np.zeros_like(Ts)    # dummy field with zeros

    # Atmospheric forcings
    dtT  = interp_z_time(input['z'][t0:t1,iloc,:], grid.z, input['dtT_dyn' ][t0:t1,iloc,:])
    dtu  = interp_z_time(input['z'][t0:t1,iloc,:], grid.z, input['dtu_dyn' ][t0:t1,iloc,:])
    dtv  = interp_z_time(input['z'][t0:t1,iloc,:], grid.z, input['dtv_dyn' ][t0:t1,iloc,:])
    dtqv = interp_z_time(input['z'][t0:t1,iloc,:], grid.z, input['dtqv_dyn'][t0:t1,iloc,:])
    dtql = interp_z_time(input['z'][t0:t1,iloc,:], grid.z, input['dtql_dyn'][t0:t1,iloc,:])
    duma = np.zeros_like(dtT)   # dummy field with zeros

    # Conversions atmospheric forcings
    p_tmp = interp_z_time(input['z'][t0:t1,iloc,:], grid.z, input['p'][t0:t1,iloc,:])
    exner = (p_tmp / constants['p0'])**(constants['rd']/constants['cp'])
    # Conversion from dtT->dtth incomplete (pressure contribution missing). Same holds for
    # conversion dtth->dtthl; differentating the theta_l eq. results in another pressure tendency term
    dtth  = dtT / exner
    dtthl = dtth - constants['lv'] / (constants['cp'] * exner) * dtql

    # Write to ls_flux.inp
    output_sfc = odict([('time', time_sec), ('wthl_s', dums), ('wqt_s', dums), \
                        ('p_s', ps), ('T_s', Ts), ('qt_s', qs)])
    output_ls  = odict([('time', time_sec), ('z', grid.z), ('ug', duma), ('vg', duma), \
                        ('dqtdt', dtqv), ('dthldt', dtthl), ('dudt', dtu), ('dvdt', dtv)])
    write_forcings('ls_flux.inp.001', output_sfc, output_ls, 'DOWA Harmonie testbed')

    # Dummy forcings for the microphysics scalars
    write_dummy_forcings('ls_fluxsv.inp.001', 2, grid.z)

    # Also create non-time dependent input file (lscale.inp)
    zero = np.zeros_like(grid.z)
    output_ls  = odict([('height',grid.z), ('ug',zero), ('vg',zero), ('wfls',zero), \
                        ('dqtdxls',zero), ('dqtdyls',zero), ('dqtdtls',zero), ('dthldt',zero)])
    write_profiles('lscale.inp.001', output_ls, 'DOWA Harmonie testbed')
