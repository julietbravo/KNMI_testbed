import numpy as np
import datetime
import sys
import os

# Add src directory to Python path, and import DALES specific tools
sys.path.append(os.path.abspath('{}/../src/'.format(os.path.dirname(os.path.abspath(__file__)))))
from DALES_tools import *

def get_file_list(path, starttime, endtime):
    """
    Get list of required DDH NetCDF files to force
    LES from `starttime` to `endtime`
    """

    # For now limited to runs starting at a complete hour, to prevent having
    # to interpolate the inital conditions
    if starttime.minute != 0 or starttime.second != 0:
        raise RuntimeError('Can only create forcings starting at a complete hour!')

    # If experiment starts at start of cycle (t=0,3,6,..) we also need the previous cycle..
    if starttime.hour % 3 == 0:
        starttime -= datetime.timedelta(hours=3)

    # Number of cycles to convert
    n_cycles = int((endtime-starttime).total_seconds() / 3600. / 3.) + 1

    # Create list of cycles to convert
    files   = []
    success = True
    for t in range(n_cycles):
        date    = starttime + t * datetime.timedelta(hours=3)
        in_file = '{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/LES_forcing_{1:04d}{2:02d}{3:02d}{4:02d}.nc'.\
            format(path, date.year, date.month, date.day, date.hour)

        files.append(in_file)

        # Check if file actually exists..
        if not os.path.exists(in_file):
            print('ERROR: Can not find input file {}'.format(in_file))
            success = False

    if not success:
        raise RuntimeError('One or more required files could not be found...')
    else:
        return files


def get_start_end_indices(start, end, time):
    """
    Get indices in `time` that correspond to the requested `start` and `end` times
    """

    t0 = np.abs(np.datetime64(start) - time).argmin()
    t1 = np.abs(np.datetime64(end)   - time).argmin() + 1

    print('Using {} (index {}) to {} (index {})'.format(time[t0], t0, time[t1-1], t1-1))

    return t0, t1


def interpz(z_input, z_output, variable):
    """
    Interpolate (linear) `variable` from input grid (`z_input`) to output grid (`z_output`)
    """
    return np.interp(z_output, z_input, variable)


def interpz_time(z_input, z_output, variable):
    """
    Interpolate time varying `variable` from input grid (`z_input`) to output grid (`z_output`)
    """
    data = np.zeros((variable.shape[0], z_output.size))
    for t in range(variable.shape[0]):
        data[t,:] = interpz(z_input[t,:], z_output, variable[t,:])
    return data


def create_initial_profiles(nc_data, grid, t0, t1, iloc):
    """
    Interpolate Harmonie data onto LES grid,
    and create the `prof.inp` and `scalar.inp` files
    """

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
    write_profiles('prof.inp.001', output, grid.kmax, '')

    # Initial scalar profiles (for microphysics) are zero
    zero = np.zeros(grid.kmax)
    output = odict([('z (m)', grid.z), ('qr (kg kg-1)', zero), ('nr (kg kg-1)', zero)])
    write_profiles('scalar.inp.001', output, grid.kmax, '')


def create_ls_forcings(nc_data, grid, t0, t1, iloc):
    """
    Create all the (partially time dependent) large scale forcings
    """

    # Time in seconds since start of run
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
    write_forcings('ls_flux.inp.001', output_sfc, output_ls, '')

    # Dummy forcings for the microphysics scalars
    write_dummy_forcings('ls_fluxsv.inp.001', 2, grid.z)

    # Also create non-time dependent file (lscale.inp), required by DALES (why?)
    zero = np.zeros_like(grid.z)
    output_ls  = odict([('height', grid.z), ('ug', zero), ('vg', zero), ('wfls', zero), \
                        ('dqtdxls', zero), ('dqtdyls', zero), ('dqtdtls', zero), ('dthldt', zero)])
    write_profiles('lscale.inp.001', output_ls, grid.kmax, '')


def create_nudging_profiles(nc_data, grid, t0, t1, iloc):
    """
    Create the nudging profiles
    """

    # Time in seconds since start of run
    time_sec    = (nc_data.time[t0:t1  ] - nc_data.time[t0]).values / 1e9

    # Vertical profiles
    T  = interpz_time( nc_data['z'][t0:t1, iloc, :], grid.z, nc_data['T' ][t0:t1, iloc, :] )
    u  = interpz_time( nc_data['z'][t0:t1, iloc, :], grid.z, nc_data['u' ][t0:t1, iloc, :] )
    v  = interpz_time( nc_data['z'][t0:t1, iloc, :], grid.z, nc_data['v' ][t0:t1, iloc, :] )
    p  = interpz_time( nc_data['z'][t0:t1, iloc, :], grid.z, nc_data['p' ][t0:t1, iloc, :] )
    qt = interpz_time( nc_data['z'][t0:t1, iloc, :], grid.z, nc_data['q' ][t0:t1, iloc, :] )
    ql = interpz_time( nc_data['z'][t0:t1, iloc, :], grid.z, nc_data['ql'][t0:t1, iloc, :] )
    zero = np.zeros_like(T)
        
    # Nudging factor (0-1) with height; is multiplied with nudging time from namelist
    nudgefac = np.ones_like(T)

    # Conversions from Harmonie -> DALES
    exner  = (p / constants['p0'])**(constants['rd']/constants['cp'])
    theta  = T / exner
    thetal = theta - constants['lv'] / (constants['cp'] * exner) * ql

    # Write to nudge.inp
    output = odict([('z (m)', grid.z), ('factor (-)', nudgefac), ('u (m s-1)', u), ('v (m s-1)', v),\
                    ('w (m s-1)', zero), ('thl (K)', thetal), ('qt (kg kg-1)', qt)])

    write_time_profiles('nudge.inp.001', time_sec, output, grid.kmax, '')
