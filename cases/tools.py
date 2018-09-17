import numpy as np
import datetime
import sys
import os

# Add src directory to Python path, and import DALES specific tools
sys.path.append(os.path.abspath('{}/../src/'.format(os.path.dirname(os.path.abspath(__file__)))))
from DALES_tools import *

def find_next_dividable_number(number, n):
    if number % n == 0:
        return number
    else:
        return ((number // n)+1)*n


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


def create_initial_profiles(nc_data, grid, t0, t1, iloc, docstring, expnr=1):
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
    output = odict([('z (m)',        grid.z),
                    ('thl (K)',      thetal),
                    ('qt (kg kg-1)', qt),
                    ('u (m s-1)',    u),
                    ('v (m s-1)',    v),
                    ('tke (m2 s-2)', tke)])

    write_profiles('prof.inp.{0:03d}'.format(expnr), output, grid.kmax, docstring)

    # Initial scalar profiles (for microphysics) are zero
    zero = np.zeros(grid.kmax)
    output = odict([('z (m)',        grid.z),
                    ('qr (kg kg-1)', zero),
                    ('nr (kg kg-1)', zero)])

    write_profiles('scalar.inp.{0:03d}'.format(expnr), output, grid.kmax, docstring)


def create_ls_forcings(nc_data, grid, t0, t1, iloc, docstring, n_accumulate=1, expnr=1, harmonie_rad=False):
    """
    Create all the (partially time dependent) large scale forcings
    """

    # Conversion of Harmonie prognostic variables (T) to LES (thl)
    exner       = (nc_data['p'][:,iloc,:] / constants['p0'])**(constants['rd']/constants['cp'])
    dtthl_dyn_f = nc_data['dtT_dyn'][:,iloc,:] / exner - constants['lv'] / (constants['cp'] * exner) * nc_data['dtql_dyn'][:,iloc,:]
    dtthl_rad_f = nc_data['dtT_rad'][:,iloc,:] / exner - constants['lv'] / (constants['cp'] * exner) * nc_data['dtql_dyn'][:,iloc,:]

    if (n_accumulate > 1):
        # Aaargg.. exceptions, exceptions...
        if (n_accumulate == 2):
            pad = 0     # "Valid" time of accumulated tendencies falls exactly at t==0
        else:
            pad = 1     # First valid time of "    " is after t=0; add padding in front to fix t==0 later...

        # Pick end time (t1) such that t1-t0 is dividable by the number of accumulation steps
        n  = find_next_dividable_number(t1-t0, n_accumulate)
        nt = int(n / n_accumulate)
        nz = nc_data.dims['level']
        
        t1 = t0 + n

        z         = np.zeros((nt+pad, nz))
        dtthl_dyn = np.zeros((nt+pad, nz))
        dtthl_rad = np.zeros((nt+pad, nz))
        dtu_dyn   = np.zeros((nt+pad, nz))
        dtv_dyn   = np.zeros((nt+pad, nz))
        dtq_dyn   = np.zeros((nt+pad, nz))

        # Average `n_accumulate` time steps
        z        [pad:,:] = nc_data['z'       ][t0:t1, iloc, :].values.reshape((-1, n_accumulate, nc_data.dims['level'])).mean(axis=1)
        dtthl_dyn[pad:,:] = dtthl_dyn_f        [t0:t1,       :].values.reshape((-1, n_accumulate, nc_data.dims['level'])).mean(axis=1)
        dtthl_rad[pad:,:] = dtthl_rad_f        [t0:t1,       :].values.reshape((-1, n_accumulate, nc_data.dims['level'])).mean(axis=1)
        dtu_dyn  [pad:,:] = nc_data['dtu_dyn' ][t0:t1, iloc, :].values.reshape((-1, n_accumulate, nc_data.dims['level'])).mean(axis=1)
        dtv_dyn  [pad:,:] = nc_data['dtv_dyn' ][t0:t1, iloc, :].values.reshape((-1, n_accumulate, nc_data.dims['level'])).mean(axis=1)
        dtq_dyn  [pad:,:] = nc_data['dtqv_dyn'][t0:t1, iloc, :].values.reshape((-1, n_accumulate, nc_data.dims['level'])).mean(axis=1)
    else:
        z         = nc_data['z'       ][t0:t1, iloc, :].values
        dtthl_dyn = dtthl_dyn_f        [t0:t1,       :].values
        dtthl_rad = dtthl_rad_f        [t0:t1,       :].values
        dtu_dyn   = nc_data['dtu_dyn' ][t0:t1, iloc, :].values
        dtv_dyn   = nc_data['dtv_dyn' ][t0:t1, iloc, :].values
        dtq_dyn   = nc_data['dtqv_dyn'][t0:t1, iloc, :].values

    # Time in seconds since start of run
    time_sec = (nc_data.time[t0:t1  ] - nc_data.time[t0]).values / 1e9

    # Forcings are accumulated, so "valid" time is half a dt earlier..
    dt = time_sec[1] - time_sec[0]      # This is a bit inaccurate....
    time_sec_ls = time_sec - 0.5 * dt

    # Calculate valid time of accumulated tendencies
    if (n_accumulate > 1):
        if (pad == 0):
            time_sec_ls = time_sec_ls.reshape((-1, n_accumulate)).mean(axis=1)
        else:
            tmp = time_sec_ls.copy()
            time_sec_ls = np.zeros(nt+1)
            time_sec_ls[1:] = tmp.reshape((-1, n_accumulate)).mean(axis=1)
    else:
        time_sec_ls[0] = 0

    # Fix t==0 forcings (if necessary)
    if (n_accumulate == 1 or n_accumulate > 2):
        dtthl_dyn[0,:] = 0.5 * (dtthl_dyn_f        [t0-1,     :] + dtthl_dyn_f        [t0,     :])
        dtthl_rad[0,:] = 0.5 * (dtthl_rad_f        [t0-1,     :] + dtthl_rad_f        [t0,     :])
        dtu_dyn  [0,:] = 0.5 * (nc_data['dtu_dyn' ][t0-1,iloc,:] + nc_data['dtu_dyn' ][t0,iloc,:])
        dtv_dyn  [0,:] = 0.5 * (nc_data['dtv_dyn' ][t0-1,iloc,:] + nc_data['dtv_dyn' ][t0,iloc,:])
        dtq_dyn  [0,:] = 0.5 * (nc_data['dtqv_dyn'][t0-1,iloc,:] + nc_data['dtqv_dyn'][t0,iloc,:])

    # Interpolate onto LES grid
    dtthl_dyn = interpz_time(z, grid.z, dtthl_dyn)
    dtthl_rad = interpz_time(z, grid.z, dtthl_rad)
    dtu_dyn   = interpz_time(z, grid.z, dtu_dyn  )
    dtv_dyn   = interpz_time(z, grid.z, dtv_dyn  )
    dtq_dyn   = interpz_time(z, grid.z, dtq_dyn  )
    zero_a    = np.zeros_like(dtthl_dyn)

    if (harmonie_rad):
        print('Adding radiative tendency from Harmonie..')
        dtthl_dyn += dtthl_rad

    # Surface forcings
    time_sec_sfc = time_sec[::n_accumulate]
    ps = nc_data['p_s'][t0:t1:n_accumulate, iloc].values
    Ts = nc_data['T_s'][t0:t1:n_accumulate, iloc].values
    qs = nc_data['q_s'][t0:t1:n_accumulate, iloc].values
    zero_s = np.zeros_like(Ts)

    # Write to ls_flux.inp.expnr
    output_sfc = odict([('time',   time_sec_sfc),
                        ('p_s',    ps          ),
                        ('T_s',    Ts          ),
                        ('qt_s',   qs          )])

    output_ls  = odict([('time',   time_sec_ls),
                        ('z',      grid.z     ),
                        ('dqtdt',  dtq_dyn    ),
                        ('dthldt', dtthl_dyn  ),
                        ('dudt',   dtu_dyn    ),
                        ('dvdt',   dtv_dyn    )])

    write_forcings('ls_flux.inp.{0:03d}'.format(expnr), output_sfc, output_ls, docstring)


    # Dummy forcings for the microphysics scalars
    write_dummy_forcings('ls_fluxsv.inp.{0:03d}'.format(expnr), 2, grid.z, docstring)

    # Also create non-time dependent file (lscale.inp), required by DALES (why?)
    zero = np.zeros_like(grid.z)
    output_ls2  = odict([('height', grid.z), ('ug', zero), ('vg', zero), ('wfls', zero), \
                         ('dqtdxls', zero), ('dqtdyls', zero), ('dqtdtls', zero), ('dthldt', zero)])
    write_profiles('lscale.inp.{0:03d}'.format(expnr), output_ls2, grid.kmax, docstring)

    # Return data from debugging....
    return output_sfc, output_ls



def create_nudging_profiles(nc_data, grid, nudgefactor, t0, t1, iloc, docstring, interval=1, expnr=1):
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
    nudgefac = np.zeros_like(T)
    nudgefac[:,:] = nudgefactor

    # Conversions from Harmonie -> DALES
    exner  = (p / constants['p0'])**(constants['rd']/constants['cp'])
    theta  = T / exner
    thetal = theta - constants['lv'] / (constants['cp'] * exner) * ql

    # Write to nudge.inp
    output = odict([('z (m)',        grid.z), 
                    ('factor (-)',   nudgefac[::interval,:]), 
                    ('u (m s-1)',    u       [::interval,:]),
                    ('v (m s-1)',    v       [::interval,:]),
                    ('w (m s-1)',    zero    [::interval,:]), 
                    ('thl (K)',      thetal  [::interval,:]),
                    ('qt (kg kg-1)', qt      [::interval,:])])

    write_time_profiles('nudge.inp.{0:03d}'.format(expnr), time_sec[::interval], output, grid.kmax, docstring)


def create_backrad(nc_data, t0, iloc, expnr=1):
    """
    Create the background profiles for RRTMG
    """

    print('Saving backrad.inp.{0:03d}.nc'.format(expnr))
        
    nc_file = nc4.Dataset('backrad.inp.{0:03d}.nc'.format(expnr), 'w')
    dims = nc_file.createDimension('lev', nc_data.dims['level'])
    
    p = nc_file.createVariable('lev', 'f4', ('lev'))
    T = nc_file.createVariable('T',   'f4', ('lev'))
    q = nc_file.createVariable('q',   'f4', ('lev'))
    
    p[:] = nc_data['p'][t0, iloc, :]
    T[:] = nc_data['T'][t0, iloc, :]
    q[:] = nc_data['q'][t0, iloc, :]
    
    nc_file.close()
