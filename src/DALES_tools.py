"""
  Some Python tools for handling DALES
  Bart van Stratum (KNMI)
"""

import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as pl
from collections import OrderedDict as odict

# ---------------------------
# "Private" help functions
# ---------------------------

def _get_or_default(dict, name, shape, default_value):
    if name in dict:
        return dict[name]
    else:
        print('No input profile for \"{}\", defaulting values at zero'.format(name))
        return default_value * np.ones(shape)


def _int_or_float_or_str(value):
    """ Helper function: convert a string to int/float/str """
    try:
        if 'true' in value:
            return True
        elif 'false' in value:
            return False
        elif ('.' in value):
            return float(value)
        else:
            return int(float(value))
    except:
        return value.rstrip()


def _convert_value(value):
    """ Helper function: convert namelist value or list """
    if ',' in value:
        value = value.split(',')
        return [_int_or_float_or_str(val) for val in value]
    else:
        return _int_or_float_or_str(value)


# ---------------------------
# DALES constants (from modglobal.f90)
# ---------------------------
constants = dict(rd = 287.04,
                 cp = 1004.,
                 lv = 2.53e6,
                 p0 = 1e5)

# ---------------------------
# Function to read DALES in/output
# ---------------------------

class Read_namelist:
    def __init__(self, namelist_file):
        self.namelist_file = namelist_file

        self.groups = {}   # Dictionary holding all the data
        with open(namelist_file) as f:
            for line in f:
                lstrip = line.strip()
                if (len(lstrip) > 0 and lstrip[0] != "#"):
                    if lstrip[0] == '&':
                        curr_group_name = lstrip[1:].lower()
                        self.groups[curr_group_name] = {}
                    elif ("=" in line):
                        var_name = lstrip.split('=')[0].strip()
                        value = _convert_value(lstrip.split('=')[1])
                        self.groups[curr_group_name][var_name] = value

    def __getitem__(self, name):
        if name in self.groups.keys():
            return self.groups[name]
        else:
            raise RuntimeError('Can\'t find group \"{}\" in .ini file'.format(name))

    def __repr__(self):
        return 'Available groups in {}:\n{}'.format(self.namelist_file, ', '.join(self.groups.keys()))



# ---------------------------
# Function to write the DALES input
# ---------------------------

def write_profiles(file_name, variables, docstring=''):
    """
    Write the prof.inp.xxx input profiles for DALES
    """

    print('Saving {}'.format(file_name))

    f = open(file_name, 'w')

    # Write header (description file)
    if docstring is '':
        f.write('DALES large-scale forcings\n')
    else:
        f.write('{}\n'.format(docstring))

    # Write header (column names)
    for var in variables.keys():
        f.write('{0:^17s} '.format(var))
    f.write('\n')

    # Number of vertical levels
    nlev = list(variables.items())[0][1].size   # yikes

    # Write data
    for k in range(nlev):
        for var in variables.keys():
            f.write('{0:+1.10E} '.format(variables[var][k]))
        f.write('\n')

    f.close()


def write_dummy_forcings(file_name, n_scalars, z):
    """
    Write dummy forcings
    """

    print('Saving {}'.format(file_name))

    f = open(file_name, 'w')

    f.write('\n\n')

    # Surface fluxes (zero)
    f.write('{0:^15s} '.format('time'))
    for i in range(n_scalars):
        f.write('{0:>10s}{1:<8d}'.format('sv', i+1))
    f.write('\n')

    for time in [0,1e6]:
        f.write('{0:+1.10E} '.format(time))
        for i in range(n_scalars):
            f.write('{0:+1.10E} '.format(0))
        f.write('\n')

    # Atmospheric forcings
    f.write('\n')

    for time in [0,1e6]:
        f.write('# {0:+1.10E}\n'.format(time))
        for k in range(z.size):
            f.write('{0:+1.10E} '.format(z[k]))
            for i in range(n_scalars):
                f.write('{0:+1.10E} '.format(0))
            f.write('\n')

    f.close()


def write_forcings(file_name, timedep_sfc, timedep_atm, docstring=''):
    """
    Write the ls_flux.inp.xxx files
    """

    print('Saving {}'.format(file_name))


    f = open(file_name, 'w')

    # Always write something; DALES expects three line header
    if docstring is '':
        f.write('DALES time dependent input\n')
    else:
        f.write('{}\n'.format(docstring))

    # Write surface variables
    f.write('{0:^15s} {1:^15s} {2:^15s} {3:^15s} {4:^15s} {5:^15s}\n'\
        .format('time', 'wthl_s', 'wqt_s', 'T_s', 'qt_s', 'p_s'))
    f.write('{0:^15s} {1:^15s} {2:^15s} {3:^15s} {4:^15s} {5:^15s}\n'\
        .format('(s)', '(K m s-1)', '(kg kg-1 m s-1)', '(K)', '(kg kg-1)', '(Pa)'))

    if timedep_sfc is None:
        # Write a large initial time, so DALES will disable the surface timedep
        f.write('{0:+1.8E} {1:+1.8E} {2:+1.8E} {3:+1.8E} {4:+1.8E} {5:+1.8E}\n'\
            .format(1e16, -1, -1, -1, -1, -1))
    else:
        nt    = timedep_sfc['time'].size
        time  = _get_or_default(timedep_sfc, 'time',  nt, 0)
        wthls = _get_or_default(timedep_sfc, 'wthl_s',nt, 0)
        wqts  = _get_or_default(timedep_sfc, 'wqt_s', nt, 0)
        Ts    = _get_or_default(timedep_sfc, 'T_s',   nt, 0)
        qts   = _get_or_default(timedep_sfc, 'qt_s',  nt, 0)
        ps    = _get_or_default(timedep_sfc, 'p_s',   nt, 0)

        for t in range(nt):
            f.write('{0:+1.8E} {1:+1.8E} {2:+1.8E} {3:+1.8E} {4:+1.8E} {5:+1.8E}\n'\
                .format(time[t], wthls[t], wqts[t], Ts[t], qts[t], ps[t]))

    if timedep_atm is not None:
        time = timedep_atm['time']
        z    = timedep_atm['z']
        nt   = time.size
        nlev = z.size

        ug   = _get_or_default(timedep_atm, 'ug',     [nt,nlev], 0)
        vg   = _get_or_default(timedep_atm, 'vg',     [nt,nlev], 0)
        wls  = _get_or_default(timedep_atm, 'wls',    [nt,nlev], 0)
        dxq  = _get_or_default(timedep_atm, 'dqtdx',  [nt,nlev], 0)
        dyq  = _get_or_default(timedep_atm, 'dqtdy',  [nt,nlev], 0)
        dtq  = _get_or_default(timedep_atm, 'dqtdt',  [nt,nlev], 0)
        dtth = _get_or_default(timedep_atm, 'dthldt', [nt,nlev], 0)
        dtu  = _get_or_default(timedep_atm, 'dudt',   [nt,nlev], 0)
        dtv  = _get_or_default(timedep_atm, 'dvdt',   [nt,nlev], 0)

        # Write atmospheric data
        for t in range(nt):
            f.write('\n')
            # Write header:
            f.write('{0:^19s} {1:^19s} {2:^19s} {3:^19s} {4:^19s} {5:^19s} {6:^19s} {7:^19s} {8:^19s} {9:^19s}\n'\
                .format('z (m)', 'u_g (m s-1)', 'v_g (m s-1)', 'w_ls (m s-1)',
                        'dqtdx (kg kg-1 m-1)', 'dqtdy (kg kg m-1)', 'dqtdt (kg kg-1 s-1)', 'dthldt (K s-1)', 'dudt (m s-2)', 'dvdt (m s-2)'))

            # Write current time:
            f.write('# {0:1.8E}\n'.format(time[t]))

            # Write profiles:
            for k in range(nlev):
                f.write('{0:+1.12E} {1:+1.12E} {2:+1.12E} {3:+1.12E} {4:+1.12E} {5:+1.12E} {6:+1.12E} {7:+1.12E} {8:+1.12E} {9:+1.12E}\n'\
                    .format(z[k], ug[t,k], vg[t,k], wls[t,k], dxq[t,k], dyq[t,k], dtq[t,k], dtth[t,k], dtu[t,k], dtv[t,k]))

    f.close()

# ---------------------------
# Vertical grids
# ---------------------------
class Grid:
    def __init__(self, kmax, dz0):
        self.kmax = kmax
        self.dz0  = dz0

        self.z = np.zeros(kmax)
        self.dz = np.zeros(kmax)
        self.zsize = None

    def plot(self):
        pl.figure()
        pl.title('zsize = {0:.1f} m'.format(self.zsize), loc='left')
        pl.plot(self.dz, self.z, '-x')
        pl.xlabel('dz (m)')
        pl.ylabel('z (m)')

class Grid_equidist(Grid):
    def __init__(self, kmax, dz0):
        Grid.__init__(self, kmax, dz0)

        self.zsize = kmax * dz0
        self.z[:]  = np.arange(dz0/2, self.zsize, dz0)
        self.dz[:] = dz0

class Grid_stretched(Grid):
    def __init__(self, kmax, dz0, nloc1, nbuf1, dz1):
        Grid.__init__(self, kmax, dz0)

        dn         = 1./kmax
        n          = np.linspace(dn, 1.-dn, kmax)
        nloc1     *= dn
        nbuf1     *= dn
        dzdn1      = dz0/dn
        dzdn2      = dz1/dn

        dzdn       = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + np.tanh((n-nloc1)/nbuf1))
        self.dz[:] = dzdn*dn

        stretch    = np.zeros(self.dz.size)

        self.z[0]  = 0.5*self.dz[0]
        stretch[0] = 1.

        for k in range(1, self.kmax):
              self.z [k] = self.z[k-1] + 0.5*(self.dz[k-1]+self.dz[k])
              stretch[k] = self.dz[k]/self.dz[k-1]

        self.zsize = self.z[kmax-1] + 0.5*self.dz[kmax-1]

class Grid_linear_stretched(Grid):
    def __init__(self, kmax, dz0, alpha):
        Grid.__init__(self, kmax, dz0)

        self.dz[:] = dz0 * (1 + alpha)**np.arange(kmax)
        zh         = np.zeros(kmax+1)
        zh[1:]     = np.cumsum(self.dz)
        self.z[:]  = 0.5 * (zh[1:] + zh[:-1])
        self.zsize = zh[-1]


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
        data[t,:] = interp_z(z_input[t,:], z_output, variable[t,:])
    return data


if __name__ == '__main__':
    #
    # Only executed (for debugging/examples) if script is called directly
    #

    pl.close('all')

    if False:
        """
        Demo of the different vertical grids
        """
        ktot = 64
        dz0  = 25

        equidist  = Grid_equidist(ktot, dz0)
        linear    = Grid_linear_stretched(ktot, dz0, 0.005)
        stretched = Grid_stretched(ktot, dz0, 40, 10, 40)

        pl.figure()
        pl.plot(equidist.dz, equidist.z, '-x', label='equidistant')
        pl.plot(linear.dz, linear.z, '-x', label='linear')
        pl.plot(stretched.dz, stretched.z, '-x', label='stretched')
        pl.legend()
        pl.xlabel('dz (m)')
        pl.ylabel('z (m)')

    if False:
        """
        Read DALES namelist
        """
        nl = Read_namelist('namoptions.001')

        print(nl)
        print('variables in domain:', nl['domain'])
        print('itot:', nl['domain']['itot'])

    if False:
        """
        Write initial vertical profiles
        Profiles which aren't specified are initialized at zero.
        """
        nlev  = 64
        zsize = 3200
        dz    = zsize/nlev
        z     = np.arange(0.5*dz, zsize, dz)

        th  = 290 + z * 0.006
        qt  = np.zeros(nlev)
        u   = np.ones(nlev) * 1
        v   = np.ones(nlev) * 2
        tke = np.ones(nlev) * 2

        data = odict({'z':z, 'thl':th, 'qt':qt, 'u':u, 'v':v, 'tke':tke})
        write_profiles('prof.inp.001', data, docstring='Example of initial DALES profiles')

    if False:
        """
        Write the time dependent surface variables,
        and large scale forcings.
        """
        nlev  = 64
        zsize = 3200
        dz    = zsize/nlev
        z     = np.arange(0.5*dz, zsize, dz)

        time     = np.arange(0,7200.01,300)
        thls     = np.ones(time.size)*300
        data_sfc = odict({'time': time, 'thl_s':thls})

        time_ls  = np.arange(0,7200.01, 1800)
        ug       = np.ones((time_ls.size, nlev))*5
        vg       = np.ones((time_ls.size, nlev))*-5
        data_ls  = odict({'time':time_ls, 'z':z, 'u_g':ug, 'v_g':vg})

        write_forcings('ls_flux.inp.001', data_sfc, data_ls, 'Example of DALES forcings')

























#class Netcdf_input:
#    def __init__(self, file_name, nlev):
#
#        # Create new file
#        self.f = nc4.Dataset(file_name, 'w')
#
#        # Create default dimensions
#        self.f.createDimension('z_f', nlev)    # full level height (m)
#        self.f.createDimension('z_h', nlev+1)  # half level height (m)
#        self.f.createDimension('t_s', None)    # time surface (s)
#        self.f.createDimension('t_a', None)    # time atmosphere (s)
#
#    def set_global_attribute(self, name, value):
#        self.f.setncattr(name, value)
#
#    def add_variable(self, data, dtype, dims, attrs):
#        # Create new variable
#        var = self.f.createVariable(attrs['name'], dtype, dims)
#
#        # Set the NetCDF variable attributes (if present)
#        if attrs is not None:
#            var.setncatts(attrs)
#
#        # Write the data
#        var[:] = data
#
#    def close(self):
#        self.f.close()
#
#f = Netcdf_input('input.nc', 64)
#f.set_global_attribute('title', 'DALES initial and time dependent input data')
#f.set_global_attribute('institution', 'KNMI')
#f.set_global_attribute('source', 'Harmonie 40h1.2tg2 DOWA reanalysis')
#
#z = np.arange(25,3200,50)
#attrs = {'name': 'zf', 'standard_name':'full_level_height', 'units':'m'}
#f.add_variable(z, np.double, ('z_f'), attrs)
#
#time = np.arange(0,7200.01,1800)
#attrs = {'name': 't_a', 'standard_name':'time_atmosphere', 'units':'s'}
#f.add_variable(time, np.double, ('t_a'), attrs)
#
#thl = np.ones((5,64))*300
#attrs = {'name': 'thl', 'standard_name':'liquid_water_potential_temperature', 'units':'K'}
#f.add_variable(thl, np.double, ('t_a','z_f'), attrs)
#
#f.close()
