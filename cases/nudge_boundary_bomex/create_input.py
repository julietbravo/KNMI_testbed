import matplotlib.pyplot as pl
import numpy as np
import xarray as xr
from scipy import interpolate
from shutil import copyfile

pl.close('all')

expnr = 1
kmax  = 75
zsize = 3000
dz    = zsize / kmax

# Uniform vertical grid:
z     = np.linspace(0.5*dz, zsize-0.5*dz, kmax)

if zsize > 3000:
    print('Warning: domain higher than the original BOMEX case!')

# Input variables LES:
thl   = np.zeros_like(z)
qt    = np.zeros_like(z)
u     = np.zeros_like(z)
v     = np.zeros_like(z)
ug    = np.zeros_like(z)
vg    = np.zeros_like(z)
wls   = np.zeros_like(z)
thlls = np.zeros_like(z)
qtls  = np.zeros_like(z)
tke   = np.zeros_like(z)

for k in range(kmax):
    # Liquid water potential temperature
    if(z[k] <= 520.):
        thl[k] = 298.7
    elif(z[k] <= 1480):
        thl[k] = 298.7 + (z[k]-520.)*(302.4-298.7)/(1480.-520.)
    elif(z[k] <= 2000):
        thl[k] = 302.4 + (z[k]-1480.)*(308.2-302.4)/(2000.-1480.)
    else:
        thl[k] = 308.2 + (z[k]-2000.)*(311.85-308.2)/(3000.-2000.)

    # Specific humidity
    if(z[k] <= 520.):
        qt[k] = 1e-3*(17.0 + z[k]*(16.3-17.0)/520.)
    elif(z[k] <= 1480):
        qt[k] = 1.e-3*(16.3 + (z[k]-520.)*(10.7-16.3)/(1480.-520.))
    elif(z[k] <= 2000):
        qt[k] = 1.e-3*(10.7 + (z[k]-1480.)*(4.2-10.7)/(2000.-1480.))
    else:
        qt[k] = 1.e-3*(4.2 + (z[k]-2000.)*(3.-4.2)/(3000.-2000.))


    # u-wind component
    if(z[k] <= 700.):
        u[k] = -8.75
    else:
        u[k] = -8.75 + (z[k]-700.)*(-4.61+8.75)/(3000.-700.)

    # ug-wind component
    ug[k] = -10. + 1.8e-3*z[k]

    # Large scale vertical velocity
    if(z[k] <= 1500):
        wls[k] = z[k]*(-0.65)/1500.
    elif(z[k] <= 2100):
        wls[k] = -0.65 + (z[k]-1500)*(0.65)/(2100.-1500.)

    # Large scale temperature tendency
    if(z[k] <= 1500):
        thlls[k] = (-2.)
    else:
        thlls[k] = (-2.) + (z[k]-1500)*(2.)/(3000.-1500.)

    # Large scale moisture tendency
    if(z[k] <= 300):
        qtls[k] = -1.2
    elif(z[k] <= 500):
        qtls[k] = -1.2 + (z[k]-300)*(1.2)/(500.-300)

    # Initial SGS-TKE
    tke[k] = 1.-z[k]/3000

tke  [tke  <0] = 0.
thlls[thlls>0] = 0.

# Reverse wind direction..
u  = -u
ug = -ug

# Normalize profiles to SI
wls    /= 100.   # from cm/s to m/s
thlls  /= 86400. # from K/d to K/s
qtls   *= 1.e-8  # 1/s

# Write initial profiles to prof.inp
f = open('prof.inp.{0:03d}'.format(expnr), 'w')
f.write('BOMEX, Siebesma et al. (2003)\n')
f.write('{0:^17s} {1:^17s} {2:^17s} {3:^17s} {4:^17s} {5:^17s}\n'.format('zf','thl','qt','u','v','tke'))
for k in range(kmax):
    f.write('{0:+1.10E} {1:+1.10E} {2:+1.10E} {3:+1.10E} {4:+1.10E} {5:+1.10E}\n'.format(z[k], thl[k], qt[k], u[k], v[k], tke[k]))
f.close()

# Write large-scale forcings to lscale.inp
f = open('lscale.inp.{0:03d}'.format(expnr), 'w')
f.write('Large-scale forcings BOMEX, Siebesma et al. (2003)\n')
f.write('{0:^17s} {1:^17s} {2:^17s} {3:^17s} {4:^17s} {5:^17s} {6:^17s} {7:^17s}\n'.format('zf','ug','vg','wfls','dqtdx','dqtdy','dqtdt','dthldt'))
for k in range(kmax):
    f.write('{0:+1.10E} {1:+1.10E} {2:+1.10E} {3:+1.10E} {4:+1.10E} {5:+1.10E} {6:+1.10E} {7:+1.10E}\n'.format(z[k], ug[k], vg[k], wls[k], 0., 0., qtls[k], thlls[k]))
f.close()

# Nudging profiles, to reference BOMEX case
f = xr.open_dataset('original_bomex/profiles.001.nc')

# Interpolate u-component wind to (potentially) new model grid
t  = np.abs(f.time - (24*3600)).argmin()
ip = interpolate.interp1d(f['zt'][:], f['u'][t,:], fill_value='extrapolate')
u  = ip(z)

f = open('nudge.inp.{0:03d}'.format(expnr), 'w')
f.write('BOMEX, Siebesma et al. (2003)\n')
f.write('{0:^17s} {1:^17s} {2:^17s} {3:^17s} {4:^17s} {5:^17s} {6:^17s}\n'.format('zf','f','u','v','w','thl','qt'))
for t in [0,1e9]:
    f.write('# {}\n'.format(t))
    for k in range(kmax):
        f.write('{0:+1.10E} {1:+1.10E} {2:+1.10E} {3:+1.10E} {4:+1.10E} {5:+1.10E} {6:+1.10E}\n'.format(z[k], 1., u[k], 0., 0., 0., 0.))
f.close()

# Copy to other experiment numbers
files = ['prof.inp', 'lscale.inp', 'nudge.inp']
for exp in [2,3,4]:
    for file in files:
        copyfile('{0:}.001'.format(file), '{0:}.{1:03d}'.format(file, exp))
