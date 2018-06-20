"""
Read / process all the required DDH output
Bart van Stratum (KNMI)
"""

import numpy as np
import datetime
import netCDF4 as nc4

import read_DDH as ddh

# Constants, stored in dict, not to mix up Harmonie (ch) & DALES (cd) constants
# Harmonie constants are from arpifs/setup/sucst.F90
ch = dict(grav=9.80665, T0=273.16, Rd=287.06)
cd = dict(p0=1.e5, Rd=287.04, cp=1004.)

class Read_DDH_files:
    def __init__(self, path, t_end, step, dt=60, quiet=False):
        """
        Read / process all DDH files
        """

        # Read first DDH file to get some settings
        f = ddh.DDH_LFA('{0:}DHFDLHARM+{1:04d}'.format(path,1))

        self.nlev = f.attributes['doc']['nlev']       # Number of vertical levels
        self.ndom = f.attributes['doc']['ndom']       # Number of output domains
        self.nt   = int(t_end/step)+1                 # Number of output time steps

        # Create empty arrays to store the individual DDH data
        self.time = np.zeros(self.nt)
        self.datetime = []

        dim3d = (self.nt, self.ndom, self.nlev)
        dim2d = (self.nt, self.ndom)

        # Atmospheric quantities
        # ----------------------
        self.cp   = np.zeros(dim3d)  # Specific heat at const pressure (J kg-1 K-1)
        self.p    = np.zeros(dim3d)  # Pressure (Pa)
        self.ph   = np.zeros(dim3d)  # Half level pressure (Pa)
        self.dp   = np.zeros(dim3d)  # Pressure difference (Pa)
        self.z    = np.zeros(dim3d)  # Geopotential height (m)

        self.u    = np.zeros(dim3d)  # u-component wind (m s-1)
        self.v    = np.zeros(dim3d)  # v-component wind (m s-1)
        self.T    = np.zeros(dim3d)  # Absolute temperature (K)

        self.qv   = np.zeros(dim3d)  # Specific humidity (vapor)   (kg kg-1)
        self.ql   = np.zeros(dim3d)  # Specific humidity (liquid)  (kg kg-1)
        self.qi   = np.zeros(dim3d)  # Specific humidity (ice)     (kg kg-1)
        self.qr   = np.zeros(dim3d)  # Specific humidity (rain)    (kg kg-1)
        self.qs   = np.zeros(dim3d)  # Specific humidity (snow)    (kg kg-1)
        self.qg   = np.zeros(dim3d)  # Specific humidity (graupel) (kg kg-1)

        # Surface quantities
        # ----------------------
        self.H    = np.zeros(dim2d)  # Surface sensible heat flux (W m-2)
        self.LE   = np.zeros(dim2d)  # Surface latent heat flux (W m-2)
        self.Tsk  = np.zeros(dim2d)  # Surface temperature (K)
        self.qsk  = np.zeros(dim2d)  # Surface specific humidity (kg kg-1)

        # Physics, dynamics and total tendencies
        # Units all in "... s-1"
        self.dtu_phy = np.zeros(dim3d)
        self.dtv_phy = np.zeros(dim3d)
        self.dtT_phy = np.zeros(dim3d)

        self.dtu_dyn = np.zeros(dim3d)
        self.dtv_dyn = np.zeros(dim3d)
        self.dtT_dyn = np.zeros(dim3d)

        self.dtu_tot = np.zeros(dim3d)
        self.dtv_tot = np.zeros(dim3d)
        self.dtT_tot = np.zeros(dim3d)

        # Specific humidity tendencies
        self.dtqv_tot = np.zeros(dim3d)
        self.dtql_tot = np.zeros(dim3d)
        self.dtqi_tot = np.zeros(dim3d)
        self.dtqr_tot = np.zeros(dim3d)
        self.dtqs_tot = np.zeros(dim3d)
        self.dtqg_tot = np.zeros(dim3d)

        self.dtqv_dyn = np.zeros(dim3d)
        self.dtql_dyn = np.zeros(dim3d)
        self.dtqi_dyn = np.zeros(dim3d)
        self.dtqr_dyn = np.zeros(dim3d)
        self.dtqs_dyn = np.zeros(dim3d)
        self.dtqg_dyn = np.zeros(dim3d)

        self.dtqv_phy = np.zeros(dim3d)
        self.dtql_phy = np.zeros(dim3d)
        self.dtqi_phy = np.zeros(dim3d)
        self.dtqr_phy = np.zeros(dim3d)
        self.dtqs_phy = np.zeros(dim3d)
        self.dtqg_phy = np.zeros(dim3d)

        # Read all files
        for tt in range(0, t_end+1, step):
            t = int(tt/step)

            if not quiet:
                print('Reading DDH file #{0:3d} (index {1:<3d})'.format(tt,t))

            if tt == 0:
                # Exception for t==0, the instantaneous variables are
                # hidden in the first output file (t=1)

                f = ddh.DDH_LFA('{0:}DHFDLHARM+{1:04d}'.format(path,1))

                # Subtract one minute..
                self.datetime.append(f.attributes['datetime']['forecast_date']-datetime.timedelta(minutes=1))

                self.cp[t,:,:] = f.read_variable('VCP0') * ch['grav']
                self.p [t,:,:] = f.read_variable('VPF0') * ch['grav']
                self.ph[t,:,:] = f.read_variable('VPH0') * ch['grav']
                self.z [t,:,:] = f.read_variable('VZF0')

                # Non-accumulated variables
                self.dp[t,:,:] = f.read_variable('VPP0')

                self.u [t,:,:] = f.read_variable('VUU0') / self.dp[t,:,:]
                self.v [t,:,:] = f.read_variable('VVV0') / self.dp[t,:,:]
                self.T [t,:,:] = f.read_variable('VCT0') / self.dp[t,:,:]/self.cp[t,:,:]

                self.qv[t,:,:] = f.read_variable('VQV0') / self.dp[t,:,:]
                self.ql[t,:,:] = f.read_variable('VQL0') / self.dp[t,:,:]
                self.qi[t,:,:] = f.read_variable('VQI0') / self.dp[t,:,:]
                self.qr[t,:,:] = f.read_variable('VQR0') / self.dp[t,:,:]
                self.qs[t,:,:] = f.read_variable('VQS0') / self.dp[t,:,:]
                self.qg[t,:,:] = f.read_variable('VQG0') / self.dp[t,:,:]

                self.H[t,:]    = f.read_variable('VSHF0')
                self.LE[t,:]   = f.read_variable('VLHF0')
                self.Tsk[t,:]  = f.read_variable('VTSK0')
                self.qsk[t,:]  = f.read_variable('VQSK0')

                # There are no tendencies for t == 0...!
            else:
                f = ddh.DDH_LFA('{0:}DHFDLHARM+{1:04d}'.format(path,tt))

                self.datetime.append(f.attributes['datetime']['forecast_date'])

                self.cp[t,:,:] = f.read_variable('VCP1') * ch['grav']
                self.p [t,:,:] = f.read_variable('VPF1') * ch['grav']
                self.ph[t,:,:] = f.read_variable('VPH1') * ch['grav']
                self.z [t,:,:] = f.read_variable('VZF1')

                # Non-accumulated variables
                self.dp[t,:,:] = f.read_variable('VPP1')

                self.u [t,:,:] = f.read_variable('VUU1') / self.dp[t,:,:]
                self.v [t,:,:] = f.read_variable('VVV1') / self.dp[t,:,:]
                self.T [t,:,:] = f.read_variable('VCT1') / self.dp[t,:,:]/self.cp[t,:,:]

                self.qv[t,:,:] = f.read_variable('VQV1') / self.dp[t,:,:]
                self.ql[t,:,:] = f.read_variable('VQL1') / self.dp[t,:,:]
                self.qi[t,:,:] = f.read_variable('VQI1') / self.dp[t,:,:]
                self.qr[t,:,:] = f.read_variable('VQR1') / self.dp[t,:,:]
                self.qs[t,:,:] = f.read_variable('VQS1') / self.dp[t,:,:]
                self.qg[t,:,:] = f.read_variable('VQG1') / self.dp[t,:,:]

                self.H[t,:]    = f.read_variable('VSHF1')
                self.LE[t,:]   = f.read_variable('VLHF1')
                self.Tsk[t,:]  = f.read_variable('VTSK1')
                self.qsk[t,:]  = f.read_variable('VQSK1')

                # Accumulated tendencies/variables/..
                self.dtu_phy[t,:,:] = f.read_variable('TUUPHY9')
                self.dtv_phy[t,:,:] = f.read_variable('TVVPHY9')
                self.dtT_phy[t,:,:] = f.read_variable('TCTPHY9')

                self.dtu_dyn[t,:,:] = f.read_variable('TUUDYN9')
                self.dtv_dyn[t,:,:] = f.read_variable('TVVDYN9')
                self.dtT_dyn[t,:,:] = f.read_variable('TCTDYN9')

                self.dtu_tot[t,:,:] = f.read_variable('TUUTOT9')
                self.dtv_tot[t,:,:] = f.read_variable('TVVTOT9')
                self.dtT_tot[t,:,:] = f.read_variable('TCTTOT9')

                # Specific humidity tendencies
                self.dtqv_tot[t,:,:] = f.read_variable('TQVTOT9')
                self.dtql_tot[t,:,:] = f.read_variable('TQLTOT9')
                self.dtqi_tot[t,:,:] = f.read_variable('TQITOT9')
                self.dtqr_tot[t,:,:] = f.read_variable('TQRTOT9')
                self.dtqs_tot[t,:,:] = f.read_variable('TQSTOT9')
                self.dtqg_tot[t,:,:] = f.read_variable('TQGTOT9')

                self.dtqv_dyn[t,:,:] = f.read_variable('TQVDYN9')
                self.dtql_dyn[t,:,:] = f.read_variable('TQLDYN9')
                self.dtqi_dyn[t,:,:] = f.read_variable('TQIDYN9')
                self.dtqr_dyn[t,:,:] = f.read_variable('TQRDYN9')
                self.dtqs_dyn[t,:,:] = f.read_variable('TQSDYN9')
                self.dtqg_dyn[t,:,:] = f.read_variable('TQGDYN9')

                self.dtqv_phy[t,:,:] = f.read_variable('TQVPHY9')
                self.dtql_phy[t,:,:] = f.read_variable('TQLPHY9')
                self.dtqi_phy[t,:,:] = f.read_variable('TQIPHY9')
                self.dtqr_phy[t,:,:] = f.read_variable('TQRPHY9')
                self.dtqs_phy[t,:,:] = f.read_variable('TQSPHY9')
                self.dtqg_phy[t,:,:] = f.read_variable('TQGPHY9')

            # Manually calculate time; DDH can't handle times < 1hour
            self.time[t] = tt/60.

        # From Python list to Numpy array..
        self.datetime    = np.array(self.datetime)
        self.hours_since = np.array([(time-datetime.datetime(2010,1,1)).total_seconds()/3600. for time in self.datetime])

        # De-accumulate the tendencies
        self.deaccumulate(self.dtu_phy, step*dt)
        self.deaccumulate(self.dtv_phy, step*dt)
        self.deaccumulate(self.dtT_phy, step*dt)

        self.deaccumulate(self.dtu_dyn, step*dt)
        self.deaccumulate(self.dtv_dyn, step*dt)
        self.deaccumulate(self.dtT_dyn, step*dt)

        self.deaccumulate(self.dtu_tot, step*dt)
        self.deaccumulate(self.dtv_tot, step*dt)
        self.deaccumulate(self.dtT_tot, step*dt)

        self.deaccumulate(self.dtqv_tot, step*dt)
        self.deaccumulate(self.dtql_tot, step*dt)
        self.deaccumulate(self.dtqi_tot, step*dt)
        self.deaccumulate(self.dtqr_tot, step*dt)
        self.deaccumulate(self.dtqs_tot, step*dt)
        self.deaccumulate(self.dtqg_tot, step*dt)

        self.deaccumulate(self.dtqv_dyn, step*dt)
        self.deaccumulate(self.dtql_dyn, step*dt)
        self.deaccumulate(self.dtqi_dyn, step*dt)
        self.deaccumulate(self.dtqr_dyn, step*dt)
        self.deaccumulate(self.dtqs_dyn, step*dt)
        self.deaccumulate(self.dtqg_dyn, step*dt)

        self.deaccumulate(self.dtqv_phy, step*dt)
        self.deaccumulate(self.dtql_phy, step*dt)
        self.deaccumulate(self.dtqi_phy, step*dt)
        self.deaccumulate(self.dtqr_phy, step*dt)
        self.deaccumulate(self.dtqs_phy, step*dt)
        self.deaccumulate(self.dtqg_phy, step*dt)

        # Sum of moisture and moisture tendencies
        self.q = self.qv + self.ql + self.qi + self.qr + self.qs + self.qg

        self.dtq_dyn = self.dtqv_dyn + self.dtql_dyn + self.dtqi_dyn +\
                       self.dtqr_dyn + self.dtqs_dyn + self.dtqg_dyn
        self.dtq_phy = self.dtqv_phy + self.dtql_phy + self.dtqi_phy +\
                       self.dtqr_phy + self.dtqs_phy + self.dtqg_phy
        self.dtq_tot = self.dtqv_tot + self.dtql_tot + self.dtqi_tot +\
                       self.dtqr_tot + self.dtqs_tot + self.dtqg_tot

        # Derived quantities
        self.exner  = (self.p / 1e5)**(ch['Rd'] / self.cp[t,:]) # Exner
        self.th     = self.T / self.exner                       # Potential temperature

        # Check...: sum of dyn+phys
        self.dtu_sum = self.dtu_phy + self.dtu_dyn
        self.dtv_sum = self.dtv_phy + self.dtv_dyn
        self.dtT_sum = self.dtT_phy + self.dtT_dyn
        self.dtq_sum = self.dtq_phy + self.dtq_dyn

        # Check...: offline tendency
        self.dtu_off  = self.calc_tendency(self.u,  step*dt)
        self.dtv_off  = self.calc_tendency(self.v,  step*dt)
        self.dtT_off  = self.calc_tendency(self.T,  step*dt)
        self.dtth_off = self.calc_tendency(self.th, step*dt)
        self.dtq_off  = self.calc_tendency(self.q , step*dt)

        # Check...: potential temperature budget
        #dtexneri  = self.calc_tendency(1./self.exner, step*dt)

        #self.dtth_phy = self.dtT_phy / self.exneri      # INCOMPLETE
        #self.dtth_dyn = self.dtT_dyn / self.exneri      # INCOMPLETE
        #self.dtth_tot = self.dtT_tot / self.exneri + self.T * dtexneri

        #self.dtth_tot_T  = self.exneri * self.dtT_tot
        #self.dtth_tot_pi = self.T * dtexneri


    def calc_tendency(self, array, dt):
        tend = np.zeros_like(array)
        tend[1:,:] = (array[1:,:] - array[:-1,:]) / dt
        return tend


    def deaccumulate(self, array, dt):
        array[1:,:] = (array[1:,:] - array[:-1,:]) / dt


    def to_netcdf(self, file_name):

        def add_variable(file, name, type, dims, ncatts, data):
            v = file.createVariable(name, type, dims, fill_value=nc4.default_fillvals['f4'])
            v.setncatts(ncatts)
            if dims[-1] in ['z', 'zh']:
                v[:] = data[:,::-1]
            else:
                v[:] = data[:]

        # Harmonie's moisture types and more descriptive names
        qtypes = {'qv':'vapor', 'ql':'liquid', 'qi':'ice', 'qr':'rain', 'qs':'snow', 'qg':'graupel'}

        # Create new NetCDF file
        f = nc4.Dataset(file_name, 'w')

        # Set some global attributes
        f.setncattr('title', 'Large-scale forcings LES')
        f.setncattr('institution', 'KNMI')
        f.setncattr('source', 'Harmonie 40h1.2tg2 DOWA reanalysis')

        # Create dimensions
        f.createDimension('time',   self.nt)
        f.createDimension('z',      self.nlev)
        f.createDimension('domain', self.ndom)

        # Dimensions in NetCDF file
        dim3d = ('time', 'domain', 'z')
        dim2d = ('time', 'domain')
        dim1d = ('time')

        # Output data type
        dtype = 'f4'

        # Create spatial/time variables
        add_variable(f, 'time', dtype, dim1d, {'units': 'hours since 2010-01-01 00:00:00', 'long_name': 'time'}, self.hours_since)
        add_variable(f, 'zg',   dtype, dim3d, {'units': 'm',  'long_name': 'Full level geopotential height'}, self.z)
        add_variable(f, 'p',    dtype, dim3d, {'units': 'Pa', 'long_name': 'Full level hydrostatic pressure'}, self.p)

        # Model variables
        add_variable(f, 'T',    dtype, dim3d, {'units': 'K',       'long_name': 'Absolute temperature'}, self.T)
        add_variable(f, 'u',    dtype, dim3d, {'units': 'm s-1',   'long_name': 'Zonal wind'}, self.u)
        add_variable(f, 'v',    dtype, dim3d, {'units': 'm s-1',   'long_name': 'Meridional wind'}, self.v)
        add_variable(f, 'q',    dtype, dim3d, {'units': 'kg kg-1', 'long_name': 'Total specific humidity'}, self.q)

        # Surface variables
        add_variable(f, 'Tsk',  dtype, dim2d, {'units': 'K', 'long_name': 'Absolute surface temperature'}, self.Tsk)
        add_variable(f, 'qsk',  dtype, dim2d, {'units': 'kg kg-1', 'long_name': 'Surface specific humidity'}, self.qsk)
        add_variable(f, 'ps',   dtype, dim2d, {'units': 'Pa', 'long_name': 'Surface pressure'}, self.ph[:,:,-1])

        for qtype,qname in qtypes.items():
            add_variable(f, qtype, dtype, dim3d, {'units': 'kg kg-1', 'long_name': 'Specific humidity ({})'.format(qname)}, getattr(self, qtype))

        # Tendencies
        add_variable(f, 'dtT_phy', dtype, dim3d, {'units': 'K s-1',  'long_name': 'Physics temperature tendency'},  self.dtT_phy)
        add_variable(f, 'dtT_dyn', dtype, dim3d, {'units': 'K s-1',  'long_name': 'Dynamics temperature tendency'}, self.dtT_dyn)
        add_variable(f, 'dtT_tot', dtype, dim3d, {'units': 'K s-1',  'long_name': 'Total temperature tendency'},    self.dtT_tot)

        add_variable(f, 'dtu_phy', dtype, dim3d, {'units': 'm s-2',  'long_name': 'Physics zonal wind tendency'},  self.dtu_phy)
        add_variable(f, 'dtu_dyn', dtype, dim3d, {'units': 'm s-2',  'long_name': 'Dynamics zonal wind tendency'}, self.dtu_dyn)
        add_variable(f, 'dtu_tot', dtype, dim3d, {'units': 'm s-2',  'long_name': 'Total zonal wind tendency'},    self.dtu_tot)

        add_variable(f, 'dtv_phy', dtype, dim3d, {'units': 'm s-2',  'long_name': 'Physics meridional wind tendency'},  self.dtv_phy)
        add_variable(f, 'dtv_dyn', dtype, dim3d, {'units': 'm s-2',  'long_name': 'Dynamics meridional wind tendency'}, self.dtv_dyn)
        add_variable(f, 'dtv_tot', dtype, dim3d, {'units': 'm s-2',  'long_name': 'Total meridional wind tendency'},    self.dtv_tot)

        add_variable(f, 'dtq_phy', dtype, dim3d, {'units': 'kg kg-1 s-1',  'long_name': 'Physics total specific humidity tendency'},  self.dtq_phy)
        add_variable(f, 'dtq_dyn', dtype, dim3d, {'units': 'kg kg-1 s-1',  'long_name': 'Dynamics total specific humidity tendency'}, self.dtq_dyn)
        add_variable(f, 'dtq_tot', dtype, dim3d, {'units': 'kg kg-1 s-1',  'long_name': 'Total total specific humidity tendency'},    self.dtq_tot)

        for qtype,qname in qtypes.items():
            add_variable(f, 'dt{}_phy'.format(qtype),  dtype, dim3d,\
                {'units': 'kg kg-1 s-1', 'long_name': 'Physics specific humidity ({}) tendency'.format(qname)},  getattr(self, 'dt{}_phy'.format(qtype)))
            add_variable(f, 'dt{}_dyn'.format(qtype),  dtype, dim3d,\
                {'units': 'kg kg-1 s-1', 'long_name': 'Dynamics specific humidity ({}) tendency'.format(qname)}, getattr(self, 'dt{}_dyn'.format(qtype)))
            add_variable(f, 'dt{}_tot'.format(qtype),  dtype, dim3d,\
                {'units': 'kg kg-1 s-1', 'long_name': 'Total specific humidity ({}) tendency'.format(qname)},    getattr(self, 'dt{}_tot'.format(qtype)))

        f.close()


if (__name__ == '__main__'):
    import matplotlib.pyplot as pl
    import matplotlib.gridspec as gridspec

    from DDH_domains import create_info_dict

    pl.close('all')

    dt    = 60      # model time step (s)
    t_end = 180     # final file to read (-)
    step  = 10       # interval to read (-)

    if (True):
        # ---------------------
        # Convert DDH to NetCDF
        # ---------------------
        year  = 2010
        month = 2
        day   = 28
        cycle = 6

        data_root = '/nobackup/users/stratum/DOWA/LES_forcing'
        data_path = '{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/'.format(data_root, year, month, day, cycle)

        data = Read_DDH_files(data_path, t_end, step)
        data.to_netcdf('test.nc')






    if (False):
        # ---------------------
        # Budget potential temperature (temperature/pressure)
        # ---------------------
        year  = 2010
        month = 2
        day   = 28

        data_root = '/nobackup/users/stratum/DOWA/LES_forcing'
        data_path1 = '{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/'.format(data_root, year, month, day, 12)
        data_path2 = '{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/'.format(data_root, year, month, day, 15)

        if 'data1' not in locals():
            data1 = Read_DDH_files(data_path1, t_end, step)
            data2 = Read_DDH_files(data_path2, t_end, step)

        k = -1
        n = 2

        # --------------------
        # T + p, plus budget theta
        # --------------------

        pl.figure()
        gs = gridspec.GridSpec(2, 2, height_ratios=[1,2])

        pl.subplot(gs[0,0])
        pl.plot(data1.time[::n]+12, data1.p[::n,k], 'k-')
        pl.plot(data2.time[::n]+15, data2.p[::n,k], 'k-')

        pl.subplot(gs[0,1])
        pl.plot(data1.time[::n]+12, data1.T[::n,k], 'k-')
        pl.plot(data2.time[::n]+15, data2.T[::n,k], 'k-')

        pl.subplot(gs[1,:])
        pl.plot(data1.time[::n]+12, data1.dtth_tot[::n,k]*3600.,   color='C0', label='total')
        pl.plot(data1.time[::n]+12, data1.dtth_tot_T[::n,k]*3600., color='C1', label='dtT')
        pl.plot(data1.time[::n]+12, data1.dtth_tot_e[::n,k]*3600., color='C2', label='dtp')

        pl.plot(data2.time[::n]+15, data2.dtth_tot[::n,k]*3600.)
        pl.plot(data2.time[::n]+15, data2.dtth_tot_T[::n,k]*3600.)
        pl.plot(data2.time[::n]+15, data2.dtth_tot_e[::n,k]*3600.)
        pl.legend()

        # --------------------
        # High/low bound tendencies T
        # --------------------
        pl.figure()
        pl.plot(data2.time[::]+12, data2.dtth_tot[::,k]*3600., label='t=0,1,2,3,... min')
        pl.plot(data2.time[::2]+12, data2.dtth_tot[::2,k]*3600., label='t=0,2,4,6,... min')
        pl.plot(data2.time[1::2]+12, data2.dtth_tot[1::2,k]*3600., label='t=1,3,5,7,... min')
        pl.legend()










    #forecast_dir = '/nobackup/users/stratum/HARMONIE/experiments/DOWA_40h12tg2_fERA5_small/20160705_00/forecast/'
    #forecast_dir = '/nobackup/users/stratum/HARMONIE/experiments/DOWA_40h12tg2_fERA5_small/20170112_12/forecast/'

    #if 'data' not in locals():
    #    forecast_dir = 'data/20100228_09_nodyntend/'
    #    data = Read_DDH_files(forecast_dir, t_end, step, tendencies=False)
    #    data.to_netcdf('20100228_09_nodyntend.nc')

    #for cycle in range(9, 16, 3):
    #    forecast_dir = 'data/20100228_00/{0:02d}/'.format(cycle)
    #    data = Read_DDH_files(forecast_dir, t_end, step)
    #    data.to_netcdf('20100228_{0:02d}.nc'.format(cycle))

    ##forecast_dir = 'data/20100228_00/06_test/'
    ##if 'data' not in locals():
    ##    data = Read_DDH_files(forecast_dir, t_end, step)

    #cc = pl.cm.jet(np.linspace(0,1,data.time.size))

    #t = 29
    #offs = -1

    if (False):
        pl.figure()
        pl.plot(data.time, data.T[:,-2])

    if (False):
        # -------------------------
        # Budgets q
        # -------------------------

        pl.figure(figsize=(8,5))
        gs = gridspec.GridSpec(2, 2, height_ratios=[3,1])

        # -------- potential temperature -------------
        ax=pl.subplot(gs[0,0])
        pl.semilogy(data.dtq_phy[t,:     ]*3600, data.z[t,:], label='phy')
        pl.semilogy(data.dtq_dyn[t,:     ]*3600, data.z[t,:], label='dyn')
        pl.semilogy(data.dtq_sum[t,:     ]*3600, data.z[t,:], label='sum')
        pl.semilogy(data.dtq_off[t+offs,:]*3600, data.z[t,:], 'x', label='offline')
        pl.legend()
        pl.xlabel(r'$\partial_t q$ (kg kg$^{-1}$ h$^{-1}$)')
        pl.ylabel('z (m)')

        ax=pl.subplot(gs[1,0])
        pl.semilogy(data.dtq_off[t+offs,:]-data.dtq_sum[t,:], data.z[t,:])
        pl.legend()
        pl.xlabel(r'$\epsilon$ (kg kg$^{-1}$ s$^{-1}$)')
        pl.ylabel('z (m)')

        gs = gridspec.GridSpec(3, 2, height_ratios=[1,1,1])
        ax=pl.subplot(gs[0,1])
        pl.title('physics', loc='left')
        pl.semilogy(data.dtqv_phy[t,:]*3600, data.z[t,:], label='qv')
        pl.semilogy(data.dtql_phy[t,:]*3600, data.z[t,:], label='ql')
        pl.semilogy(data.dtqi_phy[t,:]*3600, data.z[t,:], label='qi')
        pl.semilogy(data.dtqr_phy[t,:]*3600, data.z[t,:], label='qr')
        pl.semilogy(data.dtqs_phy[t,:]*3600, data.z[t,:], label='qs')
        pl.semilogy(data.dtqg_phy[t,:]*3600, data.z[t,:], label='qg')
        pl.xlabel(r'$\partial_t q$ (kg kg$^{-1}$ h$^{-1}$)')
        pl.legend()

        ax=pl.subplot(gs[1,1])
        pl.title('dynamics', loc='left')
        pl.semilogy(data.dtqv_dyn[t,:]*3600, data.z[t,:], label='qv')
        pl.semilogy(data.dtql_dyn[t,:]*3600, data.z[t,:], label='ql')
        pl.semilogy(data.dtqi_dyn[t,:]*3600, data.z[t,:], label='qi')
        pl.semilogy(data.dtqr_dyn[t,:]*3600, data.z[t,:], label='qr')
        pl.semilogy(data.dtqs_dyn[t,:]*3600, data.z[t,:], label='qs')
        pl.semilogy(data.dtqg_dyn[t,:]*3600, data.z[t,:], label='qg')
        pl.xlabel(r'$\partial_t q$ (kg kg$^{-1}$ h$^{-1}$)')
        pl.legend()

        ax=pl.subplot(gs[2,1])
        pl.semilogy(data.qv[t,:], data.z[t,:], label='qv')
        pl.semilogy(data.ql[t,:], data.z[t,:], label='ql')
        pl.semilogy(data.qi[t,:], data.z[t,:], label='qi')
        pl.semilogy(data.qr[t,:], data.z[t,:], label='qr')
        pl.semilogy(data.qs[t,:], data.z[t,:], label='qs')
        pl.semilogy(data.qg[t,:], data.z[t,:], label='qg')
        pl.xlabel(r'$q$ (kg kg$^{-1}$)')
        pl.legend()

    if (False):
        # -------------------------
        # Budgets theta
        # -------------------------

        pl.figure(figsize=(10,8))
        gs = gridspec.GridSpec(3, 1, height_ratios=[3,1,1])

        # -------- potential temperature -------------
        ax=pl.subplot(gs[0,0])
        pl.semilogy(data.dtth_phy[t,:     ]*3600, data.z[t,:], label='phy')
        pl.semilogy(data.dtth_dyn[t,:     ]*3600, data.z[t,:], label='dyn')
        pl.semilogy(data.dtth_sum[t,:     ]*3600, data.z[t,:], label='sum')
        pl.semilogy(data.dtth_off[t+offs,:]*3600, data.z[t,:], 'x', label='offline')
        pl.legend()
        pl.xlabel(r'$\partial_t \theta$ (K h$^{-1}$)')
        pl.ylabel('z (m)')

        ax=pl.subplot(gs[1,0])
        pl.semilogy(data.dtth_off[t+offs,:]-data.dtth_sum[t,:], data.z[t,:], label='error sum(d+p)')
        pl.semilogy(data.dtth_off[t+offs,:]-data.dtth_tot[t,:], data.z[t,:], label='error total')
        pl.legend()
        pl.xlabel(r'$\epsilon$ (K s$^{-1}$)')
        pl.ylabel('z (m)')

        ax=pl.subplot(gs[2,0])
        pl.semilogy(data.dtth_tot[t,:]*3600.,   data.z[t,:], label='total')
        pl.semilogy(data.dtth_tot_T[t,:]*3600., data.z[t,:], label=r'$\pi^{-1} \partial T/\partial t$')
        pl.semilogy(data.dtth_tot_e[t,:]*3600., data.z[t,:], label=r'$T \partial \pi^{-1}/\partial t$')
        pl.legend()
        pl.xlabel(r'$\partial_t \theta$ (K h$^{-1}$)')
        pl.ylabel('z (m)')

    if (False):
        # -------------------------
        # Budgets u,v,T
        # -------------------------

        pl.figure(figsize=(8,5))
        gs = gridspec.GridSpec(2, 3, height_ratios=[3,1])

        # -------- temperature -------------
        ax=pl.subplot(gs[0,0])
        pl.semilogy(data.dtT_phy[t,:     ]*3600, data.z[t,:], label='phy')
        pl.semilogy(data.dtT_dyn[t,:     ]*3600, data.z[t,:], label='dyn')
        pl.semilogy(data.dtT_sum[t,:     ]*3600, data.z[t,:], label='tot')
        pl.semilogy(data.dtT_off[t+offs,:]*3600, data.z[t,:], 'x', label='offl')
        pl.legend()
        pl.xlabel(r'$\partial_t T$ (K h$^{-1}$)')
        pl.ylabel('z (m)')

        ax=pl.subplot(gs[1,0])
        pl.semilogy(data.dtT_off[t+offs,:]-data.dtT_sum[t,:], data.z[t,:], label='error')
        pl.xlabel(r'$\epsilon$ (K s$^{-1}$)')
        pl.ylabel('z (m)')

        # -------- u-component -------------
        ax=pl.subplot(gs[0,1])
        pl.semilogy(data.dtu_phy[t,:     ]*3600, data.z[t,:], label='phy')
        pl.semilogy(data.dtu_dyn[t,:     ]*3600, data.z[t,:], label='dyn')
        pl.semilogy(data.dtu_sum[t,:     ]*3600, data.z[t,:], label='tot')
        pl.semilogy(data.dtu_off[t+offs,:]*3600, data.z[t,:], 'x', label='offl')
        pl.legend()
        pl.xlabel(r'$\partial_t u$ (m s$^{-1}$ h$^{-1}$)')

        ax=pl.subplot(gs[1,1])
        pl.semilogy(data.dtu_off[t+offs,:]-data.dtu_sum[t,:], data.z[t,:], label='error')
        pl.xlabel(r'$\epsilon$ (m s$^{-1}$ s$^{-1}$)')

        # -------- v-component -------------
        ax=pl.subplot(gs[0,2])
        pl.semilogy(data.dtv_phy[t,:     ]*3600, data.z[t,:], label='phy')
        pl.semilogy(data.dtv_dyn[t,:     ]*3600, data.z[t,:], label='dyn')
        pl.semilogy(data.dtv_sum[t,:     ]*3600, data.z[t,:], label='tot')
        pl.semilogy(data.dtv_off[t+offs,:]*3600, data.z[t,:], 'x', label='offl')
        pl.legend()
        pl.xlabel(r'$\partial_t v$ (m s$^{-1}$ h$^{-1}$)')

        ax=pl.subplot(gs[1,2])
        pl.semilogy(data.dtv_off[t+offs,:]-data.dtv_sum[t,:], data.z[t,:], label='error')
        pl.xlabel(r'$\epsilon$ (m s$^{-1}$ s$^{-1}$)')

        pl.tight_layout()


    #pl.figure()
    #pl.subplot(131)
    #pl.title('phys', loc='left')
    #for t in range(data.time.size):
    #    pl.semilogy(data.dtT_phy[t,:]*3600., data.z[t,:], color=cc[t], label=str(data.time[t]*60))
    #pl.legend()

    #pl.subplot(132)
    #pl.title('dyn', loc='left')
    #for t in range(data.time.size):
    #    pl.semilogy(data.dtT_dyn[t,:]*3600., data.z[t,:], color=cc[t], label=str(data.time[t]*60))
    #pl.legend()

    #pl.subplot(133)
    #pl.title('total', loc='left')
    #for t in range(data.time.size):
    #    pl.semilogy(data.dtT_tot[t,:]*3600., data.z[t,:], color=cc[t], label=str(data.time[t]*60))
    #pl.legend()





