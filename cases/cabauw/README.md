## HARMONIE LES testbed
This readme file briefly describes the steps needed to setup and run DALES on ECMWF, with the dynamic forcings from HARMONIE-AROME. 

----
### 1. DALES
----
As described in the DOWA/KNMI technical report (https://www.dutchoffshorewindatlas.nl/publications/reports/2019/12/10/knmi-report---downscaling-harmonie-arome-with-large-eddy-simulation), the dynamic forcings for LES consist of the total dynamic tendencies of the liquid water potential temperature, total specific humidity, and both horizontal wind components. The default DALES version can not handle large-scale/dynamic tendencies for momentum, so a slightly modified version of DALES is required, which can be obtained on the ECMWF login nodes using `git` (use `module load git` if required):

    git clone https://github.com/julietbravo/dales.git
    cd dales
    git checkout to4.2_knmi_testbed
    
At some point, these modifications will be merged in the main release of DALES. 

DALES needs to be compiled on the `cca` or `ccb` compute nodes, so first copy the code to the home directory at `cca/ccb` (both compute nodes share the same home directory):

    cd ..
    scp -r dales cca:~/

Next, login at `cca` or `ccb` and go the the `dales` directory:

    ssh cca
    cd dales
    
We will compile DALES inside a build directory, this will allow you to have multiple compiled versions side by side (e.g. a separate release and debug version):

    mkdir build_release
    cd build_release
    
The `to4.2_knmi_testbed` branch is configured for the ECMWF environment & compilers, if the environment variable `SYST` is set to `ECMWF-intel`:

    export SYST=ECMWF-intel

We will compile DALES with the Intel compilers, so switch the compiler environment to Intel (this will automatically load the required MPI libraries):

    prgenvswitchto intel
    
And finally load the `NetCDF4` library:

    module load netcdf4
    
Next, generate the `makefile` with `cmake`, and compile the code:

    cmake ..
    make -j 4
    
Without errors, this should produce the `dales4` executable in the `build_release/src` directory.

----
### 2. Dynamic tendencies HARMONIE & ERA5 data
----

The dynamic tendencies are available (2016 + 2017) for 12 locations, as shown in the figure below. 

![](https://i.stack.imgur.com/EwFEVl.png)

The data is stored in the ECMWF tape archive under `ec:/nkl/harmonie/DOWA/DOWA_40h12tg2_fERA5/LES_forcing/`, as compressed monthly files. Copy the required files to some location at `/scratch`:

    cd $SCRATCH
    mkdir LES_forcing
    cd LES_forcing
    ecp ec:/nkl/harmonie/DOWA/DOWA_40h12tg2_fERA5/LES_forcing/LES_forcings_201608.tar.gz .
    tar -zxvf LES_forcings_201608.tar.gz

The land surface is initialised from ERA5 data, which can easily be retrieved from MARS with a simple script (here called `request.mars`):

    retrieve,
    class=ea,
    expver=1,
    stream=oper,
    date=2016-08-01/to/2016-08-31,
    area=52.971/3.927/50.971/5.927,
    grid=0.3/0.3,
    levtype=sfc,
    type=an,
    time=0/to/23/by/1,
    param=39.128/40.128/41.128/42.128/139.128/170.128/183.128/236.128/43.128,
    target="soil_201608.grib"

which is passed to MARS as:

    mars request.mars
    
The resulting GRIB file can be converted to NetCDF with `grib_to_netcdf`:

    grib_to_netcdf soil_201608.grib -o soil_201608.nc

That should provide all the input files required to run the DALES testbed!

----
### 3. Testbed setup
----

The testbed code is also hosted at github, and therefore has to be downloaded from the login nodes. Exit cca, clone the git repository, copy the files to cca/ccb, and login to cca/ccb again:

    exit
    cd ~
    git clone https://github.com/julietbravo/KNMI_testbed.git
    scp -r KNMI_testbed cca:~/
    ssh cca
    
The Cabauw case is located in the `cases/cabauw` subdirectory:

    cd KNMI_testbed/cases/cabauw
    
The main script used to drive the testbed is `create_input.py`. Modify the script so the paths point to the LES forcings (`path`), ERA5 soil data (`path_e5`), and output directory for the experiments (`path_out`). The variable `iloc` defines which location is used:
    
    # Single column:
    0  = FINO1,        lat=54.01, lon=6.59
    1  = Goeree,       lat=51.93, lon=3.67
    2  = Europlatform, lat=52.00, lon=3.27
    3  = K13,          lat=53.22, lon=3.22
    4  = HKZ,          lat=52.30, lon=4.10
    5  = P11B,         lat=52.36, lon=3.34
    6  = F3-FB-1,      lat=54.85, lon=4.70
    7  = Cabauw,       lat=51.97, lon=4.90
    8  = Loobos,       lat=52.17, lon=5.74
    9  = Lutjewad,     lat=53.40, lon=6.35
    10 = Schiphol,     lat=52.31, lon=4.76
    11 = Rotterdam,    lat=51.91, lon=4.47
    
    # 10x10 & 30x30 km mean:
    12 / 24 = FINO1
    13 / 25 = Goeree
    14 / 26 = Europlatform
    15 / 27 = K13
    16 / 28 = HKZ
    17 / 29 = P11B
    18 / 30 = F3-FB-1
    19 / 31 = Cabauw
    20 / 32 = Loobos
    21 / 33 = Lutjewad
    22 / 34 = Schiphol
    23 / 35 = Rotterdam

Before running the script, the `dales4` executable should be copied to the `cases/cabauw` directory. Running the `create_input.py` script should now create the full case setup in the `path_out` directory, including the `run.PBS` and `post.PBS` scripts. The submission can be automated by enabling the `auto_submit` switch in `create_input.py`.

