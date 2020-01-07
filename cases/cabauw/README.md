## HARMONIE LES testbed
This readme file briefly describes the steps needed to setup and run DALES on ECMWF, with the dynamic forcings from HARMONIE-AROME. 

----
### 1. DALES
----
As described in the DOWA/KNMI technical report (TO-DO: ADD LINK), the dynamic forcings for LES consist of the total dynamic tendencies of the liquid water potential temperature, total specific humidity, and both horizontal wind components. The default DALES version can not handle large-scale/dynamic tendencies for momentum, so a slightly modified version of DALES is required, which can be obtained using `git` (use `module load git` if required):

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
### 2. Dynamic tendencies HARMONIE
----


![](https://i.stack.imgur.com/EwFEVl.png)

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

----
### 3. ERA5 soil temperature/moisture
----

blabla
