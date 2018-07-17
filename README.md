###KNMI LES testbed

#### Summary
Over the last few years, much progress has been made on running LES for real-life conditions (e.g. [Neggers et al. (2012)](https://doi.org/10.1175/BAMS-D-11-00162.1 "Neggers et al. (2012)"), [Schalkwijk et al. (2015)](https://doi.org/10.1175/MWR-D-14-00293.1 "Schalkwijk & Jonker (2015)"), [Heinze et al. (2017)](http://www.atmos-chem-phys.net/17/7083/2017/acp-17-7083-2017.pdf "Heinze et al. (2017)")). The idea -- summarized probably too briefly -- is to let LES explicitly resolve the small scale physics, while providing the (unresolved) larger-scale dynamics as an external forcing in the LES model.

The dynamic forcings are typically derived offline from the 3D output of models like IFS/ RACMO/ COSMO, by individualy estimating/reconstructing the contributions of large-scale advection, subsidence & pressure gradients (see Eq. 2 and 5 of [Schalkwijk et al. (2015)](https://doi.org/10.1175/MWR-D-14-00293.1 "Schalkwijk & Jonker (2015)")).

In this testbed setup we follow a different approach, by extracting the dynamic tendencies (combination of the aforementioned processes) directly from Harmonie......... (TO-DO: write more...)

#### Getting started

This testbed setup requires specific Harmonie output, which (for now) has been produced as part of the DOWA (Dutch Offshore Wind Atlas) Harmonie reanalysis, for the period of January 2016 to (including) December 2017.

Throughout the code we use a file structure which is similar to that of the Harmonie archive, where files for each three-hourly cycle are stored in a directory:

    BASE_PATH/yyyy/mm/dd/hh/

    Depending on ......, each directory should contain either
    - *[unlikely]* Individual `DHFDLHARM+XXXX` [LFA](http://www.umr-cnrm.fr/gmapdoc/spip.php?article44 "LFA") files, where `XXXX` denotes the Harmonie output time step. The Python script `src/convert_DDH_to_NetCDF.py` is capable of converting the LFA files to NetCDF.
        - ***NOTE:*** this step requires a compilation of the LFA library in `external/LFA`. The LFA code is pure Fortran without any dependencies, so on most systems, a `./install` from that directory does the trick. 
        - *[more likely]* Post-processed NetCDF files (`LES_forcing_yyyymmddhh.nc`).........
Harmonie LES testbed
