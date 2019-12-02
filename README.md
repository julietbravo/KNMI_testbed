### KNMI LES testbed
Code used for the KNMI/DOWA LES work

#### File structure
- `cases`:
    * `cabauw_aug2016`: Two week test period over Cabauw, used in the DOWA LES report.
    * `nudge_boundary_bomex`: Idealised inflow experiments based on BOMEX case.
    * `nudge_boundary_HARMONIE`: Experiments with nudging LES at lateral boundaries to HARMONIE.
    * `nudge_spectral`: ....
- `external/LFA`: LFA source code from Meteo France.
- `misc`: Various scripts to help working with DALES.
- `src`: Mostly code used to convert the DDH/LFA files to NetCDF.


#### Getting started

This testbed setup requires specific Harmonie output, which (for now) has been produced as part of the DOWA (Dutch Offshore Wind Atlas) Harmonie reanalysis, for the period of January 2016 to (including) December 2017, for several locations:

![](https://i.stack.imgur.com/EwFEVl.png)

For 2018 the forcings are available as well, on a ~550x550 km domain around the Netherlands. 
