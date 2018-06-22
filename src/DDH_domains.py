import numpy as np
import matplotlib.pyplot as pl

import matplotlib.patches as mpatches

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from spatial_tools import *

pl.close('all')

def create_namelist(locations, sizes):
    """
    Create Perl code required by harmonie_namelists.pm
    For averaging domains, DDH requires either the four
    corner coordinates, or the NW and SE corners
    """

    i = 1
    for plane, area_size in enumerate(sizes):

        # Write input for harmonie_namelists.pm
        # BE CAREFULL: all output locations are assigned to the same plane (BDEDDH(2,X)=1).
        # If one or more output domains overlap, they need to be given a unique plane.
        # See DDH documentation for details.
        for loc in locations:

            if area_size == -1:
                print('\'BDEDDH(1,{0:})\' => \'4\',         # {1:}, single column' .format(i, loc['name']))
                print('\'BDEDDH(2,{0:})\' => \'{1:}\',         # Plane'            .format(i, plane+1))
                print('\'BDEDDH(3,{0:})\' => \'{1:.2f}\',      # Central longitude'.format(i,loc['lon']))
                print('\'BDEDDH(4,{0:})\' => \'{1:.2f}\',     # Central latitude'  .format(i,loc['lat']))

            else:
                diff_lon = dlon(area_size/2, loc['lat'])
                diff_lat = dlat(area_size/2)

                print('\'BDEDDH(1,{0:})\' => \'3\',         # {1:}, {2:.0f}x{2:.0f}km'.format(i, loc['name'], area_size/1000.))
                print('\'BDEDDH(2,{0:})\' => \'{1:}\',         # Plane'               .format(i, plane+1))
                print('\'BDEDDH(3,{0:})\' => \'{1:.2f}\',      # West longitude'      .format(i,loc['lon']-diff_lon))
                print('\'BDEDDH(4,{0:})\' => \'{1:.2f}\',     # North latitude'       .format(i,loc['lat']+diff_lat))
                print('\'BDEDDH(5,{0:})\' => \'{1:.2f}\',      # East longitude'      .format(i,loc['lon']+diff_lon))
                print('\'BDEDDH(6,{0:})\' => \'{1:.2f}\',     # South latitude'       .format(i,loc['lat']-diff_lat))

            i += 1

def plot_locations(locations, size):

    pl.figure()

    # Map projection
    proj = ccrs.LambertConformal(central_longitude=4.9, central_latitude=51.967)

    # Create single axes object in Lambert projection
    ax=pl.axes(projection=proj)

    # Draw coast lines, resolution = {'10m','50m','110m'}
    ax.coastlines(resolution='50m', linewidth=0.8, color='black')

    # Load country geometries and lakes (for the IJselmeer)
    # from www.naturalearthdata.com and add to axes
    countries = cfeature.NaturalEarthFeature(
            category='cultural', name='admin_0_boundary_lines_land', scale='10m', facecolor='none')
    ax.add_feature(countries, edgecolor='black', linewidth=0.8)

    lakes = cfeature.NaturalEarthFeature(
            category='physical', name='lakes', scale='10m', facecolor='none')
    ax.add_feature(lakes, edgecolor='black', linewidth=0.8)

    # Set spatial extent of map, and add grid lines
    ax.set_extent([2, 7.5, 50.5, 55], ccrs.PlateCarree())
    #ax.gridlines()

    # Plot some random data
    for i,loc in enumerate(locations):
        color = 'C{}'.format(i+1)

        diff_lon = np.abs(dlon(size/2, loc['lat']))
        diff_lat = np.abs(dlat(size/2))

        ax.add_patch(mpatches.Rectangle(xy=[loc['lon']-diff_lon/2., loc['lat']-diff_lat/2.],
                                        width=2*diff_lon, height=2*diff_lat,
                                        edgecolor=color, facecolor=color, alpha=0.6,
                                        transform=ccrs.PlateCarree(), label=loc['name']))
    pl.legend()


if __name__ == '__main__':

    if (False):
        # Check influence averaging domain:
        locations = [dict(name='Cabauw',   lat=51.97, lon=4.90),
                     dict(name='IJmuiden', lat=52.85, lon=3.44)]

        sizes = [-1, 10000, 30000]

        create_namelist(locations, sizes)

        plot_locations(locations, sizes)

    if (True):
        locations = [dict(name='IJmuiden',     lat=52.85, lon=3.44),
                     dict(name='FINO1',        lat=54.01, lon=6.59),
                     dict(name='OWEZ',         lat=52.61, lon=4.39),
                     dict(name='Goeree',       lat=51.93, lon=3.67),
                     dict(name='Europlatform', lat=52.00, lon=3.27),
                     dict(name='HKZ',          lat=52.30, lon=4.10),
                     dict(name='Cabauw',       lat=51.97, lon=4.90),
                     dict(name='Schiphol',     lat=52.31, lon=4.76)]


                     #dict(name='Joyce',    lat=50.91, lon=6.41)]

        sizes = [25000]

        #create_namelist(locations, sizes)

        plot_locations(locations, 25000)


