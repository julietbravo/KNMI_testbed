import numpy as np

r_earth = 6367.47 * 1000    # Radius earth (m)

def dx(lonW, lonE, lat):
    """ Distance between longitudes in spherical coordinates """
    return r_earth * np.cos(np.deg2rad(lat)) * np.deg2rad(lonE - lonW)

def dy(latS, latN):
    """ Distance between latitudes in spherical coordinates """
    return r_earth * np.deg2rad(latN-latS)

def dlon(dx, lat):
    """ East-west distance in degrees for given dx """
    return np.rad2deg(dx / (r_earth * np.cos(np.deg2rad(lat))))

def dlat(dy):
    """ North-south distance in degrees for given dy """
    return np.rad2deg(dy / r_earth)


if __name__ == '__main__':

    locations = [
                    dict(name='Cabauw',   lat=51.97, lon=4.90),
                    dict(name='IJmuiden', lat=52.85, lon=3.44),
                    dict(name='Joyce',    lat=50.91, lon=6.41)
                ]

    area_size = 25000   # m

    # Write input for harmonie_namelists.pm
    # BE CAREFULL: all output locations are assigned to the same plane (BDEDDH(2,X)=1).
    # If one or more output domains overlap, they need to be given a unique plane.
    # See DDH documentation for details.
    for i,loc in enumerate(locations):

        diff_lon = dlon(area_size/2, loc['lat'])
        diff_lat = dlat(area_size/2)

        print('\'BDEDDH(1,{0:})\' => \'4\',        # {1:} : Area average'  .format(i+1, loc['name']))
        print('\'BDEDDH(2,{0:})\' => \'1\',        # Plane'         .format(i+1, loc['name']))
        print('\'BDEDDH(3,{0:})\' => \'{1:.2f}\',  # West longitude'.format(i+1,loc['lon']-diff_lon))
        print('\'BDEDDH(4,{0:})\' => \'{1:.2f}\',  # South latitude'.format(i+1,loc['lat']-diff_lat))
        print('\'BDEDDH(5,{0:})\' => \'{1:.2f}\',  # East longitude'.format(i+1,loc['lon']+diff_lon))
        print('\'BDEDDH(6,{0:})\' => \'{1:.2f}\',  # North latitude'.format(i+1,loc['lat']+diff_lat))
