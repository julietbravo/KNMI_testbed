import datetime
import urllib.request
import os

path = 'http://baltink.knmi.nl/~baltink/database/database/'

start = datetime.datetime(2016,8,4)
end   = datetime.datetime(2016,8,19)
dt    = datetime.timedelta(hours=24)

date = start
while date < end:

    # LiDAR images
    # ------------
    files = ['bsct','bsct_cb','cb','cb_low','cb_mid']
    for f in files:
        fname  = '{0:04d}{1:02d}{2:02d}_cabauw_ld40_{3}.png'.format(date.year, date.month, date.day, f)
        folder = '{0}/LD40/{1:04d}/{2:02d}/IMAGES/'.format(path, date.year, date.month)
        
        if not os.path.exists('figures/{}'.format(fname)):
            try:
                urllib.request.urlretrieve('{}/{}'.format(folder, fname), 'figures/{}'.format(fname))
            except:
                print('Cant find {}'.format(fname))

    # Cloud radar
    # ------------
    fname  = '{0:04d}{1:02d}{2:02d}_35GHz_CT75_dBZe.png'.format(date.year, date.month, date.day)
    folder = '{0}/PDN2/{1:04d}/{2:02d}/IMAGES/'.format(path, date.year, date.month)
    
    if not os.path.exists('figures/{}'.format(fname)):
        try:
            urllib.request.urlretrieve('{}/{}'.format(folder, fname), 'figures/{}'.format(fname))
        except:
            print('Cant find {}'.format(fname))

    # HATPRO
    # ------------
    fname  = '{0:04d}{1:02d}{2:02d}_hatpro_level1.png'.format(date.year, date.month, date.day)
    folder = '{0}/HATPRO/{1:04d}/{2:02d}/IMAGES/'.format(path, date.year, date.month)
    
    if not os.path.exists('figures/{}'.format(fname)):
        try:
            urllib.request.urlretrieve('{}/{}'.format(folder, fname), 'figures/{}'.format(fname))
        except:
            print('Cant find {}'.format(fname))

    # Nubiscope
    # ---------
    fname  = '{0:04d}{1:02d}{2:02d}_nubiscope_cloud.png'.format(date.year, date.month, date.day)
    folder = '{0}/NUBISCOPE/{1:04d}/{2:02d}/IMAGES/'.format(path, date.year, date.month)
    
    if not os.path.exists('figures/{}'.format(fname)):
        try:
            urllib.request.urlretrieve('{}/{}'.format(folder, fname), 'figures/{}'.format(fname))
        except:
            print('Cant find {}'.format(fname))


    # Video
    # ---------
    fname  = 'axis_{0:04d}{1:02d}{2:02d}.mp4'.format(date.year, date.month, date.day)
    folder = 'http://knmi-cbsql-w01p.knmi.nl/video/{1:04d}/{2:02d}/axis/'.format(path, date.year, date.month)
    
    if not os.path.exists('figures/{}'.format(fname)):
        try:
            urllib.request.urlretrieve('{}/{}'.format(folder, fname), 'figures/{}'.format(fname))
        except:
            print('Cant find {}'.format(fname))


#http://knmi-cbsql-w01p.knmi.nl/video/2016/08/axis/axis_20160804.mp4



    date += dt
