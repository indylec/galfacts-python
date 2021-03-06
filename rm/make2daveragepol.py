#!/Users/leclercq/miniconda/bin/python

import numpy as np
#import theil_sen as ts
#import matplotlib
#matplotlib.use('MacOSX')
#import matplotlib.pylab as plt
from astropy.io import fits
import sys
from numpy import pi

infile=sys.argv[1]
field=sys.argv[3]

qin=fits.open(sys.argv[1])
uin=fits.open(sys.argv[2])

polfile=field+'_polarised_intensity.fits'

header_cube=qin[0].header

map_header=header_cube.copy()
try:
    map_header.remove('ctype3')
    map_header.remove('crval3')
    map_header.remove('crpix3')
    map_header.remove('cdelt3')
    map_header.remove('crota3')

except ValueError:
    print "the header is already corrected for 2D image"

map_header['OBJECT']='GALFACTS_{0} Polarised intensity map'.format(field)

try:
    qdata=qin[0].data[0,:,:]
except IndexError:
    qdata=qin[0].data
try:
    udata=uin[0].data[0,:,:]
except IndexError:
    udata=uin[0].data

polint=np.sqrt(qdata**2+udata**2)

fits.writeto(polfile,polint,header=map_header)



