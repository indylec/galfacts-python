import numpy as np
from astropy.io import fits
import sys

infile=sys.argv[1]

header=fits.getheader(infile)

xdim=header['NAXIS1']
ydim=header['NAXIS2']

xval=header['CRVAL1']
xpix=header['CRPIX1']
xdel=header['CDELT1']

yval=header['CRVAL2']
ypix=header['CRPIX2']
ydel=header['CDELT2']

xfirst=(0.0-(xpix-1))*xdel+xval
xlast=(xdim-(xpix-1))*xdel+xval
xwidth=np.pi*np.abs(xlast-xfirst)/180.
xwidth_deg=np.abs(xlast-xfirst)
print "x-width (degrees)",xwidth_deg


yfirst=(0.0-(ypix-1))*ydel+yval
ylast=(ydim-(ypix-1))*ydel+yval
ywidth=np.pi*np.abs(ylast-yfirst)/180.
ywidth_deg=np.abs(ylast-yfirst)
if 90.0-yfirst>90.0-ylast:
    theta_1=np.pi*(90.0-yfirst)/180.
    theta_2=np.pi*(90-ylast)/180.
else:
    theta_1=np.pi*(90.0-ylast)/180.
    theta_2=np.pi*(90-yfirst)/180.
print "y -width (degrees)",ywidth_deg

surface=xwidth*(-np.cos(theta_1)+np.cos(theta_2))

print "Surface area is",surface,"sterad"
