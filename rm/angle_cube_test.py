import numpy as np
import angle_cube as ac
from astropy.io import fits
import sys

qfile=sys.argv[1]
ufile=sys.argv[2]
rmfile=sys.argv[3]

print "Opening q and u cubes..."
q=fits.getdata(qfile)
print "... Q done ..."
u=fits.getdata(ufile)
print "... U done!"

quhead= fits.getheader(qfile)

nu_size = quhead['NAXIS3']
dnu = quhead['CDELT3']
nuref = quhead['CRVAl3']

print "Getting rid of NaNs..."

q[np.isnan(q)]=0.0
print"...Q done..."

u[np.isnan(u)]=0.0

print "...U done!"

q=q.astype(np.float64)
u=u.astype(np.float64)

c2=299792458.**2

dl2 = c2/(dnu**2)

nu = np.arange(nu_size)*dnu+nuref

l2 = 0.5 * c2 * ((nu - 0.5 * dnu) ** -2 + (nu + 0.5 * dnu) ** -2)
l2 = np.flipud(l2)

#angle=0.5*np.arctan2(u,q)

#fits.writeto('angle_test.fits',angle)

rm_map=fits.getdata(rmfile)
rm_map=rm_map.astype(np.float64)

print "Making target angle cube..."
angle_corrected,target=ac.make_angle_cube(q,u,rm_map,l2)
print "...done."

#fits.writeto('angle_corrected_test_1.fits', angle_corrected)
fits.writeto('target_test_1.fits', target)

