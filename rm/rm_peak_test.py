import numpy as np
from astropy.io import fits
import itertools as it
import matplotlib.pyplot as plt
import sys

rmfile=sys.argv[1]
phi=np.arange(100)*10-500


print "Opening RM cube..."
rmcube=fits.getdata(rmfile)
print "...done."
range_no=21
fit_range=np.arange(range_no)-range_no/2
rm0=np.empty((100,100))
#print "phi:",params.phi

print "Beginning peak fits..."
for y,x in it.product(range(100),range(100)):
    print "Finding peak for pixel ({0},{1})".format(x,y)
    temp_los=rmcube[:,y,x]
    #print temp_los

    temp_peak_arg=np.nanargmax(temp_los)
    #print temp_peak_arg

    if (temp_peak_arg>=range_no/2 and temp_peak_arg<100-range_no/2):
        temp_range=fit_range+temp_peak_arg
        #print temp_peak_arg, temp_range
        temp_coeffs=np.polyfit(phi[temp_range],temp_los[temp_range],2)
        rm0[y,x]=-temp_coeffs[1]/(2.*temp_coeffs[0])
    else:
        rm0[y,x]=phi[np.nanargmax(temp_los)]
print "...done."

fits.writeto('rm0_test1.fits',rm0)
