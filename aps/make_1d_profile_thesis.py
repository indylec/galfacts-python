#Make the 1-D profile plot of the noise and the data for thesis
#NOT ACTUALLY USED, JUST SAVE THE PLOT in FINAL_RESULTS.PY

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.convolution import convolve,convolve_fft, Box1DKernel,Box2DKernel
import aplpy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
from mpl_toolkits.axes_grid1 import make_axes_locatable
pi=np.pi
import sys
from matplotlib import rc

rc('font',family='serif')
rc('text', usetex=True)


####### MAIN ######

data=np.load('')

#Set the azimuthal profile parameters
rell=5000.
rdell=100.
rpix=rell*48./1000.
rdpix=rdell*48./1000.

#get azimuth angles and sort
angle=np.ravel(phi_ell[(ell_r>rell)&(ell_r<rell+rdell)]) 
sort=np.argsort(angle)
anglesort=angle[sort]

#define boxcar kernel
box50=Box1DKernel(50)

#take azimuthal profiles and smooth with boxcar
ee_circle=np.ravel(ee_scaled[(ell_r>rell)&(ell_r<rell+rdell)])
ee_circle_sort=ee_circle[sort]
ee_circle_smooth=convolve_fft(np.abs(ee_circle_sort),box50)

def noisefit_ee(x,a,b):
        return np.sqrt(a**2*np.cos((150.*np.pi/180.)+2*x)**2+b**2*np.sin((150.*np.pi/180.)+2*x)**2)


plt.plot(anglesort,np.log10(np.abs(ee_noise_circle[sort])),'g-',alpha=0.5)
#plt.plot(np.abs(bb_noise_circle),alpha=0.5)
plt.plot(anglesort,np.log10(np.abs(ee_circle_sort)),'b-',alpha=0.5)
#plt.plot(angle[sort],np.log10(np.abs(ee_circle-ee_noise_circle)),'b-',alpha=0.4)#noisefit_ee(angle[sort],*p0),'k-',alpha=0.5)
plt.plot(anglesort,np.log10(ee_circle_smooth),'b')
plt.plot(anglesort,np.log10(ee_noise_circle_smooth),'g')
plt.plot(anglesort,np.log10(noisefit_ee(anglesort,a,b)),'r-')#*fit2[0])),'r-')
#plt.plot(anglesort,np.log10(ee_diff_fit_smoothed),'k-',alpha=0.8)
plt.show()
