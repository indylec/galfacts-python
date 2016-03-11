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
from scipy import stats
import sys
import os

#from matplotlib import rc

#rc('font',family='serif')
#rc('text', usetex=True)

latlon_out='plfits.txt'

for i in os.listdir(os.getcwd()):

    if i.endswith(".npz"):

        print "Working on "+ i

        field=i.split("_")[0]
        chunk=i.split("_")[1]

        data=np.load(i)
        ee1dcorr=data['ee1dcorr']
        bb1dcorr=data['bb1dcorr']
        final_bins=data['bins']
        gal_coords=data['gal_coords']
        eb1d=(ee1dcorr+bb1dcorr)/2.

        
        print gal_coords.size
        print gal_coords.shape
        print str(gal_coords)
        gal_coords=str(gal_coords)

        print gal_coords.split(" ")[-2]
        l=gal_coords.split(" ")[-2]
        l=float(l[1:-1])
        print gal_coords.split(" ")[-1]
        b=gal_coords.split(" ")[-1]
        b=float(b[:-2])
        
        slope,offset,c,d,e=stats.linregress(np.log10(final_bins[32:60])-np.log10(final_bins[45]),np.log10(eb1d[32:60]))

        row=field+' '+chunk[-1]+' '+' {0:.3f} '.format(l) +' {0:.3f} '.format(b) +' {0:.3f} '.format(slope) +' {0:.3f} '.format(offset) +'\n'

        if os.path.isfile(latlon_out):
            f=open(latlon_out, 'a')
            f.write(row)
            f.close()
        else:
            f=open(latlon_out,'w')
            f.write('#Field  Chunk  l   b   slope   offset \n')
            f.write(row)
            f.close()
