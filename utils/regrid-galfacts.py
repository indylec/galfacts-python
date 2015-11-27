
import sys
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from scipy import ndimage
#from scipy import interpolate as ip
from astropy.io import fits

#take command line input of fits file name
input_file=sys.argv[1]
nside=int(sys.argv[2])


def regrid (nside, inmap, channel):
    #regrids the given channel of the galfacts map into HEALPix format using interpolation

    npix=hp.pixelfunc.nside2npix(nside)

    hdulist=fits.open(inmap)
    header0=hdulist[0].header
    tt=hdulist[0].data

    naxis1=header0['naxis1']
    naxis2=header0['naxis2']
    crval1=header0['crval1']
    crpix1=header0['crpix1']
    cdelt1=header0['cdelt1']
    crval2=header0['crval2']
    crpix2=header0['crpix2']
    cdelt2=header0['cdelt2']
    
    ttonechan=tt[channel,:,:]

    #print 'Regridding channel number {!s}.'.format(channel)

    
    print "Getting healpix angles"
    #get HEALPix grid in 2D angles
    theta,phi=hp.pixelfunc.pix2ang(nside,np.arange(npix))
    dec=(np.pi*0.5-theta)*180.0/np.pi
    ra=phi*360.0/(np.pi*2.0)

    # print 'cdelt1:',cdelt1
 
    #get pixel values to interpolate at
    print "getting corresponding pixels"
    rapix = (ra - crval1)/cdelt1+crpix1-1
    
    decpix= (dec-crval2)/cdelt2 + crpix2-1

    #print 'first 100 values of rapix:',rapix[0:100]
    #print 'first 100 values of decpix:',decpix[0:100]
    #print 'ttonechan.shape is', ttonechan.shape

    print "Interpolating:"
    #interpolate the tt map onto the healpix grid!
    
    temptt=ndimage.map_coordinates(ttonechan,[decpix,rapix],order=3,cval=-1.6375e+30,prefilter=False)

    print "...done!"

    

    print 'temptt.shape is', temptt.shape
    #print 'first 500 values of temptt:',temptt[0:500]

    return temptt


map=regrid(nside,input_file,0)

#hp.mollview(map)
#plt.show()

hp.write_map('galfacts_I_healpix.fits',map)
    

   

   

    

    

    


    
    

    
    
    
