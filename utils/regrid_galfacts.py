import sys
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from scipy import ndimage
from scipy import interpolate as ip
from astropy.io import fits

#take command line input of fits file name
input_file=sys.argv[1]


def getrdvec(header):
    naxis1=header['naxis1']
    naxis2=header['naxis2']
    crval1=header['crval1']
    crpix1=header['crpix1']
    cdelt1=header['cdelt1']
    crval2=header['crval2']
    crpix2=header['crpix2']
    cdelt2=header['cdelt2']
    
    ra=(arange(naxis1) + 1 - crpix1)*cdelt1 + crval1
    dec=(arange(naxis2) + 1 - crpix2)*cdelt2 + crval2

    #Create 2d theta & phi arrays (in colatitude and longitude, radians)

    theta=np.pi/2.0-(dec*np.pi/180.0)
    phi=ra*2.*np.pi/360.0+2.*np.pi

    theta=np.tile(np.reshape(theta,(naxis2,1)),(1,naxis1))
    phi=np.tile(phi,(naxis2,1))

    rdvec=hp.pixelfunc.ang2vec(theta,phi)


    return rdvec

    

    
    

def regrid (nside, inmap):
    #regrids the given channel of the GMIMS map into HEALPix format using interpolation

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
    
    ttonechan=tt

    #print 'Regridding channel number {!s}.'.format(channel)

     #Pad the ends for better interpolation with 10 pixels on either side
    ttonechan=np.concatenate((ttonechan[:,-12:-2],ttonechan,ttonechan[:,1:11]),axis=1)
    ttlen=ttonechan.shape[1]
    print 'Padded array is {!s} pixels wide'.format(ttlen)

     #Define HEALPix directions
    x,y,z=hp.pixelfunc.pix2vec(nside,np.arange(npix),nest=True)
    rdvecs=np.column_stack((x,y,z))

    #get HEALPix grid in 2D angles
    theta,phi=hp.pixelfunc.vec2ang(rdvecs)
    dec=(np.pi*0.5-theta)*180.0/np.pi
    ra=phi*360.0/(np.pi*2.0)-360

    print 'cdelt1:',cdelt1

    #get pixel values to interpolate at
    rapix = (ra - crval1)/cdelt1+crpix1-1+10
    if np.amin(rapix) < 10:
        rapix[np.where(rapix<10)] += ttlen -1
    decpix= (dec-crval2)/cdelt2 + crpix2-1

    #print 'first 100 values of rapix:',rapix[0:100]
    #print 'first 100 values of decpix:',decpix[0:100]
    print 'ttonechan.shape is', ttonechan.shape

    print "Interpolating..."
    #interpolate the tt map onto the healpix grid!
    coords=np.vstack((rapix,decpix)).T
    x=np.arange(ttlen)
    y=np.arange(naxis2)

    #    x2d=np.tile(x,(naxis2,1))
    #    y2d=np.tile(np.reshape(y,(naxis2,1)),(1,ttlen))
    #
    #xycoords=np.vstack((x2d.ravel(),y2d.ravel())).T

    #temptt=ip.griddata(xycoords,ttonechan.ravel(),coords,method='cubic')
   
    temptt=ndimage.map_coordinates(ttonechan,[decpix,rapix],order=3,cval=-1.6375e+30,prefilter=False)

    #f=ip.interp2d(x,y,ttonechan,kind='cubic',copy=True,bounds_error=False,fill_value=-1.6375e+30)
    #temptt=f(decpix,rapix)
    

    print "...done!"

    

    print 'temptt.shape is', temptt.shape
    #print 'first 500 values of temptt:',temptt[0:500]

    return temptt


map=regrid(2048,input_file)
hp.mollview(map,nest=True)
plt.show()

    

   

   

    

    

    


    
    

    
    
    
