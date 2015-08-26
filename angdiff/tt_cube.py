import numpy as np
import theil_sen_uncertain as ts
#import matplotlib
#matplotlib.use('MacOSX')
import matplotlib.pylab as plt
from astropy.io import fits
from astropy.io.fits import Column
import sys
from numpy import pi

print 'Reading in data...'

###Open Q or U maps from both GALFACTS and DRAO - input is galfacts q/u, wolleben q/u, field name and number of subregions per axis 
gal,galhdr=fits.getdata(sys.argv[1],header=True)
wol,wolhdr=fits.getdata(sys.argv[2],header=True)
field=sys.argv[3]
no_divs=int(sys.argv[4])


#get the sizes of both maps, make sure dimensions match

if gal.shape[1] != wol.shape[0] or gal.shape[2] != wol.shape[1]:
    sys.exit("Error: Both maps must have the same x and y dimensions")

#trim the size of the maps to be a multiple of 16 either side
y_subdim=gal.shape[1]/no_divs
x_subdim=gal.shape[2]/no_divs

x_crop=gal.shape[2]%no_divs
y_crop=gal.shape[1]%no_divs

#Make sure we lop off the right number of pixels
if x_crop % 2:
    gal=gal[:,:,x_crop/2+1:-x_crop/2]
    wol=wol[:,x_crop/2+1:-x_crop/2]
else:
    gal=gal[:,:,x_crop/2:-x_crop/2]
    wol=wol[:,x_crop/2:-x_crop/2]
    
if y_crop % 2:
    gal=gal[:,y_crop/2+1:-y_crop/2,:]
    wol=wol[y_crop/2+1:-y_crop/2,:]
    
else:
    gal=gal[:,y_crop/2:-y_crop/2,:]
    wol=wol[y_crop/2:-y_crop/2,:]


#create output products

out_cube=np.empty(gal.shape)
tab_chan=np.empty(no_divs**2*gal.shape[0])
tab_x=np.empty(no_divs**2*gal.shape[0])
tab_y=np.empty(no_divs**2*gal.shape[0])
tab_offset=np.empty(no_divs**2*gal.shape[0])
tab_err_lo=np.empty(no_divs**2*gal.shape[0])
tab_err_hi=np.empty(no_divs**2*gal.shape[0])

#loop over planes of cube

for i in range(gal.shape[0]):

    tab_chan[i*no_divs**2:(i+1)*no_divs**2]=i

    #loop over subregions of map (no_divs**2)
    for j in range (no_divs**2):
        jx=j/no_divs
        jy=j%no_divs

        tab_x[i*no_divs**2+j]=jx
        tab_y[i*no_divs**2+j]=jy

        #prepare data for tt-plot
        gal_temp=np.ravel(gal[i,jy*y_subdim:(jy+1)*y_subdim,jx*x_subdim:(jx+1)*x_subdim])
        wol_temp=np.ravel(wol[jy*y_subdim:(jy+1)*y_subdim,jx*x_subdim:(jx+1)*x_subdim])

        good=np.arange(np.shape(gal_temp)[0])

        good_gal=good[~np.isnan(gal_temp)]
        good_wol=good[~np.isnan(wol_temp)]

        good_both=np.intersect1d(good_gal,good_wol)

        #calculate t-t plot, find offset, store in empty cube, store other parameters in table
        print "processing channel",i
        print "x,y",jx,jy
        print "region",i*no_divs**2+j,"out of",gal.shape[0]*no_divs**2
        params=ts.theil_sen(wol_temp[good_both],gal_temp[good_both])

        offset_temp=params[1]
        err_lo_temp=params[4]
        err_hi_temp=params[5]

        out_cube[i,jy*y_subdim:(jy+1)*y_subdim,jx*x_subdim:(jx+1)*x_subdim]=offset_temp
        tab_offset[i*no_divs**2+j]=offset_temp
        tab_err_lo[i*no_divs**2+j]=err_lo_temp
        tab_err_hi[i*no_divs**2+j]=err_hi_temp
        
        
#save final fits file with cube and table

c1=Column(name='freq_channel',format='3I',array=tab_chan)
c2=Column(name='region_x',format='J',array=tab_x)
c3=Column(name='region_y',format='J',array=tab_y)
c4=Column(name='offset',format='E',unit='K',array=tab_offset)
c5=Column(name='err_lo',format='E',unit='K',array=tab_err_lo)
c6=Column(name='err_hi',format='E',unit='K',array=tab_err_hi)

imhdu=fits.PrimaryHDU(out_cube)
tabhdu=fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6])

out_hdulist=fits.HDUList([imhdu,tabhdu])

out_hdulist.writeto('offset_cube.fits')




