import numpy as np
from astropy.io import fits
from astropy import wcs
import sys

#Read catalogue and replace found sources of given shape and size with background level provided. Only use sources w peak brightness greater than 4 sigma in Q/U

table_in=sys.argv[1]
q_in=sys.argv[2]
#u_in=sys.argv[3]
rms_factor=float(sys.argv[3])

tableHDU=fits.open(table_in)
tbhead=tableHDU[1].header
tbdata=tableHDU[1].data

#filter relevant rows

qtbdata=tbdata[np.abs(tbdata['Peak_2'])>rms_factor*np.abs(tbdata['rms_2'])]
#utbdata=tbdata[tbdata['Peak_3']>rms_factor*tbdata['rms_3']]

print 'Total number of sources:',len(qtbdata)
#read in q image, get dimensions, coordinate conversion factors etc.

qHDU=fits.open(q_in)
qim=qHDU[0].data
qim=qim[0,:,:]
qhead=qHDU[0].header

xdim=qhead['NAXIS1']
ydim=qhead['NAXIS2']

w_sky_q=wcs.WCS(qhead,naxis=2)


#set up grid with dimensions of map
y,x=np.mgrid[0:ydim,0:xdim]

#start for loop


#source-by-source, use info to determine which pixels are replaced (gaussian overlay with semimajor & semiminor axes as provided in catalogue - need to convert from arcsec, RA, DEC into pixels)

count=1
for tabline in qtbdata:

    print 'Subtracting source number',count
    #convert arcsec to pixels - one arcsec = 1/60. pixel
#get pixel coords of source center: convert from RA,DEC: use astropy WCS functionality
    temp_pixel=w_sky_q.wcs_world2pix([[tabline['RA_2a'],tabline['DEC_2a']]],0)
    xc_temp,yc_temp=temp_pixel[0][0],temp_pixel[0][1]
    a_temp=tabline['a_2']/60.
    b_temp=tabline['b_2']/60.
    theta_temp=tabline['pa_2']*np.pi/180.

#this equation tells you which pixels are contained within an ellipse with a,b,theta and centered on xc and yc (all in pixels, theta in rad)
    inside=np.where((((x-xc_temp)*np.cos(theta_temp)-(y-yc_temp)*np.sin(theta_temp))**2/a_temp**2+((x-xc_temp)*np.sin(theta_temp)+(y-yc_temp)*np.cos(theta_temp))**2/b_temp**2)<1)

#replace all pixels given by "inside" with bg value

    qim[inside]=tabline['bkg_2']

    count+=1

#output new map



newhdu=fits.PrimaryHDU(qim,header=qhead)
newhdu.writeto('S1_q_ssub_4rms_origin.fits')

#inspect by hand, compute new power spectra
