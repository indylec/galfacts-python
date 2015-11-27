import numpy as np
from astropy.io import fits
from astropy import wcs
import sys
import inpaint

#Read catalogue and replace found sources of given shape and size with NANs before running inpainting algorithm. This is done three times with subsets of gradually lower peak brightness and slightly different inpainting parameters. Runs single average map, given a Topcat-matched I/Q/U Aegean catalogue

table_in=sys.argv[1]
im_in=sys.argv[2]
field=sys.argv[3]
stokes=sys.argv[4]


tableHDU=fits.open(table_in)
tbhead=tableHDU[1].header
tbdata=tableHDU[1].data

#read in q image, get dimensions, coordinate conversion factors etc.

HDU=fits.open(im_in)
im=HDU[0].data
im=im[0,:,:]
head=HDU[0].header

xdim=head['NAXIS1']
ydim=head['NAXIS2']

w_sky=wcs.WCS(head,naxis=2)

#set up grid with dimensions of map
y_big,x_big=np.mgrid[0:ydim,0:xdim]

#select Stokes
if stokes == 'I':
    rms='rms_1'
    peak='Peak_1'
    ra='RA_1a'
    dec='DEC_1a'

elif stokes == 'Q':
    rms='rms_2'
    peak='Peak_2'
    ra='RA_2a'
    dec='DEC_2a'

elif stokes == 'U':
    rms='rms_3'
    peak='Peak_3'
    ra='RA_3a'
    dec='DEC_3a'
    

#filter relevant rows
tbdata4=tbdata[np.abs(tbdata[peak])>4.0*np.abs(tbdata[rms])]
tbdata3=tbdata[np.logical_and(4.0*np.abs(tbdata[rms])>np.abs(tbdata[peak]),np.abs(tbdata[peak])>3.0*np.abs(tbdata[rms]))]
tbdata2=tbdata[np.logical_and(3.0*np.abs(tbdata[rms])>np.abs(tbdata[peak]),np.abs(tbdata[peak])>2.0*np.abs(tbdata[rms]))]
tbdata1=tbdata[np.logical_and(2.0*np.abs(tbdata[rms])>np.abs(tbdata[peak]),np.abs(tbdata[peak])>1.0*np.abs(tbdata[rms]))]

#start for loop


#source-by-source, use info to determine which pixels are replaced (gaussian overlay with semimajor & semiminor axes as provided in catalogue - need to convert from arcsec, RA, DEC into pixels)
print 'Total number of sources to blank this iteration:',len(tbdata4)
count=1
for tabline in tbdata4:

    print 'Blanking source number',count
    #convert arcsec to pixels - one arcsec = 1/60. pixel
#get pixel coords of source center: convert from RA,DEC: use astropy WCS functionality
    temp_pixel=w_sky.wcs_world2pix([[tabline['RA_2a'],tabline['DEC_2a']]],0)
    xc_temp,yc_temp=temp_pixel[0][0],temp_pixel[0][1]
    #increase source ellipse size by 50%
    a_temp=tabline['a_2']/60.*1.5
    b_temp=tabline['b_2']/60.*1.5
    theta_temp=tabline['pa_2']*np.pi/180.

    #select a subregion on the grid centered on your ellipse to speed things up
    box_x1=xc_temp-2*a_temp
    if box_x1<0.0:
        box_x1=0.0
    box_y1=yc_temp-2*a_temp
    if box_y1<0.0:
        box_y1=0.0
    box_x2=box=xc_temp+2*a_temp
    if box_x2>xdim:
        box_x2=xdim
    box_y2=box=yc_temp+2*a_temp
    if box_y2>ydim:
        box_y2=ydim

    #print box_x1,box_x2,box_y1,box_y2

    y=y_big[box_y1:box_y2,box_x1:box_x2]
    x=x_big[box_y1:box_y2,box_x1:box_x2]
#this equation tells you which pixels are contained within an ellipse with a,b,theta and centered on xc and yc (all in pixels, theta in rad)
    inside=np.where((((x-xc_temp)*np.cos(theta_temp)-(y-yc_temp)*np.sin(theta_temp))**2/a_temp**2+((x-xc_temp)*np.sin(theta_temp)+(y-yc_temp)*np.cos(theta_temp))**2/b_temp**2)<1)
    #print y[inside]
    #print x[inside]

#replace all pixels given by "inside" with nan

    im[y[inside],x[inside]]=np.nan

    count+=1
#Now apply inpainting

im=np.float64(im)

print 'Inpainting...'
im=inpaint.replace_nans(im,10,1E-6,2,method='localmean')
print '...done.'

#Second round of blanking and inpainting with less-bright sources, use actual source sizes 
print 'Total number of sources to blank this iteration:',len(tbdata3)
count=1
for tabline in tbdata3:

    print 'Blanking source number',count
    #convert arcsec to pixels - one arcsec = 1/60. pixel
#get pixel coords of source center: convert from RA,DEC: use astropy WCS functionality
    temp_pixel=w_sky.wcs_world2pix([[tabline['RA_2a'],tabline['DEC_2a']]],0)
    xc_temp,yc_temp=temp_pixel[0][0],temp_pixel[0][1]
    #increase source ellipse size by 50%
    a_temp=tabline['a_2']/60.*1.25
    b_temp=tabline['b_2']/60.*1.25
    theta_temp=tabline['pa_2']*np.pi/180.

    #select a subregion on the grid centered on your ellipse to speed things up
    box_x1=xc_temp-2*a_temp
    if box_x1<0.0:
        box_x1=0.0
    box_y1=yc_temp-2*a_temp
    if box_y1<0.0:
        box_y1=0.0
    box_x2=box=xc_temp+2*a_temp
    if box_x2>xdim:
        box_x2=xdim
    box_y2=box=yc_temp+2*a_temp
    if box_y2>ydim:
        box_y2=ydim

    y=y_big[box_y1:box_y2,box_x1:box_x2]
    x=x_big[box_y1:box_y2,box_x1:box_x2]
#this equation tells you which pixels are contained within an ellipse with a,b,theta and centered on xc and yc (all in pixels, theta in rad)
    inside=np.where((((x-xc_temp)*np.cos(theta_temp)-(y-yc_temp)*np.sin(theta_temp))**2/a_temp**2+((x-xc_temp)*np.sin(theta_temp)+(y-yc_temp)*np.cos(theta_temp))**2/b_temp**2)<1)

#replace all pixels given by "inside" with bg value

    im[y[inside],x[inside]]=np.nan
    count+=1
#Now apply inpainting

im=np.float64(im)

print 'Inpainting...'
im=inpaint.replace_nans(im,10,1E-7,2,method='localmean')
print '...done.'

#Third round of blanking/inpainting
print 'Total number of sources to blank this iteration:',len(tbdata2)
count=1
for tabline in tbdata2:

    print 'Blanking source number',count
    #convert arcsec to pixels - one arcsec = 1/60. pixel
#get pixel coords of source center: convert from RA,DEC: use astropy WCS functionality
    temp_pixel=w_sky.wcs_world2pix([[tabline['RA_2a'],tabline['DEC_2a']]],0)
    xc_temp,yc_temp=temp_pixel[0][0],temp_pixel[0][1]
    #increase source ellipse size by 50%
    a_temp=tabline['a_2']/60.
    b_temp=tabline['b_2']/60.
    theta_temp=tabline['pa_2']*np.pi/180.

    #select a subregion on the grid centered on your ellipse to speed things up
    box_x1=xc_temp-2*a_temp
    if box_x1<0.0:
        box_x1=0.0
    box_y1=yc_temp-2*a_temp
    if box_y1<0.0:
        box_y1=0.0
    box_x2=box=xc_temp+2*a_temp
    if box_x2>xdim:
        box_x2=xdim
    box_y2=box=yc_temp+2*a_temp
    if box_y2>ydim:
        box_y2=ydim

    y=y_big[box_y1:box_y2,box_x1:box_x2]
    x=x_big[box_y1:box_y2,box_x1:box_x2]
#this equation tells you which pixels are contained within an ellipse with a,b,theta and centered on xc and yc (all in pixels, theta in rad)
    inside=np.where((((x-xc_temp)*np.cos(theta_temp)-(y-yc_temp)*np.sin(theta_temp))**2/a_temp**2+((x-xc_temp)*np.sin(theta_temp)+(y-yc_temp)*np.cos(theta_temp))**2/b_temp**2)<1)

#replace all pixels given by "inside" with bg value

    im[y[inside],x[inside]]=np.nan

    count+=1
#Now apply inpainting

im=np.float64(im)

print 'Inpainting...'
im_filled=inpaint.replace_nans(im,10,1E-8,2,method='localmean')
print '...done.'

#output new map
newqhdu=fits.PrimaryHDU(im_filled,header=head)
newqhdu.writeto(field+'_'+stokes+'_inpainted.fits')
