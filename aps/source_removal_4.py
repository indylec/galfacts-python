import numpy as np
from astropy.io import fits
from astropy import wcs
import sys
import inpaint

#Read catalogue and replace found sources of given shape and size with NANs before running inpainting algorithm. This is done three times with subsets of gradually lower peak brightness and slightly different inpainting parameters. Runs on Q and U average maps, given a Topcat-matched I/Q/U Aegean catalogue

field=sys.argv[1]
table_in=sys.argv[2]
q_in=sys.argv[3]
u_in=sys.argv[4]
i_in=sys.argv[5]

tableHDU=fits.open(table_in)
tbhead=tableHDU[1].header
tbdata=tableHDU[1].data

#read in q image, get dimensions, coordinate conversion factors etc.

qHDU=fits.open(q_in)
qim=qHDU[0].data
qim=qim[0,:,:]
qhead=qHDU[0].header

uHDU=fits.open(u_in)
uim=uHDU[0].data
uim=uim[0,:,:]
uhead=uHDU[0].header

iHDU=fits.open(u_in)
iim=iHDU[0].data
iim=iim[0,:,:]
ihead=iHDU[0].header


xdim=qhead['NAXIS1']
ydim=qhead['NAXIS2']

w_sky=wcs.WCS(qhead,naxis=2)

#set up grid with dimensions of map
y_big,x_big=np.mgrid[0:ydim,0:xdim]

#Work on Q first

#filter relevant rows
qtbdata4=tbdata[np.abs(tbdata['Peak_2'])>4.0*np.abs(tbdata['rms_2'])]
qtbdata3=tbdata[np.logical_and(4.0*np.abs(tbdata['rms_2'])>np.abs(tbdata['Peak_2']),np.abs(tbdata['Peak_2'])>3.0*np.abs(tbdata['rms_2']))]
qtbdata2=tbdata[np.logical_and(3.0*np.abs(tbdata['rms_2'])>np.abs(tbdata['Peak_2']),np.abs(tbdata['Peak_2'])>2.0*np.abs(tbdata['rms_2']))]
qtbdata1=tbdata[np.logical_and(2.0*np.abs(tbdata['rms_2'])>np.abs(tbdata['Peak_2']),np.abs(tbdata['Peak_2'])>1.0*np.abs(tbdata['rms_2']))]

#start for loop


#source-by-source, use info to determine which pixels are replaced (gaussian overlay with semimajor & semiminor axes as provided in catalogue - need to convert from arcsec, RA, DEC into pixels)
print 'Total number of sources to blank this iteration:',len(qtbdata4)
count=1
for tabline in qtbdata4:

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

    qim[y[inside],x[inside]]=np.nan

    count+=1
#Now apply inpainting

qim=np.float64(qim)

print 'Inpainting...'
qim=inpaint.replace_nans(qim,10,1E-6,2,method='localmean')
print '...done.'

#Second round of blanking and inpainting with less-bright sources, use actual source sizes 
print 'Total number of sources to blank this iteration:',len(qtbdata3)
count=1
for tabline in qtbdata3:

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

    qim[y[inside],x[inside]]=np.nan
    count+=1
#Now apply inpainting

qim=np.float64(qim)

print 'Inpainting...'
qim=inpaint.replace_nans(qim,10,1E-7,2,method='localmean')
print '...done.'

#Third round of blanking/inpainting
print 'Total number of sources to blank this iteration:',len(qtbdata2)
count=1
for tabline in qtbdata2:

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

    qim[y[inside],x[inside]]=np.nan

    count+=1
#Now apply inpainting

qim=np.float64(qim)

print 'Inpainting...'
qim_filled=inpaint.replace_nans(qim,10,1E-8,2,method='localmean')
print '...done.'

#output new map
newqhdu=fits.PrimaryHDU(qim_filled,header=qhead)
newqhdu.writeto(field+'/'+field+'_Q_inpaint_final.fits')

#Work on U now
#filter relevant rows
utbdata4=tbdata[np.abs(tbdata['Peak_3'])>4.0*np.abs(tbdata['rms_3'])]
utbdata3=tbdata[np.logical_and(4.0*np.abs(tbdata['rms_3'])>np.abs(tbdata['Peak_3']),np.abs(tbdata['Peak_3'])>3.0*np.abs(tbdata['rms_3']))]
utbdata2=tbdata[np.logical_and(3.0*np.abs(tbdata['rms_3'])>np.abs(tbdata['Peak_3']),np.abs(tbdata['Peak_3'])>2.0*np.abs(tbdata['rms_3']))]
utbdata1=tbdata[np.logical_and(2.0*np.abs(tbdata['rms_3'])>np.abs(tbdata['Peak_3']),np.abs(tbdata['Peak_3'])>1.0*np.abs(tbdata['rms_3']))]

#start for loop


#source-by-source, use info to determine which pixels are replaced (gaussian overlay with semimajor & semiminor axes as provided in catalogue - need to convert from arcsec, RA, DEC into pixels)
print 'Total number of sources to blank this iteration:',len(utbdata4)
count=1
for tabline in utbdata4:

    print 'Blanking source number',count
    #convert arcsec to pixels - one arcsec = 1/60. pixel
#get pixel coords of source center: convert from RA,DEC: use astropy WCS functionality
    temp_pixel=w_sky.wcs_world2pix([[tabline['RA_3a'],tabline['DEC_3a']]],0)
    xc_temp,yc_temp=temp_pixel[0][0],temp_pixel[0][1]
    #increase source ellipse size by 50%
    a_temp=tabline['a_3']/60.*1.5
    b_temp=tabline['b_3']/60.*1.5
    theta_temp=tabline['pa_3']*np.pi/180.

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

    uim[y[inside],x[inside]]=np.nan

    count+=1
#Now apply inpainting

uim=np.float64(uim)

print 'Inpainting...'
uim=inpaint.replace_nans(uim,10,1E-6,2,method='localmean')
print '...done.'

#Second round of blanking and inpainting with less-bright sources, use actual source sizes 
print 'Total number of sources to blank this iteration:',len(utbdata3)
count=1
for tabline in utbdata3:

    print 'Blanking source number',count
    #convert arcsec to pixels - one arcsec = 1/60. pixel
#get pixel coords of source center: convert from RA,DEC: use astropy WCS functionality
    temp_pixel=w_sky.wcs_world2pix([[tabline['RA_3a'],tabline['DEC_3a']]],0)
    xc_temp,yc_temp=temp_pixel[0][0],temp_pixel[0][1]
    #increase source ellipse size by 50%
    a_temp=tabline['a_3']/60.*1.25
    b_temp=tabline['b_3']/60.*1.25
    theta_temp=tabline['pa_3']*np.pi/180.

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

    uim[y[inside],x[inside]]=np.nan
    count+=1
#Now apply inpainting

uim=np.float64(uim)

print 'Inpainting...'
uim=inpaint.replace_nans(uim,10,1E-7,2,method='localmean')
print '...done.'

#Third round of blanking/inpainting
print 'Total number of sources to blank this iteration:',len(utbdata2)
count=1
for tabline in utbdata2:

    print 'Blanking source number',count
    #convert arcsec to pixels - one arcsec = 1/60. pixel
#get pixel coords of source center: convert from RA,DEC: use astropy WCS functionality
    temp_pixel=w_sky.wcs_world2pix([[tabline['RA_3a'],tabline['DEC_3a']]],0)
    xc_temp,yc_temp=temp_pixel[0][0],temp_pixel[0][1]
    #increase source ellipse size by 50%
    a_temp=tabline['a_3']/60.
    b_temp=tabline['b_3']/60.
    theta_temp=tabline['pa_3']*np.pi/180.

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

    uim[y[inside],x[inside]]=np.nan

    count+=1
#Now apply inpainting

uim=np.float64(uim)

print 'Inpainting...'
uim_filled=inpaint.replace_nans(uim,10,1E-8,2,method='localmean')
print '...done.'
#output new map
newuhdu=fits.PrimaryHDU(uim_filled,header=uhead)
newuhdu.writeto(field+'/'+field+'_U_inpaint_final.fits')

#Work on I now
#filter relevant rows
itbdata4=tbdata[np.abs(tbdata['Peak_1'])>4.0*np.abs(tbdata['rms_1'])]
itbdata3=tbdata[np.logical_and(4.0*np.abs(tbdata['rms_1'])>np.abs(tbdata['Peak_1']),np.abs(tbdata['Peak_1'])>3.0*np.abs(tbdata['rms_1']))]
itbdata2=tbdata[np.logical_and(3.0*np.abs(tbdata['rms_1'])>np.abs(tbdata['Peak_1']),np.abs(tbdata['Peak_1'])>2.0*np.abs(tbdata['rms_1']))]
itbdata1=tbdata[np.logical_and(2.0*np.abs(tbdata['rms_1'])>np.abs(tbdata['Peak_1']),np.abs(tbdata['Peak_1'])>1.0*np.abs(tbdata['rms_1']))]

#start for loop


#source-by-source, use info to determine which pixels are replaced (gaussian overlay with semimajor & semiminor axes as provided in catalogue - need to convert from arcsec, RA, DEC into pixels)
print 'Total number of sources to blank this iteration:',len(itbdata4)
count=1
for tabline in itbdata4:

    print 'Blanking source number',count
    #convert arcsec to pixels - one arcsec = 1/60. pixel
#get pixel coords of source center: convert from RA,DEC: use astropy WCS functionality
    temp_pixel=w_sky.wcs_world2pix([[tabline['RA_1a'],tabline['DEC_1a']]],0)
    xc_temp,yc_temp=temp_pixel[0][0],temp_pixel[0][1]
    #increase source ellipse size by 50%
    a_temp=tabline['a_1']/60.*1.5
    b_temp=tabline['b_1']/60.*1.5
    theta_temp=tabline['pa_1']*np.pi/180.

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

    iim[y[inside],x[inside]]=np.nan

    count+=1
#Now apply inpainting

iim=np.float64(iim)

print 'Inpainting...'
iim=inpaint.replace_nans(iim,10,1E-6,2,method='localmean')
print '...done.'

#Second round of blanking and inpainting with less-bright sources, use actual source sizes 
print 'Total number of sources to blank this iteration:',len(itbdata3)
count=1
for tabline in itbdata3:

    print 'Blanking source number',count
    #convert arcsec to pixels - one arcsec = 1/60. pixel
#get pixel coords of source center: convert from RA,DEC: use astropy WCS functionality
    temp_pixel=w_sky.wcs_world2pix([[tabline['RA_1a'],tabline['DEC_1a']]],0)
    xc_temp,yc_temp=temp_pixel[0][0],temp_pixel[0][1]
    #increase source ellipse size by 50%
    a_temp=tabline['a_1']/60.*1.25
    b_temp=tabline['b_1']/60.*1.25
    theta_temp=tabline['pa_1']*np.pi/180.

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

    iim[y[inside],x[inside]]=np.nan
    count+=1
#Now apply inpainting

iim=np.float64(iim)

print 'Inpainting...'
iim=inpaint.replace_nans(iim,10,1E-7,2,method='localmean')
print '...done.'

#Third round of blanking/inpainting
print 'Total number of sources to blank this iteration:',len(itbdata2)
count=1
for tabline in itbdata2:

    print 'Blanking source number',count
    #convert arcsec to pixels - one arcsec = 1/60. pixel
#get pixel coords of source center: convert from RA,DEC: use astropy WCS functionality
    temp_pixel=w_sky.wcs_world2pix([[tabline['RA_1a'],tabline['DEC_1a']]],0)
    xc_temp,yc_temp=temp_pixel[0][0],temp_pixel[0][1]
    #increase source ellipse size by 50%
    a_temp=tabline['a_1']/60.
    b_temp=tabline['b_1']/60.
    theta_temp=tabline['pa_1']*np.pi/180.

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

    iim[y[inside],x[inside]]=np.nan

    count+=1
#Now apply inpainting

iim=np.float64(iim)

print 'Inpainting...'
iim_filled=inpaint.replace_nans(iim,10,1E-8,2,method='localmean')
print '...done.'
#output new map
newuhdu=fits.PrimaryHDU(iim_filled,header=ihead)
newuhdu.writeto(field+'/'+field+'_I_inpaint_final.fits')


#inspect by hand, compute new power spectra
