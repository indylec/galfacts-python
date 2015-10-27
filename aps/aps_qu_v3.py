#/Users/leclercq/miniconda/bin/python

import numpy as np 
import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pylab as plt
import aplpy
import gaussbeam
import sys
from numpy import pi
from astropy.io import fits
from scipy import stats
pi=np.pi

q_in=sys.argv[1]
u_in=sys.argv[2]
pol_in=sys.argv[3]
field=sys.argv[4]
width=int(sys.argv[5])



def apodize(na,nb,radius):

    ni=int(na*radius)
    dni=na-ni
    nj=int(radius*nb)
    dnj=nb-nj

    tap1d_x=np.ones(na)
    tap1d_y=np.ones(nb)

    tap1d_x[0:dni]= (np.cos(3.*pi/2.+pi/2.*(np.arange(dni,dtype=np.float32)/(dni-1))))
    tap1d_x[na-dni:]= (np.cos(0.+pi/2.*(np.arange(dni,dtype=np.float32)/(dni-1))))
    tap1d_y[0:dnj]= (np.cos(3.*pi/2.+pi/2.*(np.arange(dnj,dtype=np.float32)/(dnj-1))))
    tap1d_x[nb-dnj:]= (np.cos(0.+pi/2.*(np.arange(dnj,dtype=np.float32)/(dnj-1))))

    taper=np.empty((nb,na))
    for i in range (nb):
        taper[i,:]=tap1d_x
    for i in range (na):
        taper[:,i]=taper[i,:]*tap1d_y

    return taper


#define pixel step and pixel area in rad
pixstep=pi/(60.*180.)
pixarea=pixstep**2
area1024=1024.*1024.*pixarea
area_width=width*width*pixarea
#width of chunk in degrees
width_deg=width/60.


# first step is to get images and import header and data into numpy arrays
q_hdu=fits.open(q_in) 
q_im=q_hdu[0].data

if len(q_im.shape) == 3:
    q_im=q_im[0,:,:]

q_head=q_hdu[0].header

xw=q_head['NAXIS1']
yw=q_head['NAXIS2']

u_hdu=fits.open(u_in) 
u_im=u_hdu[0].data

if len(u_im.shape) == 3:
    u_im=q_im[0,:,:]

u_head=u_hdu[0].header

#Assign zero to blanks

q_nan=np.where(np.isnan(q_im))
q_im[q_nan]=0.0

u_nan=np.where(np.isnan(u_im))
u_im[u_nan]=0.0

#Divide images into chunks

nochunks=int(xw/width)
xcrop=int((xw-width*nochunks)/2.)
ycrop=int((yw-width)/2.)

pad=(1024.-width)/2.

clims_pix_x=np.arange(nochunks)*width+xcrop
clims_pix_y=np.ones(nochunks)*ycrop

center_pix_x=np.arange(nochunks)*width+(xcrop)+width/2.
center_pix_y=np.ones(nochunks)*(yw/2)

#Make taper
ft_taper=apodize(width,width,0.98)
   
#get fft spatial frequencies
freq_1d=np.fft.fftshift(np.fft.fftfreq(1024,pixstep))

ell_1d= 2*pi*freq_1d

ellx,elly=np.meshgrid(ell_1d,ell_1d)

#define radii of constant ell
ell_r=np.sqrt((ellx)**2+(elly)**2)
ell_ref=ell_r
ell_max=np.max(ell_ref)

#set up grid of phis
phi_ell=np.arctan2(elly,ellx)

#Make ell bins, log-spaced, and bin ells
bins=np.logspace(np.log10(10.0),np.log10(ell_max),100).astype(np.uint64)
ell_scale=bins*(bins+1)/2.*pi
print "ell bins:", bins
ell_hist=np.histogram(ell_r,bins)[0]

#axis stuff for plots
bins_center=np.zeros((bins.size)-1)

for i in range((bins.size)-1):
    bins_center[i]=bins[i]+(bins[i+1]-bins[i])/2.

bins_axis=bins_center

#use square of beam, divide out from ee, bb
ft_beam=gaussbeam.makeFTgaussian(1024,fwhm=3.5)
ft_beam_sq=np.abs(ft_beam)**2

#Iterate over chunks

for i in range (nochunks):

    print "Working on chunk",i
    
    q_chunk=q_im[clims_pix_y[i]:clims_pix_y[i]+width,clims_pix_x[i]:clims_pix_x[i]+width]
    u_chunk=u_im[clims_pix_y[i]:clims_pix_y[i]+width,clims_pix_x[i]:clims_pix_x[i]+width]

#remove mean

    q_chunk=u_chunk-np.mean(q_chunk)
    u_chunk=u_chunk-np.mean(u_chunk)

#Do FTs
#taper, pad
    q_chunk=np.pad(q_chunk*ft_taper,(pad,),mode='constant')
    u_chunk=np.pad(u_chunk*ft_taper,(pad,),mode='constant')
   
#ft
    qft_final=np.fft.fftshift(np.fft.fft2(q_chunk))
    uft_final=np.fft.fftshift(np.fft.fft2(u_chunk))
    
#Calculate E and B mode functions

    emode=qft_final*np.cos(2.*phi_ell)+uft_final*np.sin(2.*phi_ell)
    bmode=-qft_final*np.sin(2.*phi_ell)+uft_final*np.cos(2.*phi_ell)

#compute correlations (EE,BB)
    ee=np.abs(emode)**2
    bb=np.abs(bmode)**2

#(divide out beam &) account for size of array:
    ee=ee/1024**2
    bb=bb/1024**2

    ee_scaled=ee/ft_beam
    bb_scaled=bb/ft_beam

#Bin the C_l for E and B to calculate radial average
    
    ee_hist=np.histogram(ell_r,bins,weights=ee_scaled)[0]
    ee_average=np.zeros(ee_hist.size)
    nonzero_ee=np.where(ee_hist!=0)
    ee_average[nonzero_ee]=area_width*ee_hist[nonzero_ee]/ell_hist[nonzero_ee]

    #print ee_average[0:30]

    bb_hist=np.histogram(ell_r,bins,weights=bb_scaled)[0]
    bb_average=np.zeros(bb_hist.size)
    nonzero_bb=np.where(bb_hist!=0)
    bb_average[nonzero_bb]=area_width*bb_hist[nonzero_bb]/ell_hist[nonzero_bb]


    #print bins_axis[30:70]-bins_axis[50]
    #print np.log10(bins_axis[30:70])-np.log10(bins_axis[50])    
    
#Fit the power-law linearly, for intermediate ells, subtract center value
    slope,offset,c,d,e=stats.linregress(np.log10(bins_axis[30:70])-np.log10(bins_axis[50]),np.log10((ee_average[30:70]+bb_average[30:70])/2))

    print slope,offset
    pwr_label="power law, index="+str(slope)[:6]
    

#plot plot plot

    fig=plt.figure(figsize=(12,6))
    title='APS of GALFACTS 3.1.2 field '+field
    fig.suptitle(title,size='medium')

#Plot the chunk on the right hand sde of the figure
    f1 = aplpy.FITSFigure(pol_in, figure=fig, subplot=[0.55,0.1,0.4,0.8])
    f1.tick_labels.set_font(size='x-small')
    f1.axis_labels.set_font(size='small')
    f1.show_colorscale(cmap='afmhot')
    x_coord,y_coord=f1.pixel2world(center_pix_x[i],center_pix_y[i])
    f1.recenter(x_coord,y_coord,width=width_deg,height=width_deg)
    f1.add_colorbar()
    f1.colorbar.set_font(size='x-small')
    f1.colorbar.set_axis_label_text('K')
    f1.axis_labels.hide()

#Plot the APS on the left hand side of the figure
    ax=fig.add_axes([0.065,0.1,0.4,0.8])

    ax.set_xlabel('$\ell$',fontsize='medium' )
    ee_lin,= ax.plot(bins_axis[nonzero_ee],ee_average[nonzero_ee],'r-',alpha=0.5)
    ee_mark,= ax.plot(bins_axis[nonzero_ee],ee_average[nonzero_ee],'r|',markersize=8)
    bb_lin,=ax.plot(bins_axis[nonzero_bb],bb_average[nonzero_bb],'b-',alpha=0.5)
    bb_mark,=ax.plot(bins_axis[nonzero_bb],bb_average[nonzero_bb],'b|',markersize=8)


    power_law,= ax.plot(bins_axis[30:70],10**(slope*(np.log10(bins_axis[30:70])-np.log10(bins_axis[50]))+offset),'k-',linewidth=0.8)


    beam_cut =ax.axvline(x=180/(3.5/60.),color='k',linestyle='dotted',alpha=0.7)

    #ax.set_ylim(1E-6,10)
    ax.set_xlim(10,20000.) 

    ax.set_ylabel('$C_{\ell}[K^2]$',fontsize='medium')
    ax.tick_params(labelsize='small')
    ax.legend([(ee_lin,ee_mark),(bb_lin,bb_mark),power_law,beam_cut],["EE","BB",pwr_label,"beamwidth scale"],fontsize='small')
    ax.set_xscale('log') 
    ax.set_yscale('log')


<<<<<<< HEAD
    fig.savefig("/Users/leclercq/galfacts/aps/plots/v3.1/"+field+"_apsv3.1_dqa3.2_c"+str(i)+".pdf",dpi=100)
=======
    fig.savefig("/Users/leclercq/galfacts/aps/plots/v4/"+field+"_apsv4_dqa3.2_c"+str(i)+".pdf",dpi=100)
>>>>>>> bfab1cee915774ee0026ff6e68d072827b950814
#plt.show()
