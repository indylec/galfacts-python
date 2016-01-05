#/Users/leclercq/miniconda/bin/python

import numpy as np 
import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import aplpy
import gaussbeam
import sys
from numpy import pi
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from scipy import stats
import os.path
pi=np.pi

q_in=sys.argv[1]
u_in=sys.argv[2]
#pol_in=sys.argv[3]
#latlon_out=sys.argv[4]
Field=sys.argv[3]
width=int(sys.argv[4])



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
    u_im=u_im[0,:,:]

u_head=u_hdu[0].header

#get WCS object from pol. map

#w=WCS(pol_in)

#Assign zero to blanks

q_nan=np.where(np.isnan(q_im))
q_im[q_nan]=0.0

u_nan=np.where(np.isnan(u_im))
u_im[u_nan]=0.0

#Divide images into chunks

nochunks=int(xw/width)
xcrop=int((xw-width*nochunks)/2.)
ycrop=int((yw-width)/2.)

pad=int((1024.-width)/2.)

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
print ell_r[512,:]
ell_ref=ell_r
ell_max=np.max(ell_ref)

#set up grid of phis
phi_ell=np.arctan2(elly,ellx)

#Make ell bins, log-spaced, and bin ells
bins=np.logspace(np.log10(10.0),np.log10(ell_max),100).astype(np.uint64)
ell_scale=bins*(bins+1)/2.*pi
#print "ell bins:", bins
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

    q_chunk=q_chunk-np.mean(q_chunk)
    u_chunk=u_chunk-np.mean(u_chunk)

#Do FTs
#taper, pad
    q_chunk=np.pad(q_chunk*ft_taper,(pad,),mode='constant')
    u_chunk=np.pad(u_chunk*ft_taper,(pad,),mode='constant')
   
#ft
    qft_final=np.fft.fftshift(np.fft.fft2(q_chunk))
    uft_final=np.fft.fftshift(np.fft.fft2(u_chunk))

#output 2d image

    ## plt.imshow(np.log(np.abs(qft_final)))
    ## #plt.show()
    ## plt.savefig('/Users/leclercq/galfacts/inpainting/'+field+'_2dft_'+str(i)+'.png',dpi=100)
    
    ## #plt.imshow(np.log(np.abs(uft_final)))
    ## #plt.savefig(field+'_2dft_u_'+str(i)+'.png',dpi=100)
    
    
    
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


    #fig=plt.figure()

    gs1 = gridspec.GridSpec(2,2)#, wspace=0 ,hspace=0, )
    

    ax1=plt.subplot(gs1[0,0])
    ax2=plt.subplot(gs1[0,1])
    ax3=plt.subplot(gs1[1,0])
    ax4=plt.subplot(gs1[1,1])

    from matplotlib.patches import Circle

    d=Circle((256,256),radius=24)
    e=Circle((256,256),radius=48)
    f=Circle((256,256),radius=95)
    g=Circle((256,256),radius=142)
    h=Circle((256,256),radius=190)

    d.set_facecolor("none")
    e.set_facecolor("none")
    f.set_facecolor("none")
    g.set_facecolor("none")
    h.set_facecolor("none")

    d.set_alpha(0.5)
    e.set_alpha(0.5)
    f.set_alpha(0.5)
    g.set_alpha(0.5)
    h.set_alpha(0.5)

    circ1=str(int(ell_r[512,512+24]))
    circ2=str(int(ell_r[512,512+48]))
    circ3=str(int(ell_r[512,512+95]))
    circ4=str(int(ell_r[512,512+142]))
    circ5=str(int(ell_r[512,512+190]))

    d1=Circle((256,256),radius=24)
    e1=Circle((256,256),radius=48)
    f1=Circle((256,256),radius=95)
    g1=Circle((256,256),radius=142)
    h1=Circle((256,256),radius=190)

    d1.set_facecolor("none")
    e1.set_facecolor("none")
    f1.set_facecolor("none")
    g1.set_facecolor("none")
    h1.set_facecolor("none")

    d1.set_alpha(0.5)
    e1.set_alpha(0.5)
    f1.set_alpha(0.5)
    g1.set_alpha(0.5)
    h1.set_alpha(0.5)

    d2=Circle((256,256),radius=24)
    e2=Circle((256,256),radius=48)
    f2=Circle((256,256),radius=95)
    g2=Circle((256,256),radius=142)
    h2=Circle((256,256),radius=190)

    d2.set_facecolor("none")
    e2.set_facecolor("none")
    f2.set_facecolor("none")
    g2.set_facecolor("none")
    h2.set_facecolor("none")

    d2.set_alpha(0.5)
    e2.set_alpha(0.5)
    f2.set_alpha(0.5)
    g2.set_alpha(0.5)
    h2.set_alpha(0.5)

    d3=Circle((256,256),radius=24)
    e3=Circle((256,256),radius=48)
    f3=Circle((256,256),radius=95)
    g3=Circle((256,256),radius=142)
    h3=Circle((256,256),radius=190)

    d3.set_facecolor("none")
    e3.set_facecolor("none")
    f3.set_facecolor("none")
    g3.set_facecolor("none")
    h3.set_facecolor("none")

    d3.set_alpha(0.5)
    e3.set_alpha(0.5)
    f3.set_alpha(0.5)
    g3.set_alpha(0.5)
    h3.set_alpha(0.5)

    im1=ax1.imshow(np.log(np.abs(qft_final))[256:768,256:768])
    cbar1=plt.colorbar(mappable=im1, ax=ax1)
    cbar1.set_label('$\mathrm{log}_{10} \; C_{\ell} [\mathrm{K}^2]$',size=12)
    cbar1.ax.tick_params(labelsize=7)

    ax1.text(12,32,'circles at $\ell$ = '+circ1+', '+circ2+', '+circ3+', '+circ4+', '+circ5,fontsize=7)
    ax1.text(12,450,'Q')

    ax1.add_artist(d)
    ax1.add_artist(e)
    ax1.add_artist(f)
    ax1.add_artist(g)
    ax1.add_artist(h)
    

    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.set_xticks([])
    ax1.set_yticks([])
    #ax1.xaxis.set_ticks([128.,256.,384.,512.,640.,768.,896.])
    #ax1.xaxis.set_ticklabels(str(int(ell_r[512,128])),str(int(ell_r[512,256])),str(int(ell_r[512,384])),str(int(ell_r[512,512])),str(int(ell_r[512,640])),str(int(ell_r[512,768])),str(int(ell_r[512,896])))
    #ax1.yaxis.set_ticks([128.,256.,384.,512.,640.,768.,896.])
    #ax1.yaxis.set_ticklabels(str(int(ell_r[512,128])),str(int(ell_r[512,256])),str(int(ell_r[512,384])),str(int(ell_r[512,512])),str(int(ell_r[512,640])),str(int(ell_r[512,768])),str(int(ell_r[512,896])))

    
    im2=ax2.imshow(np.log(np.abs(uft_final))[256:768,256:768])
    cbar2=plt.colorbar(mappable=im2, ax=ax2)
    cbar2.set_label('$\mathrm{log}_{10} \; C_{\ell} [\mathrm{K}^2]$',size=12)
    cbar2.ax.tick_params(labelsize=7)

    ax2.text(12,32,'circles at $\ell$ = '+circ1+', '+circ2+', '+circ3+', '+circ4+', '+circ5,fontsize=7)
    ax2.text(12,450,'U')

    ax2.add_artist(d1)
    ax2.add_artist(e1)
    ax2.add_artist(f1)
    ax2.add_artist(g1)
    ax2.add_artist(h1)

    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.set_xticks([])
    ax2.set_yticks([])

    #ax2.xaxis.set_ticks([128.,256.,384.,512.,640.,768.,896.])
    #ax2.xaxis.set_ticklabels(str(int(ell_r[512,128])),str(int(ell_r[512,256])),str(int(ell_r[512,384])),str(int(ell_r[512,512])),str(int(ell_r[512,640])),str(int(ell_r[512,768])),str(int(ell_r[512,896])))
    #ax2.yaxis.set_ticks([128.,256.,384.,512.,640.,768.,896.])
    #ax2.yaxis.set_ticklabels(str(int(ell_r[512,128])),str(int(ell_r[512,256])),str(int(ell_r[512,384])),str(int(ell_r[512,512])),str(int(ell_r[512,640])),str(int(ell_r[512,768])),str(int(ell_r[512,896])))


    im3=ax3.imshow(np.log(np.abs(ee_scaled))[256:768,256:768])
    cbar3=plt.colorbar(mappable=im3, ax=ax3)
    cbar3.set_label('$\mathrm{log}_{10} \; C_{\ell} [\mathrm{K}^2]$',size=12)
    cbar3.ax.tick_params(labelsize=7)
    
    ax3.text(12,32,'circles at $\ell$ = '+circ1+', '+circ2+', '+circ3+', '+circ4+', '+circ5,fontsize=7)
    ax3.text(12,450,'EE')

    ax3.add_artist(d2)
    ax3.add_artist(e2)
    ax3.add_artist(f2)
    ax3.add_artist(g2)
    ax3.add_artist(h2)

    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
    ax3.set_xticks([])
    ax3.set_yticks([])

    #ax3.xaxis.set_ticks([128.,256.,384.,512.,640.,768.,896.])
    #ax3.xaxis.set_ticklabels(str(int(ell_r[512,128])),str(int(ell_r[512,256])),str(int(ell_r[512,384])),str(int(ell_r[512,512])),str(int(ell_r[512,640])),str(int(ell_r[512,768])),str(int(ell_r[512,896])))
    #ax3.yaxis.set_ticks([128.,256.,384.,512.,640.,768.,896.])
    #ax3.yaxis.set_ticklabels(str(int(ell_r[512,128])),str(int(ell_r[512,256])),str(int(ell_r[512,384])),str(int(ell_r[512,512])),str(int(ell_r[512,640])),str(int(ell_r[512,768])),str(int(ell_r[512,896])))


    im4=ax4.imshow(np.log(np.abs(bb_scaled))[256:768,256:768])
    cbar4=plt.colorbar(mappable=im4, ax=ax4)
    cbar4.set_label('$\mathrm{log}_{10} \; C_{\ell} [\mathrm{K}^2]$',size=12)
    cbar4.ax.tick_params(labelsize=7)
    
    ax4.text(12,32,'circles at $\ell$ = '+circ1+', '+circ2+', '+circ3+', '+circ4+', '+circ5,fontsize=7)
    ax4.text(12,450,'BB')

    ax4.add_artist(d3)
    ax4.add_artist(e3)
    ax4.add_artist(f3)
    ax4.add_artist(g3)
    ax4.add_artist(h3)

    ax4.set_xticklabels([])
    ax4.set_yticklabels([])
    ax4.set_xticks([])
    ax4.set_yticks([])

    #ax4.xaxis.set_ticks([128.,256.,384.,512.,640.,768.,896.])
    #ax4.xaxis.set_ticklabels(str(int(ell_r[512,128])),str(int(ell_r[512,256])),str(int(ell_r[512,384])),str(int(ell_r[512,512])),str(int(ell_r[512,640])),str(int(ell_r[512,768])),str(int(ell_r[512,896])))
    #ax4.xaxis.set_ticks([128.,256.,384.,512.,640.,768.,896.])
    #ax4.xaxis.set_ticklabels(str(int(ell_r[512,128])),str(int(ell_r[512,256])),str(int(ell_r[512,384])),str(int(ell_r[512,512])),str(int(ell_r[512,640])),str(int(ell_r[512,768])),str(int(ell_r[512,896])))

    field=Field.lower() 

    plt.savefig(field+str(i)+'_2d_noise_ft_zoom.png',dpi=300)

## #Bin the C_l for E and B to calculate radial average
    
##     ee_hist=np.histogram(ell_r,bins,weights=ee_scaled)[0]
##     ee_average=np.zeros(ee_hist.size)
##     nonzero_ee=np.where(ee_hist!=0)
##     ee_average[nonzero_ee]=area_width*ee_hist[nonzero_ee]/ell_hist[nonzero_ee]

##     #print ee_average[0:30]

##     bb_hist=np.histogram(ell_r,bins,weights=bb_scaled)[0]
##     bb_average=np.zeros(bb_hist.size)
##     nonzero_bb=np.where(bb_hist!=0)
##     bb_average[nonzero_bb]=area_width*bb_hist[nonzero_bb]/ell_hist[nonzero_bb]


##     #print bins_axis[30:70]-bins_axis[50]
##     #print np.log10(bins_axis[30:70])-np.log10(bins_axis[50])    
    
## #Fit the power-law linearly, for intermediate ells, subtract center value
##     slope,offset,c,d,e=stats.linregress(np.log10(bins_axis[30:70])-np.log10(bins_axis[50]),np.log10((ee_average[30:70]+bb_average[30:70])/2))

##     print slope,offset
##     pwr_label="power law, index="+str(slope)[:6]

## #Store the center pixel galactic coordinates and APS slope value

##     lon,lat=w.all_pix2world(center_pix_x[i],center_pix_y[i],0)

##     icrs_coords=SkyCoord(lon,lat,frame='icrs',unit='deg')

##     gal_coords=icrs_coords.galactic

##     row=field+'{0:2d} '.format(i) +gal_coords.to_string('decimal')+' {0:.3f} '.format(slope) +'\n'

##     if os.path.isfile(latlon_out):
##         f=open(latlon_out, 'a')
##         f.write(row)
##         f.close()
##     else:
##         f=open(latlon_out,'w')
##         f.write('#Field  Chunk  l   b   slope\n')
##         f.write(row)
##         f.close()

    
## #plot plot plot

##     fig=plt.figure(figsize=(12,6))
##     title='APS of GALFACTS 3.1.2 field '+field
##     fig.suptitle(title,size='medium')

## #Plot the chunk on the right hand sde of the figure
##     f1 = aplpy.FITSFigure(pol_in, figure=fig, subplot=[0.55,0.1,0.4,0.8])
##     f1.tick_labels.set_font(size='x-small')
##     f1.axis_labels.set_font(size='small')
##     f1.show_colorscale(cmap='afmhot')
##     x_coord,y_coord=f1.pixel2world(center_pix_x[i],center_pix_y[i])
##     f1.recenter(x_coord,y_coord,width=width_deg,height=width_deg)
##     f1.add_colorbar()
##     f1.colorbar.set_font(size='x-small')
##     f1.colorbar.set_axis_label_text('K')
##     f1.axis_labels.hide()

## #Plot the APS on the left hand side of the figure
##     ax=fig.add_axes([0.065,0.1,0.4,0.8])

##     ax.set_xlabel('$\ell$',fontsize='medium' )
##     ee_lin,= ax.plot(bins_axis[nonzero_ee],ee_average[nonzero_ee],'r-',alpha=0.5)
##     ee_mark,= ax.plot(bins_axis[nonzero_ee],ee_average[nonzero_ee],'r|',markersize=8)
##     bb_lin,=ax.plot(bins_axis[nonzero_bb],bb_average[nonzero_bb],'b-',alpha=0.5)
##     bb_mark,=ax.plot(bins_axis[nonzero_bb],bb_average[nonzero_bb],'b|',markersize=8)


##     power_law,= ax.plot(bins_axis[30:70],10**(slope*(np.log10(bins_axis[30:70])-np.log10(bins_axis[50]))+offset),'k-',linewidth=0.8)


##     beam_cut =ax.axvline(x=180/(3.5/60.),color='k',linestyle='dotted',alpha=0.7)

##     #ax.set_ylim(1E-6,10)
##     ax.set_xlim(10,20000.) 

##     ax.set_ylabel('$C_{\ell}[K^2]$',fontsize='medium')
##     ax.tick_params(labelsize='small')
##     ax.legend([(ee_lin,ee_mark),(bb_lin,bb_mark),power_law,beam_cut],["EE","BB",pwr_label,"beamwidth scale"],fontsize='small')
##     ax.set_xscale('log') 
##     ax.set_yscale('log')



##     fig.savefig("/Users/leclercq/galfacts/aps/plots/v4.2/"+field+"_apsv4.2_dqa3.1.2_c"+str(i)+".pdf",dpi=100)

## #plt.show()
