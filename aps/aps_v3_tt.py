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

map_in=sys.argv[1]
field=sys.argv[2]
width=int(sys.argv[3])


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
map_hdu=fits.open(map_in) 
map_im=map_hdu[0].data
if len(map_im.shape)==3:
    map_im=map_im[0,:,:]
    
map_head=map_hdu[0].header

xw=map_head['NAXIS1']
yw=map_head['NAXIS2']

#Assign zero to blanks

map_nan=np.where(np.isnan(map_im))
map_im[map_nan]=0.0

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

    map_chunk=map_im[clims_pix_y[i]:clims_pix_y[i]+width,clims_pix_x[i]:clims_pix_x[i]+width]

    #Do FTs
#taper, pad
    map_chunk=np.pad(map_chunk*ft_taper,(pad,),mode='constant')

    ft_final=np.fft.fftshift(np.fft.fft2(map_chunk))

    tt=np.abs(ft_final)**2
    tt=tt/1024**2
    tt_scaled=tt/ft_beam

    tt_hist=np.histogram(ell_r,bins,weights=tt_scaled)[0]
    tt_average=np.zeros(tt_hist.size)
    nonzero_tt=np.where(tt_hist!=0)
    tt_average[nonzero_tt]=area_width*tt_hist[nonzero_tt]/ell_hist[nonzero_tt]

    #plot plot plot

    fig=plt.figure(figsize=(12,6))
    title='APS of GALFACTS 3.1.2 field '+field
    fig.suptitle(title,size='medium')

#Plot the chunk on the right hand sde of the figure
    f1 = aplpy.FITSFigure(map_in, figure=fig, subplot=[0.55,0.1,0.4,0.8])
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
    tt_lin, =ax.plot(bins_axis[nonzero_tt],tt_average[nonzero_tt],'k-',alpha=0.5)
    tt_mark,=ax.plot(bins_axis[nonzero_tt],tt_average[nonzero_tt],'k|',markersize=8)
    beam_cut=ax.axvline(x=180/(3.5/60.),color='k',linestyle='dotted',alpha=0.7)

    ax.set_xlim(10,20000.)  
    
    ax.set_ylabel('$C_{\ell}[K^2]$',fontsize='medium')
    ax.tick_params(labelsize='small')
    ax.legend([(tt_lin,tt_mark),beam_cut],["TT","beamwidth scale"],fontsize='small')
    ax.set_xscale('log') 
    ax.set_yscale('log')


    fig.savefig("/Users/leclercq/galfacts/aps/plots/v4/"+field+"_aps_sources_dqa3.2_c"+str(i)+".pdf",dpi=100)
