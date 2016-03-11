#/Users/leclercq/miniconda/

###PLOTTING CODE; USE THIS TO PLOT THE VARIOUS COMBINATIONS OF MAPS, 1D APS and 2D APS WITH OR WITHOUT THE CUT

#This is the code that caluclates the APS of a given Q and U field. It divides the field up into chunks of a given width, pads to the nearest power of 2, and does the 2D FT. The output displays: 1D aps for QQ, UU, 2D QQ and UU, Q, U and P. We then mask out the bow-tie cut, corresponding to the striping noise, and do the radial averaging to try and figureout what, if any, scaling relation there is between the QQ and UU components.

import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
from mpl_toolkits.axes_grid1 import make_axes_locatable
import aplpy
import gaussbeam
import sys
from numpy import pi
from numpy import ma
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
#from scipy import stats
from scipy import optimize as sciopt
import os.path
pi=np.pi

from matplotlib import rc

rc('font',family='serif')
rc('text', usetex=True)

################ MODULES #############

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


# first step is to get images and import header and data into numpy arrays
def import_qu(q_in,u_in, qn_in, un_in):
    q_hdu=fits.open(q_in) 
    q_im=q_hdu[0].data

    if len(q_im.shape) == 3:
        q_im=q_im[0,:,:]

    q_head=q_hdu[0].header

    u_hdu=fits.open(u_in) 
    u_im=u_hdu[0].data

    if len(u_im.shape) == 3:
        u_im=u_im[0,:,:]

    qn_hdu=fits.open(qn_in) 
    qn_im=qn_hdu[0].data

    if len(qn_im.shape) == 3:
        qn_im=qn_im[0,:,:]

    qn_head=qn_hdu[0].header

    un_hdu=fits.open(un_in) 
    un_im=un_hdu[0].data

    if len(un_im.shape) == 3:
        un_im=un_im[0,:,:]

    un_head=un_hdu[0].header
    
    xw=q_head['NAXIS1']
    yw=q_head['NAXIS2']

    xnw=qn_head['NAXIS1']
    ynw=qn_head['NAXIS2']

   #get WCS object from pol. map

    w=WCS(pol_in)

    #Assign zero to blanks

    q_nan=np.where(np.isnan(q_im))
    q_im[q_nan]=0.0

    u_nan=np.where(np.isnan(u_im))
    u_im[u_nan]=0.0

    return q_im, u_im, qn_im, un_im, xw, yw, xnw, ynw, w

#Divide images into chunks
def make_chunks(width, xw, yw, xnw, ynw):
    nochunks=int(xw/width)
    xcrop=int((xw-width*nochunks)/2.)
    ycrop=int((yw-width)/2.)

    pad=int((1024.-width)/2.)

    clims_pix_x=np.arange(nochunks)*width+xcrop
    clims_pix_y=np.ones(nochunks)*ycrop

    center_pix_x=np.arange(nochunks)*width+(xcrop)+width/2.
    center_pix_y=np.ones(nochunks)*(yw/2)

    nochunks_n=int(xnw/width)
    xcrop_n=int((xnw-width*nochunks_n)/2.)
    ycrop_n=int((ynw-width)/2.)

    clims_pix_xn=np.arange(nochunks_n)*width+xcrop_n
    clims_pix_yn=np.ones(nochunks_n)*ycrop_n

    center_pix_xn=np.arange(nochunks_n)*width+(xcrop_n)+width/2.
    center_pix_yn=np.ones(nochunks_n)*(ynw/2)

    return clims_pix_x,clims_pix_y,center_pix_x,center_pix_y,clims_pix_xn,clims_pix_yn,center_pix_xn,center_pix_yn,pad,nochunks, nochunks_n


#prepare for FFTs
def fft_prep(width, pixstep):
    #Make taper
    ft_taper=apodize(width,width,0.98)

    #get fft spatial frequencies
    freq_1d=np.fft.fftshift(np.fft.fftfreq(1024,pixstep))

    ell_1d= 2*pi*freq_1d

    ellx,elly=np.meshgrid(ell_1d,ell_1d)

    #define radii of constant ell
    ell_r=np.sqrt((ellx)**2+(elly)**2)

    xm=ym=np.arange(1024)
    Xm,Ym=np.meshgrid(xm,ym)
    Ym=Ym.astype('f')
    Xm=Xm.astype('f')

    xyrad=np.sqrt((Xm-512)**2+(Ym-512)**2)
    
    ell_cut=ma.masked_where((xyrad>48) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),ell_r)
    
    ell_max=np.max(ell_cut)

    #set up grid of phis
    phi_ell=np.arctan2(elly,ellx)
    phi_cut=ma.masked_where((xyrad>48) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),phi_ell)

    #Make ell bins, log-spaced, and bin ells
    bins=np.logspace(np.log10(10.0),np.log10(ell_max),100).astype(np.uint64)
    ell_scale=bins*(bins+1)/2.*pi
    #print "ell bins:", bins
    ell_hist=np.histogram(ell_r,bins)[0]
    ell_hist_cut=np.histogram(ell_cut.compressed(),bins)[0]

    #axis stuff for plots
    bins_center=np.zeros((bins.size)-1)

    for i in range((bins.size)-1):
        bins_center[i]=bins[i]+(bins[i+1]-bins[i])/2.

    bins_axis=bins_center
    

#use square of beam, divide out from ee, bb
    ft_beam=gaussbeam.makeFTgaussian(1024,fwhm=3.5)
    #ft_beam_sq=np.abs(ft_beam)**2

    return ell_r, ell_hist, ell_cut, ell_hist_cut,  phi_ell, phi_cut, bins_axis, ft_beam, ft_taper,phi_ell,bins



def chunk_fts(q_im,u_im,clims_pix_y,clims_pix_x,width,ft_taper,pad,i,area_width,beam):
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
    
#Calculate E and B mode functions

    #emode=qft_final*np.cos(2.*phi_ell)+uft_final*np.sin(2.*phi_ell)
    #bmode=-qft_final*np.sin(2.*phi_ell)+uft_final*np.cos(2.*phi_ell)

#Calculate QQ and UU functions

#Divide out beam
    if beam == True:
        qft_final=qft_final/ft_beam
        uft_final=uft_final/ft_beam


#compute correlations (EE,BB)
    qq=qft_final*np.conj(qft_final)
    uu=uft_final*np.conj(uft_final)
    

#account for size of array:
    qq_scaled=qq*area_width/1024**2
    uu_scaled=uu*area_width/1024**2

#mask out the striping features

    xm=ym=np.arange(1024)
    Xm,Ym=np.meshgrid(xm,ym)
    Ym=Ym.astype('f')
    Xm=Xm.astype('f')

    xyrad=np.sqrt((Xm-512)**2+(Ym-512)**2)
    
    qq_scaled_ma=ma.masked_where((xyrad>48) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),qq_scaled)
    uu_scaled_ma=ma.masked_where((xyrad>48) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),uu_scaled)

    

    return qq_scaled, uu_scaled, qq_scaled_ma, uu_scaled_ma

    
#Bin the C_l for E and B to calculate radial average
def binning(ee_scaled, bb_scaled, ell_r, ell_hist, bins, bins_axis):
    print np.shape(ee_scaled.compressed()), np.shape(ell_r.compressed())   
    ee_hist=np.histogram(ell_r.compressed(),bins,weights=ee_scaled.compressed())[0]
    ee_average=np.zeros(ee_hist.size)
    nonzero_ee=np.where(ee_hist!=0)
    #print ee_hist.shape, ell_hist.shape,ell_hist,nonzero_ee
    ee_average[nonzero_ee]=ee_hist[nonzero_ee]/ell_hist[nonzero_ee]

    #print ee_average[0:30]

    bb_hist=np.histogram(ell_r.compressed(),bins,weights=bb_scaled.compressed())[0]
    bb_average=np.zeros(bb_hist.size)
    nonzero_bb=np.where(bb_hist!=0)
    bb_average[nonzero_bb]=bb_hist[nonzero_bb]/ell_hist[nonzero_bb]

    nz=np.intersect1d(nonzero_ee,nonzero_bb)

    #print "There are are ",nz.shape,"nonzero elements with the same index"

    eb=(ee_average[nz]+bb_average[nz])/2.

    ez=ee_average[nz]
    bz=bb_average[nz]

    bins_z=bins_axis[nz]
    

    qu_ratio=np.mean(ez[49:59]/bz[49:59])

    #print "Q/U ratio around ell=",bins_z[55],"is ",qu_ratio

    return ez, bz, nz, bins_z


#Figure out the appropriate noise scaling
def noise_scale(qq_scaled_ma, uu_scaled_ma, qqnoise_scaled_ma, qz, uz, qnz, bins_z):

    #below_beam=np.where(bins_z>180/(3.5/60.))
    #print below_beam

    q_noise_ratio=np.mean(qz[60:71]/qnz[60:71])
    print q_noise_ratio
    u_noise_ratio=np.mean(uz[60:71]/qnz[60:71])
    print u_noise_ratio

    qu_ratio=q_noise_ratio/u_noise_ratio

    q_scaled_noise=qnz*q_noise_ratio
    u_scaled_noise=qnz*u_noise_ratio
    
    qq_average_nonoise=qz-q_scaled_noise

    qq_nonoise=qq_scaled_ma-(qqnoise_scaled_ma*q_noise_ratio)

    uu_average_nonoise=uz-u_scaled_noise

    uu_nonoise=uu_scaled_ma-(qqnoise_scaled_ma*u_noise_ratio)

    return qq_nonoise,uu_nonoise,qq_average_nonoise,uu_average_nonoise,q_scaled_noise,u_scaled_noise, qu_ratio

#plot plot plot
def plot_all(q_in, u_in, pol_in, w, center_pix_x, center_pix_y, width_deg, bins_z, qz, uz, qz_nonoise, uz_nonoise, qnz, q_scaled_noise, u_scaled_noise, qq_scaled_ma, qq_nonoise, chunk, field, beam, qu_ratio):

    fig=plt.figure(figsize=(4,10))
    #title='APS of GALFACTS 3.1.2 field '+field
    #fig.suptitle(title,size='medium')
 
    #3x1 subplots: P, QQ, QQ-cut
    #                  P , Q , U 

#Plot the images on the bottom row
    #f1 = aplpy.FITSFigure(pol_in, figure=fig, subplot=[0.05,0.05,0.30,0.45])
    f1 = aplpy.FITSFigure(q_in, figure=fig, subplot=[0.05,0.71,0.95,0.28])
    f1.tick_labels.set_font(size='large')
    #f1.axis_labels.set_font(size='small')
    #f1.axis_labels.hide()
    #f1.tick_labels.hide()
    f1.show_colorscale(cmap='afmhot', vmin=-0.05, vmax=0.10)
    x_coord,y_coord=f1.pixel2world(center_pix_x[chunk],center_pix_y[chunk])
    f1.recenter(x_coord,y_coord,width=width_deg,height=width_deg)
    f1.add_colorbar()
    f1.colorbar.set_font(size='medium')
    f1.colorbar.set_axis_label_text('K')
    f1.axis_labels.set_xtext("RA")
    f1.axis_labels.set_ytext("Dec")
    #f1.axis_labels.hide_x()
    #f1.tick_labels.hide_x()
    f1.tick_labels.set_xformat('hh:mm')
    f1.tick_labels.set_yformat('dd')
    #f1.axis_labels.hide_x()
    #f1.add_label(0.1,0.1,'P',relative=True, color='white', size='18', weight='bold')
    f1.add_label(0.1,0.1,'Q',relative=True, color='black', size='18', weight='bold')

    ## f2 = aplpy.FITSFigure(q_in, figure=fig, subplot=[0.33,0.05,0.30,0.45])
    ## #f2.tick_labels.set_font(size='x-small')
    ## #f2.axis_labels.set_font(size='small')
    ## f2.axis_labels.hide()
    ## f2.tick_labels.hide()
    ## f2.show_colorscale(cmap='afmhot')
    ## x_coord,y_coord=f2.pixel2world(center_pix_x[chunk],center_pix_y[chunk])
    ## f2.recenter(x_coord,y_coord,width=width_deg,height=width_deg)
    ## f2.add_colorbar()
    ## f2.colorbar.set_font(size='x-small')
    ## f2.colorbar.set_axis_label_text('K')
    ## f2.axis_labels.hide_x()
    ## f2.add_label(0.1,0.1,'Q',relative=True, color='white', size='18', weight='bold')
    

    ## f3 = aplpy.FITSFigure(u_in, figure=fig, subplot=[0.62,0.05,0.30,0.45])
    ## f3.tick_labels.set_font(size='x-small')
    ## #f3.axis_labels.set_font(size='small')
    ## f3.axis_labels.hide()
    ## f3.tick_labels.hide()
    ## f3.show_colorscale(cmap='afmhot')
    ## x_coord,y_coord=f3.pixel2world(center_pix_x[chunk],center_pix_y[chunk])
    ## f3.recenter(x_coord,y_coord,width=width_deg,height=width_deg)
    ## f3.add_colorbar()
    ## f3.colorbar.set_font(size='x-small')
    ## f3.colorbar.set_axis_label_text('K')
    ## f3.add_label(0.1,0.1,'U',relative=True, color='white', size='18', weight='bold')

    

#Plot the APS top left
    ## ax1=fig.add_axes([0.05,0.55,0.35,0.40])
    ## #ax.set_autoscale_on(False)

    ## ax1.set_xlabel('$\ell$',fontsize='medium' )
    ## #qz_lin,= ax1.plot(bins_z,qz,'r-',alpha=0.4)
    ## #qz_mark,= ax1.plot(bins_z,qz,'ro',markersize=3)
    ## #uz_lin,=ax1.plot(bins_z,uz,'b-',alpha=0.4)
    ## #uz_mark,=ax1.plot(bins_z,uz,'bo',markersize=3)
    ## qz_nonoise_lin,=ax1.plot(bins_z,qz_nonoise,'r-',alpha=0.4)
    ## qz_nonoise_mark,=ax1.plot(bins_z,qz_nonoise,'r^',markersize=3)
    ## uz_nonoise_lin,=ax1.plot(bins_z,uz_nonoise,'b-',alpha=0.4)
    ## uz_nonoise_mark,=ax1.plot(bins_z,uz_nonoise,'b^',markersize=3)
    ## q_scaled_noise_lin,=ax1.plot(bins_z,q_scaled_noise,'g-', alpha=0.4)
    ## #q_scaled_noise_mark,=ax1.plot(bins_z,q_scaled_noise,'gd', markersize=3)
    ## u_scaled_noise_lin,=ax1.plot(bins_z,u_scaled_noise,'g-', alpha=0.4)
    ## #u_scaled_noise_mark,=ax1.plot(bins_z,u_scaled_noise,'gD', markersize=3)
    
    ## #eefitlin = ax.plot(bins_axis[nonzero_ee],apsfit_ee(bins_axis[nonzero_ee],pfit_ee[0],pfit_ee[1],pfit_ee[2]), 'r-')
    ## #bbfitlin = ax.plot(bins_axis[nonzero_bb],apsfit_bb(bins_axis[nonzero_bb],pfit_bb[0],pfit_bb[1],pfit_bb[2]), 'b-')
    ## #ebfitlin = ax.plot(bins_axis[nz],apsfit_eb(bins_axis[nz],pfit_eb[0],pfit_eb[1],pfit_eb[2]), 'k-')

    ## #power_law,= ax.plot(bins_axis[30:70],10**(slope*(np.log10(bins_axis[30:70])-np.log10(bins_axis[50]))+offset),'k-',linewidth=0.8)


    ## beam_cut = ax1.axvline(x=180/(3.5/60.),color='k',linestyle='dashed',alpha=0.8)

    ## ax1.axvline(x=506.,color='k',linestyle='dotted',alpha=0.6)
    ## ax1.axvline(x=1012.,color='k',linestyle='dotted',alpha=0.6)
    ## ax1.axvline(x=2003.,color='k',linestyle='dotted',alpha=0.6)
    ## ax1.axvline(x=2995.,color='k',linestyle='dotted',alpha=0.6)

    #axis limits to ensure consistent plots across fields
    #N2 ylims
    #ax.set_ylim(1E-5,10)

    #S2 ylims
    s2_ymin=1E-7
    s2_ymax=1E-2
    ## ax1.set_ylim(s2_ymin,s2_ymax)
    
    ## #ax1.set_xlim(90.,4007.) 

    ## ax1.set_ylabel('$C_{\ell}[K^2]$',fontsize='medium')
    ## ax1.tick_params(labelsize='small')
    ## #ax1.legend([(qz_mark,qz_lin),(uz_mark,uz_lin),(qz_nonoise_mark,qz_nonoise_lin),(uz_nonoise_mark,uz_nonoise_lin),q_scaled_noise_lin,u_scaled_noise_lin,beam_cut],["QQ","UU","QQ-noise","UU-noise","QQ noise power","UU noise power","beamwidth scale"],fontsize='medium',loc=0)
    ## ax1.legend([(qz_nonoise_mark,qz_nonoise_lin),(uz_nonoise_mark,uz_nonoise_lin),beam_cut],["QQ-noise","UU-noise","beamwidth scale"],fontsize='medium',loc=0)
    ## ax1.set_xscale('log') 
    ## ax1.set_yscale('log')
    ## ax1.text(30, 1E-6, '$Q/U$ noise ratio = {:.2f} '.format(qu_ratio))

    #Plot 2D QQ top middle
    ax2=fig.add_axes([0.05,0.38,0.95,0.28])

    from matplotlib.patches import Circle

    ## d=Circle((256,256),radius=24)
    ## e=Circle((256,256),radius=48)
    ## f=Circle((256,256),radius=95)
    ## g=Circle((256,256),radius=142)
    ## h=Circle((256,256),radius=190)

    ## d.set_facecolor("none")
    ## e.set_facecolor("none")
    ## f.set_facecolor("none")
    ## g.set_facecolor("none")
    ## h.set_facecolor("none")

    ## d.set_alpha(0.5)
    ## e.set_alpha(0.5)
    ## f.set_alpha(0.5)
    ## g.set_alpha(0.5)
    ## h.set_alpha(0.5)

    ## circ1=str(int(ell_r[512,512+24]))
    ## circ2=str(int(ell_r[512,512+48]))
    ## circ3=str(int(ell_r[512,512+95]))
    ## circ4=str(int(ell_r[512,512+142]))
    ## circ5=str(int(ell_r[512,512+190]))

    im1=ax2.imshow(np.log10(np.abs(qq_scaled_ma))[256:768,256:768], clim=(np.log10(s2_ymin),np.log10(s2_ymax)))
    plt.gca()
    locs2=np.asarray([-190, -94, 0, 94, 189])
    labs2=[4000, 2000, 0, 2000, 4000]
    #print ell_r[xlocs2], ell_r[ylocs2]
    #plt.xticks(locs2+256,ell_r[512,locs2+512].astype(int), size='small')
    #plt.yticks(locs2+256,ell_r[locs2+512,512].astype(int),size='small')
    plt.xticks(locs2+256,labs2, size='large')
    plt.yticks(locs2+256,labs2,size='large')
    plt.xlabel('$\ell_x$', size=16, weight='bold')
    plt.ylabel('$\ell_y$', size=16, rotation='horizontal', weight='bold')
    divider1 = make_axes_locatable(ax2)
    cax1 = divider1.append_axes("right", size="5%", pad=0.05)
    cbar1=plt.colorbar(mappable=im1, cax=cax1)
    cbar1.set_label('$\mathrm{log}_{10} \; QQ [\mathrm{K}^2]$',size=16)
    cbar1.ax.tick_params(labelsize=12)

    #ax2.text(12,32,'circles at $\ell$ = '+circ1+', '+circ2+', '+circ3+', '+circ4+', '+circ5,fontsize=7)
    txt2=ax2.text(51,461,'QQ', size=18,color='white', weight='bold')
    #plt.setp(txt2,path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    ## ax2.add_artist(d)
    ## ax2.add_artist(e)
    ## ax2.add_artist(f)
    ## ax2.add_artist(g)
    ## ax2.add_artist(h)
    
    #ax2.set_xticklabels([])
    #ax2.set_yticklabels([])
    #ax2.set_xticks([])
    #ax2.set_yticks([])

    #Plot BB top right

    ax3=fig.add_axes([0.05,0.05,0.95,0.28])

    ## d1=Circle((256,256),radius=24)
    ## e1=Circle((256,256),radius=48)
    ## f1=Circle((256,256),radius=95)
    ## g1=Circle((256,256),radius=142)
    ## h1=Circle((256,256),radius=190)

    ## d1.set_facecolor("none")
    ## e1.set_facecolor("none")
    ## f1.set_facecolor("none")
    ## g1.set_facecolor("none")
    ## h1.set_facecolor("none")

    ## d1.set_alpha(0.5)
    ## e1.set_alpha(0.5)
    ## f1.set_alpha(0.5)
    ## g1.set_alpha(0.5)
    ## h1.set_alpha(0.5)

    im2=ax3.imshow(np.log10(np.abs(qq_nonoise))[256:768,256:768], clim=(np.log10(s2_ymin),np.log10(s2_ymax)))
    plt.gca()
    #locs2=np.asarray([-256, -142, -95, -48, 0, 48, 95, 142,256])
    #print ell_r[xlocs2], ell_r[ylocs2]
    plt.xticks(locs2+256,labs2, size='large')
    plt.yticks(locs2+256,labs2,size='large')
    plt.xlabel('$\ell_x$', size=16, rotation='horizontal', weight='bold')
    plt.ylabel('$\ell_y$', size=16, rotation='horizontal', weight='bold')
    divider2 = make_axes_locatable(ax3)
    cax2 = divider2.append_axes("right", size="5%", pad=0.05)
    cbar2=plt.colorbar(mappable=im2, cax=cax2)
    cbar2.set_label('$\mathrm{log}_{10} \; QQ [\mathrm{K}^2]$',size=16)
    cbar2.ax.tick_params(labelsize=12)

    #ax3.text(12,32,'circles at $\ell$ = '+circ1+', '+circ2+', '+circ3+', '+circ4+', '+circ5,fontsize=7)
    txt3=ax3.text(51,461,'QQ w/ cut', size=18, color='white', weight='bold')
    #plt.setp(txt3,path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])

    ## ax3.add_artist(d1)
    ## ax3.add_artist(e1)
    ## ax3.add_artist(f1)
    ## ax3.add_artist(g1)
    ## ax3.add_artist(h1)

    ## ax3.set_xticklabels([])
    ## ax3.set_yticklabels([])
    ## ax3.set_xticks([])
    ## ax3.set_yticks([])

    


    if beam == True:
        fig.savefig("/Users/leclercq/thesis/aps/figures/images/stripenoise_"+field+"_c"+str(chunk)+".pdf",dpi=200, bbox_inches='tight')
    else:
        fig.savefig("/Users/leclercq/thesis/aps/figures/images/stripenoise_"+field+"_c"+str(chunk)+".pdf",dpi=200, bbox_inches='tight')


##############    START OF MAIN PROGRAM    ################

#input arguments
q_in=sys.argv[1]
u_in=sys.argv[2]
qnoise_in=sys.argv[3]
unoise_in=sys.argv[4]
pol_in=sys.argv[5]
#latlon_out=sys.argv[4]
field=sys.argv[6]
width=int(sys.argv[7])
beam=int(sys.argv[8]) #this needs to be either 1 or 0; setting it to one turns beam removal on
#chunk=int(sys.argv[9]) 
noise_chunk=int(sys.argv[9])

#setting up some constants:
#define pixel step and pixel area in rad
pixstep=pi/(60.*180.)
pixarea=pixstep**2
area1024=1024.*1024.*pixarea
area_width=width*width*pixarea
#width of chunk in degrees
width_deg=width/60.


q_im, u_im, qnoise, unoise, xw, yw, xnw, ynw, w = import_qu(q_in, u_in, qnoise_in, unoise_in)

clims_pix_x,clims_pix_y,center_pix_x,center_pix_y,clims_pix_xn,clims_pix_yn,center_pix_xn,center_pix_yn,pad,nochunks, nochunks_n = make_chunks(width, xw, yw, xnw, ynw)

ell_r, ell_hist, ell_cut, ell_hist_cut, phi_ell, phi_cut, bins_axis, ft_beam, ft_taper, phi_ell, bins = fft_prep(width, pixstep)

#Work on selected noise chunk

qqnoise_scaled, uunoise_scaled, qqnoise_scaled_ma, uunoise_scaled_ma = chunk_fts(qnoise, unoise, clims_pix_yn, clims_pix_xn, width, ft_taper, pad, noise_chunk, area_width, beam)

qnz,unz,nnz,bins_nz = binning(qqnoise_scaled_ma, uunoise_scaled_ma, ell_cut, ell_hist_cut, bins, bins_axis)

#Loop over map chunks

#for i in range (nochunks):
i=0

qq_scaled, uu_scaled, qq_scaled_ma, uu_scaled_ma = chunk_fts(q_im, u_im, clims_pix_y, clims_pix_x, width, ft_taper, pad, i, area_width, beam)

qz,uz,nz,bins_z = binning(qq_scaled_ma, uu_scaled_ma, ell_cut, ell_hist_cut, bins, bins_axis)

#Scale noise APS to map level beyond beam scale, subtract from masked 2D FFT, 

qq_nonoise,uu_nonoise,qz_nonoise,uz_nonoise,q_scaled_noise,u_scaled_noise, qu_ratio=noise_scale(qq_scaled_ma, uu_scaled_ma, qqnoise_scaled_ma, qz, uz, qnz, bins_z)

#Plot

plot_all(q_in, u_in, pol_in, w, center_pix_x, center_pix_y, width_deg, bins_z, qz, uz, qz_nonoise, uz_nonoise, qnz, q_scaled_noise, u_scaled_noise, qq_scaled, qq_scaled_ma, i, field, beam, qu_ratio)



#Scale noise APS to map level beyond beam scale, subtract from masked 2D FFT, 



    




#print bins_axis[30:70]-bins_axis[50]
    #print np.log10(bins_axis[30:70])-np.log10(bins_axis[50])    
    
#Fit the power-law linearly, for intermediate ells, subtract center value
    #slope,offset,c,d,e=stats.linregress(np.log10(bins_axis[30:70])-np.log10(bins_axis[50]),np.log10((ee_average[30:70]+bb_average[30:70])/2))

    #Fit modified power law accounting for effect of beam and noise
    ## print "Pivot point is ell=",bins_axis[50]
    
    ## def apsfit_ee(x,a,w,n):
    ##     return ee_average[50]*w*(x/bins_axis[50])**a + n
    ## def apsfit_bb(x,a,w,n):
    ##     return bb_average[50]*w*(x/bins_axis[50])**a + n
    ## def apsfit_eb(x,a,w,n):
    ##     return eb[50]*w*(x/bins_axis[50])**a + n

    ## pfit_ee,pcov_ee=sciopt.curve_fit(apsfit_ee,bins_axis[nonzero_ee],ee_average[nonzero_ee])
    ## pfit_bb,pcov_bb=sciopt.curve_fit(apsfit_bb,bins_axis[nonzero_bb],bb_average[nonzero_bb])
    ## pfit_eb,pcov_eb=sciopt.curve_fit(apsfit_eb,bins_axis[nz],eb)

    ## print "ee fit params (index, weight, noise):", pfit_ee
    ## print "bb fit params (index, weight, noise):", pfit_bb
    ## print "eb fit params (index, weight, noise):", pfit_eb

    #print slope,offset
    #pwr_label="power law, index="+str(slope)[:6]

#Store the center pixel galactic coordinates and APS slope value

    ## lon,lat=w.all_pix2world(center_pix_x[i],center_pix_y[i],0)

    ## icrs_coords=SkyCoord(lon,lat,frame='icrs',unit='deg')

    ## gal_coords=icrs_coords.galactic

    ## row=field+'{0:2d} '.format(i) +gal_coords.to_string('decimal')+' {0:.3f} '.format(slope) +'\n'

    ## if os.path.isfile(latlon_out):
    ##     f=open(latlon_out, 'a')
    ##     f.write(row)
    ##     f.close()
    ## else:
    ##     f=open(latlon_out,'w')
    ##     f.write('#Field  Chunk  l   b   slope\n')
    ##     f.write(row)
    ##     f.close()
