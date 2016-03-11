#Output (using np.save) all the necessary results for a given field (for each chunk in that field):
# image center pixels, galactic coords at center, slope of power-law fit, 2d EE, 2d BB, 2d QQ, 2d UU, EEnoise, BB noise, 1D EE raw, 1D EE corrected, 1D BB raw, 1D BB corrected, 1D EE noise, ell_r, ell_cut, phi_ell, phi_cut, bin_centers

import numpy as np 
#import matplotlib
#matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Circle
import aplpy
import gaussbeam
import sys
from numpy import pi
from numpy import ma
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.convolution import convolve,convolve_fft, Box1DKernel,Box2DKernel
from scipy import stats
from scipy import optimize as sciopt
import os.path
pi=np.pi
from matplotlib import rc

#rc('font',family='serif')
#rc('text', usetex=True)

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

    w=WCS(q_in)

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

    #xm=ym=np.arange(1024)
    #Xm,Ym=np.meshgrid(xm,ym)
    #Ym=Ym.astype('f')
    #Xm=Xm.astype('f')

    #xyrad=np.sqrt((Xm-512)**2+(Ym-512)**2)
    
    ell_cut=ma.masked_where((ell_r>500) & (np.abs((elly)/(ellx))<np.abs(np.tan(np.pi/12.))),ell_r)
    
    ell_max=np.max(ell_cut)
    print ell_max

    #set up grid of phis
    phi_ell=np.arctan2(elly,ellx)
    phi_cut=ma.masked_where((ell_r>500) & (np.abs((elly)/(ellx))<np.abs(np.tan(np.pi/12.))),phi_ell)

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


def chunk_fts(q_im,u_im,clims_pix_y,clims_pix_x,width,ft_taper,pad,i,area_width,beam,alfa=0.):
    #print "Working on chunk",i
    
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

#Divide out beam
    if beam == True:
        qft_final=qft_final*ft_beam
        uft_final=uft_final*ft_beam
    
#Calculate E and B mode functions
    emode=qft_final*np.cos(2.*phi_ell+alfa)+uft_final*np.sin(2.*phi_ell+alfa)
    bmode=-qft_final*np.sin(2.*phi_ell+alfa)+uft_final*np.cos(2.*phi_ell+alfa)

#Calculate QQ and UU functions

#compute correlations (EE,BB)
    qq=qft_final*np.conj(qft_final)
    uu=uft_final*np.conj(uft_final)
    ee=emode*np.conj(emode)
    bb=bmode*np.conj(bmode)

#account for size of array:
    qq_scaled=qq*area_width/1024**2
    uu_scaled=uu*area_width/1024**2
    ee_scaled=ee*area_width/1024**2
    bb_scaled=bb*area_width/1024**2

    #mask out the striping features

    xm=ym=np.arange(1024)
    Xm,Ym=np.meshgrid(xm,ym)
    Ym=Ym.astype('f')
    Xm=Xm.astype('f')

    xyrad=np.sqrt((Xm-512)**2+(Ym-512)**2)
    
    #qq_scaled_ma=ma.masked_where((xyrad>48) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),qq_scaled)
    #uu_scaled_ma=ma.masked_where((xyrad>48) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),uu_scaled)

    ee_scaled_ma=ma.masked_where((ell_r>500) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),ee_scaled)
    bb_scaled_ma=ma.masked_where((ell_r>500) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),bb_scaled)
    

    return ee_scaled, bb_scaled, ee_scaled_ma, bb_scaled_ma, qq_scaled, uu_scaled

#Bin the C_l for E and B to calculate radial average
def binning(ee_scaled, bb_scaled, ell_r, ell_hist, bins, bins_axis):
    #print np.shape(ee_scaled.compressed()), np.shape(ell_r.compressed())   
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


################### MAIN #####################

q_in=sys.argv[1]
u_in=sys.argv[2]
qnoise_in=sys.argv[3]
unoise_in=sys.argv[4]
noise_chunk=int(sys.argv[5])
field=sys.argv[6]
#im_chunk=int(sys.argv[6])

pixstep=pi/(60.*180.)
width=900
pixarea=pixstep**2
area1024=1024.*1024.*pixarea
area_width=width*width*pixarea
#width of chunk in degrees
width_deg=width/60.


#get q u images and noise
q_im, u_im, qnoise, unoise, xw, yw, xnw, ynw, w = import_qu(q_in, u_in, qnoise_in, unoise_in)
print "Images read in"

#divide up into chunks
clims_pix_x,clims_pix_y,center_pix_x,center_pix_y,clims_pix_xn,clims_pix_yn,center_pix_xn,center_pix_yn,pad,nochunks, nochunks_n = make_chunks(width, xw, yw, xnw, ynw)
print "Divided into chunks"

#prep the ell-space, histogram bins etc
ell_r, ell_hist, ell_cut, ell_hist_cut, phi_ell, phi_cut, bins_axis, ft_beam, ft_taper, phi_ell, bins = fft_prep(width, pixstep)
print "FT prep done"


#Set the azimuthal profile parameters
rell=5000.
rdell=100.
rpix=rell*48./1000.
rdpix=rdell*48./1000.

#set theta
if field=="S1":
    theta=90.
else:
    theta=150.

#loop over chunks
for i in range (3,4):
    print "Working on chunk", i

    #get the map fts
    ee_scaled, bb_scaled, ee_scaled_ma, bb_scaled_ma, qq_scaled, uu_scaled = chunk_fts(q_im,u_im, clims_pix_y, clims_pix_x, width, ft_taper, pad, i, area_width, beam=False)

    #get azimuth angles and sort
    angle=np.ravel(phi_ell[(ell_r>rell)&(ell_r<rell+rdell)]) 
    sort=np.argsort(angle)
    anglesort=angle[sort]

    #define boxcar kernel
    box50=Box1DKernel(50)

    #take azimuthal profiles and smooth with boxcar
    ee_circle=np.ravel(ee_scaled[(ell_r>rell)&(ell_r<rell+rdell)])
    ee_circle_sort=ee_circle[sort]
    ee_circle_smooth=convolve_fft(np.abs(ee_circle_sort),box50)


    bb_circle=np.ravel(bb_scaled[(ell_r>rell)&(ell_r<rell+rdell)])
    bb_circle_sort=bb_circle[sort]
    bb_circle_smooth=convolve_fft(np.abs(bb_circle_sort),box50)

    #find local minima to get noise shape for one half-period
    if field=='S1':
        lim1=np.argmin(ee_circle_smooth)-1000
        lim2=np.argmin(ee_circle_smooth)+1000
    else:
        lim1=np.argmin(ee_circle_smooth[:2000])
        lim2=np.argmin(ee_circle_smooth[2000:4000])+2000

    #Fit ee sin function to one half-period of noisy ee map

    def noisefit_ee(x,a,b):
        return np.sqrt(a**2*np.cos((theta*np.pi/180.)+2*x)**2+b**2*np.sin((theta*np.pi/180.)+2*x)**2)

    guess_a=2.*np.nanmean(qnoise)/1024.
    guess_b=np.nanmean(unoise)/1024.
    guess_k=1.
    p0=[guess_a, guess_b]

    fit=sciopt.curve_fit(noisefit_ee,anglesort[lim1:lim2],ee_circle_smooth[lim1:lim2],p0=p0)
    fit2=sciopt.curve_fit(noisefit_ee,anglesort[lim1:lim2:],ee_circle_smooth[lim1:lim2],p0=fit[0])

    print fit2[0]

    #print fit[0]
    #set the noise level to differentiate Q & U
    a,b=fit2[0]
    qnoise_corrected=qnoise*a
    unoise_corrected=unoise*b

    #Get the noise Emodes and bmodes
    eenoise_scaled, bbnoise_scaled,eenoise_scaled_ma,bbnoise_scaled_ma, qq_noise, uu_noise = chunk_fts(qnoise_corrected, unoise_corrected, clims_pix_yn, clims_pix_xn, width, ft_taper, pad, noise_chunk, area_width, beam=False, alfa=theta*np.pi/180.)

    #take a radial profile and smooth
    ee_noise_circle=np.ravel(eenoise_scaled[(ell_r>rell)&(ell_r<rell+rdell)])
    ee_noise_circle_smooth=convolve_fft(np.abs(ee_noise_circle[sort]),box50)

    #scale the profile to the fitted function
    ee_scale=np.amax(noisefit_ee(anglesort,a,b))/np.amax(np.sqrt(ee_noise_circle_smooth))

    #use the scale for the 2d map and 1D smoothed profile
    eenoise_scaled_ma=np.sqrt(np.abs(eenoise_scaled_ma))*ee_scale
    eenoise_scaled=np.sqrt(np.abs(eenoise_scaled))*ee_scale
    ee_noise_circle=np.sqrt(ee_noise_circle)*ee_scale
    ee_noise_circle_smooth=np.sqrt(ee_noise_circle_smooth)*ee_scale

    bbnoise_scaled_ma=np.sqrt(np.abs(bbnoise_scaled_ma))*ee_scale
    bbnoise_scaled=np.sqrt(np.abs(bbnoise_scaled))*ee_scale

    anglesort=anglesort*180./np.pi
    anglesort_rad=anglesort*np.pi/180.


    fig=plt.figure(figsize=(10,6))
    ax=fig.add_subplot(111)
    noiselin,=ax.plot(anglesort,np.abs(ee_noise_circle[sort]),'g-',alpha=0.4)
    #ax.plot(np.abs(bb_noise_circle),alpha=0.5)
    datlin,=ax.plot(anglesort,np.abs(ee_circle_sort),'b-',alpha=0.3)
    #ax.plot(angle[sort],np.log10(np.abs(ee_circle-ee_noise_circle)),'b-',alpha=0.4)#noisefit_ee(angle[sort],*p0),'k-',alpha=0.5)
    datsm,=ax.plot(anglesort,ee_circle_smooth,'b')
    noisesm,=ax.plot(anglesort,ee_noise_circle_smooth,'g')
    fitsm,=ax.plot(anglesort,noisefit_ee(anglesort_rad,a,b),'r-')#*fit2[0])),'r-')
    #ax.plot(anglesort,np.log10(ee_diff_fit_smoothed),'k-',alpha=0.8)
    ax.set_yscale('log')
    ax.set_xlabel('$\phi \; [^{\circ}]$',fontsize='large')
    ax.set_ylabel('$C_{\ell} [\mathrm{K}^2]$',fontsize='large')
    ax.legend([datlin, datsm, noiselin, noisesm, fitsm], ['$EE$ profile', 'Smoothed $EE$', '$EE_{\mathrm{noise}}$ profile', 'Smoothed $EE_{\mathrm{noise}}$', 'Noise model'], fontsize='medium', loc=4)
    #plt.show()

    fig.savefig("/Users/leclercq/thesis/aps/figures/images/"+field+"_c"+str(i)+"_azimuth_profile.pdf",dpi=150, bbox_inches='tight')
    #bin the 2D map and the 2D noise

    exit()

    eenz, bbnz,nnz,bins_nz = binning(eenoise_scaled_ma, bbnoise_scaled_ma, ell_cut, ell_hist_cut, bins, bins_axis)

    eez,bbz,nz,bins_z=binning(ee_scaled_ma,bb_scaled_ma,ell_cut,ell_hist_cut, bins, bins_axis)

    #print eenz.size, eez.size

    #subtract the binned noise from the binned map

    final_bins=np.intersect1d(bins_nz,bins_z)
    #print final_bins.size


    ee1d=eez-eenz
    bb1d=bbz-bbnz
    #print ee1d
    #print bb1d

    eb1d=(ee1d+bb1d)/2.

    #print final_bins[30:62]
    #print np.log10(eb1d[30:62])

    #Fit the power law between ell=1000 and ell = 3000

    #print np.log10(final_bins[30:62])-np.log10(final_bins[46])

    slope,offset,c,d,e=stats.linregress(np.log10(final_bins[30:62])-np.log10(final_bins[46]),np.log10(eb1d[30:62]))

    #Put the center pixel and associated coordinates in handy arrays



    cpix=[center_pix_x[i],center_pix_y[i]]
    lon,lat=w.all_pix2world(cpix[0],cpix[1],0)
    icrs_coords=SkyCoord(lon,lat,frame='icrs',unit='deg')
    gal_coords=icrs_coords.galactic

    #Store all the things!

    outfile='/Users/leclercq/galfacts/aps/final_results/data/data_500/'+field+'_c'+str(i)+'_full_results_500'

    np.savez_compressed(outfile, cpix=cpix, gal_coords=gal_coords, slope=slope, ee=ee_scaled, bb=bb_scaled, qq=qq_scaled, uu=uu_scaled, een=eenoise_scaled, bbn=bbnoise_scaled, ee1draw=eez, ee1dnoise=eenz, ee1dcorr=ee1d, bb1draw=bbz, bb1dnoise=bbnz, bb1dcorr=bb1d, ell_r=ell_r, phi_ell=phi_ell, bins=final_bins, fit=fit2[0])

    print "Data written to "+outfile

