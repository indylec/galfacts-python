#This is code to test the E-mode and B-mode noise simulation


import numpy as np 
import matplotlib
matplotlib.use('MacOSX')
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
#from scipy import stats
from scipy import optimize as sciopt
import os.path
pi=np.pi

################ MODULES #############

def cumsum_sma(array, period):
    ret = np.cumsum(array, dtype=float)
    ret[period:] = ret[period:] - ret[:-period]
    return ret[period - 1:] / period



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

    xm=ym=np.arange(1024)
    Xm,Ym=np.meshgrid(xm,ym)
    Ym=Ym.astype('f')
    Xm=Xm.astype('f')

    xyrad=np.sqrt((Xm-512)**2+(Ym-512)**2)
    
    ell_cut=ma.masked_where((xyrad>48) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),ell_r)
    
    ell_max=np.max(ell_cut)
    print ell_max

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


def chunk_fts(q_im,u_im,clims_pix_y,clims_pix_x,width,ft_taper,pad,i,area_width,beam,alfa=0.):
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

#Divide out beam
    if beam == True:
        qft_final=qft_final*ft_beam
        uft_final=uft_final*ft_beam
    
#Calculate E and B mode functions
    emode=qft_final*np.cos(2.*phi_ell+alfa)+uft_final*np.sin(2.*phi_ell+alfa)
    bmode=-qft_final*np.sin(2.*phi_ell+alfa)+uft_final*np.cos(2.*phi_ell+alfa)

#Calculate QQ and UU functions

#compute correlations (EE,BB)
    #qq=qft_final*np.conj(qft_final)
    #uu=uft_final*np.conj(uft_final)
    ee=emode*np.conj(emode)
    bb=bmode*np.conj(bmode)

#account for size of array:
    #qq_scaled=qq*area_width/1024**2
    #uu_scaled=uu*area_width/1024**2
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

    ee_scaled_ma=ma.masked_where((xyrad>48) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),ee_scaled)
    bb_scaled_ma=ma.masked_where((xyrad>48) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),bb_scaled)
    

    return ee_scaled, bb_scaled, ee_scaled_ma, bb_scaled_ma


################### MAIN #####################

q_in=sys.argv[1]
u_in=sys.argv[2]
qnoise_in=sys.argv[3]
unoise_in=sys.argv[4]
noise_chunk=int(sys.argv[5])
im_chunk=int(sys.argv[6])
beam=int(sys.argv[7])


#setting up some constants:
#define pixel step and pixel area in rad
pixstep=pi/(60.*180.)
width=900
pixarea=pixstep**2
area1024=1024.*1024.*pixarea
area_width=width*width*pixarea
#width of chunk in degrees
width_deg=width/60.
#various alfa angles
alfa_45=45.
alfa_45_rad=alfa_45*np.pi/180.
alfa_60=60.
alfa_60_rad=alfa_60*np.pi/180.
alfa_15=15.
alfa_15_rad=alfa_15*np.pi/180.

#get q u images and noise
q_im, u_im, qnoise, unoise, xw, yw, xnw, ynw, w = import_qu(q_in, u_in, qnoise_in, unoise_in)

#print np.nanstd(qnoise)
#print np.nanstd(unoise)

#divide up into chunks
clims_pix_x,clims_pix_y,center_pix_x,center_pix_y,clims_pix_xn,clims_pix_yn,center_pix_xn,center_pix_yn,pad,nochunks, nochunks_n = make_chunks(width, xw, yw, xnw, ynw)

#prep the ell-space, histogram bins etc
ell_r, ell_hist, ell_cut, ell_hist_cut, phi_ell, phi_cut, bins_axis, ft_beam, ft_taper, phi_ell, bins = fft_prep(width, pixstep)

#get the map fts
ee_scaled,bb_scaled, ee_scaled_ma,bb_scaled_ma = chunk_fts(q_im,u_im, clims_pix_y, clims_pix_x, width, ft_taper, pad, im_chunk, area_width, beam=False)

#exit()

rell=5000.
rdell=100.
rpix=rell*48./1000.
rdpix=rdell*48./1000.
cmin=-9
cmax=-4

rellmax=np.amax(ell_r[256:768,256:768])


#get azimuth angles and sort
angle=np.ravel(phi_ell[(ell_r>rell)&(ell_r<rell+rdell)]) 
sort=np.argsort(angle)
anglesort=angle[sort]


#bb_noise_circle=np.ravel(bbnoise_scaled[(ell_r>180/(3.5/60.)) & (ell_r<180/(3.5/60.)+500)])
#define boxcar kernels
box50=Box1DKernel(50)
box2d50=Box2DKernel(50)
box2d10=Box2DKernel(10)

#take azimuthal profile and smooth with boxcar
ee_circle=np.ravel(ee_scaled[(ell_r>rell)&(ell_r<rell+rdell)])
ee_circle_smooth=convolve_fft(np.abs(ee_circle[sort]),box50)

#smooth 2d map - unnecessary

#ee_scaled_sm=convolve(np.abs(ee_scaled),box2d10)
#ee_scaled_sm_ma=ma.masked_where(ell_cut,ee_scaled)

#find local minima to get noise shape for one half-period
lim1=np.argmin(ee_circle_smooth[:2000])
lim2=np.argmin(ee_circle_smooth[2000:4000])+2000

#anglesort_sm=anglesort[49:]

#Fit sin function to one half-period of noisy map

def noisefit_ee(x,a,b):
    return np.sqrt(a**2*np.cos((150.*np.pi/180.)+2*x)**2+b**2*np.sin((150.*np.pi/180.)+2*x)**2)

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
qnoise*=a
unoise*=b

#Get the noise Emodes and bmodes
eenoise_scaled, bbnoise_scaled,eenoise_scaled_ma,bbnoise_scaled_ma = chunk_fts(qnoise, unoise, clims_pix_yn, clims_pix_xn, width, ft_taper, pad, noise_chunk, area_width, beam=False, alfa=150.*np.pi/180.)

#take a radial profile and smooth
ee_noise_circle=np.ravel(eenoise_scaled[(ell_r>rell)&(ell_r<rell+rdell)])
ee_noise_circle_smooth=convolve_fft(np.abs(ee_noise_circle[sort]),box50)

#scale the profile to the fitted function
ee_scale=np.amax(noisefit_ee(anglesort,a,b))/np.amax(np.sqrt(ee_noise_circle_smooth))
#ee_scale=k
#use the scale for the 2d map and smoothed profile
eenoise_scaled_ma=np.sqrt(np.abs(eenoise_scaled_ma))*ee_scale
eenoise_scaled=np.sqrt(np.abs(eenoise_scaled))*ee_scale
ee_noise_circle=np.sqrt(ee_noise_circle)*ee_scale
ee_noise_circle_smooth=np.sqrt(ee_noise_circle_smooth)*ee_scale

#smooth the 2d noise - unneccessary
#eenoise_scaled_sm=convolve(np.abs(eenoise_scaled),box2d10)

ee_circle_sort=ee_circle[sort]

#construct 2d noise fit
eenoise_2dfit=noisefit_ee(phi_ell,a,b)

#2D noise removal - unnecessary
#ee_nonoise=ee_scaled_ma
#ee_nonoise-=eenoise_scaled_sm

#ee_nonoise_sm=ee_scaled_sm
#ee_nonoise_sm-=eenoise_scaled_sm

#ee_nonoise_fit=ee_scaled_ma
#ee_nonoise_fit-=eenoise_2dfit

#ee_nonoise_fit_sm=ee_scaled_sm
#ee_nonoise_fit_sm-=eenoise_2dfit

#get 1d differences: map - fit
ee_diff_fit=ee_circle_sort-noisefit_ee(anglesort,a,b)
ee_diff_fit_smoothed=np.abs(convolve_fft(ee_diff_fit,box50))

#map - smoothed noise
ee_diff_noise=ee_circle_sort-ee_noise_circle_smooth
ee_diff_noise_smoothed=np.abs(convolve_fft(ee_diff_noise,box50))

#smoothed map - smoothed noise
ee_smooth_diff=np.abs(ee_circle_smooth-ee_noise_circle_smooth)



#ee_diff_circle=ee_circle-ee_noise_circle

#print noise_circle

plt.plot(anglesort,np.log10(np.abs(ee_noise_circle[sort])),'g-',alpha=0.5)
#plt.plot(np.abs(bb_noise_circle),alpha=0.5)
plt.plot(anglesort,np.log10(np.abs(ee_circle_sort)),'b-',alpha=0.5)
#plt.plot(angle[sort],np.log10(np.abs(ee_circle-ee_noise_circle)),'b-',alpha=0.4)#noisefit_ee(angle[sort],*p0),'k-',alpha=0.5)
plt.plot(anglesort,np.log10(ee_circle_smooth),'b')
plt.plot(anglesort,np.log10(ee_noise_circle_smooth),'g')
plt.plot(anglesort,np.log10(noisefit_ee(anglesort,a,b)),'r-')#*fit2[0])),'r-')
plt.plot(anglesort,np.log10(ee_diff_fit_smoothed),'k-',alpha=0.8)
#plt.plot(np.log10(ee_diff_noise_smoothed),'k:',alpha=0.8)
#plt.plot(np.log10(ee_smooth_diff),'k-.',alpha=0.8)
#plt.plot(angle[sort],np.log10(noisefit_ee(angle[sort],*p0)),'k-')
plt.xlabel('Azimuthal angle (rad)')
plt.ylabel('log10 C_ell (K^2)')
#plt.savefig('/Users/leclercq/Desktop/test_figs/ee_bb_noise_data_ell'+str(rell)+'.pdf')
#plt.show()
#plt.clf()

######### Plot the images #####

f1, ax1=plt.subplots(1,1)
f2, ax2=plt.subplots(1,1)
#f3, ax3=plt.subplots(1,1)
#f4, ax4=plt.subplots(1,1)
#f5, ax5=plt.subplots(1,1)
#f6, ax6=plt.subplots(1,1)

#f,(ax1, ax2)=plt.subplots(1,2)
#f,(ax1, ax2, ax3)=plt.subplots(1,3)
#f,(ax1, ax2, ax3, ax4)=plt.subplots(1,4)
#f,(ax1, ax2, ax3, ax4, ax5)=plt.subplots(1,5)



d1=Circle((256,256),radius=rpix,lw=0.5)
d1.set_facecolor("none")
e1=Circle((256,256),radius=rpix+rdpix,lw=0.5)
e1.set_facecolor("none")
im1=ax1.imshow(np.log10(np.abs(ee_scaled[256:768,256:768])),clim=(cmin,cmax))
plt.colorbar(mappable=im1, ax=ax1)
ax1.add_artist(d1)
ax1.add_artist(e1)
ax1.set_title('E-mode from map')

ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.set_xticks([])
ax1.set_yticks([])

d2=Circle((256,256),radius=rpix,lw=0.5)
d2.set_facecolor("none")
e2=Circle((256,256),radius=rpix+rdpix,lw=0.5)
e2.set_facecolor("none")
im2=ax2.imshow(np.log10(np.abs(eenoise_scaled[256:768,256:768])),clim=(cmin,cmax))
plt.colorbar(mappable=im2, ax=ax2)
ax2.add_artist(d2)
ax2.add_artist(e2)
ax2.set_title('Noise scaled to fit')

ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.set_xticks([])
ax2.set_yticks([])

## d3=Circle((256,256),radius=rpix,lw=0.5)
## d3.set_facecolor("none")
## e3=Circle((256,256),radius=rpix+rdpix,lw=0.5)
## e3.set_facecolor("none")
## im3=ax3.imshow(np.log10(np.abs(eenoise_2dfit[256:768,256:768])),clim=(cmin,cmax))
## plt.colorbar(mappable=im3, ax=ax3)
## ax3.add_artist(d3)
## ax3.add_artist(e3)
## ax3.set_title('Noise fit from map')

## ax3.set_xticklabels([])
## ax3.set_yticklabels([])
## ax3.set_xticks([])
## ax3.set_yticks([])

## d4=Circle((256,256),radius=rpix,lw=0.5)
## d4.set_facecolor("none")
## e4=Circle((256,256),radius=rpix+rdpix,lw=0.5)
## e4.set_facecolor("none")
## im4=ax4.imshow(np.log10(np.abs(ee_nonoise_fit_sm[256:768,256:768])),clim=(cmin,cmax))
## plt.colorbar(mappable=im4, ax=ax4)
## ax4.add_artist(d4)
## ax4.add_artist(e4)
## ax4.set_title('Smoothed map - fit')

## ax4.set_xticklabels([])
## ax4.set_yticklabels([])
## ax4.set_xticks([])
## ax4.set_yticks([])

## d5=Circle((256,256),radius=rpix,lw=0.5)
## d5.set_facecolor("none")
## e5=Circle((256,256),radius=rpix+rdpix,lw=0.5)
## e5.set_facecolor("none")
## im5=ax5.imshow(np.log10(np.abs(eenoise_scaled_sm[256:768,256:768])),clim=(cmin,cmax))
## plt.colorbar(mappable=im5, ax=ax5)
## ax5.add_artist(d5)
## ax5.add_artist(e5)
## ax5.set_title('2D noise fit from map')

## ax5.set_xticklabels([])
## ax5.set_yticklabels([])
## ax5.set_xticks([])
## ax5.set_yticks([])

## d6=Circle((256,256),radius=rpix,lw=0.5)
## d6.set_facecolor("none")
## e6=Circle((256,256),radius=rpix+rdpix,lw=0.5)
## e6.set_facecolor("none")
## im6=ax6.imshow(np.log10(np.abs(ee_nonoise_sm[256:768,256:768])),clim=(cmin,cmax))
## plt.colorbar(mappable=im6, ax=ax6)
## ax6.add_artist(d6)
## ax6.add_artist(e6)
## ax6.set_title('Map-2D fit')

## ax6.set_xticklabels([])
## ax6.set_yticklabels([])
## ax6.set_xticks([])
## ax6.set_yticks([])

#plt.savefig("/Users/leclercq/Desktop/test_figs/een_bbn_noise_test1.pdf",dpi=300,bbox_inches="tight")

plt.show()

plt.clf()



####################  ALL THE FT STUFF IS BELOW HERE  ####################

#een_ft=np.fft.fft(ee_noise_circle[sort]**2,norm="ortho")
#ee_ft=np.fft.fft(ee_circle[sort]**2,norm="ortho")
#ang_freq=np.fft.fftfreq(angle[sort].size,d=2*np.pi/angle[sort].size)

#ee_ft_filt=ee_ft[(np.abs(ang_freq)<123.)]
#een_ft_filt=een_ft[(np.abs(ang_freq)<123.)]

#een_ft_filt_scaled=(ee_ft_filt/een_ft_filt
#plt.plot(ang_freq,np.abs(ee_ft))
#plt.plot(ang_freq,np.abs(een_ft))
#plt.plot(ang_freq[(np.abs(ang_freq)<123.)],ee_ft_filt)
#plt.plot(ang_freq[(np.abs(ang_freq)<123.)],een_ft_filt)
#plt.xlabel("angular frequency (rad^-1)")
#plt.show()
#plt.clf()

#print np.amax(np.abs(een_ft))
#print np.abs(een_ft)[np.argsort(np.abs(een_ft))[-2]]
#print np.abs(een_ft)[np.argsort(np.abs(een_ft))[-3]]


#model_ft=np.zeros(een_ft.size,dtype='complex')
#model_ft[np.argsort(np.abs(een_ft))[-2]]=ee_ft[np.argsort(np.abs(een_ft))[-2]]
#model_ft[np.argsort(np.abs(een_ft))[-3]]=ee_ft[np.argsort(np.abs(een_ft))[-3]]
#model_ft[np.argsort(np.abs(een_ft))[-4]]=ee_ft[np.argsort(np.abs(een_ft))[-4]]
#model_ft[np.argsort(np.abs(een_ft))[-5]]=ee_ft[np.argsort(np.abs(een_ft))[-5]]



#print ee_circle_smooth.size
#print ee_circle.size

#model=np.fft.ifft(model_ft,norm="ortho")
#posmodel=model.real
#posmodel[(posmodel<0)]=0.
#plt.plot(angle[sort],posmodel,'r-')
#plt.plot(angle[sort],np.abs(ee_noise_circle[sort])**2,'g-',alpha=0.5)
## #plt.plot(np.abs(bb_noise_circle),alpha=0.5)
#plt.plot(np.log10(np.abs(ee_circle[sort])**2),'b-',alpha=0.5)
##plt.plot(angle[sort],np.abs(ee_circle[sort])**2-posmodel,'r-')
#plt.plot(np.log10(ee_circle_smooth),'k',alpha=0.5)
## #plt.plot(angle[sort],np.log10(np.abs(ee_circle-ee_noise_circle)),'b-',alpha=0.4)#noisefit_ee(angle[sort],*p0),'k-',alpha=0.5)
## #plt.plot(angle[sort],np.log10(noisefit_ee(angle[sort],*fit[0])),'r-')
## #plt.plot(angle[sort],np.log10(noisefit_ee(angle[sort],*p0)),'k-')
## plt.xlabel('Azimuthal angle (rad)')
## plt.ylabel('log10 C_ell (K^2)')
#plt.show()
#plt.clf()


#plt.plot(ang_freq,ee_ft)
#plt.show()
#plt.clf()
#een_filt=np.fft.fft(een_ft_filt,norm="ortho")
#plt.plot(een_filt)
#plt.show()


#exit()



