import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import sys

x_pix=int(sys.argv[1])
y_pix=int(sys.argv[2])

x_pix_im=x_pix-4124
y_pix_im=y_pix-800




q_in='/Users/leclercq/galfacts/aps/rm/S1_Q_poster.fits'

u_in='/Users/leclercq/galfacts/aps/rm/S1_U_poster.fits'

phi_in='/Users/leclercq/galfacts/aps/rm/S1_phi_poster.fits'

weights_in='/Users/leclercq/galfacts/aps/rm/S1_weights.txt'

synrm_in='/Users/leclercq/galfacts/aps/rm/S1_4_syn_rm.fits'

rmtf_in='/Users/leclercq/galfacts/aps/rm/GALFACTS_S1_dicube.txt'

rm_in='/Users/leclercq/galfacts/3.1.2/fits_files/rm/GALFACTS_S1_RM.fits'

angle_in='/Users/leclercq/galfacts/3.1.2/fits_files/rm/GALFACTS_S1_angle.fits'

weights=np.loadtxt(weights_in)
goodch=np.nonzero(weights)[0]

phi_range,rmtf_re,rmtf_im=np.loadtxt(rmtf_in,unpack=True)
rmtf_abs=np.sqrt(rmtf_re**2+rmtf_im**2)

print "Getting fits data"

rm=fits.getdata(rm_in)
cuberm=rm[y_pix,x_pix]
print cuberm

angle0=fits.getdata(angle_in)
cubeangle=angle0[y_pix,x_pix]
print cubeangle

synrm=fits.getdata(synrm_in)
synrm=synrm[800:,0:895]
synrm=synrm[y_pix_im,x_pix_im]

qcube=fits.getdata(q_in)
qlos=qcube[:,y_pix_im,x_pix_im]
del qcube

ucube=fits.getdata(u_in)
ulos=ucube[:,y_pix_im,x_pix_im]
del ucube

phicube=fits.getdata(phi_in)
phi_los=phicube[:,y_pix_im,x_pix_im]
del phicube

print "done. Making l2, angles..."

nu_size=376.
dnu=-420000.
nuref=1524717952.

c2=299792458.**2
    
dl2 = c2/(dnu**2)

nu = np.arange(nu_size)*dnu+nuref

l2 = 0.5 * c2 * ((nu - 0.5 * dnu) ** -2 + (nu + 0.5 * dnu) ** -2)
l2 = np.flipud(l2)

cchan=nu_size/2

angle=0.5*np.arctan2(ulos,qlos)
cangle=angle[cchan]
target=synrm*(l2-l2[cchan])+cangle
npi=np.around((target-angle)/np.pi)
rm_angle=angle+np.pi*npi

l2g=l2[goodch]
qg=qlos[goodch]
ug=ulos[goodch]
rmangg=rm_angle[goodch]
#print rmangg

print"done. Plotting!"

fig=plt.figure(figsize=(20,14.88))

gs1 = gridspec.GridSpec(2,2)
gs1.update(wspace=0.,hspace=0.05)

ax1a=fig.add_subplot(gs1[0,:])
ax2=fig.add_subplot(gs1[1,0])
ax3=fig.add_subplot(gs1[1,1])

qline=ax1a.plot(l2g,qg,'b-',alpha=0.7,linewidth=5)
uline=ax1a.plot(l2g,ug,'g-',alpha=0.7,linewidth=5)
ax1a.set_ylabel('$\mathrm{Q},\mathrm{U} [\mathrm{K}]$',fontsize=35)
ax1a.set_xlabel('$\lambda^2$', fontsize=40)
ax1a.xaxis.tick_top()
ax1a.xaxis.set_label_position('top')
ax1a.tick_params(labelsize=22)
ax1a.xaxis.set_major_locator(MaxNLocator(prune='both'))


ax1b=ax1a.twinx()
angline=ax1b.plot(l2g,rmangg,'r-',alpha=0.7,linewidth=5)
fitline=ax1b.plot(l2g,cuberm*l2g+cubeangle,'k-',linewidth=5)
ax1b.set_ylabel('$\Psi [\mathrm{rad}]$', fontsize=40)
ax1b.tick_params(labelsize=22)

#ax1a.legend((qline,uline,angline),("Q","U","$\Psi$"),fontsize='15')

ax2.plot(phi_range,phi_los,'k-',linewidth=5)
ax2.set_ylabel('$|\mathrm{F}(\phi)|[\mathrm{K}]$',fontsize=35)
ax2.set_xlabel('$\phi[\mathrm{rad}/\mathrm{m}^2]$', fontsize=40)
ax2.tick_params(labelsize=22)
ax2.xaxis.set_major_locator(MaxNLocator(prune='both'))


ax3.plot(phi_range,rmtf_abs,'b-',linewidth=5)
ax3.set_ylabel('$|\mathrm{R}(\phi)|$',fontsize=35)
ax3.set_xlabel('$\phi[\mathrm{rad}/\mathrm{m}^2]$', fontsize=40)
ax3.yaxis.set_label_position('right')
ax3.yaxis.tick_right()
ax3.tick_params(labelsize=22)
ax3.xaxis.set_major_locator(MaxNLocator(prune='both'))

fig.savefig("/Users/leclercq/galfacts/aps/rm/S1_"+str(x_pix)+"_"+str(y_pix)+"v2_plot.png",dpi=300,bbox_inches="tight")





