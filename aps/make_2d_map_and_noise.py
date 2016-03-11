#Make 2D map and noise illustration for thesis


import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.convolution import convolve,convolve_fft, Box1DKernel,Box2DKernel
import aplpy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
from mpl_toolkits.axes_grid1 import make_axes_locatable
pi=np.pi
import sys
import os
from matplotlib import rc

rc('font',family='serif')
rc('text', usetex=True)

####### MAIN ######

#qin1='/Users/leclercq/galfacts/inpainting/S1_Q_inpainted.fits'
#uin1='/Users/leclercq/galfacts/inpainting/S1_U_inpainted.fits'

#qin2='/Users/leclercq/galfacts/inpainting/N2_Q_inpainted.fits'
#uin2='/Users/leclercq/galfacts/inpainting/N2_U_inpainted.fits'
width_deg=900./60.

for i in os.listdir(os.getcwd()):

    #if i.endswith(".npz"):

    if "N2_c3" in i:
        print "Working on "+ i

        field=i.split("_")[0]
        chunk=i.split("_")[1]

        data=np.load(i)

        #pin1='/Users/leclercq/galfacts/inpainting/'+field+'_P_inpainted.fits'

        cpix1=data['cpix']
        eemap=data['ee']
        eenoise=data['een']
        bbmap=data['bb']
        bbnoise=data['bbn']
        
        #Plot! 2x1 figure, map on the left, noise on the right

        qq_ymin=1E-9
        qq_ymax=1E-3

        fig=plt.figure(figsize=(8,6))

        ax2=fig.add_subplot(221)
        im1=ax2.imshow(np.log10(np.abs(eemap))[256:768,256:768], clim=(np.log10(qq_ymin),np.log10(qq_ymax)))
        plt.gca()
        locs2=np.asarray([-190, 0, 189])
        labs2=[4000, 0, 4000]
        #print ell_r[xlocs2], ell_r[ylocs2]
        #plt.xticks(locs2+256,ell_r[512,locs2+512].astype(int), size='small')
        #plt.yticks(locs2+256,ell_r[locs2+512,512].astype(int),size='small')
        plt.xticks(locs2+256,labs2, size='medium')
        plt.yticks(locs2+256,labs2,size='medium')
        plt.xlabel('$\ell_x$', size=16, weight='bold')
        plt.ylabel('$\ell_y$', size=16, rotation='horizontal', weight='bold')
        divider1 = make_axes_locatable(ax2)
        cax1 = divider1.append_axes("right", size="5%", pad=0.05)
        cbar1=plt.colorbar(mappable=im1, cax=cax1)
        cbar1.set_label('$\mathrm{log}_{10} \; EE \; [\mathrm{K}^2]$',size=14)
        cbar1.ax.tick_params(labelsize=12)
        #txt2=ax2.text(51,461,'QQ', size=18,color='white', weight='bold')


        ax3=fig.add_subplot(222)
        im1=ax3.imshow(np.log10(np.abs(eenoise))[256:768,256:768], clim=(np.log10(qq_ymin),np.log10(qq_ymax)))
        plt.gca()
        locs3=np.asarray([-190, 0, 189])
        labs3=[4000, 0, 4000]
        #print ell_r[xlocs3], ell_r[ylocs3]
        #plt.xticks(locs3+256,ell_r[512,locs3+512].astype(int), size='small')
        #plt.yticks(locs3+256,ell_r[locs3+512,512].astype(int),size='small')
        plt.xticks(locs3+256,labs3, size='medium')
        plt.yticks(locs3+256,labs3,size='medium')
        plt.xlabel('$\ell_x$', size=16, weight='bold')
        plt.ylabel('$\ell_y$', size=16, rotation='horizontal', weight='bold')
        divider1 = make_axes_locatable(ax3)
        cax1 = divider1.append_axes("right", size="5%", pad=0.05)
        cbar1=plt.colorbar(mappable=im1, cax=cax1)
        cbar1.set_label('$\mathrm{log}_{10} \; EE_{\mathrm{noise}} \; [\mathrm{K}^2]$',size=14)
        cbar1.ax.tick_params(labelsize=12)
        #txt3=ax3.text(51,461,'UU', size=18,color='white', weight='bold')

        ax4=fig.add_subplot(223)
        im1=ax4.imshow(np.log10(np.abs(bbmap))[256:768,256:768], clim=(np.log10(qq_ymin),np.log10(qq_ymax)))
        plt.gca()
        locs2=np.asarray([-190, 0, 189])
        labs2=[4000, 0, 4000]
        #print ell_r[xlocs2], ell_r[ylocs2]
        #plt.xticks(locs2+256,ell_r[512,locs2+512].astype(int), size='small')
        #plt.yticks(locs2+256,ell_r[locs2+512,512].astype(int),size='small')
        plt.xticks(locs2+256,labs2, size='medium')
        plt.yticks(locs2+256,labs2,size='medium')
        plt.xlabel('$\ell_x$', size=16, weight='bold')
        plt.ylabel('$\ell_y$', size=16, rotation='horizontal', weight='bold')
        divider1 = make_axes_locatable(ax4)
        cax1 = divider1.append_axes("right", size="5%", pad=0.05)
        cbar1=plt.colorbar(mappable=im1, cax=cax1)
        cbar1.set_label('$\mathrm{log}_{10} \; BB \; [\mathrm{K}^2]$',size=14)
        cbar1.ax.tick_params(labelsize=12)
        #txt2=ax2.text(51,461,'QQ', size=18,color='white', weight='bold')


        ax5=fig.add_subplot(224)
        im1=ax5.imshow(np.log10(np.abs(bbnoise))[256:768,256:768], clim=(np.log10(qq_ymin),np.log10(qq_ymax)))
        plt.gca()
        locs3=np.asarray([-190, 0, 189])
        labs3=[4000, 0, 4000]
        #print ell_r[xlocs3], ell_r[ylocs3]
        #plt.xticks(locs3+256,ell_r[512,locs3+512].astype(int), size='small')
        #plt.yticks(locs3+256,ell_r[locs3+512,512].astype(int),size='small')
        plt.xticks(locs3+256,labs3, size='medium')
        plt.yticks(locs3+256,labs3,size='medium')
        plt.xlabel('$\ell_x$', size=16, weight='bold')
        plt.ylabel('$\ell_y$', size=16, rotation='horizontal', weight='bold')
        divider1 = make_axes_locatable(ax5)
        cax1 = divider1.append_axes("right", size="5%", pad=0.05)
        cbar1=plt.colorbar(mappable=im1, cax=cax1)
        cbar1.set_label('$\mathrm{log}_{10} \; BB_{\mathrm{noise}} \; [\mathrm{K}^2]$',size=14)
        cbar1.ax.tick_params(labelsize=12)
        #txt3=ax3.text(51,461,'UU', size=18,color='white', weight='bold')

        plt.tight_layout()
        #plt.show()

        #outfile=

        fig.savefig("/Users/leclercq/galfacts/aps/final_results/plots/"+field+"_"+chunk+"2d_map_and_noise.pdf",dpi=150, bbox_inches='tight')
