#Make 2D map and noise illustration for thesis


import numpy as np
from numpy import ma
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

    if i.endswith(".npz"):

        print "Working on "+ i

        field=i.split("_")[0]
        chunk=i.split("_")[1]

        data=np.load(i)

        p_in='/Users/leclercq/galfacts/inpainting/'+field+'_P_inpainted.fits'

        cpix1=data['cpix']
        eemap=data['ee']
        eenoise=data['een']
        bbmap=data['bb']
        bbnoise=data['bbn']

        ell_r=data['ell_r']

        xm=ym=np.arange(1024)
        Xm,Ym=np.meshgrid(xm,ym)
        Ym=Ym.astype('f')
        Xm=Xm.astype('f')


        eemap_ma=ma.masked_where((ell_r>1000) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),eemap)

        bbmap_ma=ma.masked_where((ell_r>1000) & (np.abs((Ym-512)/(Xm-512))<np.abs(np.tan(np.pi/12.))),bbmap)
        
        #Plot! 2x1 figure, map on the left, noise on the right

        qq_ymin=1E-9
        qq_ymax=1E-3

        fig=plt.figure(figsize=(18,8))

        f1 = aplpy.FITSFigure(p_in, figure=fig, subplot=[0.05,0.05,0.26,0.9])
        f1.tick_labels.set_font(size='large')
        #f1.axis_labels.set_font(size='small')
        #f1.axis_labels.hide()
        #f1.tick_labels.hide()
        f1.show_colorscale(cmap='afmhot', vmin=-0.05, vmax=0.10)
        x_coord,y_coord=f1.pixel2world(cpix1[0],cpix1[1])
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
        f1.add_label(0.1,0.1,'P',relative=True, color='white', size='18', weight='bold')
        ax2=fig.add_axes([0.38,0.05,0.26,0.9])
        im1=ax2.imshow(np.log10(np.abs(eemap_ma))[256:768,256:768], clim=(np.log10(qq_ymin),np.log10(qq_ymax)))
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


        
        ax4=fig.add_axes([0.71,0.05,0.26,0.9])
        im1=ax4.imshow(np.log10(np.abs(bbmap_ma))[256:768,256:768], clim=(np.log10(qq_ymin),np.log10(qq_ymax)))
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


        
        #plt.tight_layout()
        #plt.show()

        #outfile=

        fig.savefig("/Users/leclercq/galfacts/aps/final_results/plots/"+field+"_"+chunk+"_2d_pol_bowtie.pdf",dpi=150, bbox_inches='tight')
