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
from scipy import stats
import sys
import os

#from matplotlib import rc

#rc('font',family='serif')
#rc('text', usetex=True)
field=sys.argv[1]
chunk=sys.argv[2]

width_deg=900./60.


for i in os.listdir(os.getcwd()):

    #if i.endswith(".npz"):
    if field+"_c"+chunk in i:

        print "Working on "+ i

        field=i.split("_")[0]
        chunk=i.split("_")[1]

        data=np.load(i)

        pol_in='/Users/leclercq/galfacts/inpainting/'+field+'_P_inpainted.fits'

        
        ee1dcorr=data['ee1dcorr']
        bb1dcorr=data['bb1dcorr']
        final_bins=data['bins']

        cpix1=data['cpix']
        
        final_bins=final_bins[:64]
        ee1dcorr=ee1dcorr[:64]
        bb1dcorr=bb1dcorr[:64]
        

        #print final_bins
        #exit()


        eb1d=(ee1dcorr+bb1dcorr)/2.
        slope,offset,c,d,e=stats.linregress(np.log10(final_bins[32:60])-np.log10(final_bins[45]),np.log10(eb1d[32:60]))

        fig=plt.figure(figsize=(18,6))

        ax=fig.add_axes([0.05,0.05,0.6,0.9])
        #ax.set_autoscale_on(False)

        ax.set_xlabel('$\ell$',fontsize='medium' )

        #bbcorr_lin,= ax.plot(final_bins,bb1dcorr,'b-',alpha=0.7, linewidth=1.5)
        #eecorr_lin,=ax.plot(final_bins,ee1dcorr,'r-',alpha=0.7, linewidth=1.5)

        eb1d_mark,=ax.plot(final_bins[:64],eb1d[:64],color='purple', marker='d',markersize=7,linestyle='')
        #eb1d_lin,= ax.plot(final_bins,eb1d,color='purple',alpha=0.9, linewidth=1.5)

        power_law,= ax.plot(final_bins[32:60],10**(slope*(np.log10(final_bins[32:60])-np.log10(final_bins[45]))+offset),'k-',linewidth=3)
        
        beam_cut =ax.axvline(x=180/(3.5/60.),color='k',linestyle='dashed',alpha=0.8)

        ymin=1E-7
        ymax=10.
        ax.set_ylim(ymin,ymax)
        ax.set_xlim(10,8000)

        ax.set_ylabel('$C_{\ell}[K^2]$',fontsize='large')
        ax.tick_params(labelsize='medium')
        
        ax.legend([eb1d_mark, power_law, beam_cut],["(E+B)/2", "Power-law fit", "beamwidth scale"],fontsize='medium',loc=0)
        ax.set_xscale('log') 
        ax.set_yscale('log')
        ax.text(15, 5E-6, r'$\beta$ = {:.2f} '.format(-slope),size='large')

        f1 = aplpy.FITSFigure(pol_in, figure=fig, subplot=[0.68,0.05,0.30,0.9])
        #f1.tick_labels.set_font(size='x-small')
        #f1.axis_labels.set_font(size='small')
        f1.axis_labels.hide()
        f1.tick_labels.hide()
        f1.show_colorscale(cmap='afmhot')
        x_coord,y_coord=f1.pixel2world(cpix1[0],cpix1[1])
        f1.recenter(x_coord,y_coord,width=width_deg,height=width_deg)
        f1.add_colorbar()
        f1.colorbar.set_font(size='medium')
        f1.colorbar.set_axis_label_text('K')
        f1.colorbar.set_axis_label_font(size='medium')
        f1.axis_labels.hide_x()
        f1.add_label(0.1,0.1,'P',relative=True, color='white', size='18', weight='bold')

        #plt.tight_layout()

        outfile=field+"_"+chunk+"_eb_results.pdf"

        #plt.show()
        fig.savefig("/Users/leclercq/galfacts/aps/final_results/plots/"+outfile,dpi=150, bbox_inches='tight')

    else:
        continue


    
