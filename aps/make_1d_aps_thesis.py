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

for i in os.listdir(os.getcwd()):

    if i.endswith(".npz"):

        print "Working on "+ i

        field=i.split("_")[0]
        chunk=i.split("_")[1]

        data=np.load(i)

        ee1draw=data['ee1draw']
        ee1dnoise=data['ee1dnoise']
        ee1dcorr=data['ee1dcorr']
        bb1draw=data['bb1draw']
        bb1dnoise=data['bb1dnoise']
        bb1dcorr=data['bb1dcorr']
        final_bins=data['bins']
        

        #print final_bins
        #exit()


        eb1d=(ee1dcorr+bb1dcorr)/2.
        slope,offset,c,d,e=stats.linregress(np.log10(final_bins[32:60])-np.log10(final_bins[45]),np.log10(eb1d[32:60]))

        fig=plt.figure(figsize=(11,6))

        ax=fig.add_subplot(111)
        #ax.set_autoscale_on(False)

        ax.set_xlabel('$\ell$',fontsize='medium' )

        eeraw_lin,= ax.plot(final_bins,ee1draw,'r-',alpha=0.4)
        #eeraw_mark,= ax.plot(final_bins,ee1draw,'ro',markersize=3)
        #bbraw_lin,=ax.plot(final_bins,bb1draw,'b-',alpha=0.4)
        #bbraw_mark,=ax.plot(final_bins,bb1draw,'bo',markersize=3)

        eecorr_lin,=ax.plot(final_bins,ee1dcorr,'r-',alpha=0.4)
        #eecorr_mark,=ax.plot(final_bins,ee1dcorr,'r^',markersize=3)
        #bbcorr_lin,=ax.plot(final_bins,bb1dcorr,'b-',alpha=0.4)
        #bbcorr_mark,=ax.plot(final_bins,bb1dcorr,'b^',markersize=3)

        #eb1d_mark,=ax.plot(final_bins[:64],eb1d[:64],'b^',markersize=5)

        #ee_noise_lin,=ax.plot(final_bins,ee1dnoise,'g:', alpha=0.6)
        #bb_noise_lin,=ax.plot(final_bins,bb1dnoise,'g--', alpha=0.6)

        #power_law,= ax.plot(final_bins[32:60],10**(slope*(np.log10(final_bins[32:60])-np.log10(final_bins[45]))+offset),'k-',linewidth=0.8)

        beam_cut =ax.axvline(x=180/(3.5/60.),color='k',linestyle='dashed',alpha=0.8)

        ymin=1E-9
        ymax=2.
        ax.set_ylim(ymin,ymax)
        ax.set_xlim(10,10000)

        ax.set_ylabel('$C_{\ell}[K^2]$',fontsize='medium')
        ax.tick_params(labelsize='small')
        #ax.legend([(eeraw_mark,eeraw_lin),(bbraw_mark,bbraw_lin),(eecorr_mark,eecorr_lin),(bbcorr_mark,bbcorr_lin),ee_noise_lin,bb_noise_lin, power_law, beam_cut],["EE w/noise","BB w/noise","EE","BB","EE noise","BB noise","EB power-law fit","beamwidth scale"],fontsize='medium',loc=0)
        #ax.legend([(eecorr_mark,eecorr_lin),(bbcorr_mark,bbcorr_lin), power_law, beam_cut],["EE","BB","EB power-law fit","beamwidth scale"],fontsize='medium',loc=0)
        #ax.legend([eb1d_mark, power_law, beam_cut],["EB","EB power-law fit","beamwidth scale"],fontsize='medium',loc=0)
        ax.legend([eeraw_lin, eecorr_lin, beam_cut],["raw EE","corrected EE","beamwidth scale"],fontsize='medium',loc=0)
        ax.set_xscale('log') 
        ax.set_yscale('log')
        ax.text(9500, 1E-4, r'$\beta$ = {:.2f} '.format(-slope))

        plt.tight_layout()

        outfile=field+"_"+chunk+"_final_aps_ee_plot.pdf"

        #plt.show()
        fig.savefig("/Users/leclercq/galfacts/aps/final_results/plots/"+outfile,dpi=150, bbox_inches='tight')

    else:
        continue


    
