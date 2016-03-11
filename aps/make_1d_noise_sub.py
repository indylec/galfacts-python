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

    #if i.endswith(".npz"):
    if "N2_c3" in i:

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

        fig=plt.figure(figsize=(10,8))

        ax=fig.add_subplot(211)
        #ax.set_autoscale_on(False)

        ax.set_xlabel('$\ell$',fontsize='medium' )

        eeraw_lin,= ax.plot(final_bins,ee1draw,'r--',alpha=0.8, linewidth=2)
        eecorr_lin,=ax.plot(final_bins,ee1dcorr,'r-',alpha=0.8, linewidth=2)
        eenoise_lin,= ax.plot(final_bins,ee1dnoise,'r:',alpha=0.8, linewidth=2)

        beam_cut =ax.axvline(x=180/(3.5/60.),color='k',linestyle='dashed',alpha=0.8)

        ymin=1E-9
        ymax=2.
        ax.set_ylim(ymin,ymax)
        ax.set_xlim(10,30000)

        ax.set_ylabel('$C_{\ell}[K^2]$',fontsize='medium')
        ax.tick_params(labelsize='small')
        
        ax.legend([eeraw_lin, eecorr_lin, eenoise_lin, beam_cut],["raw EE","corrected EE","EE noise", "beamwidth scale"],fontsize='medium',loc=0)
        ax.set_xscale('log') 
        ax.set_yscale('log')
        #ax.text(9500, 1E-4, r'$\beta$ = {:.2f} '.format(-slope))


        ax1=fig.add_subplot(212)
        #ax1.set_autoscale_on(False)

        ax1.set_xlabel('$\ell$',fontsize='medium' )

        bbraw_lin,= ax1.plot(final_bins,bb1draw,'b--',alpha=0.4, linewidth=2)
        bbcorr_lin,=ax1.plot(final_bins,bb1dcorr,'b-',alpha=0.4, linewidth=2)
        bbnoise_lin,= ax1.plot(final_bins,bb1dnoise,'b:',alpha=0.8, linewidth=2)
        beam_cut =ax1.axvline(x=180/(3.5/60.),color='k',linestyle='dashed',alpha=0.8)

        ymin=1E-9
        ymax=2.
        ax1.set_ylim(ymin,ymax)
        ax1.set_xlim(10,30000)

        ax1.set_ylabel('$C_{\ell}[K^2]$',fontsize='medium')
        ax1.tick_params(labelsize='small')
        
        ax1.legend([bbraw_lin, bbcorr_lin, bbnoise_lin, beam_cut],["raw BB","corrected BB", "BB noise", "beamwidth scale"],fontsize='medium',loc=0)
        ax1.set_xscale('log') 
        ax1.set_yscale('log')
        #ax1.text(9500, 1E-4, r'$\beta$ = {:.2f} '.format(-slope))

        plt.tight_layout()

        outfile=field+"_"+chunk+"_1d_noise_sub.pdf"

        #plt.show()
        fig.savefig("/Users/leclercq/galfacts/aps/final_results/plots/"+outfile,dpi=150, bbox_inches='tight')

    else:
        continue


    
