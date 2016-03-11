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


for i in os.listdir(os.getcwd()):

    #if i.endswith(".npz"):
    if field+"_c"+chunk in i:

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

        final_bins=final_bins[:64]
        ee1dcorr=ee1dcorr[:64]
        bb1dcorr=bb1dcorr[:64]
        

        print final_bins
        #exit()


        eb1d=(ee1dcorr+bb1dcorr)/2.
        slope,offset,c,d,e=stats.linregress(np.log10(final_bins[32:60])-np.log10(final_bins[45]),np.log10(eb1d[32:60]))

        fig=plt.figure(figsize=(10,8))

        ax=fig.add_subplot(211)
        #ax.set_autoscale_on(False)

        ax.set_xlabel('$\ell$',fontsize='medium' )

        bbcorr_lin,= ax.plot(final_bins,bb1dcorr,'b-',alpha=0.7, linewidth=1.5)
        eecorr_lin,=ax.plot(final_bins,ee1dcorr,'r-',alpha=0.7, linewidth=1.5)
        eb1d_lin,= ax.plot(final_bins,eb1d,color='purple',alpha=0.9, linewidth=1.5)

        beam_cut =ax.axvline(x=180/(3.5/60.),color='k',linestyle='dashed',alpha=0.8)

        ymin=1E-7
        ymax=2.
        ax.set_ylim(ymin,ymax)
        ax.set_xlim(10,8000)

        ax.set_ylabel('$C_{\ell}[K^2]$',fontsize='medium')
        ax.tick_params(labelsize='small')
        
        ax.legend([eecorr_lin, bbcorr_lin, eb1d_lin, beam_cut],["E-mode","B-mode","(E+B)/2", "beamwidth scale"],fontsize='medium',loc=0)
        ax.set_xscale('log') 
        ax.set_yscale('log')
        #ax.text(9500, 1E-4, r'$\beta$ = {:.2f} '.format(-slope))

        plt.tight_layout()

        outfile=field+"_"+chunk+"_eb_example.pdf"

        #plt.show()
        fig.savefig("/Users/leclercq/galfacts/aps/final_results/plots/"+outfile,dpi=150, bbox_inches='tight')

    else:
        continue


    
