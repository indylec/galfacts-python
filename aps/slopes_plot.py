import numpy as np
import healpy as hp
from astropy.io import ascii
import sys
import matplotlib.pyplot as plt
from matplotlib import cm

from matplotlib import rc

#rc('font',family='serif')
#rc('text', usetex=True)

tabfile=sys.argv[1]
#nside=int(sys.argv[2])

tab=ascii.read(tabfile)

glat=np.asarray(tab['b'])
glon=np.asarray(tab['l'])
slopes=np.asarray(tab['slope'])
fieldnos=np.empty(len(slopes))
fieldnos[0:5]=1
fieldnos[5:11]=2
fieldnos[11:17]=3
fieldnos[17:23]=4
fieldnos[23:29]=5
fieldnos[29:35]=6
fieldnos[35:41]=7


fig=plt.figure(figsize=(18,9))
ax1=fig.add_subplot(111)

cax=ax1.scatter(glat,slopes,marker='D',facecolors='none',edgecolors='k',s=50, linewidth='2') #,c=fieldnos
#ylims=ax1.get_ylim()
ax1.axvline(-20.,color='k',linestyle='dotted')
ax1.axvline(20.,color='k',linestyle='dotted')
#cbar = fig.colorbar(cax, ticks=[1, 2, 3, 4, 5, 6, 7])
#cbar.ax.set_yticklabels(['S1', 'S2', 'S3', 'S4', 'N2', 'N3', 'N4'])
ax1.tick_params(labelsize='large')

plt.xlabel('$\\rm b (^{\\circ})$', fontsize=30)
plt.ylabel('$-\\beta$', fontsize=30)

#plt.plot()
fig.savefig("/Users/leclercq/galfacts/aps/plots/slopes_plot_poster.png",dpi=300, bbox_inches='tight')
