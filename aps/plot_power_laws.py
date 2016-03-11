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
tab=ascii.read(tabfile)

glat=np.asarray(tab['b'])
glon=np.asarray(tab['l'])
slopes=np.asarray(tab['slope'])
offsets=np.asarray(tab['offset'])
fields=tab['Field']

fig=plt.figure(figsize=(18,9))
ax1=fig.add_subplot(111)

## colours=np.empty(len(slopes),dtype='string')
## colours[np.where(fields=='S1')]='b'
## colours[np.where(fields=='S2')]='g'
## colours[np.where(fields=='S3')]='r'
## colours[np.where(fields=='S4')]='c'
## colours[np.where(fields=='N2')]='m'
## colours[np.where(fields=='N3')]='darkorange'
## colours[np.where(fields=='N4')]='greenyellow'

x=np.arange(10000)

for i in range(0,6):
    ax1.plot(x[453:4849],10**(slopes[i]*(np.log10(x[453:4849])-np.log10(x[1482]))+offsets[i]),color='b')
for i in range(6,12):
    ax1.plot(x[453:4849],10**(slopes[i]*(np.log10(x[453:4849])-np.log10(x[1482]))+offsets[i]),color='g')
for i in range(12,18):
    ax1.plot(x[453:4849],10**(slopes[i]*(np.log10(x[453:4849])-np.log10(x[1482]))+offsets[i]),color='r')
for i in range(18,23):
    ax1.plot(x[453:4849],10**(slopes[i]*(np.log10(x[453:4849])-np.log10(x[1482]))+offsets[i]),color='c')
for i in range(23,29):
    ax1.plot(x[453:4849],10**(slopes[i]*(np.log10(x[453:4849])-np.log10(x[1482]))+offsets[i]),color='m')
for i in range(29,35):
    ax1.plot(x[453:4849],10**(slopes[i]*(np.log10(x[453:4849])-np.log10(x[1482]))+offsets[i]),color='darkorange')
for i in range(35,41):
    ax1.plot(x[453:4849],10**(slopes[i]*(np.log10(x[453:4849])-np.log10(x[1482]))+offsets[i]),color='greenyellow')

ax1.set_xlabel('$\ell$',fontsize='medium' )
ax1.set_xlim(300,6000)
ax1.set_xscale('log') 
ax1.set_yscale('log')

import matplotlib.lines as mlines
s1=mlines.Line2D([], [], color='b', label='S1')
s2=mlines.Line2D([], [], color='g', label='S2')
s3=mlines.Line2D([], [], color='r', label='S3')
s4=mlines.Line2D([], [], color='c', label='S4')
n2=mlines.Line2D([], [], color='m', label='N2')
n3=mlines.Line2D([], [], color='darkorange', label='N3')
n4=mlines.Line2D([], [], color='greenyellow', label='N4')

ax1.legend(handles=[s1, s2, s3, s4, n2, n3, n4], loc=0)

plt.show()
