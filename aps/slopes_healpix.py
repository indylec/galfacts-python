import numpy as np
import healpy as hp
from astropy.io import ascii
import sys
import matplotlib.pyplot as plt
from matplotlib import cm

tabfile=sys.argv[1]
nside=int(sys.argv[2])

tab=ascii.read(tabfile)

glat=np.asarray(tab['b'])
glon=np.asarray(tab['l'])
slopes=np.asarray(tab['slope'])

map=np.ones(hp.nside2npix(nside))*hp.UNSEEN

cmap=cm.jet
cmap.set_under("white")

for i in range(len(slopes)):
    theta=(90.0-glat[i])*np.pi/180.
    if glon[i] < 0:
        phi= (360.+glon[i])*np.pi/180.
    else: phi = glon[i]*np.pi/180.
    temp_pix=hp.ang2pix(nside,theta,phi)
    print temp_pix
    map[temp_pix]=slopes[i]

hp.mollview(map, cmap=cmap)

plt.show()

