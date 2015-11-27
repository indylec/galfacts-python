#!/Users/leclercq/miniconda/bin/python

import numpy as np
from astropy.io import fits
import sys
from numpy import pi
import glob

field=sys.argv[1]

qglob=field+'*_Q_*.fits'
uglob=field+'*_U_*.fits'

qfile=glob.glob(qglob)
#print qfile[0]
#sys.exit('')
ufile=glob.glob(uglob)

qin=fits.open(qfile[0])
uin=fits.open(ufile[0])

polfile='/Users/leclercq/galfacts/inpainting/'+field+'_P_inpainted.fits'

header_cube=qin[0].header

map_header=header_cube.copy()
map_header.remove('ctype3')
map_header.remove('crval3')
map_header.remove('crpix3')
map_header.remove('cdelt3')
map_header.remove('crota3')
map_header['OBJECT']='GALFACTS_{0} Polarised intensity map'.format(field)

if len(qin[0].data.shape) == 3:
    qdata=qin[0].data[0,:,:]
else:
    qdata=qin[0].data
if len(uin[0].data.shape) == 3:
    udata=uin[0].data[0,:,:]
else:
    udata=uin[0].data



polint=np.sqrt(qdata**2+udata**2)

fits.writeto(polfile,polint,header=map_header)



