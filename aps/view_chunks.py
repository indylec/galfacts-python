import numpy as np
from astropy.io import fits
import sys
from numpy import pi

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as mpl

import aplpy

infile=sys.argv[1]
field=sys.argv[3]
width=int(sys.argv[4])

width_deg=width*(1/60.)

qin=fits.open(sys.argv[1])
uin=fits.open(sys.argv[2])

polfile=field+'_polarised_intensity.fits'

header_cube=qin[0].header

map_header=header_cube.copy()
map_header.remove('ctype3')
map_header.remove('crval3')
map_header.remove('crpix3')
map_header.remove('cdelt3')
map_header.remove('crota3')
map_header['OBJECT']='GALFACTS_{0} Polarised intensity map'.format(field)

xw=map_header['NAXIS1']
yw=map_header['NAXIS2']

qdata=qin[0].data[0,:,:]
udata=uin[0].data[0,:,:]

polint=np.sqrt(qdata**2+udata**2)

fits.writeto(polfile,polint,header=map_header)

fig = mpl.figure(figsize=(14,7))

f1 = aplpy.FITSFigure(polfile, figure=fig)

f1.tick_labels.set_font(size='x-small')
f1.axis_labels.set_font(size='small')
f1.show_colorscale(cmap='afmhot')
f1.add_colorbar()
f1.colorbar.set_axis_label_text('brightness [K]')


nochunks=int(xw/width)
crop=(xw-width*nochunks)/2.
clims_pix_x=np.arange(nochunks)*width+(crop)+width/2.
clims_pix_y=np.ones(nochunks)*(yw/2)

clims_world_x, clims_world_y=f1.pixel2world(clims_pix_x,clims_pix_y)

for i in range(nochunks):
    f1.show_rectangles(clims_world_x[i],clims_world_y[i], width_deg, width_deg)

fig_file=field+'pol_withchunks.pdf'

fig.savefig(fig_file, bbox_inches='tight')
