import numpy as np
from astropy.io import fits
import sys
from numpy import pi

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as mpl

import aplpy

infile=sys.argv[1]
field=sys.argv[2]
width=sys.argv[3]

width_deg=width*(1/60.)

polin=fits.open(sys.argv[1])
header_cube=polin[0].header

xw=header_cube['NAXIS1']
yw=header_cube['NAXIS2']

fig = mpl.figure(figsize=(14,7))

f1 = aplpy.FITSFigure(infile, figure=fig)

f1.tick_labels.set_font(size='x-small')
f1.axis_labels.set_font(size='small')
f1.show_colorscale(cmap='afmhot')
f1.add_colorbar()
f1.colorbar.set_axis_label_text('brightness [K]')

nochunks=int(xw/width)
crop=(xw-width*nochunks)/2.
clims_pix_x=np.arange(nochunks)*width+(crop)+width/2.
clims_pix_y=np.ones(nochunks)*(yw/2)

for i in range(nochunks):
    f1.show_rectangles(clims_world_x[i],clims_world_y[i], width_deg, width_deg)

fig_file=field+'pol_withchunks.pdf'

fig.savefig(fig_file, bbox_inches='tight')
