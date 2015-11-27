#/Users/leclercq/miniconda/bin/python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import aplpy
from astropy.io import fits
import sys

field=sys.argv[1]

dpi = 300
figsize_inch = 20, 13
fig = plt.figure(figsize=figsize_inch, dpi=dpi)

rmfile='GALFACTS_'+field+'_RM.fits'
polfile=field+'_polarised_intensity.fits'
chisqfile='GALFACTS_'+field+'_chisq.fits'

try:
    fits.delval(rmfile,'CROTA3')    
except KeyError:
    print "The RM fits header has already been fixed"

try:
    fits.delval(chisqfile,'CROTA3')
except KeyError:
    print "The chisq fits header has already been fixed"

try:
    fits.delval(polfile,'CROTA3')
except KeyError:
    print "The pol fits header has already been fixed"

f1 = aplpy.FITSFigure(rmfile, figure=fig, subplot=(3,1,2))

f1.tick_labels.set_font(size='x-small')
f1.axis_labels.set_font(size='small')
f1.axis_labels.hide()
f1.ticks.hide_x()
f1.tick_labels.hide_x()
f1.tick_labels.set_xformat('dd')
f1.tick_labels.set_yformat('dd')
f1.show_colorscale(cmap='seismic',vmin=-150,vmax=150)
f1.add_colorbar()
f1.add_grid()
f1.grid.set_color('black')
#f1.colorbar.set_axis_label_text('brightness [K]')
#f1.show_rectangles(0,32,1000,1000)
#f1.show_rectangles(1000,32,1000,1000)
#f1.show_rectangles(2000,32,1000,1000)
#f1.show_rectangles(3000,32,1000,1000)
#f1.show_rectangles(4000,32,1000,1000)

f2 = aplpy.FITSFigure(chisqfile, figure=fig,subplot=(3,1,3))

f2.tick_labels.set_font(size='x-small')
f2.axis_labels.set_font(size='small')
f2.show_colorscale(cmap='Blues_r',vmin=1.0,vmax=10.0)
f2.tick_labels.set_xformat('dd')
f2.tick_labels.set_yformat('dd')
f2.add_colorbar()
f2.add_grid()
f2.grid.set_color('black')
#f2.colorbar.set_axis_label_text('RM [rad/m^2]')

f3 = aplpy.FITSFigure(polfile, figure=fig,subplot=(3,1,1))

f3.tick_labels.set_font(size='x-small')
f3.axis_labels.set_font(size='small')
f3.axis_labels.hide()
f3.ticks.hide_x()
f3.tick_labels.hide_x()
f3.tick_labels.set_xformat('dd')
f3.tick_labels.set_yformat('dd')
f3.show_colorscale(cmap='afmhot',vmin=0.0,vmax=0.35)
f3.add_colorbar()
f3.add_grid()
f3.grid.set_color('white')
#f3.colorbar.set_axis_label_text('RM [rad/m^2]')

## f1.axis_labels.hide_x()
## f1.tick_labels.hide_x()

## f2.axis_labels.hide_x()
## f2.tick_labels.hide_x()

plt.savefig(field+'_rm_chisq.png',dpi=dpi, bbox_inches='tight')

