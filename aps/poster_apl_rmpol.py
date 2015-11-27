#/Users/leclercq/miniconda/bin/python
#run from ~/galfacts

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import aplpy
from astropy.io import fits
import sys

#field=sys.argv[1]
center_pix_x=4362
center_pix_y=937
width_x=1315
width_y=274

width_x_deg=width_x/60.
width_y_deg=width_y/60.

dpi = 300
figsize_inch = 14.26, 6.22
fig = plt.figure(figsize=figsize_inch, dpi=dpi)

x1,y1=4670-4124,1040-800
x2,y2=4754-4124,932-800



rmfile='aps/rm/GALFACTS_S1_RM.fits'
polfile='inpainting/S1_P_inpainted.fits'
#chisqfile='GALFACTS_'+field+'_chisq.fits'

try:
    fits.delval(rmfile,'CROTA3')    
except KeyError:
    print "The RM fits header has already been fixed"

## try:
##     fits.delval(chisqfile,'CROTA3')
## except KeyError:
##     print "The chisq fits header has already been fixed"

try:
    fits.delval(polfile,'CROTA3')
except KeyError:
    print "The pol fits header has already been fixed"

f1 = aplpy.FITSFigure(rmfile, figure=fig, subplot=[0.1,0.57,0.9,0.38])
f1.tick_labels.set_xformat('dd')
f1.tick_labels.set_yformat('dd')
f1.tick_labels.set_font(size='large')
f1.axis_labels.set_font(size='x-large')
#f1.axis_labels.set_xpad(20)
#f1.axis_labels.set_ypad(15)
#f1.axis_labels.set_xtext("RA")
f1.axis_labels.set_ytext("Dec")
f1.axis_labels.hide_x()
f1.tick_labels.hide_x()
f1.tick_labels.set_xformat('dd')
f1.tick_labels.set_yformat('dd')
f1.show_colorscale(cmap='seismic',vmin=-130,vmax=130)
x_coord,y_coord=f1.pixel2world(center_pix_x,center_pix_y)
x1w,y1w=f1.pixel2world(x1,y1)
x2w,y2w=f1.pixel2world(x2,y2)

f1.recenter(x_coord,y_coord,width=width_x_deg,height=width_y_deg)
f1.add_colorbar()
f1.colorbar.set_box([0.95,0.58,0.02,0.36])
#f1.add_grid()
#f1.grid.set_color('black')
f1.colorbar.set_axis_label_text('$\mathrm{rad}/\mathrm{m}^2$')
f1.colorbar.set_axis_label_font(size=25)
f1.show_circles(x1w,y1w,10./60.,linewidth=4,facecolor='none',edgecolor='black')
f1.show_circles(x2w,y2w,10./60.,linewidth=4,facecolor='none',edgecolor='black')
#f1.show_rectangles(0,32,1000,1000)
#f1.show_rectangles(1000,32,1000,1000)
#f1.show_rectangles(2000,32,1000,1000)
#f1.show_rectangles(3000,32,1000,1000)
#f1.show_rectangles(4000,32,1000,1000)

f2 = aplpy.FITSFigure(polfile, figure=fig,subplot=[0.1,0.15,0.9,0.38])
f2.tick_labels.set_xformat('dd')
f2.tick_labels.set_yformat('dd')
f2.tick_labels.set_font(size='large')
f2.axis_labels.set_font(size='x-large')
#f1.axis_labels.set_xpad(20)
#f1.axis_labels.set_ypad(15)
f2.axis_labels.set_xtext("RA")
f2.axis_labels.set_ytext("Dec")
#f2.axis_labels.hide_x()
#f2.tick_labels.hide_x()
f2.tick_labels.set_xformat('dd')
f2.tick_labels.set_yformat('dd')
f2.show_colorscale(cmap='afmhot',vmax=0.25)
x_coord,y_coord=f1.pixel2world(center_pix_x,center_pix_y)
f2.recenter(x_coord,y_coord,width=width_x_deg,height=width_y_deg)
f2.add_colorbar()
f2.colorbar.set_box([0.95,0.16,0.02,0.36])
#f1.add_grid()
#f1.grid.set_color('black')
f2.colorbar.set_axis_label_text('$\mathrm{K}$')
f2.colorbar.set_axis_label_font(size=25)
f2.colorbar.set_axis_label_pad(20)

plt.savefig('poster_rm_pol_circ.png',dpi=dpi, bbox_inches='tight')

