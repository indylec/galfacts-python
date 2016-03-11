#Extract Q, U, QQ and UU from the numpy savefile and plot (two example chunks, S1 and N2)

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
import sys
from matplotlib import rc

rc('font',family='serif')
rc('text', usetex=True)

####### MAIN ######

field1='S1'
chunk1='1'
field2='N2'
chunk2='3'

#read in the relevant data

infile1='/Users/leclercq/galfacts/aps/final_results/data/data_500/'+field1+'_c'+chunk1+'_full_results_500.npz'
infile2='/Users/leclercq/galfacts/aps/final_results/data/data_500/'+field2+'_c'+chunk2+'_full_results_500.npz'

qin1='/Users/leclercq/galfacts/inpainting/S1_Q_inpainted.fits'
uin1='/Users/leclercq/galfacts/inpainting/S1_U_inpainted.fits'

qin2='/Users/leclercq/galfacts/inpainting/N2_Q_inpainted.fits'
uin2='/Users/leclercq/galfacts/inpainting/N2_U_inpainted.fits'

width_deg=900./60.

data1=np.load(infile1)
cpix1=data1['cpix']
qq1=data1['qq']
uu1=data1['uu']

data2=np.load(infile2)
cpix2=data2['cpix']
qq2=data2['qq']
uu2=data2['uu']

#Plot! 2x4 figure, S1 on the left, N2 on the right

fig=plt.figure(figsize=(9,14))

f1 = aplpy.FITSFigure(qin1, figure=fig, subplot=[0.05,0.7625,0.4,0.1875])
f1.tick_labels.set_font(size='medium')
f1.axis_labels.set_font(size='medium')
#f1.axis_labels.hide()
#f1.tick_labels.hide()
f1.show_colorscale(cmap='afmhot', vmin=-0.15, vmax=0.10)
x_coord,y_coord=f1.pixel2world(cpix1[0],cpix1[1])
f1.recenter(x_coord,y_coord,width=width_deg,height=width_deg)
f1.add_colorbar()
f1.colorbar.set_font(size='medium')
f1.colorbar.set_axis_label_text('K')
f1.axis_labels.set_xtext("RA")
f1.axis_labels.set_ytext("Dec")
#f1.axis_labels.hide_x()
#f1.tick_labels.hide_x()
f1.tick_labels.set_xformat('hh:mm')
f1.tick_labels.set_yformat('dd')
#f1.axis_labels.hide_x()
#f1.add_label(0.1,0.1,'P',relative=True, color='white', size='18', weight='bold')
f1.add_label(0.17,0.08,field1+'-R'+chunk1+' Q',relative=True, color='black', size='14', weight='bold')

f2 = aplpy.FITSFigure(uin1, figure=fig, subplot=[0.55,0.7625,0.4,0.1875])
f2.tick_labels.set_font(size='medium')
f2.axis_labels.set_font(size='medium')
#f2.axis_labels.hide()
#f2.tick_labels.hide()
f2.show_colorscale(cmap='afmhot', vmin=-0.15, vmax=0.10)
x_coord,y_coord=f2.pixel2world(cpix1[0],cpix1[1])
f2.recenter(x_coord,y_coord,width=width_deg,height=width_deg)
f2.add_colorbar()
f2.colorbar.set_font(size='medium')
f2.colorbar.set_axis_label_text('K')
f2.axis_labels.set_xtext("RA")
f2.axis_labels.set_ytext("Dec")
#f2.axis_labels.hide_x()
#f2.tick_labels.hide_x()
f2.tick_labels.set_xformat('hh:mm')
f2.tick_labels.set_yformat('dd')
#f2.axis_labels.hide_x()
#f2.add_label(0.1,0.1,'P',relative=True, color='white', size='18', weight='bold')
f2.add_label(0.17,0.1,field1+'-R'+chunk1+' U',relative=True, color='black', size='14', weight='bold')

f3 = aplpy.FITSFigure(qin2, figure=fig, subplot=[0.05,0.2875,0.4,0.1875])
f3.tick_labels.set_font(size='medium')
f3.axis_labels.set_font(size='medium')
#f3.axis_labels.hide()
#f3.tick_labels.hide()
f3.show_colorscale(cmap='afmhot', vmin=-0.05, vmax=0.10)
x_coord,y_coord=f3.pixel2world(cpix2[0],cpix2[1])
f3.recenter(x_coord,y_coord,width=width_deg,height=width_deg)
f3.add_colorbar()
f3.colorbar.set_font(size='medium')
f3.colorbar.set_axis_label_text('K')
f3.axis_labels.set_xtext("RA")
f3.axis_labels.set_ytext("Dec")
#f3.axis_labels.hide_x()
#f3.tick_labels.hide_x()
f3.tick_labels.set_xformat('hh:mm')
f3.tick_labels.set_yformat('dd')
#f3.axis_labels.hide_x()
#f3.add_label(0.1,0.1,'P',relative=True, color='white', size='18', weight='bold')
f3.add_label(0.17,0.08,field2+'-R'+chunk2+' Q',relative=True, color='black', size='14', weight='bold')

f4 = aplpy.FITSFigure(uin2, figure=fig, subplot=[0.55,0.2875,0.4,0.1875])
f4.tick_labels.set_font(size='medium')
f4.axis_labels.set_font(size='medium')
#f4.axis_labels.hide()
#f4.tick_labels.hide()
f4.show_colorscale(cmap='afmhot', vmin=-0.05, vmax=0.10)
x_coord,y_coord=f4.pixel2world(cpix2[0],cpix2[1])
f4.recenter(x_coord,y_coord,width=width_deg,height=width_deg)
f4.add_colorbar()
f4.colorbar.set_font(size='medium')
f4.colorbar.set_axis_label_text('K')
f4.axis_labels.set_xtext("RA")
f4.axis_labels.set_ytext("Dec")
#f4.axis_labels.hide_x()
#f4.tick_labels.hide_x()
f4.tick_labels.set_xformat('hh:mm')
f4.tick_labels.set_yformat('dd')
#f4.axis_labels.hide_x()
#f4.add_label(0.1,0.1,'P',relative=True, color='white', size='18', weight='bold')
f4.add_label(0.17,0.1,field2+'-R'+chunk2+' U',relative=True, color='black', size='14', weight='bold')


qq_ymin=1E-7
qq_ymax=1E-2

ax2=fig.add_axes([0.05,0.525,0.4,0.1875])
im1=ax2.imshow(np.log10(np.abs(qq1))[256:768,256:768], clim=(np.log10(qq_ymin),np.log10(qq_ymax)))
plt.gca()
locs2=np.asarray([-190, 0, 189])
labs2=[4000, 0, 4000]
#print ell_r[xlocs2], ell_r[ylocs2]
#plt.xticks(locs2+256,ell_r[512,locs2+512].astype(int), size='small')
#plt.yticks(locs2+256,ell_r[locs2+512,512].astype(int),size='small')
plt.xticks(locs2+256,labs2, size='medium')
plt.yticks(locs2+256,labs2,size='medium')
plt.xlabel('$\ell_x$', size=16, weight='bold')
plt.ylabel('$\ell_y$', size=16, rotation='horizontal', weight='bold')
divider1 = make_axes_locatable(ax2)
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar1=plt.colorbar(mappable=im1, cax=cax1)
cbar1.set_label('$\mathrm{log}_{10} \; QQ [\mathrm{K}^2]$',size=12)
cbar1.ax.tick_params(labelsize=12)
txt2=ax2.text(20,480,field1+'-R'+chunk1+' QQ', size=14,color='white', weight='bold')


ax3=fig.add_axes([0.55,0.525,0.4,0.1875])
im1=ax3.imshow(np.log10(np.abs(uu1))[256:768,256:768], clim=(np.log10(qq_ymin),np.log10(qq_ymax)))
plt.gca()
locs3=np.asarray([-190, 0, 189])
labs3=[4000, 0, 4000]
#print ell_r[xlocs3], ell_r[ylocs3]
#plt.xticks(locs3+256,ell_r[512,locs3+512].astype(int), size='small')
#plt.yticks(locs3+256,ell_r[locs3+512,512].astype(int),size='small')
plt.xticks(locs3+256,labs3, size='medium')
plt.yticks(locs3+256,labs3,size='medium')
plt.xlabel('$\ell_x$', size=16, weight='bold')
plt.ylabel('$\ell_y$', size=16, rotation='horizontal', weight='bold')
divider1 = make_axes_locatable(ax3)
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar1=plt.colorbar(mappable=im1, cax=cax1)
cbar1.set_label('$\mathrm{log}_{10} \; QQ [\mathrm{K}^2]$',size=12)
cbar1.ax.tick_params(labelsize=12)
txt3=ax3.text(20,480,field1+'-R'+chunk1+' UU',size=14, color='black', weight='bold')

ax4=fig.add_axes([0.05,0.05,0.4,0.1875])
im1=ax4.imshow(np.log10(np.abs(qq2))[256:768,256:768], clim=(np.log10(qq_ymin),np.log10(qq_ymax)))
plt.gca()
locs4=np.asarray([-190, 0, 189])
labs4=[4000, 0, 4000]
#print ell_r[xlocs4], ell_r[ylocs4]
#plt.xticks(locs4+256,ell_r[512,locs4+512].astype(int), size='small')
#plt.yticks(locs4+256,ell_r[locs4+512,512].astype(int),size='small')
plt.xticks(locs4+256,labs4, size='medium')
plt.yticks(locs4+256,labs4,size='medium')
plt.xlabel('$\ell_x$', size=16, weight='bold')
plt.ylabel('$\ell_y$', size=16, rotation='horizontal', weight='bold')
divider1 = make_axes_locatable(ax4)
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar1=plt.colorbar(mappable=im1, cax=cax1)
cbar1.set_label('$\mathrm{log}_{10} \; QQ [\mathrm{K}^2]$',size=12)
cbar1.ax.tick_params(labelsize=12)
txt4=ax4.text(20,480,field2+'-R'+chunk2+' QQ',size=14,color='white', weight='bold')

ax5=fig.add_axes([0.55,0.05,0.4,0.1875])
im1=ax5.imshow(np.log10(np.abs(uu2))[256:768,256:768], clim=(np.log10(qq_ymin),np.log10(qq_ymax)))
plt.gca()
locs5=np.asarray([-190, 0, 189])
labs5=[4000, 0, 4000]
#print ell_r[xlocs5], ell_r[ylocs5]
#plt.xticks(locs5+256,ell_r[512,locs5+512].astype(int), size='small')
#plt.yticks(locs5+256,ell_r[locs5+512,512].astype(int),size='small')
plt.xticks(locs5+256,labs5, size='medium')
plt.yticks(locs5+256,labs5,size='medium')
plt.xlabel('$\ell_x$', size=16, weight='bold')
plt.ylabel('$\ell_y$', size=16, rotation='horizontal', weight='bold')
divider1 = make_axes_locatable(ax5)
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar1=plt.colorbar(mappable=im1, cax=cax1)
cbar1.set_label('$\mathrm{log}_{10} \; QQ [\mathrm{K}^2]$',size=12)
cbar1.ax.tick_params(labelsize=12)
txt5=ax5.text(20,480,field2+'-R'+chunk2+' UU',size=14,color='black', weight='bold')

plt.show()
fig.savefig("/Users/leclercq/thesis/aps/figures/images/aps_qqnoise_illustration.pdf",dpi=150, bbox_inches='tight')
