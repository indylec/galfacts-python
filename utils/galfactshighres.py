import matplotlib
matplotlib.use("agg")
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

infile="/Users/leclercq/galfacts/3.1.2/fits_files/avg/galfacts_I_healpix.fits"

m=hp.read_map(infile)
#out_filename = input_filename.split(".")[2]


#planck_cmap = ListedColormap(np.loadtxt("/Users/leclercq/repos/galfacts-python/utils/Planck_Parchment_RGB.txt")/255.)
#planck_cmap.set_bad("gray") # color of missing pixels
#planck_cmap.set_under("white") # color of background, necessary if you want to use this colormap directly with hp.mollview(m, cmap=planck_cmap)

#cmap = planck_cmap
#out_filename += "_planck_cmap"
from matplotlib import cm
myafm=cm.afmhot
myafm.set_bad("black")
myafm.set_under("white")

cmap=myafm

out_filename="/Users/leclercq/galfacts/3.1.2/fits_files/avg/galfacts_I_hires_afm_cmap_scale1"

dpi = 300
figsize_inch = 30, 20
fig = plt.figure(figsize=figsize_inch, dpi=dpi)

print "Mollview"
# removed the colorbar, the map range is -500 / +500 microK
hp.mollview(m, fig=fig.number, xsize=figsize_inch[0]*dpi,coord=['C','G'], title="",min=3.0, max=7.0, cbar=False, cmap=cmap,notext=True)
print "Save"
plt.savefig(out_filename + ".png", dpi=dpi, bbox_inches="tight")
