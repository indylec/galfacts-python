import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt


#read in source fluxes
icat=Table.read('/Users/leclercq/galfacts/aps/source_removal/S1_sources/Averaged_I.fits')

flux=icat['S_int']
print np.amin(flux)

#read in NVSS

nvss=Table.read('/Users/leclercq/galfacts/aps/source_removal/nvss_s1_only.fits')
nvss_flux=0.001*nvss['S1.4']

#set up Katgert counts:

xr=np.arange(-1,4,0.01)
lx=10**(xr-3)

def kat_poly(x):
    a0=0.878
    a1=0.614
    a2=0.191
    a3=-0.07
    poly3=a3*x**3+a2*x**2+a1*x+a0
    return poly3

#make bins for galfacts sources and figure out fraction of sky (hard coded: use utils/surface_area.py to find this for various fields)

sterad_area=0.461937072074

#make bins where delta log (S)=0.3
logbins=np.arange(np.log10(np.amin(flux)),np.log10(np.amax(flux))+0.3,0.3)
logbins10=10**logbins

#bin the sources by flux
counts,bin_edges=np.histogram(flux,bins=logbins10)

#work out bin parameters for plotting (center, width)
width=bin_edges[1:]-bin_edges[:-1]
bin_centers=(bin_edges[:-1]+bin_edges[1:])/2.
#print bin_centers
normcounts=counts/width/sterad_area

#do the same for NVSS
logbins_nvss=np.arange(np.log10(np.amin(nvss_flux)),np.log10(np.amax(nvss_flux))+0.3,0.3)
logbins10_nvss=10**logbins_nvss

#bin the sources by flux
counts_nvss,bin_edges_nvss=np.histogram(nvss_flux,bins=logbins10_nvss)

#work out bin parameters for plotting (center, width)
width_nvss=bin_edges_nvss[1:]-bin_edges_nvss[:-1]
bin_centers_nvss=(bin_edges_nvss[:-1]+bin_edges_nvss[1:])/2.
#print bin_centers_nvss
normcounts_nvss=counts_nvss/width_nvss/sterad_area

#plot
plt.plot(lx,10**kat_poly(xr),bin_centers,normcounts*bin_centers**(5/2),bin_centers_nvss,normcounts_nvss*bin_centers_nvss**(5/2))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('S [Jy]')
plt.ylabel(r'N(S) S$^{5/2}$ [Jy$^{-1}$ sr$^{-1}$]')
plt.legend(['Katgert counts','GALFACTS S1', 'NVSS sources for S1'],loc=2)

#plt.show()
plt.savefig('/Users/leclercq/galfacts/aps/source_removal/source_counts.pdf',bbox_inches='tight')

