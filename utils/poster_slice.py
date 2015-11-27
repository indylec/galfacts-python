import numpy as np
from astropy.io import fits
import sys

q_in='/local/scratch/GALFACTS/raw_cubes_3.1.2/s1/inpainted/stokes_q/GALFACTS_S1_Q_full_inpainted.fits'

u_in='/local/scratch/GALFACTS/raw_cubes_3.1.2/s1/inpainted/stokes_u/GALFACTS_S1_U_full_inpainted.fits'

dicube_in='/local2/scratch/GALFACTS/rm_synthesis/rmsyn/s1/GALFACTS_S1_dicube_4.fits'

q_full=fits.getdata(q_in)

q_trunc=q_full[:,800:1074,4124:5019]

fits.writeto('S1_Q_poster.fits',q_trunc)

u_full=fits.getdata(u_in)

u_trunc=u_full[:,800:1074,4124:5019]

fits.writeto('S1_U_poster.fits',u_trunc)

di_full=fits.getdata(dicube_in)

di_trunc=di_full[:,800:1074,0:895]

fits.writeto('S1_phi_poster.fits',di_trunc)

