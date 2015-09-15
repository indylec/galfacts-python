import numpy as np
from astropy.io import fits
import sys

field=sys.argv[1]
q_in=sys.argv[2]
u_in=sys.argv[3]

#read in q cube & u cube, sum along frequency axis to find nans

qcube=fits.getdata(q_in)

qcube_sum=np.sum(qcube,axis=0)

qcube=0.0

ucube=fits.getdata(u_in)

ucube_sum=np.sum(ucube,axis=0)

#create union of nan boolean arrays (True means not nan) and convert to float

nan_bool=np.logical_and(~np.isnan(qcube_sum),~np.isnan(ucube_sum))

mask=nan_bool.astype(float)

#output mask as fits file

outfile=field+"_pyrmsynth_mask.fits"

fits.writeto(outfile,mask)
