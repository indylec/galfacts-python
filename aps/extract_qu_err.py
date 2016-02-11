import numpy as np
from astropy.io import fits
import sys

qin=sys.argv[1]
uin=sys.argv[2]

qerr,qhead=fits.getdata(qin,1,header=True)
uerr,uhead=fits.getdata(uin,1,header=True)


qerr_avg=np.mean(qerr,axis=0)
uerr_avg=np.mean(uerr,axis=0)


qout="S1_binned_err_Q.fits"
qavgout="S1_binned_err_avg_Q.fits"

uout="S1_binned_err_U.fits"
uavgout="S1_binned_err_avg_U.fits"

fits.writeto(qout,qerr,qhead)
fits.writeto(uout,uerr,uhead)

fits.writeto(qavgout,qerr_avg,qhead)
fits.writeto(uavgout,uerr_avg,uhead)

