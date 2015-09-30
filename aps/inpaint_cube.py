#Does inpainting on full cubes to remove nans. Run from within the raw_cubes field directory

import numpy as np
from astropy.io import fits
from astropy import wcs
import sys
import inpaint

#read in data

field=sys.argv[1]
stokes=sys.argv[2]
cube_in=sys.argv[3]

cube=fits.getdata(cube_in, ignore_missing_end=True)
cubehead=fits.getheader(cube_in, ignore_missing_end=True)

filled_cube=np.empty(cube.shape)

#loop through planes of cube to inpaint

for i in range(cube.shape[0]):
    print 'Inpainting channel',i
    filled_cube[i,:,:]=inpaint.replace_nans(np.float64(cube[i,:,:]),10,1E-6,2,method='localmean')

#write cube

outname='inpainted/stokes_'+stokes+'/GALFACTS_'+field+'_'+stokes+'_full_inpainted.fits'
cubehead['COMMENT']='All NaNs have been inpainted'

fits.writeto(outname, filled_cube, header=cubehead)



