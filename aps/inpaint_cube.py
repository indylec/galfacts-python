import numpy as np
from astropy.io import fits
from astropy import wcs
import sys
import inpaint

#read in data

field=sys.argv[1]
cube_in=sys.argv[2]

cube=fits.getdata(cube_in)
cubehead=fits.getheader(cube_in)

filled_cube=np.empty(cube.shape)

#loop through planes of cube to inpaint

for i in range cube.shape[0]:
    print 'Inpainting channel',i
    filled_cube[i,:,:]=inpaint.replace_nans(cube[i,:,:],10,1E-6,2,method='localmean')

#write cube

outname='GALFACTS_'+field+'_full_inpainted.fits'
cubehead['COMMENT']='All NaNs have been inpainted'

fits.writeto(outname, filled_cube, header=cubehead)



