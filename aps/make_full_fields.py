#input:  individual chunks of rm, angle
#output: rm map, rm_err_map for each field
#        angle map, angle_err_map for each field
#run from the rmsyn directory

import numpy as np
from astropy.io import fits
import sys
import glob

field=sys.argv[1]
avgfile=sys.argv[2]
field_low=field.lower()
field_up=field.upper()

angle_files=glob.glob(field_low+'/*[0-4]_cube_angle*')

rm_files=glob.glob(field_low+'/*[0-4]_cube_rm*')

rm_files.sort()
angle_files.sort()
print "Stitching together",angle_files,rm_files

for i in range(len(angle_files)):
    temp_ang=fits.getdata(angle_files[i])
    temp_ang_err=fits.getdata(angle_files[i],1)

    temp_rm=fits.getdata(rm_files[i])
    temp_rm_err=fits.getdata(rm_files[i],1)

    if i == 0:
        angle_map=temp_ang
        angle_err=temp_ang_err

        rm_map=temp_rm
        rm_err=temp_rm_err

    else:
        angle_map=np.hstack((angle_map,temp_ang))
        angle_err=np.hstack((angle_err,temp_ang_err))

        rm_map=np.hstack((rm_map,temp_rm))
        rm_err=np.hstack((rm_err,temp_rm_err))

ref_header=fits.getheader(avgfile)
new_header=ref_header.copy()

new_header.remove('NAXIS3')
new_header.remove('CTYPE3')
new_header.remove('CRPIX3')
new_header.remove('CRVAL3')
new_header.remove('CDELT3')

rm_head=new_header.copy()
rm_head['OBJECT']= 'GALFACTS '+field_up+' RM map'
rm_head['BUNIT']='rad/m^2'

rm_err_head=new_header.copy()
rm_err_head['OBJECT']= 'GALFACTS '+field_up+' RM error map'
rm_err_head['BUNIT']='rad/m^2'

angle_head=new_header.copy()
angle_head['OBJECT']= 'GALFACTS '+field_up+' ANGLE map'
angle_head['BUNIT']='rad'

angle_err_head=new_header.copy()
angle_err_head['OBJECT']= 'GALFACTS '+field_up+' ANGLE error map'
angle_err_head['BUNIT']='rad'

fits.writeto(field_low+'/GALFACTS_'+field_up+'_RM.fits',rm_map,rm_head)
fits.writeto(field_low+'/GALFACTS_'+field_up+'_RMerr.fits',rm_err,rm_err_head)
fits.writeto(field_low+'/GALFACTS_'+field_up+'_angle.fits',angle_map,angle_head)
fits.writeto(field_low+'/GALFACTS_'+field_up+'_angle_err.fits',angle_err,angle_err_head)



