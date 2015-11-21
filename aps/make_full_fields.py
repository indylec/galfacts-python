#input:  individual chunks of rm, angle, chisq
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

chisq_files=glob.glob(field_low+'/*[0-4]_chisq.fits')

rm_files.sort()
angle_files.sort()
chisq_files.sort()

print "Stitching together",angle_files,rm_files,chisq_files

for i in range(len(angle_files)):
    temp_ang=fits.getdata(angle_files[i])
    temp_ang_err=fits.getdata(angle_files[i],1)

    temp_rm=fits.getdata(rm_files[i])
    temp_rm_err=fits.getdata(rm_files[i],1)

    temp_chi=fits.getdata(chisq_files[i])
    

    if i == 0:
        angle_map=temp_ang
        angle_err=temp_ang_err

        rm_map=temp_rm
        rm_err=temp_rm_err

        chi_map=temp_chi
        

    else:
        angle_map=np.hstack((angle_map,temp_ang))
        angle_err=np.hstack((angle_err,temp_ang_err))

        rm_map=np.hstack((rm_map,temp_rm))
        rm_err=np.hstack((rm_err,temp_rm_err))

        chi_map=np.hstack((chi_map,temp_chi))
        

        

angle_map=np.asarray(angle_map % (2.*np.pi))

#print angle_map.shape

#print np.logical_and(0.5*np.pi<angle_map, angle_map<=1.5*np.pi).shape

#print np.asarray(np.where(np.logical_and(0.5*np.pi<angle_map, angle_map<=1.5*np.pi)))

###print np.where(np.logical_and(0.5*np.pi<angle_map, angle_map<=1.5*np.pi))



angle_map[np.where(np.logical_and(0.5*np.pi<angle_map, angle_map<=1.5*np.pi))]-=np.pi
angle_map[np.where(angle_map>1.5*np.pi)]-=2*np.pi

chi_map=chi_map/373.#This is the DOF factor, N=376 (376 channels), DOF = N-1-(no_parameters) = 2 in this case, slope and intercept

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

chi_head=new_header.copy()
chi_head['OBJECT']= 'GALFACTS '+field_up+' CHI^2 map'
chi_head['BUNIT']='rad'

fits.writeto(field_low+'/GALFACTS_'+field_up+'_RM.fits',rm_map,rm_head)
fits.writeto(field_low+'/GALFACTS_'+field_up+'_RMerr.fits',rm_err,rm_err_head)
fits.writeto(field_low+'/GALFACTS_'+field_up+'_angle.fits',angle_map,angle_head)
fits.writeto(field_low+'/GALFACTS_'+field_up+'_angle_err.fits',angle_err,angle_err_head)
fits.writeto(field_low+'/GALFACTS_'+field_up+'_chisq.fits',chi_map,chi_head)




