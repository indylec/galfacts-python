#import all the stuff
import numpy as np
from astropy.io import fits
import sys

qfile=sys.argv[1]
outdir=sys.argv[2]
field=sys.argv[3]
stokes=sys.argv[4]
weightfile=sys.argv[5]


nochunks=5
binsize=10

#read in Q
q=fits.getdata(qfile)
qhead=fits.getheader(qfile)

weights=np.loadtxt(weightfile)
bad_chans=np.where(weights==0)

for bad in bad_chans:
    q[bad,:,:]=0.


#loop over chunks

chunksize=q.shape[2]/nochunks
chunk_remain=q.shape[2]%nochunks

nobins=q.shape[0]/binsize
bin_remain=q.shape[0]%binsize


print "Working on field",field

for i in range(nochunks):
    print 'Working on chunk',i
    if chunk_remain==0:
        temp_q=q[:,:,i*chunksize:(i+1)*chunksize]
        #u=u[:,:,i*chunksize:(i+1)*chunksize]
    elif  i <=chunk_remain:
        temp_q=q[:,:,i*(chunksize+1):(i+1)*(chunksize+1)]
        #u=u[:,:,i*(chunksize+1):(i+1)*(chunksize+1)]

    elif i>chunk_remain:
        temp_q=q[:,:,i*chunksize+chunk_remain:(i+1)*chunksize+chunk_remain]
        #u=u[:,:,i*chunksize+chunk_remain:(i+1)*chunksize+chunk_remain]

#loop over channel bins

    temp_q_err=np.empty(np.shape(temp_q))

    for bin in range (nobins):
        #print 'working on bin',bin
        if bin_remain == 0:
            temp_q_err[bin*binsize:(bin+1)*binsize,:,:]=np.nanstd(temp_q[bin*binsize:(bin+1)*binsize,:,:],axis=0,keepdims=True)

        elif bin<=bin_remain:
            temp_q_err[bin*(binsize+1):(bin+1)*(binsize+1),:,:]=np.nanstd(temp_q[bin*(binsize+1):(bin+1)*(binsize+1),:,:],axis=0,keepdims=True)

        elif bin>bin_remain:
            temp_q_err[bin*binsize+bin_remain:(bin+1)*binsize+bin_remain,:,:]=np.nanstd(temp_q[bin*binsize+bin_remain:(bin+1)*binsize+bin_remain,:,:],axis=0,keepdims=True)
            

    chunkrange=np.empty(4)
    chunkrange[0]=0
    chunkrange[1]=temp_q.shape[2]-1
    chunkrange[2]=0
    chunkrange[3]=temp_q.shape[1]-1

    new_header=qhead.copy()

    cpix_ra = (chunkrange[1] - chunkrange[0]) / 2. + 1.
    cpix_dec = (chunkrange[3] - chunkrange[2]) / 2. + 1.

    cpix_ra_old = qhead['CRPIX1'] - chunkrange[0]
    cpix_dec_old = qhead['CRPIX2'] - chunkrange[2]

    crval_ra = qhead['CRVAL1'] + (cpix_ra - cpix_ra_old) * qhead['CDELT1']
    crval_dec = qhead['CRVAL2'] + (cpix_dec - cpix_dec_old) * qhead['CDELT2']

    new_header['CRPIX1'] = cpix_ra
    new_header['CRVAL1'] = crval_ra
    new_header['NAXIS1'] = chunkrange[1]-chunkrange[0]

    new_header['CRPIX2'] = cpix_dec
    new_header['CRVAL2'] = crval_dec
    new_header['NAXIS2'] = chunkrange[3]-chunkrange[2]

    temp_q=temp_q.astype(np.float32)
    temp_q_err=temp_q_err.astype(np.float32)

    qhdu=fits.PrimaryHDU(temp_q,new_header)
    qerrhdu=fits.ImageHDU(temp_q_err)

    hdulist=fits.HDUList([qhdu,qerrhdu])

    hdulist.writeto(outdir+'GALFACTS_'+field+'_'+stokes+'_chunk_'+str(i)+'.fits')
    



#write out Q chunk and Q stdev as fits file and image extension

#same with U

