import numpy as np
import sys
import os
import datetime
import argparse
from astropy.io import ascii
from astropy.io import fits
import rmsyn_dicube as di
import bottleneck as bn

class Params:
    """
    Class to store all the parameters necessary to compute the dirty cube & FPSF
    """

    def __init__(self):
        """
        Params()

        Initialize parameter class with defaults.
    
        """
        self.qfile = ''
        self.ufile = ''
        self.outfile = ''
        self.outdir = ''
        self.l2 = None
        self.l20 = 0.
        self.weights = None
        self.nphi = 100
        self.dphi = 10
        self.phi_min = -500
        temp = np.arange(self.nphi, dtype=int)
        self.phi = self.phi_min + temp * self.dphi
        self.phi_max=0.
        self.phi_scale=0.
        self.fwhm=0.
        self.field = ''
        self.range = None
        self.nohead = False
        

def get_l2_params(params):
    """
    get_l2(params)

    Gets the frequency axis from FITS cube header and returns the corresponding
    lambda-squared axis

    Inputs:
        params

    Outputs: none but saves the following in params:
        
        l2 - numpy array containing the frequencies converted to wavelength-squared

        l20 - the weighted mean of the l2 array

        phi_max - theoretical maximum Faraday depth which can be probed

        phi_scale - largest scale which can be measured (for extended components)

        
        
    """
    print "Computing l2, l20 and weights..."
    
    c2=299792458.**2

    header=fits.getheader(params.qfile)
    
    dnu=header['CDELT3']
    nchan=header['NAXIS3']
    nuref=header['CRVAL3']

    dl2 = c2/(dnu**2)

    nu = np.arange(nchan)*dnu+nuref

    l2 = 0.5 * c2 * ((nu - 0.5 * dnu) ** -2 + (nu + 0.5 * dnu) ** -2)
    params.l2 = np.flipud(l2)
        
    if params.weights != None:
        if nchan != params.weights.shape[0]:
            raise Exception ('Weights and freq.axis have different sizes')
        else:
            params.l20 = np.sum(params.weights*params.l2)/np.sum(params.weights)
    else:
        params.weights=np.ones(nchan)
        params.l20 = np.sum(params.l2)/nchan

    #Replace nans with 0.0

    #print "Checking for bad channels..."
    #nonzeroq=np.asarray(np.nonzero(np.isnan(qcube)))
    #bad_q=np.unique(nonzeroq[0])
    #print "...Q done..."

    #nonzerou=np.asarray(np.nonzero(np.isnan(ucube)))
    #bad_u=np.unique(nonzerou[0])
    #print "...U done."

    
    

    #compute resolution parameters

    params.phi_max = np.sqrt(3.)/dl2
    
    params.phi_scale = np.pi/np.amin(params.l2)

    print "phi_max is ",params.phi_max
    print "phi_scale is ",params.phi_scale

    print "...done."

            
def compute_fpsf(params):
    """
    compute_fpsf(params)

    Computes the Faraday Point Spread Function for a given lambda-squared distribution
    and phi axis. Also returns resolution parameters specifying the width, FWHM and
    maximum phi of the FPSF. Weighting of the frequency axis can optionally be
    specified.

    input:
        params

    output:
        fpsf - complex numpy array containing the FPSF, binned according to the
        provided phi distribution

        fwhm_fpsf - FWHM of the FPSF, in rad/m^2, providing a measure of the beam
        resolution in Faraday depth space. Computed from B&dB 2005.
        
    """
    fpsf=np.empty(params.phi.shape[0],dtype=complex)
    
    
    
    #if params.weights.shape[0] != params.l2.shape[0]:
    #        raise Exception ('weights array must have same length as lambda-squared\
    #            array')
    print "Making FPSF..."   

    for i in range(params.phi.shape[0]):
        fpsf[i] = np.sum(params.weights*np.exp(-2.*1j*params.phi[i]*(params.l2-params.l20)))/np.sum(params.weights)

    params.fwhm = 2.0*np.sqrt(3)/(np.amax(params.l2)-np.amin(params.l2))

    print "...done."

    return fpsf
        
        

## def make_di_cube(qcube, ucube, params):
##     """
##     make_di_cube(qcube, ucube, phi, l2, l20)

##     Creates a 'dirty' RM cube using slow discrete Fourier transform of the complex
##     polarisation.

##     input:
##         qcube - 3D numpy array (RA, dec, frequency) of Stokes Q intensity. Must have
##         same dimensions as Stokes U cube.

##         ucube - 3D numpy array (RA, dec, frequency) of Stokes U intensity. Must have
##         same dimensions as Stokes Q cube.

##         phi - numpy array of faraday depth bins

##         l2 - numpy array containing wavelength squared axis of data cubes

##         l20 - weighted average of lambda-squared

##         weights - (optional) numpy array containing the weights of each lambda-squared
##         value. The weights mus be between 0 and 1. If this is omitted, all values are
##         assumed to have weight 1.

##     output:
##         rm_di_cube - 3D complex numpy array (RA, dec, phi) containing complex
##         polarisation as a function of faraday depth.
        
##     """
##     #if params.weights == None:
##     #    params.weights = np.ones(params.l2.shape[0])
##     #else:
##     #    if params.weights.shape[0] != params.l2.shape[0]:
##     #        raise Exception ('weights array must have same length as lambda-squared\
##     #            array')
    
    
##     if qcube.shape != ucube.shape:
##         raise Exception('qcube and ucube have different dimensions! Not allowed.')

##     s=qcube.shape

##     rm_di_cube=np.empty((params.phi.shape[0],s[1],s[2]),dtype=complex)

##     temp_los=np.empty((s[0]), dtype=complex)

##     print "Making RM cube..."
##     for i in range(s[2]):
##         for j in range(s[1]):
##             print 'Synthesising pixel {0},{1}'.format(i,j)
##             temp_los=qcube[:,j,i]+1j*ucube[:,j,i]
##             #check for nans in the los
##             nans=np.where(np.isnan(temp_los))
##             #set these to zero in the los and set corresponding weights to 0
##             temp_los[nans]=0.0
##             temp_weights=params.weights
##             temp_weights[nans]=0.0
            
##             #############debugging stuff#################
##             #print temp_los
##             #print "argument",-2.*1j*params.phi[0]*(params.l2-params.l20)
##             #print "exponential",np.exp(-2.*1j*params.phi[0]*(params.l2-params.l20))
##             #print "weight sum",np.sum(temp_weights)
            
##             #print "exp times weights and los",temp_los*temp_weights*np.exp(-2.*1j*params.phi[0]*(params.l2-params.l20))
            
##             #print np.sum(temp_los*temp_weights*np.exp(-2.*1j*params.phi[0]*(params.l2-params.l20)))/np.sum(temp_weights)
##             #############################################
##             for p in range(params.phi.shape[0]):
##                 rm_di_cube[p,j,i] = np.sum(temp_los*temp_weights*\
##                                         np.exp(-2.*1j*params.phi[p]*(params.l2-params.l20)))/np.sum(temp_weights)
##             #print rm_di_cube[0,j,i],rm_di_cube[50,j,i]
##     print "...done."    
    
##     return rm_di_cube


def open_and_trim(params):

    print "Opening q and u cubes..."
    q=fits.getdata(params.qfile)
    print "... Q done ..."
    u=fits.getdata(params.ufile)
    print "... U done!"

    print "Getting rid of NaNs..."

    bn.replace(q,np.nan,0.0)

    print"...Q done..."

    bn.replace(u,np.nan,0.0)

    print "...U done!"

    if params.range != None:
        q=q[:,params.range[2]:params.range[3],params.range[0]:params.range[1]]
        u=u[:,params.range[2]:params.range[3],params.range[0]:params.range[1]]
        
    else:
        params.range=np.empty(4)
        params.range[0]=0
        params.range[1]=q.shape[2]-1
        params.range[2]=0
        params.range[3]=q.shape[1]-1

    #print np.sum(q),np.sum(u)

    return q.newbyteorder(),u.newbyteorder()


def new_header(params):

    old_header=fits.getheader(params.qfile)

    header=old_header.copy()

    date=datetime.datetime.today()

    header['NAXIS3']=params.nphi
    header['CTYPE3']=('phi', 'Faraday depth')
    header['CRPIX3']=(1,'reference pixel')
    header['CRVAL3']=(params.phi_min,'reference value')
    header['CDELT3']=(params.dphi,'phi step')

    header['OBJECT']= 'GALFACTS '+params.field+' RM cube'

    cpix_ra = (params.range[1] - params.range[0]) / 2. + 1.
    cpix_dec = (params.range[3] - params.range[2]) / 2. + 1.

    cpix_ra_old = old_header['CRPIX1'] - params.range[0]
    cpix_dec_old = old_header['CRPIX2'] - params.range[2]

    crval_ra = old_header['CRVAL1'] + (cpix_ra - cpix_ra_old) * old_header['CDELT1']
    crval_dec = old_header['CRVAL2'] + (cpix_dec - cpix_dec_old) * old_header['CDELT2']

    header['CRPIX1'] = cpix_ra
    header['CRVAL1'] = crval_ra
    header['NAXIS1'] = params.range[1]-params.range[0]

    header['CRPIX2'] = cpix_dec
    header['CRVAL2'] = crval_dec
    header['NAXIS2'] = params.range[3]-params.range[2]

    
    header.add_comment('Made with rmsyn-di.py')
    header.add_comment('on '+str(date))

    header.add_comment('FWHM of FPSF: {0:10.2f} rad/m^2'.format(params.fwhm))
    header.add_comment('Largest phi-scale: {0:10.2f} rad/m^2'.format(params.phi_scale))

    return header


def output_cube_fpsf(cube,fpsf,params):

    print "Output is in "+params.outdir

    print "Writing RM cube to "+params.outfile+".fits "
    
    if params.nohead:
        fits.writeto(params.outdir+params.outfile+".fits",np.abs(cube))
    else:
        rm_header=new_header(params)
        fits.writeto(params.outdir+params.outfile+".fits",np.abs(cube), rm_header)

        
    print "Writing FPSF to "+params.outfile+".txt"

    fpsf_out=np.zeros((params.nphi,3))
    fpsf_out[:,0]=params.phi
    fpsf_out[:,1]=fpsf.real
    fpsf_out[:,2]=fpsf.imag

    np.savetxt(params.outdir+params.outfile+".txt", fpsf_out)



def params_from_args():
    print 'Reading in parameters...'
    parser = argparse.ArgumentParser()
    parser.add_argument("q_in",help="Name of Q fits file")
    parser.add_argument("u_in",help="Name of U fits file")
    parser.add_argument("outfile",help="Output file prefix")
    parser.add_argument("outdir",help="Output directory")
    parser.add_argument("field",help="GALFACTS Field")
    parser.add_argument("-w","--weights",dest="weights_in",help="Optional weight file")
    parser.add_argument("-p", "--phi",dest="phi_in", nargs=3, type=int,help="phi axis parameters: nphi,dphi,phi_min")
    parser.add_argument("-r", "--range", dest="xyrange", nargs=4, type=int, help="x and y pixel range: xmin, xmax, ymin, ymax")
    parser.add_argument("-n","--nohead",help="output dirty cube without header", action="store_true")

    args=parser.parse_args()

    parameters=Params()
    parameters.qfile=args.q_in
    parameters.ufile=args.u_in
    parameters.outfile=args.outfile
    parameters.outdir=args.outdir
    parameters.field=args.field
    if args.phi_in != None:
        parameters.nphi=args.phi_in[0]
        parameters.dphi=args.phi_in[1]
        parameters.phi_min=args.phi_in[2]
    if args.weights_in != None:
        parameters.weights=np.loadtxt(args.weights_in)
    if args.xyrange != None:
        parameters.range=np.asarray(args.xyrange)
    if args.nohead != None:
        parameters.nohead = args.nohead

    print "...done."

    return parameters

def main():

    params=params_from_args()

    #open fits files, get headers

    q,u = open_and_trim(params)

    #get l2
    
    #params.l2,params.l20,params.phi_max,params.phi_scale=get_l2_params(params)
    get_l2_params(params)

    #get fpsf
    
    #fpsf,params.fwhm=compute_fpsf(params.l2,params.l20,params.phi,params.weights)
    fpsf=compute_fpsf(params)

    #get dirty cube
    
    #rm_cube=make_di_cube(q,u,params.phi,params.l2,params.l20,params.weights)
    #rm_cube=make_di_cube(q,u,params)
    rm_cube=di.compute_dicube(q,u,params.phi,params.l2,params.weights,params.l20)
    
    #write out dirty cube, fpsf

    output_cube_fpsf(rm_cube,fpsf,params)
    
    
if __name__ == "__main__":
    main()
