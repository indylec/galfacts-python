import numpy as np
import sys
import os
import datetime
import argparse
from astropy.io import ascii
from astropy.io import fits

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
        temp = np.arange(self.nphi)
        self.phi = self.phi_min + temp * self.dphi
        self.phi_max=0.
        self.phi_scale=0.
        self.fwhm=0.
        self.field = ''

def get_l2_params(params):
    """
    get_l2(header, weights=None)

    Gets the frequency axis from FITS cube header and returns the corresponding
    lambda-squared axis

    Inputs:
        file - one of the data cubes being used for the rm synthesis 

        weights - (optional) array containing the weights of the different frequency
        channels. Weights must be between 0 and 1, and the array must be of the same
        length as the frequency axis. If this is omitted, all values are assumed to
        have  weight 1.

    Outputs:
        
        l2 - numpy array containing the frequencies converted to wavelength-squared

        l20 - the weighted mean of the l2 array

        phi_max - theoretical maximum Faraday depth which can be probed

        phi_scale - largest scale which can be measured (for extended components)

        
        
    """
    c2=299792458.**2

    header=fits.getheader(params.qfile)
    
    dnu=header['CDELT3']
    nchan=header['NAXIS3']
    nuref=header['CRVAL3']

    dl2 = c2/dnu**2

    nu = np.arange(nchan)*dnu+nuref

    l2 = 0.5 * c2 * ((nu - 0.5 * dnu) ** -2 + (nu + 0.5 * dnu) ** -2)
    l2 = np.flipud(l2)
        
    if params.weights != None:
        if header['NAXIS3']!= params.weights.shape[0]:
            raise Exception ('Weights and freq.axis have different sizes')
        else:
            l20 = np.sum(params.weights*l2)/np.sum(params.weights)
    else:
        l20 = np.sum(l2)/nchan

    phi_max = np.sqrt(3.)/dl2
    
    phi_scale = np.pi/np.amin(l2)

    return l2, l20, phi_max, phi_scale

            
def compute_fpsf(l2, l20, phi, weights=None):
    """
    compute_fpsf(l2, phi, weights=None)

    Computes the Faraday Point Spread Function for a given lambda-squared distribution
    and phi axis. Also returns resolution parameters specifying the width, FWHM and
    maximum phi of the FPSF. Weighting of the frequency axis can optionally be
    specified.

    input:
        l2 - numpy array containing the lambda-squared values

        l20 - weighted mean of the lambda-squared values 

        phi - numpy array containing the faraday depth (phi) values

        weights - (optional) numpy array containing the weights of each lambda-squared
        value. The weights mus be between 0 and 1. If this is omitted, all values are
        assumed to have weight 1.

    output:
        fpsf - complex numpy array containing the FPSF, binned according to the
        provided phi distribution

        fwhm_fpsf - FWHM of the FPSF, in rad/m^2, providing a measure of the beam
        resolution in Faraday depth space. Computed from B&dB 2005.
        
    """
    fpsf=np.empty(phi.shape[0],dtype=complex)
    
    if weights == None:
        weights = np.ones(l2.shape[0])
    else:
        if weights.shape[0] != l2.shape[0]:
            raise Exception ('weights array must have same length as lambda-squared\
                array')
    print "Making FPSF..."   

    for i in range(phi.shape[0]):
        fpsf[i] = np.sum(weights*np.exp(-2.*1j*phi[i]*(l2-l20)))/np.sum(weights)

    fwhm_fpsf = 2.0*np.sqrt(3)/(np.amax(l2)-np.amin(l2))

    print "...done."

    return fpsf, fwhm_fpsf
        
        

def make_di_cube(qcube, ucube, phi, l2, l20, weights = None):
    """
    make_di_cube(qcube, ucube, phi, l2, l20)

    Creates a 'dirty' RM cube using slow discrete Fourier transform of the complex
    polarisation.

    input:
        qcube - 3D numpy array (RA, dec, frequency) of Stokes Q intensity. Must have
        same dimensions as Stokes U cube.

        ucube - 3D numpy array (RA, dec, frequency) of Stokes U intensity. Must have
        same dimensions as Stokes Q cube.

        phi - numpy array of faraday depth bins

        l2 - numpy array containing wavelength squared axis of data cubes

        l20 - weighted average of lambda-squared

        weights - (optional) numpy array containing the weights of each lambda-squared
        value. The weights mus be between 0 and 1. If this is omitted, all values are
        assumed to have weight 1.

    output:
        rm_di_cube - 3D complex numpy array (RA, dec, phi) containing complex
        polarisation as a function of faraday depth.
        
    """
    if weights == None:
        weights = np.ones(l2.shape[0])
    else:
        if weights.shape[0] != l2.shape[0]:
            raise Exception ('weights array must have same length as lambda-squared\
                array')
    
    if qcube.shape != ucube.shape:
        raise Exception('qcube and ucube have different dimensions! Not allowed.')

    s=qcube.shape

    rm_di_cube=np.empty((phi.shape[0],s[1],s[2]),dtype=complex)

    temp_los=np.empty((s[0]), dtype=complex)

    print "Making RM cube..."
    for i in range(s[2]):
        for j in range(s[1]):
            print 'Synthesising pixel {0},{1}'.format(i,j)
            temp_los=qcube[:,j,i]+1j*ucube[:,j,i]
            for p in range(phi.shape[0]):
                rm_di_cube[p,j,i] = np.sum(temp_los*weights*\
                                        np.exp(-2.*1j*phi[p]*(l2-l20)))/np.sum(weights)
    print "...done."    
    
    return rm_di_cube

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
    header.add_comment('Made with rmsyn-di.py')
    header.add_comment('on '+str(date))

    header.add_comment('FWHM of FPSF: {0:10.2f} rad/m^2'.format(params.fwhm))
    header.add_comment('Largest phi-scale: {0:10.2f} rad/m^2'.format(params.phi_scale))

    return header


def output_cube_fpsf(cube,fpsf,params):

    rm_header=new_header(params)

    print "Output is in "+params.outdir

    print "Writing RM cube to "+params.outfile+".fits "
    
    fits.writeto(params.outdir+params.outfile+".fits",np.abs(cube),rm_header)

    print "Writing FPSF to "+params.outfile+".txt"

    fpsf_out=np.zeros((params.nphi,3))
    fpsf_out[:,0]=params.phi
    fpsf_out[:,1]=fpsf.real
    fpsf_out[:,2]=fpsf.imag

    np.savetxt(params_outdir+params.outfile+".txt", fpsf_out)



def params_from_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("q_in",help="Name of Q fits file")
    parser.add_argument("u_in",help="Name of U fits file")
    parser.add_argument("outfile",help="Output file prefix")
    parser.add_argument("outdir",help="Output directory")
    parser.add_argument("field",help="GALFACTS Field")
    parser.add_argument("-w","--weights",dest="weights_in",help="Optional weight file")
    parser.add_argument("-p", "--phi",dest="phi_in", nargs=3, type=int,help="phi axis parameters: nphi,dphi,phi_min")

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
        parameters.weights=np.loadtxt(args.weights_in).astype(int)

    return parameters

def main():

    params=params_from_args()

    #open fits files, get headers

    q=fits.getdata(params.qfile)
    u=fits.getdata(params.ufile)

    #get l2
    
    params.l2,params.l20,params.phi_max,params.phi_scale=get_l2_params(params)

    #get fpsf
    
    fpsf,params.fwhm=compute_fpsf(params.l2,params.l20,params.phi,params.weights)

    #get dirty cube
    
    rm_cube=make_di_cube(q,u,params.phi,params.l2,params.l20,params.weights)

    #write out dirty cube, fpsf

    output_cube_fpsf(rm_cube,fpsf,params)
    
    
if __name__ == "__main__":
    main()
