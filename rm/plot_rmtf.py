import numpy as np
from astropy.io import fits
import argparse
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib import rc

rc('font',family='serif')
rc('text', usetex=True)
import sys

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
        self.nochunks=1
        self.chunkrange=None

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

    

    nu = np.arange(nchan)*dnu+nuref

    l2 = 0.5 * c2 * ((nu - 0.5 * dnu) ** -2 + (nu + 0.5 * dnu) ** -2)
    params.l2 = np.flipud(l2)

    params.dl2 = np.abs(l2[1]-l2[0])

    
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

    params.phi_max = np.sqrt(3.)/params.dl2
    
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

    params.fwhm = 3.8/(np.amax(params.l2)-np.amin(params.l2))
    
    

    print "...done."

    print "Writing RMSF to "+params.outfile+"_"+params.field+"_RMSF.txt"

    fpsf_out=np.zeros((params.nphi,3))
    fpsf_out[:,0]=params.phi
    fpsf_out[:,1]=fpsf.real
    fpsf_out[:,2]=fpsf.imag

    np.savetxt(params.outdir+params.outfile+"_"+params.field+"_RMSF.txt", fpsf_out)

    print "Writing l2/RMSF parameters to "+params.outfile+"_"+params.field+"_RMSF_params.txt"

    f=open(params.outdir+params.outfile+"_"+params.field+"_RMSF_params.txt",'a')
    f.write('######### RMSF /l2 distribution parameters ###########\n')
    f.write('Del_l2 fwhm l2_min phi_scale del_l2 phi_max\n')
    f.write(str(np.amax(params.l2)-np.amin(params.l2))+' '+str(params.fwhm)+' '+str(np.amin(params.l2))+' '+str(params.phi_scale)+' '+str(params.dl2)+' '+str(params.phi_max)+'\n')
    f.write('######################################################')
    f.close()
    

    return fpsf

def params_from_args():
    print 'Reading in parameters...'
    parser = argparse.ArgumentParser()
    parser.add_argument("q_in",help="Name of Q fits file")
    #parser.add_argument("u_in",help="Name of U fits file")
    parser.add_argument("outfile",help="Output file prefix")
    parser.add_argument("outdir",help="Output directory")
    parser.add_argument("field",help="GALFACTS Field")
    parser.add_argument("-w","--weights",dest="weights_in",help="Optional weight file")
    parser.add_argument("-p", "--phi",dest="phi_in", nargs=3, type=int,help="phi axis parameters: nphi,dphi,phi_min")
    parser.add_argument("-r", "--range", dest="xyrange", nargs=4, type=int, help="x and y pixel range: xmin, xmax, ymin, ymax")
    parser.add_argument("-n","--nohead",help="output dirty cube without header", action="store_true")
    parser.add_argument("-c","--chunks",type=int ,help="Number of chunks to split the field into for memory purposes. These will be written out as separate fits files.")

    args=parser.parse_args()

    parameters=Params()
    parameters.qfile=args.q_in
    #parameters.ufile=args.u_in
    parameters.outfile=args.outfile
    parameters.outdir=args.outdir
    parameters.field=args.field
    if args.phi_in:
        parameters.nphi=args.phi_in[0]
        parameters.dphi=args.phi_in[1]
        parameters.phi_min=args.phi_in[2]
    if args.weights_in:
        parameters.weights=np.loadtxt(args.weights_in)
    if args.xyrange:
        parameters.range=np.asarray(args.xyrange)
    if args.nohead:
        parameters.nohead = args.nohead
    if args.chunks:
        parameters.nochunks=args.chunks

    print "...done."

    return parameters

params=params_from_args()

#open fits files, get headers

get_l2_params(params)

rmtf=compute_fpsf(params)

rmtf_abs=np.sqrt(rmtf.real**2+rmtf.imag**2)

fig=plt.figure(figsize=(10,6))

ax=fig.add_subplot(111)

ax.plot(params.phi,rmtf_abs,'k-',linewidth=2)
ax.set_ylabel('$|\mathrm{R}(\phi)|$',fontsize=24)
ax.set_xlabel('$\phi[\mathrm{rad}/\mathrm{m}^2]$', fontsize=26)
#ax.yaxis.set_label_position('right')
#ax.yaxis.tick_right()
ax.tick_params(labelsize=22)
ax.xaxis.set_major_locator(MaxNLocator(prune='both'))

fig.savefig("/Users/leclercq/galfacts/aps/rm/"+params.field+"_RMSF_plot.pdf",dpi=200,bbox_inches="tight")
