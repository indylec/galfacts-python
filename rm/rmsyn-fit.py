import numpy as np
import sys
import os
import datetime
import argparse
from astropy.io import ascii
from astropy.io import fits
import bottleneck as bn
import itertools as it
import cylinfit as cfit
import angle_cube as ac

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
        self.rmfile = ''
        self.field = ''
        #self.synth = False
        #self.fit = False
        self.ra_size = 0
        self.dec_size = 0
        self.nu_size = 0
        self.dnu = 0.
        self.nuref = 0.
        self.nphi = 0
        self.dphi = 0
        self.phi_min = 0
        temp = np.arange(self.nphi, dtype=int)
        self.phi = self.phi_min + temp * self.dphi
        self.l2 = None
        self.l20 = 0.
        self.weights = None

def fit_rm_peak(params):

    print "Opening RM cube..."
    rmcube=fits.getdata(params.rmfile)
    print "...done."
    range_no=21
    fit_range=np.arange(range_no)-range_no/2
    rm0=np.empty((params.dec_size,params.ra_size))

    print "Beginning peak fits..."
    for y,x in it.product(range(params.dec_size),range(params.ra_size)):
        print "Finding peak for pixel ({0},{1})".format(x,y)
        temp_los=rmcube[:,y,x]
        temp_peak_arg=bn.nanmaxarg(temp_los)
        if (temp_peak_arg>=range_no/2 and temp_peak_arg<=params.nphi-range_no/2):
            temp_range=fit_range+temp_peak_arg
            temp_coeffs=np.polyfit(params.phi[range],temp_los[range],2)
            rm0[y,x]=-temp_coeffs[1]/(2.*temp_coeffs[0])
        elif temp_peak_arg<range_no/2:
            temp_range=np.arange(2*temp_peak_arg+1)
            temp_coeffs=np.polyfit(params.phi[temp_range],temp_los[temp_range],2)
            rm0[y,x]=-temp_coeffs[1]/(2.*temp_coeffs[0])
        else:
            diff=params.nphi-temp_peak_arg
            temp_range=np.arange(2*diff+1)-diff+temp_peak_arg
            temp_coeffs=np.polyfit(params.phi[temp_range],temp_los[temp_range],2)
            rm0[y,x]=-temp_coeffs[1]/(2.*temp_coeffs[0])
        print "...done."

    return rm0

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

    return q,u

def make_angle_cube(rm_map, params):

    qcube,ucube=open_and_trim(params)

    print "Making angle cube..."
    angle=a.make_angle_cube(qcube,ucube,rm_map,params.l2)
    print "...done."

##     cchan=params.nu_size/2

##     angle=0.5*np.atan2(ucube,qcube)

##     rm=np.repeat(rm_map.reshape(1,params.dec_size,params.ra_size),params.nphi,axis=0)

##     target=rm*np.reshape(params.l2-params.l2[cchan],(sz,1,1))

##     target += angle[cchan].reshape(1,sy,sx)

##     npi=np.around(((target-angle)/np.pi))
##     #npi[np.where(np.isnan(npi))]=0
##     angle += pi * npi

    return angle


def get_l2_l0(params):
    print "Computing l2, l20 and weights..."
    
    c2=299792458.**2
    
    dl2 = c2/(params.dnu**2)

    nu = np.arange(params.nu_size)*params.dnu+params.nuref

    l2 = 0.5 * c2 * ((nu - 0.5 * params.dnu) ** -2 + (nu + 0.5 * params.dnu) ** -2)
    params.l2 = np.flipud(l2)
        
    if params.weights != None:
        if params.nu_size != params.weights.shape[0]:
            raise Exception ('Weights and freq.axis have different sizes')
        else:
            params.l20 = np.sum(params.weights*params.l2)/np.sum(params.weights)
    else:
        params.weights=np.ones(params.nu_size)
        params.l20 = np.sum(params.l2)/params.nu_size
    print"...done."

def find_syn_angle(cube,rmmap,params):
    syn_angle=np.empty((params.dec_size,params.ra_size))
    from scipy import interpolate
    for y,x in it.product(range(params.dec_size),range(params.ra_size)):
        print "Finding zero-angle for pixel ({0},{1})".format(x,y)
        spl=interpolate.UnivariateSpline(params.l2,cube[:,y,x])
        syn_angle[y,x]=spl(params.l20)-rmmap[y,x]*params.l20

    return syn_angle

def fit_angle_cube(angle_cube,rm0,params):

    print "Fitting angle cube..."
    rm_map,angle0 = cfit.fit_cube(angle_cube,params.l2)
    print "...done."

    print "Getting synthesis zero-angles..."
    syn_angle=find_syn_angle(angle_cube,rm0,params)
    print "...done."
    
    return rm_map,angle0,syn_angle

def params_from_args():

    print 'Reading in parameters...'
    parser = argparse.ArgumentParser()
    parser.add_argument("rm_in",help="Name of RM cube")
    parser.add_argument("q_in",help="Name of Q fits file")
    parser.add_argument("u_in",help="Name of U fits file")
    parser.add_argument("outfile",help="Output file prefix")
    parser.add_argument("outdir",help="Output directory")
    parser.add_argument("field",help="GALFACTS Field")
    parser.add_argument("-w","--weights",dest="weights_in",help="Optional weight file")
    #parser.add_argument("-s","--synth",help="output rm0 and angle-0 from synthesis", action="store_true")
    #parser.add_argument("-f","--fit",help="output rm0 and angle-0 from  angle cube fit", action="store_true")
    

    args=parser.parse_args()

    parameters=Params()
    parameters.qfile=args.q_in
    parameters.ufile=args.u_in
    parameters.outfile=args.outfile
    parameters.outdir=args.outdir
    parameters.field=args.field
    parameters.rmfile=args.rm_in

    if args.weights_in != None:
        parameters.weights=np.loadtxt(args.weights_in)

    #if args.synth != None:
    #   parameters.synth=args.synth
    #if args.fit != None:
    #   parameters.fit=args.fit

    rmhead=fits.getheader(parameters.rmfile)

    parameters.ra_size=rmhead['NAXIS1']
    parameters.dec_size=rmhead['NAXIS2']
    parameters.nphi=rmhead['NAXIS3']
    parameters.phi_min=rmhead['CRVAL3']
    parameters.dphi=['CDELT3']

    quhead= fits.getheader(parameters.qfile)

    parameters.nu_size = quhead['NAXIS3']
    parameters.dnu = quhead['CDELT3']
    parameters.nuref = quhead['CRVAl3']
    
    print "...done."

    return parameters

def new_header(params,object):

    old_header=fits.getheader(params.qfile)

    header=old_header.copy()

    date=datetime.datetime.today()

    header.remove('NAXIS3')
    header.remove('CTYPE3')
    header('CRPIX3')
    header('CRVAL3')
    header('CDELT3')

    header['OBJECT']= 'GALFACTS '+params.field+' '+object

    header.add_comment('Made with rmsyn-fit.py')
    header.add_comment('on '+str(date))


    return header

def output_maps(cube_rm,cube_angle,syn_rm,syn_angle):

    print "Output is in "+params.outdir

    print "Writing Cube-derived RM map to "+params.outfile+"cube_rm.fits "
        
    cube_rm_header=new_header(params,"cube RM map")
    fits.writeto(params.outdir+params.outfile+"cube_rm.fits",cube_rm, cube_rm_header)

    print "Writing Cube-derived angle map to "+params.outfile+"cube_angle.fits "
        
    cube_angle_header=new_header(params,"cube angle map")
    fits.writeto(params.outdir+params.outfile+"cube_angle.fits",cube_angle, cube_angle_header)

    print "Writing synthesis-derived RM map to "+params.outfile+"cube_rm.fits "
        
    syn_rm_header=new_header(params,"syn RM map")
    fits.writeto(params.outdir+params.outfile+"syn_rm.fits",np.abs(cube), syn_rm_header)

    print "Writing synthesis-derived RM map to "+params.outfile+"cube_rm.fits "
        
    syn_angle_header=new_header(params,"syn angle map")
    fits.writeto(params.outdir+params.outfile+"syn_angle.fits",np.abs(cube), syn_angle_header)


def main():
    params=params_from_args()

    get_l2_l20(params)

    syn_rm=fit_rm_peak(params)

    angle=make_angle_cube(params)

    cube_rm,cube_angle,syn_angle = fit_angle_cube(angle,rm0,params)

    output_maps(cube_rm,cube_angle,syn_rm,syn_angle)


if __name__ == "__main__":
    main()
    

    
