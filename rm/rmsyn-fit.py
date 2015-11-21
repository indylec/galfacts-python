import numpy as np
import sys
import os
import datetime
import argparse
from astropy.io import ascii
from astropy.io import fits
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
        self.nphi = 100
        self.dphi = 10
        self.phi_min = -500
        temp = np.arange(self.nphi, dtype=int)
        self.phi = self.phi_min + temp * self.dphi
        self.l2 = None
        self.l20 = 0.
        self.weights = None
        self.rms_q = None
        self.rms_u = None


def fit_rm_peak(params):

    print "Opening RM cube..."
    rmcube=fits.getdata(params.rmfile)
    print rmcube.shape
    print "...done."
    range_no=21
    fit_range=np.arange(range_no)-range_no/2
    rm0=np.empty((params.dec_size,params.ra_size))
    print rm0.shape
    #print "phi:",params.phi

    print "Beginning peak fits..."
    for y,x in it.product(range(params.dec_size),range(params.ra_size)):
        #print "Finding peak for pixel ({0},{1})".format(x,y)
        temp_los=rmcube[:,y,x]
        #print temp_los
        
        temp_peak_arg=np.nanargmax(temp_los)
        #print temp_peak_arg
        
        if (temp_peak_arg>=range_no/2 and temp_peak_arg<params.nphi-range_no/2):
            temp_range=fit_range+temp_peak_arg
            #print 1,temp_range
            temp_coeffs=np.polyfit(params.phi[temp_range],temp_los[temp_range],2)
            rm0[y,x]=-temp_coeffs[1]/(2.*temp_coeffs[0])
        else:
            rm0[y,x]=params.phi[np.nanargmax(temp_los)]
    print "...done."

    return rm0

def open_and_trim(params):

    print "Opening q and u cubes..."
    print params.qfile
    qhdu=fits.open(params.qfile)

    qdata=qhdu[0]

    qerrdata=qhdu[1]

    q=qdata.data

    qerr=qerrdata.data

    print "... Q done ..."
    uhdu=fits.open(params.ufile)
    print params.ufile
    udata=uhdu[0]
    
    uerrdata=uhdu[1]

    u=udata.data

    uerr=uerrdata.data

    print "... U done!"

    print "Getting rid of NaNs..."

    q[np.isnan(q)]=0.0
    print"...Q done..."

    u[np.isnan(u)]=0.0

    print "...U done!"

    return q.astype(np.float64),u.astype(np.float64),qerr.astype(np.float64),uerr.astype(np.float64)

def make_angle_and_err(rm_map,params):

    qcube,ucube,qerr,uerr=open_and_trim(params)

    print qcube.shape,ucube.shape,qerr.shape,uerr.shape,rm_map.shape

    print "Making angle cube and errors..."
    angle,angle_err=ac.make_angle_cube(qcube,ucube,qerr,uerr,rm_map,params.l2)
    print "...done."


    return angle,angle_err


def get_l2_l20(params):
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


#def find_syn_angle(cube,syn_rm,params):
#    syn_angle=np.empty((params.dec_size,params.ra_size))
#    from scipy import interpolate
#    for y,x in it.product(range(params.dec_size),range(params.ra_size)):
#        #print "Finding zero-angle for pixel ({0},{1})".format(x,y)
#        spl=interpolate.UnivariateSpline(params.l2,cube[:,y,x])
#        syn_angle[y,x]=spl(params.l20)-syn_rm[y,x]*params.l20
#
#    return syn_angle

def fit_angle_cube(angle_cube,angle_err,rm0,params):

    good_chans=np.nonzero(params.weights)[0]
    
    print "Fitting angle cube..."

    angle0, rm_map, ang0err, rm_err, chisq = cfit.fit_cube(angle_cube[good_chans,:,:],angle_err[good_chans,:,:],params.l2[good_chans])

    print "...done."

    #print "Getting synthesis zero-angles..."
    #syn_angle=find_syn_angle(angle_cube,rm0,params)
    #print "...done."
    

    return rm_map,angle0,ang0err,rm_err,chisq


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

    if args.weights_in:
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

def new_header(params,object,unit):

    old_header=fits.getheader(params.qfile)

    header=old_header.copy()

    date=datetime.datetime.today()

    header.remove('NAXIS3')
    header.remove('CTYPE3')
    header.remove('CRPIX3')
    header.remove('CRVAL3')
    header.remove('CDELT3')

    header['OBJECT']= 'GALFACTS '+params.field+' '+object

    header['BUNIT']=unit

    header.add_comment('Made with rmsyn-fit.py')
    header.add_comment('on '+str(date))


    return header

def output_maps(cube_rm,cube_rm_err,cube_angle,cube_angle_err,syn_rm,chisq,params):

    print "Output is in "+params.outdir

    print "Writing Cube-derived RM map and error to "+params.outfile+"_cube_rm.fits "
        
    cube_rm_header=new_header(params,"cube RM map","rad/m^2")
    cube_rm_err_header=new_header(params,"cube RM errors", "rad/m^2")

    crmhdu=fits.PrimaryHDU(cube_rm,cube_rm_header)
    crmerrhdu=fits.ImageHDU(cube_rm_err,cube_rm_err_header)

    crmhdulist=fits.HDUList([crmhdu,crmerrhdu])

    crmhdulist.writeto(params.outdir+params.outfile+"_cube_rm.fits")


    print "Writing Cube-derived angle map and error to "+params.outfile+"_cube_angle.fits "
        
    cube_angle_header=new_header(params,"cube ANGLE map","rad/m^2")
    cube_angle_err_header=new_header(params,"cube ANGLE errors", "rad/m^2")

    canglehdu=fits.PrimaryHDU(cube_angle,cube_angle_header)
    cangleerrhdu=fits.ImageHDU(cube_angle_err,cube_angle_err_header)

    canglehdulist=fits.HDUList([canglehdu,cangleerrhdu])

    canglehdulist.writeto(params.outdir+params.outfile+"_cube_angle.fits")

    print "Writing synthesis-derived RM map to "+params.outfile+"_syn_rm.fits "
        
    syn_rm_header=new_header(params,"syn RM map","rad/m^2")
    fits.writeto(params.outdir+params.outfile+"_syn_rm.fits",syn_rm, syn_rm_header)

    print "Writing chi-squared  map to "+params.outfile+"_chisq.fits "
        
    chisq_header=new_header(params,"chi-squared map"," ")
    fits.writeto(params.outdir+params.outfile+"_chisq.fits",chisq, chisq_header)


def main():
    params=params_from_args()

    get_l2_l20(params)

    syn_rm=fit_rm_peak(params)

    angle,ang_err=make_angle_and_err(syn_rm,params)


    cube_rm,cube_angle,cube_angle_err,cube_rm_err,chisq = fit_angle_cube(angle,ang_err,syn_rm,params)


    output_maps(cube_rm,cube_rm_err,cube_angle,cube_angle_err,syn_rm,chisq,params)


if __name__ == "__main__":
    main()
    

    
