import numpy as np
cimport numpy as np
cimport cython

DTYPEf = np.float #this implies float64!
ctypedef np.float_t DTYPEf_t

DTYPEi = np.int
ctypedef np.int_t DTYPEi_t

cdef cython_fsum(np.ndarray[DTYPEf_t, ndim=1] y):
    cdef int N = y.shape[0]
    cdef float x = y[0]
    cdef int i
    for i in xrange(N):
        x += y[i]
    return x

#linear regression function from Numerical Recipes
#a is offset, b is slope

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False) # turn off negative indexing for entire function

cdef linfit(np.ndarray[DTYPEf_t, ndim=1] y, np.ndarray[DTYPEf_t, ndim=1] x,np.ndarray[DTYPEf_t, ndim=1] y_unc):

    if y.shape[0] != x.shape[0]:
        raise ValueError('x and y must have the same size!')

    cdef np.ndarray[DTYPEf_t,ndim=1] ryu2 = np.empty([y.shape[0]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=1] xs = np.empty([y.shape[0]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=1] ys = np.empty([y.shape[0]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=1] t = np.empty([y.shape[0]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=1] t2 = np.empty([y.shape[0]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=1] ty = np.empty([y.shape[0]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=1] yfit = np.empty([y.shape[0]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=1] chi_term = np.empty([y.shape[0]],dtype=DTYPEf)

    cdef float S, Sx, Sy, Stt, b, a,tysum
    cdef int i,j
    cdef int chan=y.shape[0]

    for i in range(chan):
        if y_unc[i]==0:
            y_unc[i]=1.
        ryu2[i]=1./y_unc[i]**2
        xs[i]=x[i]*ryu2[i]
        ys[i]=y[i]*ryu2[i]

    S     = cython_fsum(ryu2)
    Sx    = cython_fsum(xs)
    Sy    = cython_fsum(ys)

    for i in range(chan):
        t[i]=1./y_unc[i] * (x[i]-Sx/S)
        t2[i]=t[i]**2
        ty[i]=t[i]*y[i]/y_unc[i]
	
    Stt = cython_fsum(t2)
    tysum = cython_fsum(ty)    

    b = tysum/Stt
    a = (Sy - Sx * b) / S

    covab= -Sx / (S*Stt)
    sa = np.sqrt(1. / S * (1. - Sx * covab))
    sb=np.sqrt(1./Stt)

    for j in range(chan):
        yfit[j] = a + b * x[j]
        chi_term[j]=((y[j]-yfit[j]) / y_unc[j])**2

    chisq=cython_fsum(chi_term)

    return a, b, sa, sb, chisq
    

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False) # turn off negative indexing for entire function

def fit_cube(np.ndarray[DTYPEf_t, ndim=3] cube,np.ndarray[DTYPEf_t, ndim=3] angerr,np.ndarray[DTYPEf_t, ndim=1] l2):

    cdef int i,j

    cdef np.ndarray[DTYPEf_t,ndim=2] rm = np.empty([cube.shape[1],cube.shape[2]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=2] angle0 = np.empty([cube.shape[1],cube.shape[2]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=2] rmerr = np.empty([cube.shape[1],cube.shape[2]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=2] angle0err = np.empty([cube.shape[1],cube.shape[2]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=2] chi2 = np.empty([cube.shape[1],cube.shape[2]],dtype=DTYPEf)

    for i in range(cube.shape[1]):
        for j in range(cube.shape[2]):
            angle0[i,j],rm[i,j],angle0err[i,j],rmerr[i,j],chi2[i,j]=linfit(cube[:,i,j],l2,angerr[:,i,j])

    return angle0, rm, angle0err, rmerr, chi2

    
