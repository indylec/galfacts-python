#cython module containing functions to compute the fpsf and dirty cube

import numpy as np
cimport numpy as np

from libc.math cimport exp #use exp from c math library

DTYPEf = np.float #this implies float64!
ctypedef np.float_t DTYPEf_t
DTYPEi = np.int
ctypedef np.int_t DTYPEi_t
DTYPEc = np.complex128
ctypedef= np.complex128_t DTYPEc_t


def cython_sum(np.ndarray[DTYPEf_t, ndim=1] y):
    cdef int N = y.shape[0]
    cdef int x = y[0]
    cdef int i
    for i in xrange(N):
        x += y[i]
    return x

def cython_mean(np.ndarray[DTYPEf_t, ndim=1] y):
    cdef int N = y.shape[0]
    cdef int x = y[0]
    cdef int i
    for i in xrange(N):
        x += y[i]
    return x/N


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False) # don't check for negative indices (also for entire function)
def compute_fpsf( np.ndarray[DTYPEf_t, ndim=1] nu, np.ndarray[DTYPEf_t, ndim=1] phi, np.ndarray[DTYPEf_t, dim=1] weights = None):

    cdef enum:
        c=299792458
    cdef int i,j,k,iphi
    cdef unsigned int nl2 = nu.shape[0]
    cdef unsigned int nphi = phi.shape[0]
    cdef np.ndarray[DTYPEf_t, ndim=1] l2 = np.zeros(nl2, dtype=DTYPEf_t)
    cdef np.ndarray[DTYPEc_t, ndim=1] fpsf = np.zeros(nphi, dtype=DTYPEc_t)

#calculate lambda2 from nu
    for i in range (nl2):
        l2[i]=(c/nu[i])**2


#find the mean lambda2
    l0=cython_mean(l2)

#now compute FPSF
    for iphi in range nphi:
        fpsf[nphi]=
        

    
    

    
    

