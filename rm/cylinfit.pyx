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

#this function only calculates regression without uncertainties!

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False) # turn off negative indexing for entire function

cdef linfit(np.ndarray[DTYPEf_t, ndim=1] y, np.ndarray[DTYPEf_t, ndim=1] x, ):

    if y.shape[0] != x.shape[0]:
        raise ValueError('x and y must have the same size!')

    cdef np.ndarray[DTYPEf_t,ndim=1] ryu2 = np.ones([y.shape[0]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=1] t = np.empty([y.shape[0]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=1] t2 = np.empty([y.shape[0]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=1] ty = np.empty([y.shape[0]],dtype=DTYPEf)


    cdef float S, Sx, Sy, Stt, b, a,tysum
    cdef int i
    cdef int chan=y.shape[0]

    S     = cython_fsum(ryu2)
    Sx    = cython_fsum(x)
    Sy    = cython_fsum(y)

    for i in range(chan):
        t[i]=(x[i]-Sx/S)
        t2[i]=t[i]**2
        ty=t[i]*y[i]
	
    Stt = cython_fsum(t2)
    tysum = cython_fsum(ty)    

    b = tysum/Stt
    a = a = (Sy - Sx * b) / S

    return a, b

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False) # turn off negative indexing for entire function

def fit_cube(np.ndarray[DTYPEf_t, ndim=3] cube,np.ndarray[DTYPEf_t, ndim=1] l2):

    cdef int i,j

    cdef np.ndarray[DTYPEf_t,ndim=2] rm = np.empty([cube.shape[1],cube.shape[2]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=2] angle0 = np.empty([cube.shape[1],cube.shape[2]],dtype=DTYPEf)

    for i in range(cube.shape[1]):
        for j in range(cube.shape[2]):
            rm[i,j],angle0[i,j]=linfit(cube[:,i,j],l2)

    return rm, angle0

    
