import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport atan2
from libc.math cimport round #use atan2, round from c math library

cdef extern from "math.h":
    double M_PI

DTYPEf = np.float #this implies float64!
ctypedef np.float_t DTYPEf_t

DTYPEi = np.int
ctypedef np.int_t DTYPEi_t


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False) # turn off negative indexing for entire function

def make_angle_cube(np.ndarray[DTYPEf_t, ndim=3] qcube, np.ndarray[DTYPEf_t, ndim=3] ucube, np.ndarray[DTYPEf_t, ndim=2] rmmap, np.ndarray[DTYPEf_t, ndim=1] l2):

    cdef int i,j,k

    cdef int cchan = qcube.shape[0]/2

    cdef np.ndarray[DTYPEf_t,ndim=3] angle = np.empty([qcube.shape[0],qcube.shape[1],qcube.shape[2]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=3] target = np.empty([qcube.shape[0],qcube.shape[1],qcube.shape[2]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=3] npi = np.empty([qcube.shape[0],qcube.shape[1],qcube.shape[2]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEf_t,ndim=2] cangle = np.empty([qcube.shape[1],qcube.shape[2]],dtype=DTYPEf)
    

    for j in range(qcube.shape[1]):
        #print "Doing dec slice",j
        for k in range(qcube.shape[2]):
            cangle[j,k]=0.5*atan2(ucube[cchan,j,k],qcube[cchan,j,k])
            for i in range(qcube.shape[0]):
                angle[i,j,k]=0.5*atan2(ucube[i,j,k],qcube[i,j,k])
                target[i,j,k]=rmmap[j,k]*(l2[i]-l2[cchan])
                target[i,j,k]+=cangle[j,k]
                npi[i,j,k]=round((target[i,j,k]-angle[i,j,k])/M_PI)
                angle[i,j,k]+=M_PI*npi[i,j,k]

    return angle
    

    



    