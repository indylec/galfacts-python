#cython module containing functions to compute the fpsf and dirty cube

import numpy as np
cimport numpy as np
cimport cython


#from libc.math cimport exp #use exp from c math library

cdef extern from "complex.h":
    double complex cexp(double complex)

DTYPEf = np.float #this implies float64!
ctypedef np.float_t DTYPEf_t
DTYPEi = np.int
ctypedef np.int_t DTYPEi_t
DTYPEc = np.complex128
ctypedef np.complex128_t DTYPEc_t

cdef cython_fsum(np.ndarray[DTYPEf_t, ndim=1] y):
    cdef int N = y.shape[0]
    cdef float x = y[0]
    cdef int i
    for i in xrange(N):
        x += y[i]
    return x

cdef cython_csum(np.ndarray[DTYPEc_t, ndim=1] y):
    cdef int N = y.shape[0]
    cdef double complex x = y[0]
    cdef int i
    for i in xrange(N):
        x += y[i]
    return x

#cdef cython_fmean(np.ndarray[DTYPEf_t, ndim=1] y):
#    cdef int N = y.shape[0]
#    cdef float x = y[0]
#    cdef int i
#    for i in xrange(N):
#        x += y[i]
#    return x/N


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False) # turn off negative indexing for entire function

#make sure there are no remaining nans before feeding the cubes to the function

def compute_dicube(np.ndarray[DTYPEf_t, ndim=3] qcube, np.ndarray[DTYPEf_t, ndim=3] ucube, np.ndarray[DTYPEi_t, ndim=1] phi, np.ndarray[DTYPEf_t, ndim=1] l2, np.ndarray[DTYPEf_t, ndim=1] weight, float l20, ):

    cdef int i,j,k,p
    cdef int ra = qcube.shape[2]
    cdef int dec = qcube.shape[1]
    cdef double complex losk,expk,expk_2 #,weightedpk
    cdef float weightk,expk_1
    
    

    cdef np.ndarray[DTYPEc_t,ndim=3] dicube = np.empty([phi.shape[0],qcube.shape[1],qcube.shape[2]],dtype=DTYPEc)

    cdef np.ndarray[DTYPEc_t,ndim=1] temp_los = np.empty([qcube.shape[0]],dtype=DTYPEc)
    cdef np.ndarray[DTYPEc_t,ndim=1] weighted_p = np.empty([qcube.shape[0]],dtype=DTYPEc)
    cdef np.ndarray[DTYPEc_t,ndim=1] temp_exp = np.empty([qcube.shape[0]],dtype=DTYPEc)
    cdef np.ndarray[DTYPEf_t,ndim=1] temp_l2 = np.empty([qcube.shape[0]],dtype=DTYPEf)

    cdef np.ndarray[DTYPEc_t,ndim=1] sum_term = np.empty([qcube.shape[0]],dtype=DTYPEc)   
    cdef float inverse_sum_weight = 1./cython_fsum(weight)
    
    for i in range(ra):
        print 'Working on ra pixel',i
        for j in range(dec):
            temp_los=qcube[:,j,i]+ucube[:,j,i]*1j
            for p in range(phi.shape[0]):
                for k in range(qcube.shape[0]):
                    losk=temp_los[k]
                    weightk=weight[k]

                    weighted_p[k]=weightk*losk
                    #weighted_p[k]=temp_los[k]*weight[k]

                    temp_l2[k]=l2[k]-l20
                    expk_1=-2.*phi[p]*temp_l2[k]
                    expk_2=1j

                    expk=expk_1*expk_2
                    temp_exp[k]=cexp(expk)
		    #weightedpk=weighted_p[k]
		    

                    sum_term[k]=weighted_p[k]*temp_exp[k]
            dicube[p,j,i]=inverse_sum_weight*cython_csum(sum_term)


    return dicube
		
	    	
	    
	    
	    
	    


    
    