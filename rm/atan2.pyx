import numpy as np
cimport numpy as np
cimport cython

DTYPEf = np.float #this implies float64!
ctypedef np.float_t DTYPEf_t

DTYPEi = np.int
ctypedef np.int_t DTYPEi_t

def make_angle_cube(np.ndarray[DTYPEf_t, ndim=3] qcube, np.ndarray[DTYPEf_t, ndim=3] ucube)