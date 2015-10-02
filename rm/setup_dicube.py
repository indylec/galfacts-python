from distutils.core import setup,Extension
from Cython.Build import cythonize
import numpy

module1 = Extension(    name         = "rmsyn_dicube",
                        sources      = ["rmsyn_dicube.pyx"],
                        include_dirs = [numpy.get_include()],
                    )

setup(
    ext_modules = cythonize(module1)
)
