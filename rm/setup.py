from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("rmsyn_dicube.pyx")
)
