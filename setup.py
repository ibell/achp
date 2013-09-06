
from distutils.core import setup, Extension
import subprocess,shutil,os,sys

# Obtain the numpy include directory.  This logic works across numpy versions.
## import numpy
## try:
##     numpy_include = numpy.get_include()
## except AttributeError:
##     numpy_include = numpy.get_numpy_include()
    
    
sys.argv += ['build_ext','--inplace','--reswig']

if '--reswig' in sys.argv:
    import subprocess
    subprocess.check_call(['swig','-python','-c++','Compressor.i'],cwd = 'src')
    sys.argv.remove('--reswig')

numpy_include=['']
                         
commons = dict(include_dirs = ['externals/coolprop/CoolProp'],
               libraries = ['CoolPropLib_MD'],
               library_dirs = ['externals/coolprop/wrappers/StaticLibrary/VS2008'],
               )
               
Compressor_module = Extension('src._Compressor',
                           sources=['src/Compressor_wrap.cxx', 'src/Compressor.cpp'],
                           **commons
                           )

setup (name = 'ACHP',
       version = '3.0.0dev',
       author      = "Ian Bell",
       author_email='ian.h.bell@gmail.com',
       url='http://achp.sourceforge.net',
       description = """ Steady-state system model using moving boundary models """,
       ext_modules = [Compressor_module],
       )
       
sys.path.insert(0,'src')
import Compressor
dv = Compressor.vectord(10)
C = Compressor.CompressorClass()
import numpy as np
C.set_P([1.0,2.0,3.0])
C.P = dv
print C.P
