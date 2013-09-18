
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
    subprocess.check_output(['swig','-python','-c++','-I../externals/coolprop','ACHP.i'],cwd = 'src')
    sys.argv.remove('--reswig')

numpy_include=['']
                         
commons = dict(include_dirs = ['externals/coolprop'],
               libraries = ['CoolPropLib_MD'],
               library_dirs = ['externals/coolprop/wrappers/StaticLibrary/VS2008'],
               )
               
ACHP_module = Extension('src._ACHP',
                        sources=['src/ACHP_wrap.cxx', 'src/Compressor.cpp', 'src/BPHE.cpp'],
                        **commons
                        )

setup (name = 'ACHP',
       version = '0.0.1dev',
       author      = "Ian Bell",
       author_email='ian.h.bell@gmail.com',
       url='http://achp.sourceforge.net',
       description = """ Steady-state system model using moving boundary models """,
       ext_modules = [ACHP_module],
       )
