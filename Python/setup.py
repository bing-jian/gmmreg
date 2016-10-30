#!/usr/bin/env python
#coding=utf-8

"""
python setup.py build  -c mingw32
python setup.py install --skip-build
"""

from os.path import join
from distutils.core import setup, Extension
import numpy

module_pycvgmi = Extension('gmmreg._extension',
                            define_macros = [('MAJOR_VERSION', '1'),
                                             ('MINOR_VERSION', '0')],
                            sources = [join('c_extension','py_extension.c'),
                                       join('c_extension','GaussTransform.c'),
                                       join('c_extension','DistanceMatrix.c')],
                            include_dirs =[numpy.get_include()])


setup (name = 'gmmreg',
       version = '1.0',
       description = 'A Python package for robust point set registration using mixture of Gaussians.',
       author = 'Bing Jian',
       author_email = 'bing.jian@gmail.com',
       url = 'http://gmmreg.googlecode.com/',
       long_description = '''
              This is a Python package for robust point set registration using mixture of Gaussians.
              ''',
       package_dir={'gmmreg': ''},
       packages=['gmmreg'],
       ext_modules = [module_pycvgmi])
