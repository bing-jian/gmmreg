#!/usr/bin/env python
#coding=utf-8

from distutils.core import setup, Extension

pygmmreg = Extension('pygmmreg',
                     define_macros = [('MAJOR_VERSION', '1'),
                                      ('MINOR_VERSION', '0')],
                     include_dirs = ['/usr/include'],
                     libraries = ['gmmreg_api', 'vnl'],
                     library_dirs = ['/usr/lib64', './build'],
                     sources = ['pygmmreg.cpp'])

setup (name = 'pygmmreg',
              version = '1.0',
              description = 'Python wrapper of GMMREG algorithms.',
              author = 'Bing Jian',
              author_email = 'bing.jian@gmail.com',
              url = 'http://www.python.org/doc/current/ext/building.html',
              long_description = '''
              This is a Python wrapper of the C++ implementation of GMMREG algorithms.
              ''',
              ext_modules = [pygmmreg])
