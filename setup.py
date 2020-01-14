#!/usr/bin/env python
# -*- coding: utf-8 -*-
from glob import glob
from distutils.core import setup
from distutils.extension import Extension
from distutils.dist import Distribution
import os,sys,re

version = '0.1'
description = ''
long_description = ''

classifiers = [
  'Development Status :: ? - Development/Unstable',
  'Programming Language :: Python :: 2',
  'License :: OSI Approved :: GNU Lesesr General Public License v3 '
  'or  later (LGPLv3+)',
  'Topic :: Scientific/Engineering :: Astronomy',
  'Intended Audience :: Science/Research']

if os.path.exists('fdlsgm.pyx'):
  USE_CYTHON = True
  filename   = 'fdlsgm.pyx'
else:
  USE_CYTHON = False
  filename   = 'fdlsgm.cc'


if __name__ == '__main__':
  try:
    import numpy
  except ImportError:
    raise SystemExit('NumPy is not available.')

  sources      = [filename] + glob(os.path.join('src', '*.cc'))
  depends      = glob(os.path.join('src', '*.h'))
  libraries    = ['m',]
  include_dirs = [numpy.get_include(), 'src']
  compile_args = ['-std=c++11']
  extensions   = [
    Extension('fdlsgm',
              language='c++',
              sources=sources,
              libraries=libraries,
              include_dirs=include_dirs,
              depends=depends,
              extra_compile_args=compile_args)
    ]

  if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, language_level=3)

  setup(
    name='fdlsgm',
    version=version,
    description=description,
    long_description=long_description,
    license='LGPLv3',
    platforms=['POSIX',],
    classifiers=classifiers,
    url='',
    author='Ryou Ohsawa',
    author_email='ohsawa@ioa.s.u-tokyo.ac.jp',
    requires=['numpy'],
    ext_modules=extensions)
