#!/usr/bin/env python
# -*- coding: utf-8 -*-
from glob import glob
import os

from setuptools import Extension, setup


def get_extensions():
    from Cython.Build import cythonize
    import numpy

    extensions = [
        Extension(
            'fdlsgm',
            language='c++',
            sources=['fdlsgm.pyx'] + glob(os.path.join('src', '*.cc')),
            libraries=['m'],
            include_dirs=[numpy.get_include(), 'src'],
            depends=glob(os.path.join('src', '*.h')),
            extra_compile_args=['-std=c++11', '-O2'],
        )
    ]

    return cythonize(
        extensions,
        build_dir=os.path.join('build', 'cython'),
        language_level=3,
    )


setup(ext_modules=get_extensions())
