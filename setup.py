#!/usr/bin/env python
# -*- coding: utf-8 -*-
from glob import glob
import os

from setuptools import Extension, setup


def get_extensions():
    try:
        import numpy
    except ImportError as exc:
        raise SystemExit('NumPy is required to build fdlsgm.') from exc
    use_cython = os.path.exists('fdlsgm.pyx')
    source = 'fdlsgm.pyx' if use_cython else 'fdlsgm.cpp'
    extensions = [
        Extension(
            'fdlsgm',
            language='c++',
            sources=[source] + glob(os.path.join('src', '*.cc')),
            libraries=['m'],
            include_dirs=[numpy.get_include(), 'src'],
            depends=glob(os.path.join('src', '*.h')),
            extra_compile_args=['-std=c++11', '-O2'],
        )
    ]

    if use_cython:
        from Cython.Build import cythonize

        return cythonize(
            extensions,
            build_dir=os.path.join('build', 'cython'),
            language_level=3,
        )

    return extensions


setup(ext_modules=get_extensions())
