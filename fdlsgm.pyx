#!/usr/bin/env python
from libc.math cimport sin
from libcpp.vector cimport vector
from libcpp.list cimport list


cdef double f (double x):
  return sin(x * x)
