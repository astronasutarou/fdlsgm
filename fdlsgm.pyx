#!/usr/bin/env cython
# -*- coding: utf-8 -*-
from libcpp.list cimport list as clist
import numpy as np


cdef clist[int] count_number(int N):
  cdef int n = 0
  cdef clist[int] ret
  for n in range(N):
    ret.push_back(n)
  return ret

def count(int N):
  cdef clist[int] aaa = count_number(N)
  for x in aaa:
    print(x)
