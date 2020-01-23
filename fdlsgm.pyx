#!/usr/bin/env cython
# -*- coding: utf-8 -*-
from libcpp.vector cimport vector as vec
from libcpp.list cimport list as clist
from numpy cimport ndarray
import numpy as np


cdef extern from 'fdlsgm.h' namespace 'fdlsgm':
  ctypedef struct parameter:
    double argument_limit_base
    double argument_limit_element
    double distance_limit
    double gap_limit
    int size_limit

  cdef cppclass dls:
    dls(const double&, const double&, const double&,
        const double&, const double&, const double&) except+

  ctypedef struct baseline_view:
    double v0[3]
    double v1[3]
    int size
    int* elements

  cdef const vec[baseline_view] find_segments(
    const size_t& n_size, const double* pool, const size_t& size_threshold)


cdef extern from *:
  ctypedef int N360 '360'


cdef view2dict(baseline_view v):
  cdef int n = v.size
  cdef list e = [v.elements[i] for i in range(n)]
  return dict(v0=v.v0, v1=v.v1, n=n, e=e)


def solve(ndarray pool, const int limit = 6):
  if isinstance(pool,np.ndarray):
    pool = pool.flatten()
  cdef vec[baseline_view] result
  cdef const double[:] view = pool
  cdef size_t n_size = len(pool)//6
  result = find_segments(n_size=n_size, pool=&view[0], size_threshold=limit)
  baseline = list()
  for n in range(result.size()):
    baseline.append(view2dict(result[n]))
  return baseline


def simple_solver_test(
    int group=10, int frame=10, int jammer=100, double scatter=0.05):
  cdef int n
  cdef double x0,y0,z0,x1,y1,z1,x2,y2,z2,vx,vy,vz
  cdef list baseline
  from datetime import datetime
  segments = list()
  shape = (frame,3)
  for n in range(group):
    p0 = np.random.uniform(-100,100, size=(1,3))
    v0 = np.random.normal(size=(1,3))
    v = np.vstack((np.zeros((1,3)),np.cumsum(v0*np.ones((frame,3)),axis=0)))
    p = np.hstack((p0+v[0:frame],p0+v[1:frame+1]))
    p = p + np.random.normal(0.0, scatter, size=p.shape)
    segments.extend(p)
  p0 = np.random.uniform(-100,100, size=(jammer,3))
  p1 = p0 + np.random.normal(size=(jammer,3))
  p = np.hstack((p0,p1))
  segments.extend(p)

  segments = np.array(segments)

  t0 = datetime.now()
  baseline = solve(segments)
  t1 = datetime.now()
  print('## elements')
  for s in segments:
    print('{} {} {} {} {} {}'.format(*s))
  print('\n\n## baselines')
  for n,b in enumerate(baseline):
    print('{} {} {} {} {} {}'.format(*b['v0'],*b['v1']))
  print('')
  print('# found {} segments'.format(len(baseline)))
  print('# elapsed time: {}ms'.format((t1-t0).microseconds/1e3))
