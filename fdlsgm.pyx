#!/usr/bin/env cython
# -*- coding: utf-8 -*-
from libcpp.vector cimport vector as vec
from libcpp.list cimport list as clist
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

  cdef cppclass baseline:
    double x0()
    double y0()
    double z0()
    double x1()
    double y1()
    double z1()
    int size()
    void dprint()
    clist[long] elements()

  cdef cppclass accumulator[N]:
    accumulator() except+
    void insert(const dls&, const parameter&)
    void reallocate(const parameter&)
    void coalesce(const parameter&)
    int count_baseline()
    int count_segment()
    const baseline& operator[]()

  cdef const vec[baseline] find_segments(
    const vec[dls]& pool, const size_t& size_threshold)


cdef extern from *:
  ctypedef int N360 '360'


class solver(object):
  def __init__(self, limit=6):
    self.pool = list()
    self.limit = limit

  def insert(self, double x1, double y1, double z1,
             double x2, double y2, double z2):
    self.pool.append((x1,y1,z1,x2,y2,z2))

  def calculate(self):
    cdef vec[dls] pool
    cdef vec[baseline] result
    cdef clist[long] elem
    for x1,y1,z1,x2,y2,z2 in self.pool:
      pool.push_back(dls(x1,y1,z1,x2,y2,z2))
    result = find_segments(pool=pool, size_threshold=self.limit)
    baseline = list()
    bvector = list()
    for n in range(result.size()):
      elem = result[n].elements()
      baseline.append([e for e in elem])
      bvector.append(
        (result[n].x0(),result[n].y0(),result[n].z0(),
         result[n].x1(),result[n].y1(),result[n].z1()))
    return baseline, bvector


def simple_solver_test(int group=10, int frame=10, int scatter=1000):
  cdef int m, n
  cdef double x0,y0,z0,x1,y1,z1,x2,y2,z2,vx,vy,vz
  from datetime import datetime
  segments = list()
  s = solver()
  for m in range(group):
    x0 = np.random.uniform(-100,100)
    y0 = np.random.uniform(-100,100)
    z0 = np.random.uniform(-100,100)
    vx = np.random.normal()
    vy = np.random.normal()
    vz = np.random.normal()
    for n in range(frame):
      x1 = x0+vx*n+np.random.normal()/10.
      y1 = y0+vy*n+np.random.normal()/10.
      z1 = z0+vz*n+np.random.normal()/10.
      x2 = x0+vx*(n+1)+np.random.normal()/10.
      y2 = y0+vy*(n+1)+np.random.normal()/10.
      z2 = z0+vz*(n+1)+np.random.normal()/10.
      segments.append((m,x1,y1,z1,x2,y2,z2))
      s.insert(x1,y1,z1,x2,y2,z2)
  for m in range(scatter):
    x1 = np.random.uniform(-100,100)
    y1 = np.random.uniform(-100,100)
    z1 = np.random.uniform(-100,100)
    x2 = x1+np.random.normal()
    y2 = y1+np.random.normal()
    z2 = z1+np.random.normal()
    segments.append((-1,x1,y1,z1,x2,y2,z2))
    s.insert(x1,y1,z1,x2,y2,z2)

  t0 = datetime.now()
  elem,vector = s.calculate()
  t1 = datetime.now()
  print('## elements')
  for s in segments:
    print('{} {} {} {} {} {} {}'.format(*s))
  print('\n\n## baselines')
  for n,v in enumerate(vector):
    print('{} {} {} {} {} {} {}'.format(n,*v))
  print('\n# elapsed time: {}ms'.format((t1-t0).microseconds/1e3))
