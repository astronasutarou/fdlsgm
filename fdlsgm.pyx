#!/usr/bin/env cython
# -*- coding: utf-8 -*-
from libcpp.vector cimport vector as vec
from libcpp.list cimport list as clist
from numpy cimport ndarray
import numpy as np


cdef extern from 'fdlsgm.h' namespace 'fdlsgm':

  ### Container of grouping parameters.
  ctypedef struct parameter:
    double argument_limit_base
    double argument_limit_element
    double distance_limit
    double gap_limit
    int size_limit

  ### View of Baseline Class.
  ctypedef struct baseline_view:
    double v0[3]
    double v1[3]
    int size
    vec[int] elements

  ### Find baselines by Grouping Method.
  cdef const vec[baseline_view] find_segments(
    const size_t& n_size, const double* pool,
    const size_t& size_threshold, const parameter& param_insert,
    const parameter& param_reallocate, const parameter& param_merge)


cdef extern from *:
  ### Accumulator Size (fixed to 360).
  ctypedef int N360 '360'


### Convert baseline_view into a Python baseline object
cdef dump_view(baseline_view v):
  cdef int n = v.size
  cdef list e = [v.elements[i] for i in range(n)]
  return baseline(vertex0=v.v0, vertex1=v.v1, size=n, elements=e)


### set parameter values to a parameter struct
cdef void set_parameter(parameter& param, object elem):
  param.argument_limit_base = elem.argument_limit_base
  param.argument_limit_element = elem.argument_limit_element
  param.distance_limit = elem.distance_limit
  param.gap_limit = elem.gap_limit
  param.size_limit = int(elem.size_limit)


class baseline(object):
  ''' Copy of baseline_view in Python.

  Attributes:
    vertex0 (list): x-, y-, and z-coordinates of the starting point.
    vertex1 (list): x-, y-, and z-coordinates of the terminal point.
    size (int): the number of the elements composing the baseline.
    elements (list): the indexes of the elements composing the baseline.
  '''
  def __init__(self, list vertex0, list vertex1, int size, list elements):
    self.vertex0  = vertex0
    self.vertex1  = vertex1
    self.size     = size
    self.elements = elements

  def __len__(self):
    return self.size


class parameter_collection(object):
  ''' Parameter collection for procedures.

  The parameter collection for the insert, reallocation, and merge
  procedures inside the `find_segments` function. The meanings of the
  attributes are written below. Note that `size_limit` attribute is not
  used in the insert and reallocation procedures.

  Attributes:
    argument_limit_base: threshold of angle difference
    argument_limit_element: threshold of angle differnce
    distance_limit: threshold of lateral distance between segments
    gap_limit: threshold of gap between segments
    size_limit: threshold for the size of baseline
  '''
  def __init__(self, lb=None, le=None, ld=None, lg=None, ls=None):
    self.argument_limit_base = lb or 10.0*np.pi/180.0
    self.argument_limit_element = le or 20.0*np.pi/180.0
    self.distance_limit = ld or 3.0
    self.gap_limit = lg or 0.5
    self.size_limit = ls or 0


class solve_parameters(object):
  ''' Whole parameter collections for the solver.

  The instance contains the parameters for the solver. The attribute `limit`
  is the integer, which determine the minimum size of the baselines. The
  other attributes --- insert, reallocate, and merge --- are the instances
  of `parameter_collection`, which contain the sets of parameters for the
  collesponding operations.

  Attributes:
    limit (int): baselines whose size are smaller than `limit` are dropped.
    insert (list): parameter collection for the "insert" operation.
    reallocate (list): parameter collection for the "reallocate" operation.
    merge (list): parameter collection for the "merge" operation.
  '''
  def __init__(self, limit = 6,
        insert = (10.0*np.pi/180.0, 30.0*np.pi/180.0, 5.0, 0.5, 0),
        reallocate = (10.0*np.pi/180.0, 30.0*np.pi/180.0, 5.0, 0.5, 0),
        merge = (5.0*np.pi/180.0, 10.0*np.pi/180.0, 3.0, 0.5, 3)):
    self.limit = limit
    self.insert = parameter_collection(*insert)
    self.reallocate = parameter_collection(*reallocate)
    self.merge = parameter_collection(*merge)

  def __str__(self):
    return ''.join((
      '## FDLSGM parameters:\n'
      '#   Limit: {}\n'.format(self.limit),
      '#   insert    : ({:.4f}, {:.4f}, {:.1f}, {:.1f}, {:1d})\n'.format(
        self.insert.argument_limit_base,
        self.insert.argument_limit_element,
        self.insert.distance_limit,
        self.insert.gap_limit, self.insert.size_limit),
      '#   reallocate: ({:.4f}, {:.4f}, {:.1f}, {:.1f}, {:1d})\n'.format(
        self.reallocate.argument_limit_base,
        self.reallocate.argument_limit_element,
        self.reallocate.distance_limit,
        self.reallocate.gap_limit, self.reallocate.size_limit),
      '#   merge     : ({:.4f}, {:.4f}, {:.1f}, {:.1f}, {:1d})'.format(
        self.merge.argument_limit_base,
        self.merge.argument_limit_element,
        self.merge.distance_limit,
        self.merge.gap_limit, self.merge.size_limit),
      ))


def default_parameters():
  ''' Return a default parameter set. '''
  return solve_parameters()


def solve(ndarray pool, object param = None):
  ''' Find line segments from a set of elementary line segments.

  Parameters:
    pool (numpy.ndarray): N x 6 array containing the coordinates.
    param (solve_parameters): Parameter set for solver.

  Return:
    list of baseline.
  '''
  if isinstance(pool,np.ndarray):
    pool = pool.flatten()
  if param is None:
    param = default_parameters()
  if not isinstance(param, solve_parameters):
    raise ValueError('"param" should be an instance of "solve_parameters".')
  cdef vec[baseline_view] result
  cdef double[:] view = pool
  cdef size_t n_size = len(pool)//6
  cdef size_t limit = param.limit
  cdef parameter param_insert, param_reallocate, param_merge
  set_parameter(param_insert, param.insert)
  set_parameter(param_reallocate, param.reallocate)
  set_parameter(param_merge, param.merge)
  result = find_segments(
    n_size=n_size, pool=&view[0],
    size_threshold=limit, param_insert=param_insert,
    param_reallocate=param_reallocate, param_merge=param_merge)
  baseline = list()
  for n in range(result.size()):
    baseline.append(dump_view(result[n]))
  return baseline


def simple_solver_test(
    int group=30, int frame=12, int jammer=100, double scatter=0.01):
  ''' Functional test with a simple situation.

  Put `group` line segments composed of `frame` elementary segments in the
  2000 x 2000 x 2000 space. The loci, lengths, and directions of the line
  segments are randomly selected. The line segments are slightly affected
  by disturbance. The degree of the disturbance is given by `scatter`.
  Additionally, `jammer` obstacle segments are appended.

  This function calls the `solve` routine, which will search line segments
  in the space. The function prints the starting and terminal points of the
  all the elementary line segments, the starting and terminal points of the
  detected line segments, the default search parameter, the number of the
  detected line segments, and the elapsed time in milliseconds.

  Parameters:
    group (int): the number of line segments.
    frame (int): the number of elements composing a line segment.
    jammer (int): the number of obstacle segments.
    scatter (double): the scatter of line segments.
  '''
  cdef int n
  cdef double x0,y0,z0,x1,y1,z1,x2,y2,z2,vx,vy,vz
  cdef list baseline
  from datetime import datetime
  segments = list()
  shape = (frame,3)
  for n in range(group):
    p0 = np.random.uniform(-1000,1000, size=(1,3))
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
  param = default_parameters()

  t0 = datetime.now()
  baseline = solve(segments, param)
  t1 = datetime.now()
  print('## elements')
  for s in segments:
    print('{} {} {} {} {} {}'.format(*s))
  print('\n\n## baselines')
  for b in baseline:
    print('{} {} {} {} {} {}'.format(*b.vertex0,*b.vertex1))
  print('')
  print(param)
  print('#')
  print('# found {} segments'.format(len(baseline)))
  print('# elapsed time: {}ms'.format((t1-t0).total_seconds()*1e3))
