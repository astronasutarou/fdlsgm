# Fast Directed Line Segment Grouping Method

## Overview
This package will provide an algorithm and interface to group elementary line segments in terms of direction and vicinity. The algorithm used is based on the algorithm developed by Jang &amp; Hong (2002)[^JH2002]. The software efficiently finds line segments from a bunch of elemental directed line segments placed in a three-dimensional space.


## Interfaces
### C++

``` c++
const std::vector<baseline_view>
find_segments(const std::vector<dls>& pool,
              const size_t& size_threshold = 6,
              const parameter& param_insert = default_param_insert,
              const parameter& param_reallocate = default_param_reallocate,
              const parameter& param_coalesce = default_param_coalesce);
```


``` c++
const std::vector<baseline_view>
find_segments(const size_t& n_elements,
              const double* pool,
              const size_t& size_threshold = 6,
              const parameter& param_insert = default_param_insert,
              const parameter& param_reallocate = default_param_reallocate,
              const parameter& param_coalesce = default_param_coalesce);
```

### Python

``` Python
def solve(ndarray pool, object param = None):
  ''' Find line segments from a set of elementary line segments.

  Parameters:
    pool (numpy.ndarray): N x 6 array containing the coordinates.
    param (solve_parameters): Parameter set for solver.

  Return:
    list of baseline.
  '''
```

## Dependencies
The library is written in _C++11_ and do not depends on any library outside of `STL`. The Python interface is depends on `NumPy`. The library is developed on `g++` version 5.4 installed in Linux Mint 18.1 (serena) and Python 3.7.1.


## References
[^JH2002]: Jeong-Hun Jang &amp; Ki-Sang Hong, Pattern Recognition 35 (2002), 2235&ndash;2247 (doi: [10.1016/S0031-3203(01)00175-3](https://doi.org/10.1016/S0031-3203(01)00175-3 "Jand & Hong (2002)"))
