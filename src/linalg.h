/**
 * @file fdlsgm.h
 * @brief Header file for the FDLSGM library
 * @author Ryou Ohsawa
 * @year 2020
 */
#ifndef __LINALG_H_INCLUDE
#define __LINALG_H_INCLUDE

#include<cstdio>
#include<cmath>
#include<limits>
#include<array>
#include<algorithm>
#include<utility>
#include<exception>


namespace linalg {

  /** Signed index type. */
  typedef signed long int index_t;

  /** Two-Dimensional Vector. */
  template <class T> using vector2 = std::array<T, 2>;

  /** Three-Dimensional Vector. */
  template <class T> using vector3 = std::array<T, 3>;

  /** Four-Dimensional Vector. */
  template <class T> using vector4 = std::array<T, 4>;

  /** 3x3 Matrix */
  template <class T> using matrix3x3 = std::array<vector3<T>,3>;

  /** 4x4 Matrix */
  template <class T> using matrix4x4 = std::array<vector4<T>,4>;

  constexpr matrix4x4<double> zeros4x4 =
    {vector4<double>{0,0,0,0},vector4<double>{0,0,0,0},
     vector4<double>{0,0,0,0},vector4<double>{0,0,0,0}};

  constexpr matrix4x4<double> ident4x4 =
    {vector4<double>{1,0,0,0},vector4<double>{0,1,0,0},
     vector4<double>{0,0,1,0},vector4<double>{0,0,0,1}};

  const vector2<double>
  normalize_vector(const vector2<double>&);
  const vector3<double>
  normalize_vector(const vector3<double>&);

  const vector3<double>
  outer_product(const vector3<double>&, const vector3<double>&);
  const vector3<double>
  outer_product(const vector4<double>&, const vector4<double>&);

  const vector3<double>
  linalg_gauss_solve(const matrix3x3<double>& A, const vector3<double>& b);
  const vector3<double>
  linalg_gauss_solve(const matrix4x4<double>&);

  /** Calculate Eigen Vectors with the Jacobi algorithm */
  const matrix4x4<double>
  eigenvector_jacobi_4x4(const matrix4x4<double>&);
}

#endif // __LINALG_H_INCLUDE
