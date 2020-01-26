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

  /** signed index type. */
  typedef signed long int index_t;

  /** two-dimensional vector. */
  template <class T> using vector2 = std::array<T, 2>;

  /** three-dimensional vector. */
  template <class T> using vector3 = std::array<T, 3>;

  /** four-dimensional vector. */
  template <class T> using vector4 = std::array<T, 4>;

  /** 3x3 matrix. */
  template <class T> using matrix3x3 = std::array<vector3<T>,3>;

  /** 4x4 matrix. */
  template <class T> using matrix4x4 = std::array<vector4<T>,4>;

  /** alias to 4x4 zero matrix. */
  constexpr matrix4x4<double> zeros4x4 =
    {vector4<double>{0,0,0,0},vector4<double>{0,0,0,0},
     vector4<double>{0,0,0,0},vector4<double>{0,0,0,0}};

  /** alias to 4x4 identity matrix. */
  constexpr matrix4x4<double> ident4x4 =
    {vector4<double>{1,0,0,0},vector4<double>{0,1,0,0},
     vector4<double>{0,0,1,0},vector4<double>{0,0,0,1}};

  /**
   * @brief normalize two-dimensional vector.
   * @param[in] v: input two-dimensional vector.
   * @return normalized two-dimensional vector.
   */
  const vector2<double>
  normalize_vector(const vector2<double>& v);
  /**
   * @brief normalize three-dimensional vector.
   * @param[in] v: input three-dimensional vector.
   * @return normalized three-dimensional vector.
   */
  const vector3<double>
  normalize_vector(const vector3<double>& v);

  /**
   * @brief calculate outer product of three-dimensional vectors.
   * @param[in] v0: three-dimensional vector, left-hand-side.
   * @param[in] v1: three-dimensional vector, right-hand-side.
   * @return the result of outer product (v0 x v1).
   */
  const vector3<double>
  outer_product(const vector3<double>& v0, const vector3<double>& v1);
  /**
   * @brief calculate outer product of four-dimensional vectors.
   * @param[in] v0: four-dimensional vector, left-hand-side.
   * @param[in] v1: four-dimensional vector, right-hand-side.
   * @return the result of outer product (v0 x v1).
   * @note only the first three elements are used for the calculation.
   */
  const vector3<double>
  outer_product(const vector4<double>& v0, const vector4<double>& v1);

  /**
   * @brief solve the Ax = b type equation.
   * @param[in] A: 3x3 coefficient matrix.
   * @param[in] b: 3x1 vector.
   * @return the resultant three-dimensional vector.
   */
  const vector3<double>
  linalg_gauss_solve(const matrix3x3<double>& A, const vector3<double>& b);
  /**
   * @brief solve the Ax = b type equation.
   * @param[in] Ab: 4x4 coefficient matrix.
   * @return the resultant three-dimensional vector.
   * @note The fourth column is dropped. The first 3x3 part is used as the
   *       coefficient matrix A. The remaining part is used as the vector b.
   */
  const vector3<double>
  linalg_gauss_solve(const matrix4x4<double>& Ab);

  /**
   * @brief calculate eigen vectors with the Jacobi algorithm.
   * @param[in] M: input 4x4 matrix
   * @return 4x4 matrix containing four eigen vectors.
   * @note the eigen vectors are sorted in the ascending order of the
   *       corresponding singular values. The first eigen vector is the
   *       one for the least significant singular value.
   */
  const matrix4x4<double>
  eigenvector_jacobi_4x4(const matrix4x4<double>& M);
}

#endif // __LINALG_H_INCLUDE
