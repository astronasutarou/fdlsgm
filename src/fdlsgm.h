/**
 * @file fdlsgm.h
 * @brief Header file for the FDLSGM library
 * @author Ryou Ohsawa
 * @year 2020
 */
#ifndef __FDLSGM_H_INCLUDE
#define __FDLSGM_H_INCLUDE

#define _USE_MATH_DEFINES
#include<cstdio>
#include<cstdint>
#include<cstdlib>
#include<cmath>
#include<limits>
#include<array>
#include<set>
#include<map>
#include<unordered_map>
#include<algorithm>
#include<utility>
#include<exception>

#include "linalg.h"


namespace fdlsgm {
  using linalg::index_t;
  using linalg::vector2;
  using linalg::vector3;
  using linalg::vector4;
  using linalg::matrix3x3;
  using linalg::matrix4x4;
  using linalg::zeros4x4;
  using linalg::ident4x4;
  using linalg::normalize_vector;
  using linalg::outer_product;
  using linalg::eigenvector_jacobi_4x4;

  class dls;
  class baseline;

  /** Line Segment. */
  template <class T> using segment = std::array<vector3<T>, 2>;

  /** */
  typedef std::pair<size_t,dls> ndls;

  /** */
  typedef std::map<size_t,dls> dlspool;


  /**
   * @brief Directed Line Segment
   */
  class dls {
  public:
    /**
     * @brief Constuct a DLS with a pair of three-dimensional vectors.
     * @param[in] __v1: The first vertex of the DLS.
     * @param[in] __v2: The second vertex of the DLS.
     * @note The z-coordinate of the first vertex should be smaller than that
     *       of the second vertex.
     */
    dls(const vector3<double>&, const vector3<double>&);

    /**
     * @brief Constuct a DLS with six coordinate values.
     * @note Two Z-values should not be the same to define the direction.
     */
    dls(const double, const double, const double,
        const double, const double, const double);

    double x0() const; /** x-coordinate of the first vertex */
    double y0() const; /** y-coordinate of the first vertex */
    double z0() const; /** z-coordinate of the first vertex */
    double x1() const; /** x-coordinate of the second vertex */
    double y1() const; /** y-coordinate of the second vertex */
    double z1() const; /** z-coordinate of the second vertex */
    double cx() const; /** x-coordinate of the DLS center */
    double cy() const; /** y-coordinate of the DLS center */
    double cz() const; /** z-coordinate of the DLS center */
    double dx() const; /** difference in x-coordinates */
    double dy() const; /** difference in y-coordinates */
    double dz() const; /** difference in z-coordinates */

    double ex() const; /** x-element of the unit vector */
    double ey() const; /** y-element of the unit vector */
    double ez() const; /** z-element of the unit vector */

    const segment<double> vertices() const;

    /** unit vector */
    const vector3<double> unit_vector() const;

    /**
     * @brief Position Angle of the DLS in radian.
     * @note The y-axis is the origin and the anti-clockwise is positive.
     */
    double pa() const;

    /**
     * @brief Length of the DLS projected onto the XY-plane.
     */
    double radius() const;

    /**
     * @brief Length of the DLS.
     */
    double length() const;

    /** Dot product of the two DLSs. */
    double dot(const dls&) const;
    double dot(const ndls&) const;
    double dot(const baseline&) const;

    /**
     * @brief Angle between two DLSs in radius
     */
    double argument(const dls&) const;
    double argument(const ndls&) const;
    double argument(const baseline&) const;

    /** Print function for debugging. */
    void dprint() const;

  private:
    const double _x0, _y0, _z0; /** first vertex  */
    const double _x1, _y1, _z1; /** second vertex */
    const double _r;            /** radius */
    const double _l;            /** length */
    const double _pa;           /** position angle */
  };


  class baseline {
  public:
    /**
     * @brief Constuct a baseline.
     * @param[in] dls: a pair of (size_t, DLS)
     */
    baseline();
    baseline(const ndls&);
    baseline(const baseline&);

    bool append(const ndls&);

    double x0() const; /** x-coordinate of the first vertex */
    double y0() const; /** y-coordinate of the first vertex */
    double z0() const; /** z-coordinate of the first vertex */
    double x1() const; /** x-coordinate of the second vertex */
    double y1() const; /** y-coordinate of the second vertex */
    double z1() const; /** z-coordinate of the second vertex */
    double cx() const; /** x-coordinate of the baseline center */
    double cy() const; /** y-coordinate of the baseline center */
    double cz() const; /** z-coordinate of the baseline center */
    double dx() const; /** difference in x-coordinates */
    double dy() const; /** difference in y-coordinates */
    double dz() const; /** difference in z-coordinates */

    double ex() const; /** x-element of the unit vector */
    double ey() const; /** y-element of the unit vector */
    double ez() const; /** z-element of the unit vector */

    const segment<double> vertices() const;

    /** unit vector */
    const vector3<double> unit_vector() const;

    /**
     * @brief Position Angle of the DLS.
     * The y-axis is the origin and the anti-clockwise is positive.
     */
    double pa() const;

    /**
     * @brief Velocity of the DLS projected onto the XY-plane.
     */
    double radius() const;

    /** total length */
    double length() const;

    size_t size() const;

    double dot(const dls&) const;
    double dot(const ndls&) const;
    double dot(const baseline&) const;

    double argument(const dls&) const;
    double argument(const ndls&) const;
    double argument(const baseline&) const;

    double point_distance(const vector3<double>&) const;

    double lateral_distance(const dls&) const;
    double lateral_distance(const ndls&) const;
    double lateral_distance(const baseline&) const;

    double gap_length(const dls&) const;
    double gap_length(const baseline&) const;

    /** debug function */
    void dprint() const;

    friend const baseline
    merge_baseline(const baseline&, const baseline&);
  private:
    std::set<size_t> _elements;
    double _x0,_y0,_z0,_x1,_y1,_z1;
    double _pa, _r, _l, _ncx, _ncy, _ncz;
    matrix4x4<double> _f;

    double root_position(const vector3<double>&) const;
    double overlap_length(const dls&) const;

    void update_matrix(const dls&);
    void update_direction();
    void update_parameters();
  };


  class accumulator {
  public:
  private:
    size_t _nsize;
    std::multimap<size_t, baseline> _array;
  };

}

#endif // __FDLSGM_H_INCLUDE
