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
#include<queue>
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
  typedef std::pair<index_t,index_t> connector;

  /**
   * @brief Directed Line Segment
   */
  class dls {
  public:
    /**
     * @brief Constuct a DLS with a pair of three-dimensional vectors.
     * @param[in] v1: The first vertex of the DLS.
     * @param[in] v2: The second vertex of the DLS.
     * @note The z-coordinate of the first vertex should be smaller than that
     *       of the second vertex.
     */
    dls(const vector3<double>& v1, const vector3<double>& v2);

    /**
     * @brief Constuct a DLS with six coordinate values.
     * @note Two Z-values should not be the same to define the direction.
     */
    dls(const double x1, const double y1, const double z1,
        const double x2, const double y2, const double z2);

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
    double dot(const dls& dls) const;
    double dot(const baseline& baseline) const;

    /**
     * @brief Angle between two DLSs in radius
     */
    double argument(const dls& dls) const;
    double argument(const baseline& baseline) const;

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
     * @param[in] dls: a pair of (index_t, DLS)
     */
    baseline();
    baseline(const index_t& n, const dls& dls);
    baseline(const baseline& baseline);

    bool append(const index_t& n, const dls& dls);

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

    double dot(const dls& dls) const;
    double dot(const baseline& baseline) const;

    double argument(const dls& dls) const;
    double argument(const baseline& baseline) const;

    double point_distance(const vector3<double>& v) const;

    double lateral_distance(const dls& dls) const;
    double lateral_distance(const baseline& baseline) const;

    double gap_length(const dls& dls) const;
    double gap_length(const baseline& baseline) const;

    /** debug function */
    void dprint() const;

    friend const baseline
    merge_baseline(const baseline& first, const baseline& second);
  private:
    std::set<index_t> _elements;
    double _x0,_y0,_z0,_x1,_y1,_z1;
    double _pa, _r, _l, _ncx, _ncy, _ncz;
    matrix4x4<double> _f;

    double root_position(const vector3<double>& v) const;
    double overlap_length(const dls& dls) const;

    void update_matrix(const dls& dls);
    void update_direction();
    void update_parameters();
  };


  template <index_t N>
  class accumulator {
  public:

    void insert(const dls& dls, const bool& allow_create = false);
    const baseline& operator[](const index_t& n) const;
  private:
    const double tics = 2.0*M_PI/N;
    std::vector<dls> _segments;
    std::vector<baseline> _baselines;
    std::unordered_multimap<index_t, index_t> _connector;

    baseline& get(const index_t& n) const;

    void push(const dls& dls);
    void push(const baseline& baseline);

    void append(const index_t& baseline_index, const index_t& dls_index);

    std::queue<index_t> pop (const index_t& pa_index, const index_t& range);
    std::queue<index_t> pop (const double& pa, const index_t& range);

    const index_t index(const double& pa) const;
  };

  template<index_t N>
  const baseline&
  accumulator<N>::operator[](const index_t& n) const
  {
    return _baselines[n];
  }

  template<index_t N>
  void
  accumulator<N>::insert(const dls& dls,
                         const bool& allow_create)
  {
    _segments.push_back(dls);
    const index_t idx = _segments.size()-1;
  }

  template<index_t N>
  baseline&
  accumulator<N>::get(const index_t& n) const
  {
    return _baselines[n];
  }

  template<index_t N>
  void
  accumulator<N>::push(const dls& dls)
  {
    baseline bl(dls);
    push(bl);
  }
  template<index_t N>
  void
  accumulator<N>::push(const baseline& bl)
  {
    const index_t pa_idx = index(bl.pa());
    _baselines.push_back(bl);
    const index_t bl_idx = _baselines.size()-1;
    _connector.emplace(pa_idx, bl_idx);
  }

  template<index_t N>
  void
  accumulator<N>::append(const index_t& baseline_index,
                         const index_t& dls_index)
  {
    baseline& bl = get(baseline_index);
    const dls& dls = _segments[dls_index];
    bl.append(dls_index, dls);
    const index_t pa_idx = index(bl.pa());
    _connector.emplace(pa_idx, baseline_index);
  }

  template<index_t N>
  std::queue<index_t>
  accumulator<N>::pop(const index_t& pa_index, const index_t& range)
  {
    std::queue<index_t> ret;
    for (index_t i=pa_index-range; i<pa_index+range; i++) {
      auto iter = _connector.equal_range(i);
      for_each(iter.first, iter.second,
               [&ret](const connector& x){ ret.push(x.second); });
    }
    return ret;
  }
  template<index_t N>
  std::queue<index_t>
  accumulator<N>::pop(const double& pa, const index_t& range)
  {
    return pop(index(pa), range);
  }

  template<index_t N>
  const index_t
  accumulator<N>::index(const double& argument) const
  {
    return ((index_t)std::floor(argument/tics))%360;
  }
}

#endif // __FDLSGM_H_INCLUDE
