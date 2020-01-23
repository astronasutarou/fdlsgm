/**
 * @file fdlsgm.h
 * @brief header file for the FDLSGM library
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
#include<list>
#include<map>
#include<unordered_map>
#include<algorithm>
#include<utility>
#include<exception>

#include "linalg.h"


namespace fdlsgm {
  constexpr bool DEBUG_MESSAGE = false;

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

  /** directed line segment. */
  class dls;
  /** baseline. */
  class baseline;
  /** handler of DLS and baseline, N is the size of internal array. */
  template <index_t N> class accumulator;


  /**
   * @brief view of the directed line segment object.
   *
   * This struct contains two three-dimensional vectors to denote the
   * line segment. The first member `v0` is the starting point of the
   * segment. The second member `v1` is the other end of the segment.
   * Both members are given in double[3] type.
   */
  typedef struct {
    double v0[3];
    double v1[3];
  } dls_view;

  /**
   * @brief view of the baseline object.
   *
   * This struct contains two three-dimensional vectors to denote the
   * baseline. The first member `v0` is the starting point of the
   * segment. The second member `v1` is the other end of the segment.
   * Both members are given in double[3] type. The member `size` is the
   * number of the DLSs to construct the baseline. The id-numbers of the
   * DLSs are stored in `elements`.
   */
  typedef struct {
    double v0[3];
    double v1[3];
    size_t size;
    size_t* elements;
  } baseline_view;


  /**
   * @brief parameter set for find segments.
   *
   * These parameters used in accumulator::insert, accumulator::reallocate,
   * and accumulator::coalesce, accordingly in find_segments. The meaning
   * of each parameter is as follows:
   *
   *   - argment_limit_base: tolerance in argument angle in radian.
   *   - argment_limit_element: tolerance in argument angle in radian.
   *   - distance_limit: tolerance in lateral distance.
   *   - gap_limit: tolerance in gap length [0,1].
   *   - size_limit: threshold of the member size to ignore baselines.
   *
   * @note The parameter `size_limit` is used neither in accumulator::insert
   *       nor accumulator::reallocate.
   */
  typedef struct {
    double argument_limit_base;
    double argument_limit_element;
    double distance_limit;
    double gap_limit;
    size_t size_limit;
  } parameter;

  /** default parameter set for accumulator::insert. */
  constexpr parameter default_param_insert = {
    10.0*M_PI/180.0, 20.0*M_PI/180.0, 3.0, 0.5, 0 };
  /** default parameter set for accumulator::reallocate. */
  constexpr parameter default_param_reallocate = {
    10.0*M_PI/180.0, 20.0*M_PI/180.0, 3.0, 0.5, 0 };
  /** default parameter set for accumulator::coalesce. */
  constexpr parameter default_param_coalesce = {
    10.0*M_PI/180.0, 10.0*M_PI/180.0, 1.0, 0.5, 3 };

  /**
   * @brief find baseline from a set of directed line segments.
   * @param[in] pool: a set of directed line segments.
   * @param[in] param_insert: parameters for accumulator::insert.
   * @param[in] param_reallocate: parameters for accumulator::reallocate.
   * @param[in] param_coalesce: parameters for accumulator::coalesce.
   * @return a set of baselines found in the pool.
   */
  const std::vector<baseline_view>
  find_segments(const std::vector<dls>& pool,
                const size_t& size_threshold = 6,
                const parameter& param_insert = default_param_insert,
                const parameter& param_reallocate = default_param_reallocate,
                const parameter& param_coalesce = default_param_coalesce);
  /**
   * @brief find baseline from a set of directed line segments.
   * @param[in] n_elements a number of line segments in the pool.
   * @param[in] pool: a set of directed line segments.
   * @param[in] param_insert: parameters for accumulator::insert.
   * @param[in] param_reallocate: parameters for accumulator::reallocate.
   * @param[in] param_coalesce: parameters for accumulator::coalesce.
   * @return a set of baselines found in the pool.
   * @note The array pool should contain 6N double elements, where N is the
           number of the directed line segments.
   */
  const std::vector<baseline_view>
  find_segments(const size_t& n_elements,
                const double* pool,
                const size_t& size_threshold = 6,
                const parameter& param_insert = default_param_insert,
                const parameter& param_reallocate = default_param_reallocate,
                const parameter& param_coalesce = default_param_coalesce);

  /** line segment defined by two vertices. */
  template <class T> using segment = std::array<vector3<T>, 2>;

  /**
   * @brief assert a comparable value within [lo,hi]
   * @param[in] v: a number to be constrained.
   * @param[in] lo: a lower bound value.
   * @param[in] hi: an upper bound value.
   * @return a number constrained in [lo,hi].
   */
  template <typename T> inline T
  clamp(const T& v, const T& lo, const T& hi)
  {
    return ((v>lo)?(v<hi)?v:hi:lo);
  }

  class dls {
  public:
    /**
     * @brief constuct a DLS with [0,0,0] to [0,0,1].
     */
    dls();

    /**
     * @brief constuct a DLS with a pair of three-dimensional vectors.
     * @param[in] v1: the first vertex of the DLS.
     * @param[in] v2: the second vertex of the DLS.
     * @note the direction is defined from v1 to v2.
     */
    dls(const vector3<double>& v1, const vector3<double>& v2);

    /**
     * @brief constuct a DLS with six coordinate values.
     * @param[in] x1: x-coordinate of the first vertex.
     * @param[in] y1: y-coordinate of the first vertex.
     * @param[in] z1: z-coordinate of the first vertex.
     * @param[in] x2: x-coordinate of the second vertex.
     * @param[in] y2: y-coordinate of the second vertex.
     * @param[in] z2: z-coordinate of the second vertex.
     */
    dls(const double& x1, const double& y1, const double& z1,
        const double& x2, const double& y2, const double& z2);
    /**
     * @brief constuct a DLS using C-array.
     * @param[in] v: an array containing coordinates of vertecies.
     */
    dls(const double* v);

    /**
     * @brief return a view of the DLS.
     * @return an instance of the struct "dls_view".
     */
    const dls_view view() const;

    double x0() const; /** x-coordinate of the first vertex. */
    double y0() const; /** y-coordinate of the first vertex. */
    double z0() const; /** z-coordinate of the first vertex. */
    double x1() const; /** x-coordinate of the second vertex. */
    double y1() const; /** y-coordinate of the second vertex. */
    double z1() const; /** z-coordinate of the second vertex. */
    double cx() const; /** x-coordinate of the DLS center. */
    double cy() const; /** y-coordinate of the DLS center. */
    double cz() const; /** z-coordinate of the DLS center. */
    double dx() const; /** difference along with x-axis. */
    double dy() const; /** difference along with y-axis. */
    double dz() const; /** difference along with z-axis. */

    double ex() const; /** x-element of the unit vector. */
    double ey() const; /** y-element of the unit vector. */
    double ez() const; /** z-element of the unit vector. */

    /**
     * @brief vertices which define the DLS.
     * @return a pair of the vertices as segment.
     */
    const segment<double> vertices() const;

    /**
     * @brief unit vector.
     * @return the direction of the DLS as a vector of unit length.
     */
    const vector3<double> unit_vector() const;

    /**
     * @brief position angle of the DLS in radian.
     * @return position angle in radian.
     * @note the y-axis is the origin and the anti-clockwise is positive.
     */
    double pa() const;

    /**
     * @brief length of the DLS projected onto the XY-plane.
     * @return the projected length.
     */
    double radius() const;

    /**
     * @brief length of the DLS in the XYZ-space.
     * @return the length of the DLS.
     */
    double length() const;

    /**
     * @brief dot product of two DLSs.
     * @param[in] dls: another DLS.
     * @return the dot product of this DLS and another DLS.
     */
    double dot(const dls& dls) const;
    /**
     * @brief dot product of this DLS and a baseline.
     * @param[in] baseline: a baseline.
     * @return the dot product of this DLS and the baseline.
     */
    double dot(const baseline& baseline) const;

    /**
     * @brief angle between two DLSs.
     * @param[in] dls: another DLS.
     * @return the angle between this DLS and another DLS in radian.
     */
    double argument(const dls& dls) const;
    /**
     * @brief angle between this DLS and a baseline.
     * @param[in] baseline: a baseline.
     * @return the angle between this DLS and the baseline in radian.
     */
    double argument(const baseline& baseline) const;

    /** Print function for debugging. */
    void dprint() const;

  private:
    double _x0, _y0, _z0; /** first vertex  */
    double _x1, _y1, _z1; /** second vertex */
    double _r;            /** radius */
    double _l;            /** length */
    double _pa;           /** position angle */
  };


  class baseline {
  public:
    /**
     * @brief constuct a baseline without DLS.
     */
    baseline();
    /**
     * @brief constuct a baseline with a DLS.
     * @param[in] n: the ID-number of the DLS.
     * @param[in] dls: a directed line segment.
     */
    baseline(const index_t& n, const dls& dls);
    /**
     * @brief constuct a copy of a baseline.
     * @param[in] baseline: an original baseline.
     */
    baseline(const baseline& baseline);

    /**
     * @brief return a view of the baseline.
     * @return an instance of the struct "baseline_view".
     */
    const baseline_view view() const;

    /**
     * @brief append a DLS into the baseline.
     * @param[in] n: the ID-number of the DLS.
     * @param[in] dls: a directed line segment.
     * @return true if the DLS is newly registered in the baseline.
     *         false if the DLS is already registered.
     */
    bool append(const index_t& n, const dls& dls);

    double x0() const; /** x-coordinate of the first vertex. */
    double y0() const; /** y-coordinate of the first vertex. */
    double z0() const; /** z-coordinate of the first vertex. */
    double x1() const; /** x-coordinate of the second vertex. */
    double y1() const; /** y-coordinate of the second vertex. */
    double z1() const; /** z-coordinate of the second vertex. */
    double cx() const; /** x-coordinate of the baseline center. */
    double cy() const; /** y-coordinate of the baseline center. */
    double cz() const; /** z-coordinate of the baseline center. */
    double dx() const; /** difference along with x-axis. */
    double dy() const; /** difference along with y-axis. */
    double dz() const; /** difference along with z-axis. */

    double ex() const; /** x-element of the unit vector. */
    double ey() const; /** y-element of the unit vector. */
    double ez() const; /** z-element of the unit vector. */

    /**
     * @brief vertices which define the baseline.
     * @return a pair of the vertices as segment.
     */
    const segment<double> vertices() const;

    /**
     * @brief unit vector.
     * @return the direction of the baseline as a vector of unit length.
     */
    const vector3<double> unit_vector() const;

    /**
     * @brief position angle of the baseline in radian.
     * @return position angle in radian.
     * @note the y-axis is the origin and the anti-clockwise is positive.
     */
    double pa() const;

    /**
     * @brief length of the baseline projected onto the XY-plane.
     * @return the projected length.
     */
    double radius() const;

    /**
     * @brief length of the baseline in the XYZ-space.
     * @return the length of the line.
     */
    double length() const;

    /**
     * @brief number of the DLSs registered in the baseline.
     * @return the number of the DLSs.
     */
    size_t size() const;

    /**
     * @brief obtain a list of DLS ID-numbers.
     * @return a list of the ID-numbsers of the DLS members.
     */
    std::list<index_t> elements() const;

    /**
     * @brief dot product of this baseline and a DLS.
     * @param[in] dls: a directed line segment.
     * @return the dot product of this baseline and a DLS.
     */
    double dot(const dls& dls) const;
    /**
     * @brief dot product of two baselines.
     * @param[in] baseline: another baseline.
     * @return the dot product of this baseline and another baseline.
     */
    double dot(const baseline& baseline) const;

    /**
     * @brief angle between this baseline and a DLS.
     * @param[in] dls: a directed line segment.
     * @return the angle between this baseline and a DLS.
     */
    double argument(const dls& dls) const;
    /**
     * @brief angle between this DLS and a baseline.
     * @param[in] baseline: a baseline.
     * @return the angle between this DLS and the baseline in radian.
     */
    double argument(const baseline& baseline) const;

    /**
     * @brief distance between the baseline and a point.
     * @param[in] v: a vector representing coordinates of a point.
     * @return the distance between this baseline and the given position.
     */
    double point_distance(const vector3<double>& v) const;
    /** squared-distance between the baseline and a point. */
    double point_distance_squared(const vector3<double>& v) const;

    /**
     * @brief lateral distance between the baseline and a DLS.
     * @param[in] dls: a directed line segment.
     * @return the lateral disntace between this baseline and the DLS.
     * @note the lateral distance is defined as the the square-root of
     *       the squared-distance integrated along with the DLS and
     *       devided by the length of the DLS.
     */
    double lateral_distance(const dls& dls) const;
    /** squared-lateral distance between the baseline and a DLS. */
    double lateral_distance_squared(const dls& dls) const;
    /**
     * @brief lateral distance between two baselines.
     * @param[in] baseline: another baseline.
     * @return the lateral disntace between this baseline and another baseline.
     * @note the lateral distance is defined as the the square-root of
     *       the squared-distance integrated along with the baseline and
     *       devided by the length of the baseline.
     */
    double lateral_distance(const baseline& baseline) const;
    /** squared-lateral distance between two baselines. */
    double lateral_distance_squared(const baseline& baseline) const;

    /**
     * @brief gap between this baseline and a DLS.
     * @param[in] dls: a directed line segment.
     * @return a float value (>0); zero if the projection of the DLS onto
     *         this baseline is overlapped with this baseline. Otherwise,
     *         the gap-length in units of the baseline length is returned.
     */
    double gap_length(const dls& dls) const;
    /**
     * @brief gap between two baselines.
     * @param[in] baseline: another baseline.
     * @return a float value (>0); zero if the projection of the baseline
     *         on this baseline is overlapped with this baseline. Otherwise,
     *         the gap-length in units of the baseline length is returned.
     */
    double gap_length(const baseline& baseline) const;

    /**
     * @brief overlap between this baseline and a DLS.
     * @param[in] dls: a directed line segment.
     * @return a float value [0,1]; unity if the projection of the DLS
     *         onto this baseline is completely contained in this baseline.
     *         Otherwise, a fraction of the overlapped is returned.
     */
    double overlap_length(const dls& dls) const;
    /**
     * @brief overlap between this baseline and another baseline.
     * @param[in] baseline: another baseline.
     * @return a float value [0,1]; unity if the projection of the baseline
     *         onto this baseline is completely contained in this baseline.
     *         Otherwise, a fraction of the overlapped is returned.
     */
    double overlap_length(const baseline& baseline) const;

    /** debug function */
    void dprint() const;

    /**
     * @brief merge two baselines.
     * @param[in] first: the first baseline.
     * @parma[in] second: the second baseline.
     * @return merged baseline, which contains all the members contained
     *         in the two baselines.
     */
    friend const baseline
    merge_baseline(const baseline& first, const baseline& second);

  private:
    std::set<index_t> _elements; /** container of DLS. */
    double _x0,_y0,_z0;          /** first vertex */
    double _x1,_y1,_z1;          /** second vertex */
    double _r;                   /** radius */
    double _l;                   /** length */
    double _pa;                  /** position angle */
    double _ncx, _ncy, _ncz;     /** variables to calculate baseline center */
    matrix4x4<double> _f;        /** matrix for updating parameters */

    /** calculate projected root onto the baseline. */
    double root_position(const vector3<double>& v) const;

    /** update internal matrix with appended DLS. */
    void update_matrix(const dls& dls);
    /** update baseline parameters using internal matrix */
    void update_direction();
    /** update _r, _l, and _pa */
    void update_parameters();
  };


  template <index_t N>
  class accumulator {
  public:
    /**
     * @brief construct an empty accumulator.
     */
    accumulator();
    /**
     * @brief construct an empty accumulator.
     * @parma[in] reserve: reserved size of internal vector.
     */
    accumulator(const size_t& reserve);

    /**
     * @brief number of directed line segments.
     * @return the number of the inserted directed line segments.
     */
    const size_t count_segment() const;
    /**
     * @brief number of baselines.
     * @return the current number of the generated baselines.
     */
    const size_t count_baseline() const;

    /**
     * @brief insert a directed line segment.
     * @param[in] dls: a directed line segment.
     * @param[in] param: parameters for baseline matching.
     *   - argment_limit_base: tolerance in argument angle in radian.
     *   - argment_limit_element: tolerance in argument angle in radian.
     *   - distance_limit: tolerance in lateral distance.
     *   - gap_limit: tolerance in gap length [0,1].
     *   - size_limit: not used.
     */
    void insert(const dls& dls,
                const parameter& param = default_param_insert);

    /**
     * @brief assign registered DLSs into baselines until converged.
     * @param[in] param: parameters for baseline matching.
     *   - argment_limit_base: tolerance in argument angle in radian.
     *   - argment_limit_element: tolerance in argument angle in radian.
     *   - distance_limit: tolerance in lateral distance.
     *   - gap_limit: tolerance in gap length [0,1].
     *   - size_limit: not used.
     */
    void reallocate(const parameter& param = default_param_reallocate);

    /**
     * @brief merge similar baselines
     * @param[in] param: parameters for baseline matching.
     *   - argment_limit_base: tolerance in argument angle in radian.
     *   - argment_limit_element: tolerance in argument angle in radian.
     *   - distance_limit: tolerance in lateral distance.
     *   - gap_limit: tolerance in gap length [0,1].
     *   - size_limit: threshold of the member size to ignore baselines.
     */
    void coalesce(const parameter& param = default_param_coalesce);

    /** debug function */
    void dprint(const size_t& limit = 0) const;

    /**
     * @brief obtain reference to the n-th baseline.
     * @param[in] n: sequential number of the baseline.
     * @return a constant reference to the n-th baseline.
     */
    const baseline& operator[](const index_t& n) const;
  private:
    /** an alias to connect two indexes. */
    typedef std::pair<index_t,index_t> connector;
    /** width of internal position angle bins. */
    const double tics = 2.0*M_PI/N;
    /** object to map position angle -> baselines. */
    std::unordered_multimap<index_t, index_t> _connector;

    std::vector<dls> _elements;        /** container of DLSs. */
    std::vector<baseline> _baselines;  /** container of baselines. */

    /** create and push a baseline into internal container. */
    void push_baseline(const index_t& dls_index, const dls& dls);
    /** push a baseline into internal container. */
    void push_baseline(const baseline& baseline);

    /** request a list of baselines around index of position angle */
    std::list<index_t> query(const index_t& pa_index, const index_t& range);
    /** request a list of baselines around position angle */
    std::list<index_t> query(const double& pa, const index_t& range);

    /** drop a mapping specified with pa_index and baseline_index. */
    void drop(const index_t& pa_index, const index_t& baseline_index);

    /** calculate index of position angle */
    const index_t index(const double& pa) const;
  };

  template<index_t N>
  accumulator<N>::accumulator():
    accumulator(1000)
  {}
  template<index_t N>
  accumulator<N>::accumulator(const size_t& reserve)
  {
    _elements.reserve(reserve);
    _baselines.reserve(reserve);
  }


  template<index_t N>
  const size_t
  accumulator<N>::count_segment() const
  {
    return _elements.size();
  }

  template<index_t N>
  const size_t
  accumulator<N>::count_baseline() const
  {
    return _baselines.size();
  }

  template<index_t N>
  void
  accumulator<N>::insert(const dls& dls,
                         const parameter& param)
  {
    const double& arg_limit_b = param.argument_limit_base;
    const double& arg_limit_e = param.argument_limit_element;
    const double  dist_limit_sq = param.distance_limit*param.distance_limit;
    const double& gap_limit = param.gap_limit;
    const index_t range =
      clamp((index_t)std::ceil((arg_limit_b+arg_limit_e)/tics),0L,N/2);
    const double pa = dls.pa();

    _elements.push_back(dls);
    const index_t idx = _elements.size()-1;
    const auto& baseline_index = query(pa, range);
    bool appended = false;

    for (auto& n: baseline_index) {
      baseline& b = _baselines[n];
      const double arg_limit =
        arg_limit_b + arg_limit_e/std::sqrt(b.length());
      const double d = b.argument(dls);
      if (d < arg_limit) {
        const double g = b.gap_length(dls);
        if (g < gap_limit) {
          const double l = b.lateral_distance_squared(dls);
          if (l < dist_limit_sq) {
            if (DEBUG_MESSAGE) {
              printf("# match[%04lx:%04lx] (%6.4lf,%6.4f,%6.4f)\n",
                     idx,n,d/arg_limit,l/dist_limit_sq,g/gap_limit);
            }
            drop(index(b.pa()), n);
            b.append(idx, dls);
            appended = true;
            _connector.emplace(index(b.pa()), n);
          }
        }
      }
    }
    if (!appended) {
      push_baseline(idx, dls);
    }
  }

  template<index_t N>
  void
  accumulator<N>::reallocate(const parameter& param)
  {
    const index_t n_elements = count_segment();
    const double& arg_limit_b = param.argument_limit_base;
    const double& arg_limit_e = param.argument_limit_element;
    const double  dist_limit_sq = param.distance_limit*param.distance_limit;
    const double& gap_limit = param.gap_limit;
    const index_t range =
      clamp((index_t)std::ceil((arg_limit_b+arg_limit_e)/tics),0L,N/2);

    bool updated = false;
    size_t counter(0);
    while (true) {
      for (index_t idx=0; idx<n_elements; idx++) {
        const auto& dls = _elements[idx];
        const auto& baseline_index = query(dls.pa(), range);

        for (auto& n: baseline_index) {
          baseline& b = _baselines[n];
          const double arg_limit =
            arg_limit_b + arg_limit_e/std::sqrt(b.length());
          const double d = b.argument(dls);
          if (d < arg_limit) {
            const double g = b.gap_length(dls);
            if (g < gap_limit) {
              const double l = b.lateral_distance_squared(dls);
              if (l < dist_limit_sq) {
                drop(index(b.pa()), n);
                updated |= b.append(idx, dls);
                _connector.emplace(index(b.pa()), n);
              }
            }
          }
        }
      }
      if (!updated || counter > 10) break;
      updated = false; counter++;
    }
  }

  template<index_t N>
  void
  accumulator<N>::coalesce(const parameter& param)
  {
    std::set<index_t> done;
    std::vector<baseline> tmp_baseline;
    std::unordered_multimap<index_t, index_t> tmp_connector;
    const index_t n_elements = count_baseline();
    const double& arg_limit_b = param.argument_limit_base;
    const double& arg_limit_e = param.argument_limit_element;
    const double  dist_limit_sq = param.distance_limit*param.distance_limit;
    const double& gap_limit = param.gap_limit;
    const size_t& size_limit = param.size_limit;
    const index_t range =
      clamp((index_t)std::ceil((arg_limit_b+arg_limit_e)/tics),0L,N/2);

    tmp_baseline.reserve(n_elements);
    for (index_t i=0; i<n_elements; i++) {
      if (done.find(i) != done.end()) continue;
      done.insert(i);
      baseline b(_baselines[i]);
      const auto& baseline_index = query(b.pa(), range);
      for (auto& n: baseline_index) {
        const baseline& x = _baselines[n];
        if (i==n || done.find(n) != done.end()) continue;
        if (x.size() < size_limit) { done.insert(n); continue; }
        const double arg_limit =
          arg_limit_b + arg_limit_e/std::sqrt(x.length());
        const double d = b.argument(x);
        if (d < arg_limit) {
          const double g = b.gap_length(x);
          if (g < gap_limit) {
            const double l = b.lateral_distance_squared(x);
            if (l < dist_limit_sq) {
              if (DEBUG_MESSAGE) {
                printf("# merge[%04lx] (%6.4lf,%6.4f,%6.4f)\n",
                       n,d/arg_limit,l/dist_limit_sq,g/gap_limit);
              }
              drop(index(x.pa()), n);
              b = merge_baseline(b,x);
              done.insert(n);
            }
          }
        }
      }
      tmp_baseline.push_back(b);
      const index_t idx = tmp_baseline.size()-1;
      tmp_connector.emplace(b.pa(), idx);
    }
    _baselines = tmp_baseline;
    _connector = tmp_connector;
  }


  template<index_t N>
  const baseline&
  accumulator<N>::operator[](const index_t& n) const
  {
    return _baselines[n];
  }

  template<index_t N>
  void
  accumulator<N>::push_baseline(const index_t& n, const dls& dls)
  {
    baseline bl(n, dls);
    push_baseline(bl);
  }
  template<index_t N>
  void
  accumulator<N>::push_baseline(const baseline& bl)
  {
    const index_t pa_idx = index(bl.pa());
    _baselines.push_back(bl);
    const index_t bl_idx = _baselines.size()-1;
    _connector.emplace(pa_idx, bl_idx);
  }

  template<index_t N>
  std::list<index_t>
  accumulator<N>::query(const index_t& pa_index, const index_t& range)
  {
    std::list<index_t> ret;
    for (index_t i=pa_index-range; i<pa_index+range; i++) {
      auto range = _connector.equal_range(i);
      for_each(range.first, range.second,
               [&ret](const connector& x){ ret.push_back(x.second); });
    }
    return ret;
  }
  template<index_t N>
  std::list<index_t>
  accumulator<N>::query(const double& pa, const index_t& range)
  {
    return query(index(pa), range);
  }

  template<index_t N>
  void
  accumulator<N>::drop(const index_t& pa_index, const index_t& baseline_index)
  {
    auto range = _connector.equal_range(pa_index);
    for (auto& it=range.first; it!=range.second; it++) {
      if (it->second == baseline_index) {
        _connector.erase(it);
        break;
      }
    }
  }

  template<index_t N>
  const index_t
  accumulator<N>::index(const double& argument) const
  {
    return ((index_t)std::floor(argument/tics)+N)%N;
  }

  template<index_t N>
  void
  accumulator<N>::dprint(const size_t& limit) const
  {
    printf("## segments\n");
    for (auto& e: _elements) e.dprint();
    printf("## N(segment)  : %ld\n\n\n", count_segment());
    printf("## baselines\n");
    size_t cnt(0);
    size_t z = (size_t)&_baselines[0];
    for (auto& b: _baselines) {
      if (b.size()>=limit) {
        b.dprint(); cnt++;
        if (DEBUG_MESSAGE) {
          for (auto& x: _baselines) {
            if (x.size()>=limit) {
              const double d = b.argument(x);
              const double l = b.lateral_distance(x);
              const double g = b.gap_length(x);
              const double w = b.overlap_length(x);
              printf("#\t[%04lx]<->[%04lx] "
                     "(d,l,g,w) = (%12.6lf, %12.6lf, %12.6lf %12.6f)\n",
                     (size_t)&b-z,(size_t)&x-z,d,l,g,w);
            }
          }
        }
      }
    }
    printf("## N(baseline): %ld\n", count_baseline());
    printf("## N(extract) : %ld\n\n", cnt);
  }
}

#endif // __FDLSGM_H_INCLUDE
