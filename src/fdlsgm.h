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
#include<array>
#include<set>
#include<map>
#include<unordered_map>
#include<algorithm>
#include<exception>


namespace fdlsgm {
  /** Two-Dimensional Vector */
  template <class T> using vector2 = std::array<T, 2>;

  /** Three-Dimensional Vector */
  template <class T> using vector3 = std::array<T, 3>;


  /** Constrain Angle between [-180, +180] */
  double wrap_angle(const double);

  /** Calcurate One-Dimensional Angle Distance */
  double angle_separation(const double, const double);

  /**
   * @brief Directed Line Segment
   */
  class DLS {
  public:
    /**
     * @brief Constuct a DLS with a pair of three-dimensional vectors.
     * @param[in] __v1: The first vertex of the DLS.
     * @param[in] __v2: The second vertex of the DLS.
     * @note The z-coordinate of the first vertex should be smaller than that
     *       of the second vertex.
     */
    DLS(const vector3<double>&, const vector3<double>&);

    /**
     * @brief Constuct a DLS with six coordinate values.
     * @note Two Z-values should not be the same to define the direction.
     */
    DLS(const double, const double, const double,
        const double, const double, const double);

    double x1() const; /** x-coordinate of the first vertex */
    double y1() const; /** y-coordinate of the first vertex */
    double z1() const; /** z-coordinate of the first vertex */
    double x2() const; /** x-coordinate of the second vertex */
    double y2() const; /** y-coordinate of the second vertex */
    double z2() const; /** z-coordinate of the second vertex */
    double dx() const; /** difference in x-coordinates */
    double dy() const; /** difference in y-coordinates */
    double dz() const; /** difference in z-coordinates */

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
    double dot(const DLS& other) const;

    /**
     * @brief Angle between two DLSs in radius
     */
    double argument(const DLS& other) const;

    /** Print function for debugging. */
    void dprint() const;

  private:
    double _x1, _y1, _z1; /** first vertex  */
    double _x2, _y2, _z2; /** second vertex */
    double _dx, _dy, _dz; /** differences   */
    double _r;            /** radius */
    double _l;            /** length */
    double _pa;           /** position angle */
  };


  /** */
  typedef std::pair<size_t,DLS> nDLS;

  /** */
  typedef std::map<size_t,DLS> DLSpool;


  class BaseLine {
  public:
    /**
     * @brief Constuct a BaseLine.
     * @param[in] dls: a pair of (size_t, DLS)
     */
    BaseLine(const nDLS&);

    void append(const nDLS&);

    /**
     * @brief Position Angle of the DLS.
     * The y-axis is the origin and the anti-clockwise is positive.
     */
    double pa() const;

    /**
     * @brief Velocity of the DLS projected onto the XY-plane.
     */
    double radius() const;

    double angle_separation(const DLS&) const;
    double angle_separation(const nDLS&) const;
    double angle_separation(const BaseLine&) const;

    /** debug function */
    void dprint() const;
  private:
    std::set<size_t> _elements;
    double _x1,_y1,_z1,_x2,_y2,_z2,_dx,_dy,_dz;
    double _pa, _r, _l;

  };

}

#endif // __FDLSGM_H_INCLUDE
