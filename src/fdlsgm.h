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


  /**
   * @brief Constrain Angle between [-180, +180]
   */
  template <class T> T wrap_angle(T);


  /**
   * @brief Directed Line Segment
   */
  template <class T>
  class DLS {
  public:
    /**
     * @brief Constuct a DLS with a pair of three-dimensional vectors.
     * @param[in] __v1: The first vertex of the DLS.
     * @param[in] __v2: The second vertex of the DLS.
     * @note The Z-values should not be the same to define the direction.
     */
    DLS(const vector3<T>& __v1, const vector3<T>& __v2)
      : DLS(__v1[0], __v1[1], __v1[2], __v2[0], __v2[1], __v2[2]) {}
    /**
     * @brief Constuct a DLS with six coordinate values.
     * @note Two Z-values should not be the same to define the direction.
     */
    DLS(const T __x1, const T __y1, const T __z1,
        const T __x2, const T __y2, const T __z2)
    {
      if (__z1 < __z2) {
        _x1 = __x1; _y1 = __y1; _z1 = __z1;
        _x2 = __x2; _y2 = __y2; _z2 = __z2;
      } else if (__z2 < __z1) {
        _x1 = __x2; _y1 = __y2; _z1 = __z2;
        _x2 = __x1; _y2 = __y1; _z2 = __z1;
      } else {
        throw std::invalid_argument
          ("cannot define direction. z-positions should not be the same.");
      }
      _dx = _x2-_x1; _dy = _y2-_y1; _dz = _z2-_z1;
      _radial_distance = std::sqrt(_dx*_dx+_dy*_dy);
      _position_angle  = fdlsgm::wrap_angle(std::atan2(_dx, _dy));
    }

    inline T x1() const { return _x1; }
    inline T y1() const { return _y1; }
    inline T z1() const { return _z1; }
    inline T x2() const { return _x2; }
    inline T y2() const { return _y2; }
    inline T z2() const { return _z2; }
    inline T dx() const { return _dx; }
    inline T dy() const { return _dy; }
    inline T dz() const { return _dz; }

    /**
     * @brief Position Angle of the DLS.
     * The y-axis is the origin and the anti-clockwise is positive.
     */
    inline double
    position_angle() const { return _position_angle; }

    /**
     * @brief Length of the DLS projected onto the XY-plane.
     */
    inline double
    radial_distance() const { return _radial_distance; }

    /**
     * @brief Velocity of the DLS projected onto the XY-plane.
     */
    inline double
    radial_velocity() const { return _radial_distance/_dz; }

    /** Print function for debugging. */
    inline void
    dprint() const {
      printf("%8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf # ",
             _x1,_y1,_z1,_x2,_y2,_z2);
      printf("(r,v,t) = (%.2lf, %.2lf, %.2lf)\n",
             radial_distance(),
             radial_velocity(),
             position_angle()/M_PI*180.0);
    }

  private:
    T _x1, _y1, _z1; /** first vertex  */
    T _x2, _y2, _z2; /** second vertex */
    T _dx, _dy, _dz; /** differences   */
    double _radial_distance;
    double _position_angle;
  };


  template <class T> using annotatedDLS = std::pair<size_t, DLS<T> >;
  template <class T> using DLSpool = std::map<size_t, DLS<T> >;


  template <class T>
  class BaseLine {
  public:
    BaseLine(const annotatedDLS<T>& dls):
      _x1(dls.second.x1()),_y1(dls.second.y1()),_z1(dls.second.z1()),
      _x2(dls.second.x2()),_y2(dls.second.y2()),_z2(dls.second.z2()),
      _dx(dls.second.dx()),_dy(dls.second.dy()),_dz(dls.second.dz()),
      _position_angle(dls.second.position_angle()),
      _radial_velocity(dls.second.radial_velocity())
    {
      _elements.insert(dls.first);
    }

    void append(const annotatedDLS<T>& dls)
    {
      auto& n = dls.first;
      if (_elements.find(n)!=_elements.end()) {
        printf("%ld find\n", n);
      } else {
        _elements.insert(n);
      }
    }

    /**
     * @brief Position Angle of the DLS.
     * The y-axis is the origin and the anti-clockwise is positive.
     */
    inline double
    position_angle() const { return _position_angle; }

    /**
     * @brief Velocity of the DLS projected onto the XY-plane.
     */
    inline double
    radial_velocity() const { return _radial_velocity; }

    double angle_separation(const DLS<T>& dls) const;
    double angle_separation(const BaseLine<T>& bl) const;
  private:
    std::set<size_t> _elements;
    T _x1,_y1,_z1,_x2,_y2,_z2,_dx,_dy,_dz;
    double _position_angle;
    double _radial_velocity;

    inline double
    separation(double d1, double d2) const
    { return std::abs(wrap_angle(d1-d2)); }
  };

}
#endif
