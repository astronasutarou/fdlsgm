/**
 * @file fdlsgm.h
 * @brief Header file for the FDLSGM library
 * @author Ryou Ohsawa
 * @year 2020
 */
#ifndef __FDLSGM_H_INCLUDE
#define __FDLSGM_H_INCLUDE
namespace fdlsgm {
  /** Three-Dimensional Vector */
  template <class T> using vector3 = std::array<T, 3>;

  /**
   * @brief Directed Line Segment
   */
  template <class T>
  class dls {
  public:
    dls(vector3<T> __v)
      : dls(__v[0], __v[1], __v[2]) {}
    dls(const T __x, const T __y, const T __z)
      : _x(__x), _y(__y), _z(__z) {
      _r = std::sqrt(_x*_x+_y*_y);
    }
  private:
    T _x, _y, _z;
    double _r, _t;
  };
}
#endif
