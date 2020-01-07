/**
 * @file dls.cc
 * @brief Directed Line Segment Class
 * @author Ryou Ohsawa
 * @year 2020
 */
#include "fdlsgm.h"

namespace fdlsgm {


  DLS::DLS(const vector3<double>& __v1, const vector3<double>& __v2)
    : DLS(__v1[0], __v1[1], __v1[2], __v2[0], __v2[1], __v2[2]) {}

  DLS::DLS(const double __x1, const double __y1, const double __z1,
        const double __x2, const double __y2, const double __z2)
    : _x1(__x1),_y1(__y1),_z1(__z1),_x2(__x2),_y2(__y2),_z2(__z2),
      _dx(_x2-_x1),_dy(_y2-_y1),_dz(_z2-_z1)
  {
    _r  = std::sqrt(_dx*_dx+_dy*_dy);
    _l  = std::sqrt(_dx*_dx+_dy*_dy+_dz*_dz);
    _pa = fdlsgm::wrap_angle(std::atan2(_dx, _dy));
  }

  double DLS::x1() const { return _x1; }
  double DLS::y1() const { return _y1; }
  double DLS::z1() const { return _z1; }
  double DLS::x2() const { return _x2; }
  double DLS::y2() const { return _y2; }
  double DLS::z2() const { return _z2; }
  double DLS::dx() const { return _dx; }
  double DLS::dy() const { return _dy; }
  double DLS::dz() const { return _dz; }

  double
  DLS::pa() const
  {
    return _pa;
  }

  double
  DLS::radius() const
  {
    return _r;
  }

  double
  DLS::length() const
  {
    return _l;
  }

  double
  DLS::dot(const DLS& other) const
  {
    return dx()*other.dx()+dy()*other.dy()+dz()*other.dz();
  }

  double
  DLS::argument(const DLS& other) const
  {
    return std::acos(dot(other)/length()/other.length());
  }

  void
  DLS::dprint() const
  {
    printf("%8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf # ",
           _x1,_y1,_z1,_x2,_y2,_z2);
    printf("(r,l,t) = (%.2lf, %.2lf, %.2lf)\n",
           radius(), length(), pa()/M_PI*180.0);
  }

}
