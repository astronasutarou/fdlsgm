/**
 * @file misc.cc
 * @brief Miscellaneous Functions
 * @author Ryou Ohsawa
 * @year 2020
 */
#include "fdlsgm.h"


namespace fdlsgm {

  double
  wrap_angle(double angle)
  {
    if (angle > M_PI || angle < -M_PI)
      angle = angle-std::floor(angle/(M_PI*2.0))*(M_PI*2.0);
    return angle;
  }

  double
  angle_separation(const double d1, const double d2)
  {
    return std::abs(fdlsgm::wrap_angle(d1-d2));
  }

  const vector2<double>
  unit_vector(const vector2<double>& v)
  {
    double n = std::sqrt(v[0]*v[0]+v[1]*v[1]);
    return { v[0]/n, v[1]/n };
  }

  const vector3<double>
  unit_vector(const vector3<double>& v)
  {
    double n = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    return { v[0]/n, v[1]/n, v[2]/n };
  }

  const vector3<double>
  outer_product(const vector3<double>& o1, const vector3<double>& o2)
  {
    const double& x = o1[1]*o2[2] - o1[2]*o2[1];
    const double& y = o1[2]*o2[0] - o1[0]*o2[2];
    const double& z = o1[0]*o2[1] - o1[1]*o2[0];
    return { x, y, z};
  }

}
