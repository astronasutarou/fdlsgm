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

}
