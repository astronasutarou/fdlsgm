/**
 * @file misc.cc
 * @brief Miscellaneous Functions
 * @author Ryou Ohsawa
 * @year 2020
 */
#include "fdlsgm.h"

namespace fdlsgm {

  /** Wrapped Angle */
  template<>
  double wrap_angle(double angle) {
    if (angle > M_PI || angle < -M_PI)
      angle = angle-std::floor(angle/(M_PI*2.0))*(M_PI*2.0);
    return angle;
  }

}
