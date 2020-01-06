/**
 * @file baseline.cc
 * @brief Base Line Class
 * @author Ryou Ohsawa
 * @year 2020
 */
#include "fdlsgm.h"

namespace fdlsgm {

  template<>
  double
  BaseLine<float>::angle_separation(const DLS<float>& dls) const
  { return separation(position_angle(), dls.position_angle()); }
  template<>
  double
  BaseLine<double>::angle_separation(const DLS<double>& dls) const
  { return separation(position_angle(), dls.position_angle()); }

  template<>
  double
  BaseLine<float>::angle_separation(const BaseLine<float>& bl) const
  { return separation(position_angle(), bl.position_angle()); }
  template<>
  double
  BaseLine<double>::angle_separation(const BaseLine<double>& bl) const
  { return separation(position_angle(), bl.position_angle()); }

}
