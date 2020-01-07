/**
 * @file baseline.cc
 * @brief Base Line Class
 * @author Ryou Ohsawa
 * @year 2020
 */
#include "fdlsgm.h"

namespace fdlsgm {

  BaseLine::BaseLine(const nDLS& ndls):
    _x1(ndls.second.x1()),_y1(ndls.second.y1()),_z1(ndls.second.z1()),
    _x2(ndls.second.x2()),_y2(ndls.second.y2()),_z2(ndls.second.z2()),
    _dx(ndls.second.dx()),_dy(ndls.second.dy()),_dz(ndls.second.dz()),
    _pa(ndls.second.pa()), _r(ndls.second.radius()), _l(ndls.second.length())
  {
    _elements.insert(ndls.first);
  }

  void
  BaseLine::append(const nDLS& ndls)
  {
    dprint();

    auto& n = ndls.first;
    if (_elements.find(n)!=_elements.end()) {
      printf("%3ld find / ", n);
    } else {
      _elements.insert(n);
      printf("%3ld new item / ", n);
    }

    dprint();
  }

  double
  BaseLine::pa() const
  {
    return _pa;
  }

  double
  BaseLine::radius() const
  {
    return _r;
  }

  double
  BaseLine::angle_separation(const DLS& dls) const
  {
    return fdlsgm::angle_separation(pa(), dls.pa());
  }

  double
  BaseLine::angle_separation(const nDLS& ndls) const
  {
    return fdlsgm::angle_separation(pa(), ndls.second.pa());
  }

  double
  BaseLine::angle_separation(const BaseLine& bl) const
  {
    return fdlsgm::angle_separation(pa(), bl.pa());
  }

  void
  BaseLine::dprint() const
  {
    printf("[ ");
    for (auto& n: _elements)
      printf("%2ld ", n);
    printf("]\n");
  }

}
