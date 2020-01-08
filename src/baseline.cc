/**
 * @file baseline.cc
 * @brief Base Line Class
 * @author Ryou Ohsawa
 * @year 2020
 */
#include "fdlsgm.h"


namespace fdlsgm {

  baseline::baseline(const ndls& ndls):
    _x0(ndls.second.x0()),_y0(ndls.second.y0()),_z0(ndls.second.z0()),
    _x1(ndls.second.x1()),_y1(ndls.second.y1()),_z1(ndls.second.z1()),
    _pa(ndls.second.pa()), _r(ndls.second.radius()), _l(ndls.second.length()),
    _ndx(ndls.second.dx()), _ndy(ndls.second.dy()), _f(zeros4x4)
  {
    _elements.insert(ndls.first);
  }

  bool
  baseline::append(const ndls& ndls)
  {
    auto& n = ndls.first;
    auto& e = ndls.second;
    if (_elements.find(n)!=_elements.end()) {
      return false;
    } else {
      _elements.insert(n);
      if (root_position(e.vertices()[0]) < 0.0) {
        auto& v = e.vertices()[0];
        _x0 = v[0]; _y0 = v[1]; _z0 = v[2];
      }
      if (root_position(e.vertices()[1]) > 1.0) {
        auto& v = e.vertices()[1];
        _x1 = v[0]; _y1 = v[1]; _z1 = v[2];
      }
      _f[0][0] += (e.length()/2.0)*(e.x0()*e.x0()+e.x1()*e.x1());
      _f[1][1] += (e.length()/2.0)*(e.y0()*e.y0()+e.y1()*e.y1());
      _f[2][2] += (e.length()/2.0)*(e.z0()*e.z0()+e.z1()*e.z1());
      _f[3][3] += (e.length());
      _f[0][1] += (e.length()/2.0)*(e.x0()*e.y0()+e.x1()*e.y1());
      _f[1][0] += (e.length()/2.0)*(e.y0()*e.x0()+e.y1()*e.x1());
      _f[0][2] += (e.length()/2.0)*(e.x0()*e.z0()+e.x1()*e.z1());
      _f[2][0] += (e.length()/2.0)*(e.z0()*e.x0()+e.z1()*e.x1());
      _f[0][3] += (e.length()/2.0)*(e.x0()*(-1.0)+e.x1()*(-1.0));
      _f[3][0] += (e.length()/2.0)*((-1.0)*e.x0()+(-1.0)*e.x1());
      _f[1][2] += (e.length()/2.0)*(e.y0()*e.z0()+e.y1()*e.z1());
      _f[2][1] += (e.length()/2.0)*(e.z0()*e.y0()+e.z1()*e.y1());
      _f[1][3] += (e.length()/2.0)*(e.y0()*(-1.0)+e.y1()*(-1.0));
      _f[3][1] += (e.length()/2.0)*((-1.0)*e.y0()+(-1.0)*e.y1());
      _f[2][3] += (e.length()/2.0)*(e.z0()*(-1.0)+e.z1()*(-1.0));
      _f[3][2] += (e.length()/2.0)*((-1.0)*e.z0()+(-1.0)*e.z1());

      printf("\n\n[\n");
      for (size_t i=0; i<4; i++) {
        printf(" [ ");
        for (size_t j=0; j<4; j++) {
          printf("%16.8le", _f[j][i]);
        }
        printf(" ],\n");
      }
      printf("];\n");

      auto& eigen = eigenvector_jacobi_4x4(_f);

      printf("\n[\n");
      for (size_t i=0; i<4; i++) {
        printf(" [ ");
        for (size_t j=0; j<4; j++) {
          printf("%16.8le", eigen[j][i]);
        }
        printf(" ],\n");
      }
      printf("];\n\n");

      _r  = std::sqrt(dx()*dx()+dy()*dy());
      _l  = std::sqrt(dx()*dx()+dy()*dy()+dz()*dz());
      _pa = fdlsgm::wrap_angle(std::atan2(-dx(), dy()));
      _ndx += e.dx(); _ndy += e.dy();
      return true;
    }
  }

  double baseline::x0() const { return _x0; }
  double baseline::y0() const { return _y0; }
  double baseline::z0() const { return _z0; }
  double baseline::x1() const { return _x1; }
  double baseline::y1() const { return _y1; }
  double baseline::z1() const { return _z1; }
  double baseline::dx() const { return _x1-_x0; }
  double baseline::dy() const { return _y1-_y0; }
  double baseline::dz() const { return _z1-_z0; }

  double baseline::ex() const { return dx()/_l; }
  double baseline::ey() const { return dy()/_l; }
  double baseline::ez() const { return dz()/_l; }

  double baseline::pa() const { return _pa; }
  double baseline::radius() const { return _r; }
  double baseline::length() const { return _l; }

  const segment<double>
  baseline::vertices() const
  {
    return { vector3<double>{x0(), y0(), z0()},
             vector3<double>{x1(), y1(), z1()} };
  }

  const vector3<double>
  baseline::unit() const { return { ex(), ey(), ex() }; }

  template<> double
  baseline::dot(const dls& dls) const
  {
    return dx()*dls.dx()+dy()*dls.dy()+dz()*dls.dz();
  }
  template<> double
  baseline::dot(const ndls& ndls) const
  {
    return dx()*ndls.second.dx()+dy()*ndls.second.dy()+dz()*ndls.second.dz();
  }
  template<> double
  baseline::dot(const baseline& bl) const
  {
    return dx()*bl.dx()+dy()*bl.dy()+dz()*bl.dz();
  }


  template<> double
  baseline::argument(const dls& dls) const
  {
    return std::acos(dot(dls)/length()/dls.length());
  }
  template<> double
  baseline::argument(const ndls& ndls) const
  {
    return std::acos(dot(ndls.second)/length()/ndls.second.length());
  }
  template<> double
  baseline::argument(const baseline& bl) const
  {
    return std::acos(dot(bl)/length()/bl.length());
  }



  template<> double
  baseline::projected_angle(const dls& dls) const
  {
    return fdlsgm::angle_separation(pa(), dls.pa());
  }
  template<> double
  baseline::projected_angle(const ndls& ndls) const
  {
    return fdlsgm::angle_separation(pa(), ndls.second.pa());
  }
  template<> double
  baseline::projected_angle(const baseline& bl) const
  {
    return fdlsgm::angle_separation(pa(), bl.pa());
  }

  double
  baseline::point_distance(const vector3<double>& p) const
  {
    const double& Cxx =
      std::pow((p[0]-x0()),2.0)*(1.0-ex()*ex());
    const double& Cyy =
      std::pow((p[1]-y0()),2.0)*(1.0-ey()*ey());
    const double& Czz =
      std::pow((p[2]-z0()),2.0)*(1.0-ez()*ez());
    const double& Cxy =
      -2.0*(p[0]-x0())*(p[1]-y0())*ex()*ey();
    const double& Cyz =
      -2.0*(p[1]-y0())*(p[2]-z0())*ey()*ez();
    const double& Czx =
      -2.0*(p[2]-z0())*(p[0]-x0())*ez()*ex();
    return std::sqrt(Cxx+Cyy+Czz+Cxy+Cyz+Czx);
  }

  template<> double
  baseline::lateral_distance(const dls& dls) const
  {
    const double& Cxx =
      ( std::pow(dls.x1()-x0(),3.0) - std::pow(dls.x0()-x0(),3.0) )
      *(1.0-ex()*ex())/(3.0*dls.dx());
    const double& Cyy =
      ( std::pow(dls.y1()-y0(),3.0) - std::pow(dls.y0()-y0(),3.0) )
      *(1.0-ey()*ey())/(3.0*dls.dy());
    const double& Czz =
      ( std::pow(dls.z1()-z0(),3.0) - std::pow(dls.z0()-z0(),3.0) )
      *(1.0-ez()*ez())/(3.0*dls.dz());
    const double& Cxy =
      - ( 2.0*dls.dx()*dls.dy()
          + 3.0*(dls.dy()*(dls.x0()-x0())+dls.dx()*(dls.y0()-y0()))
          + 6.0*(dls.x0()-x0())*(dls.y0()-y0()) )*ex()*ey()/3.0;
    const double& Cyz =
      - ( 2.0*dls.dy()*dls.dz()
          + 3.0*(dls.dz()*(dls.y0()-y0())+dls.dy()*(dls.z0()-z0()))
          + 6.0*(dls.y0()-y0())*(dls.z0()-z0()) )*ey()*ez()/3.0;
    const double& Czx =
      - ( 2.0*dls.dz()*dls.dx()
          + 3.0*(dls.dx()*(dls.z0()-z0())+dls.dz()*(dls.x0()-x0()))
          + 6.0*(dls.z0()-z0())*(dls.x0()-x0()) )*ez()*ex()/3.0;
    return std::sqrt((Cxx+Cyy+Czz+Cxy+Cyz+Czx));
  }
  template<> double
  baseline::lateral_distance(const ndls& ndls) const
  {
    return lateral_distance(ndls.second);
  }


  void
  baseline::dprint() const
  {
    printf("# baseline: [%08lx]\n", (size_t)this);
    printf("#\tVertex 0      : (%lf %lf %lf)\n", x0(),y0(),z0());
    printf("#\tVertex 1      : (%lf %lf %lf)\n", x1(),y1(),z1());
    printf("#\tPosition Angle: %lf\n", pa()*180.0/M_PI);
    printf("#\tLength        : %lf\n", length());
    printf("#\tRadius        : %lf\n", radius());
    printf("#\tVertices      : [ ");
    for (auto& n: _elements) printf("%ld ", n);  printf("]\n");
  }


  double
  baseline::root_position(const vector3<double>& v) const
  {
    return ((v[0]-x0())*ex()+(v[1]-y0())*ey()+(v[2]-z0())*ez())/length();
  }

  double
  baseline::overlap_length(const dls& dls) const
  {
    auto& v0 = dls.vertices()[0];
    auto& v1 = dls.vertices()[1];
    double t0 = root_position(v0), t1 = root_position(v1);
    if (t1<0.0 || t0>1.0) return 0.0;
    if (t0<0.0 && t1>1.0) return 1.0;
    if (t0<0.0 && t1>0.0) return -t0;
    if (t0>0.0 && t1<1.0) return t1-t0;
    if (t0>0.0 && t1>1.0) return 1.0-t0;
    throw std::invalid_argument("something wrong");
  }

  double
  baseline::gap_length(const dls& dls) const
  {
    auto& v0 = dls.vertices()[0];
    auto& v1 = dls.vertices()[1];
    double t0 = root_position(v0), t1 = root_position(v1);
    return (t1<0.0)?-t1:(t0>1.0)?t0-1.0:0.0;
  }

}
