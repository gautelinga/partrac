#ifdef USE_DOLFIN
#ifndef __TRIANGLEFREQINTERPOL_HPP
#define __TRIANGLEFREQINTERPOL_HPP

#include "MeshInterpol.hpp"
#include "strings.hpp"
#include "FreqStamps.hpp"
#include "Triangle.hpp"
#include "cell_locate.hpp"

class TriangleFreqInterpol final
  : public MeshInterpol<Triangle>
{
public:
  TriangleFreqInterpol(const std::string& infilename);
  ~TriangleFreqInterpol() { std::cout << "Destructing TriangleFreqInterpol." << std::endl; };
  void update(const double t);
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  double get_t_min() { return dolfin_params.get<double>("t_min"); };
  double get_t_max() { return dolfin_params.get<double>("t_max"); };
protected:
  FreqStamps fs; // frequencies holder
  double alpha_t;

  //Vector3d x_min = {0., 0., 0.};
  //Vector3d x_max = {0., 0., 0.};

  Vector3d U = {0., 0., 0.};  // FIXME: 2d
  //double Uy = 0.;
  //double Uz = 0.;
  Vector3d A = {0., 0., 0.};
  double P = 0.;
  Matrix3d gradU, gradA;

  bool inside;

  //std::vector<double> Nu_, Nux_, Nuy_;
  //std::vector<double> Np_;

  std::shared_ptr<dolfin::Function> u_;
  std::shared_ptr<dolfin::Function> p_;

  /*
  std::vector<std::shared_ptr<dolfin::Function>> u__;
  std::vector<std::shared_ptr<dolfin::Function>> p__;
  std::vector<std::vector<double>> u_coefficients__;
  std::vector<std::vector<double>> p_coefficients__;
  */ 
 
  std::vector<std::vector<std::vector<double>>> u_coefficients_;
  std::vector<std::vector<std::vector<double>>> p_coefficients_;

  double omega0 = 0.;

};

#endif
#endif
