// A quadratic divergence-free velocity: the perpendicular gradient of a cubic
// stream function in 2D, the curl of a cubic potential in 3D. The coefficients
// are off every symmetry, so a wrong node order does not pass by luck.
#pragma once

#include "typedefs.hpp"

inline Vector3d poly_u(const Vector3d& p, const int dim){
  const double x = p[0], y = p[1], z = p[2];
  if (dim == 2)
    return {-0.7*x*x + 1.8*x*y + 0.75*y*y - 0.3*x + 1.0*y - 0.2,
            -1.2*x*x + 1.4*x*y - 0.9*y*y - 1.2*x + 0.3*y - 0.8, 0.};
  return {-0.8*x*x + 0.75*y*y - 1.8*z*z + 1.6*x*y + 1.8*x*z + 0.4*x + 0.4,
           1.5*x*x - 1.5*y*y + 1.2*z*z + 1.2*x*y - 0.7*x + 0.4*y + 0.4*z,
          -0.9*y*y - 0.9*z*z + 0.4*x*z + 1.4*y*z + 1.6*x - 0.8*z - 0.5};
}

// grad(i, j) is du_i/dx_j
inline Matrix3d poly_grad(const Vector3d& p, const int dim){
  const double x = p[0], y = p[1], z = p[2];
  Matrix3d g = Matrix3d::Zero();
  if (dim == 2){
    g(0, 0) = -1.4*x + 1.8*y - 0.3;   g(0, 1) = 1.8*x + 1.5*y + 1.0;
    g(1, 0) = -2.4*x + 1.4*y - 1.2;   g(1, 1) = 1.4*x - 1.8*y + 0.3;
    return g;
  }
  g(0, 0) = -1.6*x + 1.6*y + 1.8*z + 0.4;  g(0, 1) = 1.6*x + 1.5*y;       g(0, 2) = 1.8*x - 3.6*z;
  g(1, 0) = 3.0*x + 1.2*y - 0.7;           g(1, 1) = 1.2*x - 3.0*y + 0.4; g(1, 2) = 2.4*z + 0.4;
  g(2, 0) = 0.4*z + 1.6;                   g(2, 1) = -1.8*y + 1.4*z;      g(2, 2) = 0.4*x + 1.4*y - 1.8*z - 0.8;
  return g;
}
