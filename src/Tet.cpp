#ifdef USE_DOLFIN
#include "Tet.hpp"

Tet::Tet(const dolfin::Cell& cell)
{
  // std::vector<double> coords;
  // cell.get_vertex_coordinates(coords);
  // for (int i=0; i<4; ++i)
  // {
  //   xx_[i] = coords[3*i];
  //   yy_[i] = coords[3*i+1];
  //   zz_[i] = coords[3*i+2];
  // }

  std::array<double, 4> xx, yy, zz;
  for (dolfin::VertexIterator v(cell); !v.end(); ++v)
  {
    const std::size_t pos = v.pos();
    xx[pos] = v->x(0);
    yy[pos] = v->x(1);
    zz[pos] = v->x(2);
  }
  x0_ = xx[0];
  y0_ = yy[0];
  z0_ = zz[0];

  double j11 = xx[1]-xx[0];
  double j12 = yy[1]-yy[0];
  double j13 = zz[1]-zz[0];
  double j21 = xx[2]-xx[0];
  double j22 = yy[2]-yy[0];
  double j23 = zz[2]-zz[0];
  double j31 = xx[3]-xx[0];
  double j32 = yy[3]-yy[0];
  double j33 = zz[3]-zz[0];

  g2x_ = j22*j33-j23*j32;  g3x_ = j13*j32-j12*j33;  g4x_ = j12*j23-j13*j22;
  g2y_ = j23*j31-j21*j33;  g3y_ = j11*j33-j13*j31;  g4y_ = j13*j21-j11*j23;
  g2z_ = j21*j32-j22*j31;  g3z_ = j12*j31-j11*j32;  g4z_ = j11*j22-j12*j21;
  double det = j11 * g2x_ + j12 * g2y_ + j13 * g2z_;
  double d = 1.0/det;
  g2x_ *= d;  g3x_ *= d;  g4x_ *= d;
  g2y_ *= d;  g3y_ *= d;  g4y_ *= d;
  g2z_ *= d;  g3z_ *= d;  g4z_ *= d;
  g1x_ = -g2x_-g3x_-g4x_;  g1y_ = -g2y_-g3y_-g4y_;  g1z_ = -g2z_-g3z_-g4z_;
}

void Tet::xyz2bary(double x, double y, double z,
                   double &r,double &s,double &t,double &u) const
{
  double dx=x-x0_, dy=y-y0_, dz=z-z0_;
  s = g2x_*dx+g2y_*dy+g2z_*dz;
  t = g3x_*dx+g3y_*dy+g3z_*dz;
  u = g4x_*dx+g4y_*dy+g4z_*dz;
  r = 1.-s-t-u;
}

bool Tet::contains(const Vector3d& x) const 
{
  double r1, r2, r3, r4;
  xyz2bary(x[0], x[1], x[2], r1, r2, r3, r4);
  return (r1 >= 0. && r2 >= 0. && r3 >= 0. && r4 >= 0.);
}

void Tet::linearbasis(double r,
                      double s,
                      double t,
                      double u,
                      double *N) const
{
  N[0] = r;
  N[1] = s;
  N[2] = t;
  N[3] = u;
}

void Tet::linearderiv(double ,
                      double ,
                      double ,
                      double ,
                      double *Nx,
                      double *Ny,
                      double *Nz) const
{
  Nx[0] = g1x_;
  Nx[1] = g2x_;
  Nx[2] = g3x_;
  Nx[3] = g4x_;

  Ny[0] = g1y_;
  Ny[1] = g2y_;
  Ny[2] = g3y_;
  Ny[3] = g4y_;

  Nz[0] = g1z_;
  Nz[1] = g2z_;
  Nz[2] = g3z_;
  Nz[3] = g4z_;
}

void Tet::quadbasis(double r,
                    double s,
                    double t,
                    double u,
                    double *N) const
{
  N[0] = r*(2*r-1);
  N[1] = s*(2*s-1);
  N[2] = t*(2*t-1);
  N[3] = u*(2*u-1);
  N[perm_[4]] = 4*r*s;
  N[perm_[5]] = 4*s*t;
  N[perm_[6]] = 4*r*t;
  N[perm_[7]] = 4*r*u;
  N[perm_[8]] = 4*s*u;
  N[perm_[9]] = 4*t*u;
}

void Tet::quadderiv(double r,
                    double s,
                    double t,
                    double u,
                    double *Nx,
                    double *Ny,
                    double *Nz) const
{
  double a = 4.0*r-1.0;
  double b = 4.0*s-1.0;
  double c = 4.0*t-1.0;
  double d = 4.0*u-1.0;

  Nx[0] = a*g1x_;
  Nx[1] = b*g2x_;
  Nx[2] = c*g3x_;
  Nx[3] = d*g4x_;
  Nx[perm_[4]] = 4*(r*g2x_+s*g1x_);
  Nx[perm_[5]] = 4*(s*g3x_+t*g2x_);
  Nx[perm_[6]] = 4*(t*g1x_+r*g3x_);
  Nx[perm_[7]] = 4*(r*g4x_+u*g1x_);
  Nx[perm_[8]] = 4*(s*g4x_+u*g2x_);
  Nx[perm_[9]] = 4*(t*g4x_+u*g3x_);

  Ny[0] = a*g1y_;
  Ny[1] = b*g2y_;
  Ny[2] = c*g3y_;
  Ny[3] = d*g4y_;
  Ny[perm_[4]] = 4*(r*g2y_+s*g1y_);
  Ny[perm_[5]] = 4*(s*g3y_+t*g2y_);
  Ny[perm_[6]] = 4*(t*g1y_+r*g3y_);
  Ny[perm_[7]] = 4*(r*g4y_+u*g1y_);
  Ny[perm_[8]] = 4*(s*g4y_+u*g2y_);
  Ny[perm_[9]] = 4*(t*g4y_+u*g3y_);

  Nz[0] = a*g1z_;
  Nz[1] = b*g2z_;
  Nz[2] = c*g3z_;
  Nz[3] = d*g4z_;
  Nz[perm_[4]] = 4*(r*g2z_+s*g1z_);
  Nz[perm_[5]] = 4*(s*g3z_+t*g2z_);
  Nz[perm_[6]] = 4*(t*g1z_+r*g3z_);
  Nz[perm_[7]] = 4*(r*g4z_+u*g1z_);
  Nz[perm_[8]] = 4*(s*g4z_+u*g2z_);
  Nz[perm_[9]] = 4*(t*g4z_+u*g3z_);
}
#endif
