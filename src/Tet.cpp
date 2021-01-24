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

  for (dolfin::VertexIterator v(cell); !v.end(); ++v)
  {
    const std::size_t pos = v.pos();
    xx_[pos] = v->x(0);
    yy_[pos] = v->x(1);
    zz_[pos] = v->x(2);
  }

  double j11 = xx_[1]-xx_[0];
  double j12 = yy_[1]-yy_[0];
  double j13 = zz_[1]-zz_[0];
  double j21 = xx_[2]-xx_[0];
  double j22 = yy_[2]-yy_[0];
  double j23 = zz_[2]-zz_[0];
  double j31 = xx_[3]-xx_[0];
  double j32 = yy_[3]-yy_[0];
  double j33 = zz_[3]-zz_[0];

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
  double dx=x-xx_[0], dy=y-yy_[0], dz=z-zz_[0];
  s = g2x_*dx+g2y_*dy+g2z_*dz;
  t = g3x_*dx+g3y_*dy+g3z_*dz;
  u = g4x_*dx+g4y_*dy+g4z_*dz;
  r = 1.-s-t-u;
}

void Tet::linearbasis(double r,
                      double s,
                      double t,
                      double u,
                      std::array<double, 4> &N) const
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
                      std::array<double, 4> &Nx,
                      std::array<double, 4> &Ny,
                      std::array<double, 4> &Nz) const
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
                    std::array<double, 10> &N) const
{
  N[0] = r*(2*r-1);
  N[1] = s*(2*s-1);
  N[2] = t*(2*t-1);
  N[3] = u*(2*u-1);
  N[4] = 4*r*s;
  N[5] = 4*s*t;
  N[6] = 4*r*t;
  N[7] = 4*r*u;
  N[8] = 4*s*u;
  N[9] = 4*t*u;
}

void Tet::quadderiv(double r,
                    double s,
                    double t,
                    double u,
                    std::array<double, 10> &Nx,
                    std::array<double, 10> &Ny,
                    std::array<double, 10> &Nz) const
{
  double a = 4.0*r-1.0;
  double b = 4.0*s-1.0;
  double c = 4.0*t-1.0;
  double d = 4.0*u-1.0;

  Nx[0] = a*g1x_;
  Nx[1] = b*g2x_;
  Nx[2] = c*g3x_;
  Nx[3] = d*g4x_;
  Nx[4] = 4*(r*g2x_+s*g1x_);
  Nx[5] = 4*(s*g3x_+t*g2x_);
  Nx[6] = 4*(t*g1x_+r*g3x_);
  Nx[7] = 4*(r*g4x_+u*g1x_);
  Nx[8] = 4*(s*g4x_+u*g2x_);
  Nx[9] = 4*(t*g4x_+u*g3x_);


  Ny[0] = a*g1y_;
  Ny[1] = b*g2y_;
  Ny[2] = c*g3y_;
  Ny[3] = d*g4y_;
  Ny[4] = 4*(r*g2y_+s*g1y_);
  Ny[5] = 4*(s*g3y_+t*g2y_);
  Ny[6] = 4*(t*g1y_+r*g3y_);
  Ny[7] = 4*(r*g4y_+u*g1y_);
  Ny[8] = 4*(s*g4y_+u*g2y_);
  Ny[9] = 4*(t*g4y_+u*g3y_);

  Nz[0] = a*g1z_;
  Nz[1] = b*g2z_;
  Nz[2] = c*g3z_;
  Nz[3] = d*g4z_;
  Nz[4] = 4*(r*g2z_+s*g1z_);
  Nz[5] = 4*(s*g3z_+t*g2z_);
  Nz[6] = 4*(t*g1z_+r*g3z_);
  Nz[7] = 4*(r*g4z_+u*g1z_);
  Nz[8] = 4*(s*g4z_+u*g2z_);
  Nz[9] = 4*(t*g4z_+u*g3z_);
}
