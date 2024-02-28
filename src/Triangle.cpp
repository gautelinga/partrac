#ifdef USE_DOLFIN
#include "Triangle.hpp"

Triangle::Triangle(const dolfin::Cell& cell)
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
  }

  double j11 = xx_[1]-xx_[0];
  double j12 = yy_[1]-yy_[0];
  double j21 = xx_[2]-xx_[0];
  double j22 = yy_[2]-yy_[0];

  det = j11 * j22 - j12*j21;

  double d = 1.0/det;
  g2x_ = j22*d;   g3x_ = -j12*d;
  g2y_ = -j21*d;  g3y_ = j11*d;
  g1x_ = -g2x_-g3x_;  g1y_ = -g2y_-g3y_;
}

bool Triangle::contains(const Vector3d& x) const 
{
  double r1, r2, r3;
  xy2bary(x[0], x[1], r1, r2, r3);
  return (r1 >= 0. && r2 >= 0. && r3 >= 0.);
}

void Triangle::xy2bary(double x, double y,
                       double &r, double &s, double &t) const
{
  double dx=x-xx_[0], dy=y-yy_[0];
  s = g2x_*dx+g2y_*dy;
  t = g3x_*dx+g3y_*dy;
  r = 1.-s-t;
}

void Triangle::linearbasis( double r
                          , double s
                          , double t
                          //, std::array<double, 3> &N
                          , std::vector<double> &N
                          ) const
{
  N[0] = r;
  N[1] = s;
  N[2] = t;
}

std::vector<double> Triangle::dof_coords(const int index) const {
  if (index < 3){
    return {xx_[index], yy_[index]};
  }
  else if (index == perm_[3]){
    return {0.5*(xx_[0] + xx_[1]), 0.5*(yy_[0] + yy_[1])};
  }
  else if (index == perm_[4]){
    return {0.5*(xx_[1] + xx_[2]), 0.5*(yy_[1] + yy_[2])};
  }
  else if (index == perm_[5]){
    return {0.5*(xx_[0] + xx_[2]), 0.5*(yy_[0] + yy_[2])};
  }
  else {
    std::cout << "ERROR: Triangle" << std::endl;
    exit(0);
  }
};

void Triangle::linearderiv( double r
                          , double s
                          , double t
                          //, std::array<double, 3> &Nx
                          //, std::array<double, 3> &Ny
                          , std::vector<double> &Nx
                          , std::vector<double> &Ny
                          ) const
{
  Nx[0] = g1x_;
  Nx[1] = g2x_;
  Nx[2] = g3x_;

  Ny[0] = g1y_;
  Ny[1] = g2y_;
  Ny[2] = g3y_;
}

void Triangle::quadbasis( double r
                        , double s
                        , double t
                        //, std::array<double, 6> &N
                        , std::vector<double> &N
                        ) const
{
  N[0] = r*(2*r-1);
  N[1] = s*(2*s-1);
  N[2] = t*(2*t-1);
  N[perm_[3]] = 4*r*s;
  N[perm_[4]] = 4*s*t;
  N[perm_[5]] = 4*r*t;
}

void Triangle::quadderiv( double r
                        , double s
                        , double t
                        //, std::array<double, 6> &Nx
                        //, std::array<double, 6> &Ny
                        , std::vector<double> &Nx
                        , std::vector<double> &Ny
                        ) const
{
  double a = 4.0*r-1.0;
  double b = 4.0*s-1.0;
  double c = 4.0*t-1.0;

  Nx[0] = a*g1x_;
  Nx[1] = b*g2x_;
  Nx[2] = c*g3x_;
  Nx[perm_[3]] = 4*(r*g2x_+s*g1x_);
  Nx[perm_[4]] = 4*(s*g3x_+t*g2x_);
  Nx[perm_[5]] = 4*(t*g1x_+r*g3x_);

  Ny[0] = a*g1y_;
  Ny[1] = b*g2y_;
  Ny[2] = c*g3y_;
  Ny[perm_[3]] = 4*(r*g2y_+s*g1y_);
  Ny[perm_[4]] = 4*(s*g3y_+t*g2y_);
  Ny[perm_[5]] = 4*(t*g1y_+r*g3y_);
}

void Triangle::dump(){
  for (Uint i=0; i<3; ++i){
    std::cout << xx_[i] << " " << yy_[i] << std::endl;
  }
  std::cout << std::endl;
}

double Triangle::dot_grad_gi(const double vx, const double vy, const int index) const
{
  if (index == 0){
    return vx * g1x_ + vy * g1y_;
  }
  else if (index == 1){
    return vx * g2x_ + vy * g2y_;
  }
  else if (index == 2){
    return vx * g3x_ + vy * g3y_;
  }
  else {
    std::cout << "ERROR: Triangle::dot_grad_gi invalid index." << std::endl;
    exit(0);
  }
  return 0.0;
}

#endif
