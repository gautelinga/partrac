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

  std::array<double, 3> xx, yy;
  for (dolfin::VertexIterator v(cell); !v.end(); ++v)
  {
    const std::size_t pos = v.pos();
    xx[pos] = v->x(0);
    yy[pos] = v->x(1);
  }
  x0_ = xx[0];
  y0_ = yy[0];

  double j11 = xx[1]-xx[0];
  double j12 = yy[1]-yy[0];
  double j21 = xx[2]-xx[0];
  double j22 = yy[2]-yy[0];

  const double det = j11 * j22 - j12*j21;

  double d = 1.0/det;
  g2x_ = j22*d;   g3x_ = -j12*d;
  g2y_ = -j21*d;  g3y_ = j11*d;
  g1x_ = -g2x_-g3x_;  g1y_ = -g2y_-g3y_;
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
    exit(1);
  }
  return 0.0;
}

#endif
