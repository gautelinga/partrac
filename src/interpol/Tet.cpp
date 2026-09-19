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
#endif
