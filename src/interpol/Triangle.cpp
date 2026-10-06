#ifdef USE_DOLFIN
#include "Triangle.hpp"
#include <dolfin.h>

namespace {
// The cell's vertex coordinates, two per vertex, in dolfin's local order
std::array<double, 6> vertex_coords(const dolfin::Cell& cell)
{
  std::array<double, 6> c{};
  for (dolfin::VertexIterator v(cell); !v.end(); ++v)
  {
    const std::size_t pos = v.pos();
    c[2*pos]   = v->x(0);
    c[2*pos+1] = v->x(1);
  }
  return c;
}
}

Triangle::Triangle(const dolfin::Cell& cell)
{
  const std::array<double, 6> c = vertex_coords(cell);
  *this = Triangle(c.data(), c.data() + 2, c.data() + 4);
}
#endif
