#ifdef USE_DOLFIN
#include "Tet.hpp"
#include <dolfin.h>

namespace {
// The cell's vertex coordinates, three per vertex, in dolfin's local order
std::array<double, 12> vertex_coords(const dolfin::Cell& cell)
{
  std::array<double, 12> c{};
  for (dolfin::VertexIterator v(cell); !v.end(); ++v)
  {
    const std::size_t pos = v.pos();
    c[3*pos]   = v->x(0);
    c[3*pos+1] = v->x(1);
    c[3*pos+2] = v->x(2);
  }
  return c;
}
}

Tet::Tet(const dolfin::Cell& cell)
{
  const std::array<double, 12> c = vertex_coords(cell);
  *this = Tet(c.data(), c.data() + 3, c.data() + 6, c.data() + 9);
}
#endif
