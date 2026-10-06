#include "openfoam_split.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>

#include "Error.hpp"

namespace openfoam_split {

using openfoam_load::CaseData;

namespace {

using P3 = std::array<double, 3>;

// OpenFOAM's small^2, great and rootVSmall
constexpr double min_tet_quality = 1e-30;
constexpr double great = 1e15;
constexpr double root_vsmall = 1e-150;

P3 sub(const P3& a, const P3& b){ return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
P3 cross(const P3& a, const P3& b){
  return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
}
double dot(const P3& a, const P3& b){ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

// Six times the signed volume of (a, b, c, d)
double vol6(const P3& a, const P3& b, const P3& c, const P3& d){
  return dot(cross(sub(b, a), sub(c, a)), sub(d, a));
}

// tetrahedron::quality: the volume over that of the regular tet of the same circumradius
double quality(const P3& a, const P3& b, const P3& c, const P3& d){
  const P3 A = sub(b, a), B = sub(c, a), C = sub(d, a);
  const P3 ba = cross(B, A), ca = cross(C, A);
  const double lambda = dot(C, C) - dot(A, C);
  const double mu = dot(B, B) - dot(A, B);
  const double denom = dot(C, ba);
  double r;
  if (std::abs(denom) < root_vsmall){
    r = std::sqrt(3.)*great;
  }
  else {
    P3 v;
    for (int i = 0; i < 3; ++i) v[i] = (A[i] + (lambda*ba[i] - mu*ca[i])/denom)/2;
    r = std::sqrt(dot(v, v));
  }
  const double vol = vol6(a, b, c, d)/6.;
  const double rr = std::min(r, great);
  return vol/((8.0/27.0)*std::sqrt(3.)*rr*rr*rr + root_vsmall);
}

P3 point(const CaseData& c, const std::size_t p){
  return {c.points[3*p], c.points[3*p + 1], c.points[3*p + 2]};
}

P3 centre(const CaseData& c, const std::size_t i){
  return {c.cell_centres[3*i], c.cell_centres[3*i + 1], c.cell_centres[3*i + 2]};
}

// A face's points
struct FaceRef {
  const std::int32_t* p;
  std::size_t k;
  std::int32_t operator[](const std::size_t i) const { return p[i]; }
};

FaceRef face(const CaseData& c, const std::size_t f){
  return {c.face_points.data() + c.face_start[f], std::size_t(c.face_start[f + 1] - c.face_start[f])};
}

// polyMeshTetDecomposition::minQuality: the worst tet of the fan from base with apex cc
double min_quality(const CaseData& c, const P3& cc, const FaceRef& f, const bool owner_side,
                   const std::size_t base){
  double q = 1e300;
  const P3 b = point(c, f[base]);
  for (std::size_t i = 1; i < f.k - 1; ++i){
    const std::size_t a0 = (i + base) % f.k, a1 = (a0 + 1) % f.k;
    const std::int32_t pa = owner_side ? f[a0] : f[a1];
    const std::int32_t pb = owner_side ? f[a1] : f[a0];
    q = std::min(q, quality(cc, b, point(c, pa), point(c, pb)));
  }
  return q;
}

// The patch of every face, -1 inside
std::vector<std::int32_t> face_patches(const CaseData& c){
  std::vector<std::int32_t> fp(c.nfaces(), -1);
  for (std::size_t i = 0; i < c.patches.size(); ++i){
    const std::size_t f0 = std::size_t(c.patches[i].start);
    std::fill(fp.begin() + std::ptrdiff_t(f0), fp.begin() + std::ptrdiff_t(f0 + std::size_t(c.patches[i].size)),
              std::int32_t(i));
  }
  return fp;
}

// The faces of each cell, ascending
struct CellFaces {
  std::vector<std::size_t> start;
  std::vector<std::size_t> faces;
};

CellFaces cell_faces(const CaseData& c){
  const std::size_t ni = std::size_t(c.n_internal);
  CellFaces cf;
  cf.start.assign(std::size_t(c.ncells) + 1, 0);
  for (std::size_t f = 0; f < c.nfaces(); ++f){
    ++cf.start[std::size_t(c.owner[f]) + 1];
    if (f < ni) ++cf.start[std::size_t(c.neighbour[f]) + 1];
  }
  std::partial_sum(cf.start.begin(), cf.start.end(), cf.start.begin());
  cf.faces.resize(cf.start.back());
  std::vector<std::size_t> at(cf.start.begin(), cf.start.end() - 1);
  for (std::size_t f = 0; f < c.nfaces(); ++f){
    cf.faces[at[c.owner[f]]++] = f;
    if (f < ni) cf.faces[at[c.neighbour[f]]++] = f;
  }
  return cf;
}

// Union-find root
std::int32_t root(std::vector<std::int32_t>& up, std::int32_t i){
  while (up[i] != i){
    up[i] = up[up[i]];
    i = up[i];
  }
  return i;
}

// Every face's triangles in its own orientation, from its base point
struct FaceTris {
  std::vector<std::size_t> start;                // nfaces + 1, into tris
  std::vector<std::array<std::int32_t, 3>> tris;
};

FaceTris face_triangles(const CaseData& c, const std::vector<std::int32_t>& base){
  FaceTris t;
  t.start.assign(c.nfaces() + 1, 0);
  for (std::size_t f = 0; f < c.nfaces(); ++f) t.start[f + 1] = t.start[f] + face(c, f).k - 2;
  t.tris.resize(t.start.back());
  for (std::size_t f = 0; f < c.nfaces(); ++f){
    const FaceRef F = face(c, f);
    const std::size_t j = F.k == 3 ? 0 : std::size_t(base[f]);
    for (std::size_t i = 1; i < F.k - 1; ++i)
      t.tris[t.start[f] + i - 1] = {F[j], F[(j + i) % F.k], F[(j + i + 1) % F.k]};
  }
  return t;
}

using Tet = std::array<std::uint32_t, 4>;

// A tet from its apex and a face triangle, the triangle as it is or flipped
Tet tet(const std::uint32_t apex, const std::array<std::int32_t, 3>& t, const bool flip){
  return flip ? Tet{apex, std::uint32_t(t[0]), std::uint32_t(t[2]), std::uint32_t(t[1])}
              : Tet{apex, std::uint32_t(t[0]), std::uint32_t(t[1]), std::uint32_t(t[2])};
}

// The mesh points and the centres as one node set, centre of cell i at npoints + i
struct Geometry {
  const CaseData& c;
  P3 x(const std::uint32_t n) const { return n < c.npoints() ? point(c, n) : centre(c, n - c.npoints()); }
  double vol(const Tet& t) const { return vol6(x(t[0]), x(t[1]), x(t[2]), x(t[3]))/6.; }
};

// The centre fan of a cell: its owned faces' triangles, then those it neighbours, flipped
void fan(const CaseData& c, const CellFaces& cf, const FaceTris& ft, const std::size_t cell,
         const std::uint32_t apex, std::vector<Tet>& out){
  for (int side = 0; side < 2; ++side)
    for (std::size_t q = cf.start[cell]; q < cf.start[cell + 1]; ++q){
      const std::size_t f = cf.faces[q];
      const bool own = std::size_t(c.owner[f]) == cell;
      if ((side == 0) != own) continue;
      for (std::size_t k = ft.start[f]; k < ft.start[f + 1]; ++k) out.push_back(tet(apex, ft.tris[k], !own));
    }
}

// The cell's volume under its faces' triangles, from its centre
double fan_volume(const CaseData& c, const CellFaces& cf, const FaceTris& ft, const std::size_t cell){
  const Geometry g{c};
  std::vector<Tet> t;
  fan(c, cf, ft, cell, std::uint32_t(c.npoints() + cell), t);
  double v = 0.;
  for (const Tet& q : t) v += g.vol(q);
  return v;
}

// A hex: six quads on eight points
bool is_hex(const CaseData& c, const CellFaces& cf, const std::size_t cell, std::array<std::int32_t, 8>& pts){
  if (cf.start[cell + 1] - cf.start[cell] != 6) return false;
  std::array<std::int32_t, 24> p;
  std::size_t n = 0;
  for (std::size_t q = cf.start[cell]; q < cf.start[cell + 1]; ++q){
    const FaceRef F = face(c, cf.faces[q]);
    if (F.k != 4) return false;
    for (std::size_t i = 0; i < 4; ++i) p[n++] = F[i];
  }
  std::sort(p.begin(), p.end());
  if (std::unique(p.begin(), p.end()) - p.begin() != 8) return false;
  std::copy(p.begin(), p.begin() + 8, pts.begin());
  return true;
}

// Dompierre's split of a hex from a corner its three faces' diagonals meet
// at, the lowest by (key, id): the cone from it over the far faces'
// triangles, or 5 tets when no far diagonal passes the opposite corner.
// False where no corner takes all three diagonals.
bool dompierre(const CaseData& c, const CellFaces& cf, const FaceTris& ft,
               const std::vector<std::int32_t>& key, const std::size_t cell,
               const std::array<std::int32_t, 8>& P, std::vector<Tet>& out, bool& five){
  std::array<std::size_t, 6> F;
  for (std::size_t i = 0; i < 6; ++i) F[i] = cf.faces[cf.start[cell] + i];
  const auto on_face = [&](const std::size_t i, const std::int32_t p){
    const FaceRef f = face(c, F[i]);
    for (std::size_t j = 0; j < 4; ++j) if (f[j] == p) return true;
    return false;
  };
  const auto on_diagonal = [&](const std::size_t i, const std::int32_t p){
    const auto& t = ft.tris[ft.start[F[i]]];
    return t[0] == p || t[2] == p;
  };
  std::int32_t v0 = -1;
  for (const std::int32_t p : P){
    bool ok = true;
    for (std::size_t i = 0; i < 6 && ok; ++i)
      if (on_face(i, p) && !on_diagonal(i, p)) ok = false;
    if (ok && (v0 < 0 || key[p] < key[v0] || (key[p] == key[v0] && p < v0))) v0 = p;
  }
  if (v0 < 0) return false;
  std::int32_t v6 = -1;
  for (const std::int32_t p : P){
    bool near = false;
    for (std::size_t i = 0; i < 6; ++i) near = near || (on_face(i, v0) && on_face(i, p));
    if (!near) v6 = p;
  }
  std::vector<Tet> cone;
  bool through = false;
  for (std::size_t i = 0; i < 6; ++i){
    if (on_face(i, v0)) continue;
    through = through || on_diagonal(i, v6);
    const std::size_t f = F[i];
    const bool own = std::size_t(c.owner[f]) == cell;
    for (std::size_t k = 0; k < 2; ++k) cone.push_back(tet(std::uint32_t(v0), ft.tris[ft.start[f] + k], !own));
  }
  five = !through;
  if (!five){
    out.insert(out.end(), cone.begin(), cone.end());
    return true;
  }
  const std::uint32_t w6 = std::uint32_t(v6);
  std::vector<std::uint32_t> u;
  for (const Tet& t : cone){
    const bool with6 = t[1] == w6 || t[2] == w6 || t[3] == w6;
    if (!with6){
      out.push_back(t);
      continue;
    }
    for (std::size_t j = 1; j < 4; ++j) if (t[j] != w6) u.push_back(t[j]);
  }
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  if (vol6(point(c, std::size_t(v0)), point(c, u[0]), point(c, u[1]), point(c, u[2])) < 0) std::swap(u[1], u[2]);
  out.push_back({std::uint32_t(v0), u[0], u[1], u[2]});
  out.push_back({w6, u[0], u[2], u[1]});
  return true;
}

// A split under valid_fraction of the cell, or whose tets overlap, which shows
// as a sum past the volume under the faces' triangles
bool hex_valid(const CaseData& c, const CellFaces& cf, const FaceTris& ft, const std::size_t cell,
               const std::vector<Tet>& t){
  const Geometry g{c};
  const double vc = c.cell_volumes[cell];
  double tot = 0.;
  for (const Tet& q : t){
    const double v = g.vol(q);
    if (v <= valid_fraction*vc) return false;
    tot += v;
  }
  const double fv = fan_volume(c, cf, ft, cell);
  return std::abs(tot - fv) <= 1e-9*std::abs(fv);
}

// A boundary face's triangle (3D) or edge (2D), by its sorted points
template<std::size_t NF>
struct BoundaryPiece {
  std::array<std::int32_t, NF> key;
  std::int32_t patch;
};

// The boundary pieces of each cell
template<std::size_t NF>
struct CellPieces {
  std::vector<std::size_t> start;                // ncells + 1, into piece
  std::vector<BoundaryPiece<NF>> piece;
};

// The pieces visit(emit) emits as (cell, piece), twice called, grouped by cell
template<std::size_t NF, typename Visit>
CellPieces<NF> by_cell(const CaseData& c, Visit&& visit){
  CellPieces<NF> out;
  out.start.assign(std::size_t(c.ncells) + 1, 0);
  visit([&](const std::size_t cell, const BoundaryPiece<NF>&){ ++out.start[cell + 1]; });
  std::partial_sum(out.start.begin(), out.start.end(), out.start.begin());
  out.piece.resize(out.start.back());
  std::vector<std::size_t> at(out.start.begin(), out.start.end() - 1);
  visit([&](const std::size_t cell, const BoundaryPiece<NF>& q){ out.piece[at[cell]++] = q; });
  return out;
}

// The patch of every simplex facet on the boundary, each against its own
// cell's boundary pieces, and its partner across a cyclic
template<std::size_t NV>
void boundary_facets(const CellPieces<NV-1>& pieces, const CaseData& c,
                     const std::vector<std::int32_t>& master_point, SplitData& s){
  const std::size_t n = s.nsimplices();
  s.facet_patch.assign(n*NV, -1);
  s.facet_partner.assign(n*NV, -1);
  // (masters, pair, slot) of the facets on cyclics
  std::vector<std::pair<std::array<std::int32_t, NV>, std::int64_t>> cyclic;
  for (std::size_t t = 0; t < n; ++t){
    const std::size_t cell = std::size_t(s.cell_of[t]);
    const auto first = pieces.piece.begin() + std::ptrdiff_t(pieces.start[cell]);
    const auto last = pieces.piece.begin() + std::ptrdiff_t(pieces.start[cell + 1]);
    if (first == last) continue;
    for (std::size_t k = 0; k < NV; ++k){
      std::array<std::int32_t, NV-1> key;
      bool points = true;
      std::size_t m = 0;
      for (std::size_t j = 0; j < NV; ++j){
        if (j == k) continue;
        const std::int32_t p = s.node_point[s.cells[t*NV + j]];
        if (p < 0) points = false;
        key[m++] = p;
      }
      if (!points) continue;
      std::sort(key.begin(), key.end());
      const auto it = std::find_if(first, last, [&](const BoundaryPiece<NV-1>& q){ return q.key == key; });
      if (it == last) continue;
      const std::int32_t patch = it->patch;
      const auto& P = c.patches[std::size_t(patch)];
      s.facet_patch[t*NV + k] = patch;
      if (P.type != "cyclic") continue;
      std::array<std::int32_t, NV> mk{};
      for (std::size_t j = 0; j < NV - 1; ++j) mk[j] = master_point[key[j]];
      std::sort(mk.begin(), mk.end() - 1);
      mk[NV - 1] = std::min(patch, std::int32_t(P.neighbour));
      cyclic.push_back({mk, std::int64_t(t*NV + k)});
    }
  }
  std::sort(cyclic.begin(), cyclic.end());
  std::size_t unpaired = 0;
  for (std::size_t i = 0; i < cyclic.size(); ){
    std::size_t j = i + 1;
    while (j < cyclic.size() && cyclic[j].first == cyclic[i].first) ++j;
    if (j - i == 2){
      s.facet_partner[std::size_t(cyclic[i].second)] = cyclic[i + 1].second;
      s.facet_partner[std::size_t(cyclic[i + 1].second)] = cyclic[i].second;
    }
    else {
      unpaired += j - i;
    }
    i = j;
  }
  if (unpaired)
    partrac::fail(c.dir, ": ", unpaired, " cyclic facets of the split have no image across "
                  "their pair");
}

void split_3d(const CaseData& c, const int tets_per_hex, SplitData& s){
  const std::size_t np = c.npoints(), nc = std::size_t(c.ncells);
  const CellFaces cf = cell_faces(c);
  const std::vector<std::int32_t> base = face_bases(c, &s.base_failures);
  const FaceTris ft = face_triangles(c, base);
  const std::vector<std::int32_t> master = cyclic_masters(c);
  const Geometry g{c};

  // Each cell's split: the fan, or the hex's tets
  std::vector<std::vector<Tet>> hex_tets(nc);
  std::vector<char> fanned(nc, 1);
  std::vector<Tet> t;
  for (std::size_t cell = 0; cell < nc; ++cell){
    std::array<std::int32_t, 8> P;
    const bool hex = is_hex(c, cf, cell, P);
    const double vc = c.cell_volumes[cell];
    bool fan_ok = true;
    if (tets_per_hex == 12 && hex){
      t.clear();
      fan(c, cf, ft, cell, std::uint32_t(np + cell), t);
      for (const Tet& q : t) fan_ok = fan_ok && g.vol(q) > valid_fraction*vc;
    }
    if (!hex || (tets_per_hex == 12 && fan_ok)) continue;
    t.clear();
    bool five = false;
    const bool admitted = dompierre(c, cf, ft, master, cell, P, t, five);
    if (admitted && hex_valid(c, cf, ft, cell, t)){
      fanned[cell] = 0;
      hex_tets[cell] = t;
      if (five) ++s.five_tet_hexes;
      if (tets_per_hex == 12) ++s.fallback;
    }
    else if (tets_per_hex == 6){
      ++s.fallback;
    }
  }

  // Nodes: the points, then the fanned cells' centres in cell order
  std::vector<std::uint32_t> apex(nc, 0);
  std::size_t nn = np;
  for (std::size_t cell = 0; cell < nc; ++cell)
    if (fanned[cell]) apex[cell] = std::uint32_t(nn++);
  s.fan_cells = nn - np;
  s.node_x.assign(c.points.begin(), c.points.end());
  s.node_x.resize(3*nn);
  s.node_kind.assign(nn, 0);
  s.node_point.resize(nn);
  s.node_cell.assign(nn, -1);
  s.node_master.resize(nn);
  for (std::size_t p = 0; p < np; ++p){
    s.node_point[p] = std::int32_t(p);
    s.node_master[p] = std::uint32_t(master[p]);
  }
  for (std::size_t cell = 0; cell < nc; ++cell){
    if (!fanned[cell]) continue;
    const std::size_t n = apex[cell];
    for (std::size_t d = 0; d < 3; ++d) s.node_x[3*n + d] = c.cell_centres[3*cell + d];
    s.node_kind[n] = 1;
    s.node_point[n] = -1;
    s.node_cell[n] = std::int32_t(cell);
    s.node_master[n] = std::uint32_t(n);
  }

  // The tets, by cell, checked
  std::vector<char> bad(nc, 0);
  for (std::size_t cell = 0; cell < nc; ++cell){
    t.clear();
    if (fanned[cell]) fan(c, cf, ft, cell, apex[cell], t);
    else t.swap(hex_tets[cell]);
    const double vc = c.cell_volumes[cell];
    for (const Tet& q : t){
      P3 x[4];
      for (std::size_t j = 0; j < 4; ++j)
        for (std::size_t d = 0; d < 3; ++d) x[j][d] = s.node_x[3*q[j] + d];
      if (vol6(x[0], x[1], x[2], x[3])/6. <= valid_fraction*vc){
        ++s.invalid_simplices;
        bad[cell] = 1;
      }
      s.cells.insert(s.cells.end(), q.begin(), q.end());
      s.cell_of.push_back(std::int32_t(cell));
    }
  }
  s.invalid_cells = std::size_t(std::count(bad.begin(), bad.end(), 1));

  // The boundary faces' triangles, by their points
  const std::vector<std::int32_t> fp = face_patches(c);
  const auto pieces = by_cell<3>(c, [&](auto&& emit){
    for (std::size_t f = std::size_t(c.n_internal); f < c.nfaces(); ++f)
      for (std::size_t k = ft.start[f]; k < ft.start[f + 1]; ++k){
        BoundaryPiece<3> q{ft.tris[k], fp[f]};
        std::sort(q.key.begin(), q.key.end());
        emit(std::size_t(c.owner[f]), q);
      }
  });
  boundary_facets<4>(pieces, c, master, s);
}

// Twice the signed area of (a, b, c) in the plane
double area2(const std::array<double, 2>& a, const std::array<double, 2>& b, const std::array<double, 2>& c){
  return (b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0]);
}

void split_2d(const CaseData& c, const int tets_per_hex, SplitData& s){
  const std::size_t axis = std::size_t(c.empty_axis);
  s.nv = 3;
  s.gdim = 2;
  s.inplane.clear();
  for (int a = 0; a < 3; ++a) if (a != c.empty_axis) s.inplane.push_back(a);
  const std::size_t np = c.npoints(), nc = std::size_t(c.ncells), ni = std::size_t(c.n_internal);
  const std::size_t i0 = std::size_t(s.inplane[0]), i1 = std::size_t(s.inplane[1]);
  const std::vector<std::int32_t> fp = face_patches(c);
  const std::vector<std::int32_t> master = cyclic_masters(c);
  const auto in = [&](const std::size_t p){ return std::array<double, 2>{c.points[3*p + i0], c.points[3*p + i1]}; };

  // The front plane: the empty faces at the smaller coordinate
  double lo = 1e300;
  const double size = openfoam_load::box_diagonal(c);
  std::vector<std::size_t> empty_faces;
  for (std::size_t f = ni; f < c.nfaces(); ++f)
    if (c.patches[std::size_t(fp[f])].type == "empty"){
      empty_faces.push_back(f);
      lo = std::min(lo, c.bface_centres[3*(f - ni) + axis]);
    }
  std::vector<std::int64_t> front_of(nc, -1);
  for (const std::size_t f : empty_faces)
    if (std::abs(c.bface_centres[3*(f - ni) + axis] - lo) <= 1e-9*size)
      front_of[std::size_t(c.owner[f])] = std::int64_t(f);

  std::vector<char> is_front(c.nfaces(), 0);
  for (const std::int64_t f : front_of) if (f >= 0) is_front[std::size_t(f)] = 1;
  const std::vector<std::int32_t> base = face_bases(c, &s.base_failures, &is_front);

  // Each cell's front polygon, counterclockwise, and its triangles
  using Tri = std::array<std::uint32_t, 3>;
  std::vector<std::vector<Tri>> tris(nc);
  std::vector<std::vector<std::uint32_t>> fan_poly(nc);
  std::vector<char> fanned(nc, 0);
  for (std::size_t cell = 0; cell < nc; ++cell){
    if (front_of[cell] < 0) partrac::fail(c.dir, ": cell ", cell, " has no face on the front empty plane");
    const std::size_t f = std::size_t(front_of[cell]);
    const FaceRef F = face(c, f);
    const std::size_t k = F.k;
    std::vector<std::uint32_t> poly(F.p, F.p + k);
    double a = 0.;
    for (std::size_t i = 0; i < k; ++i){
      const auto p = in(poly[i]), q = in(poly[(i + 1) % k]);
      a += p[0]*q[1] - q[0]*p[1];
    }
    std::size_t j = std::size_t(base[f]);
    if (a < 0){
      std::reverse(poly.begin(), poly.end());
      j = k - 1 - j;
    }
    const double area = std::abs(a)/2;
    if (k == 3){
      tris[cell].push_back({poly[0], poly[1], poly[2]});
      continue;
    }
    if (k == 4 && tets_per_hex == 6){
      const Tri t1{poly[j], poly[(j + 1) % 4], poly[(j + 2) % 4]};
      const Tri t2{poly[j], poly[(j + 2) % 4], poly[(j + 3) % 4]};
      const bool ok = area2(in(t1[0]), in(t1[1]), in(t1[2])) > 2*valid_fraction*area
                   && area2(in(t2[0]), in(t2[1]), in(t2[2])) > 2*valid_fraction*area;
      if (ok){
        tris[cell] = {t1, t2};
        continue;
      }
      ++s.fallback;
    }
    fanned[cell] = 1;
    fan_poly[cell] = poly;
  }

  // Nodes: the front points in mesh order, then the fanned cells' centres
  std::vector<std::int64_t> node_of(np, -1);
  for (std::size_t cell = 0; cell < nc; ++cell){
    const FaceRef F = face(c, std::size_t(front_of[cell]));
    for (std::size_t i = 0; i < F.k; ++i) node_of[std::size_t(F[i])] = 0;
  }
  std::size_t nn = 0;
  for (std::size_t p = 0; p < np; ++p){
    if (node_of[p] < 0) continue;
    node_of[p] = std::int64_t(nn++);
    s.node_point.push_back(std::int32_t(p));
    s.node_cell.push_back(-1);
    s.node_kind.push_back(0);
    const auto x = in(p);
    s.node_x.insert(s.node_x.end(), x.begin(), x.end());
  }
  std::vector<std::uint32_t> apex(nc, 0);
  for (std::size_t cell = 0; cell < nc; ++cell){
    if (!fanned[cell]) continue;
    apex[cell] = std::uint32_t(nn++);
    s.node_point.push_back(-1);
    s.node_cell.push_back(std::int32_t(cell));
    s.node_kind.push_back(1);
    s.node_x.push_back(c.cell_centres[3*cell + i0]);
    s.node_x.push_back(c.cell_centres[3*cell + i1]);
  }
  s.fan_cells = std::size_t(std::count(fanned.begin(), fanned.end(), 1));
  s.node_master.resize(nn);
  for (std::size_t n = 0; n < nn; ++n)
    s.node_master[n] = s.node_point[n] >= 0 ? std::uint32_t(node_of[std::size_t(master[std::size_t(s.node_point[n])])])
                                            : std::uint32_t(n);
  const auto node = [&](const std::uint32_t p){ return std::uint32_t(node_of[p]); };

  // The triangles, by cell, checked
  std::vector<char> bad(nc, 0);
  for (std::size_t cell = 0; cell < nc; ++cell){
    std::vector<Tri> t;
    if (fanned[cell]){
      const auto& poly = fan_poly[cell];
      const std::size_t k = poly.size();
      for (std::size_t i = 0; i < k; ++i) t.push_back({apex[cell], node(poly[i]), node(poly[(i + 1) % k])});
    }
    else {
      for (const Tri& q : tris[cell]) t.push_back({node(q[0]), node(q[1]), node(q[2])});
    }
    std::vector<double> a;
    double total = 0.;
    for (const Tri& q : t){
      std::array<double, 2> x[3];
      for (std::size_t j = 0; j < 3; ++j) x[j] = {s.node_x[2*q[j]], s.node_x[2*q[j] + 1]};
      a.push_back(area2(x[0], x[1], x[2])/2);
      total += a.back();
    }
    for (std::size_t i = 0; i < t.size(); ++i){
      if (a[i] <= valid_fraction*total){
        ++s.invalid_simplices;
        bad[cell] = 1;
      }
      s.cells.insert(s.cells.end(), t[i].begin(), t[i].end());
      s.cell_of.push_back(std::int32_t(cell));
    }
  }
  s.invalid_cells = std::size_t(std::count(bad.begin(), bad.end(), 1));

  // The side faces by their two front points
  const auto pieces = by_cell<2>(c, [&](auto&& emit){
    for (std::size_t f = ni; f < c.nfaces(); ++f){
      if (c.patches[std::size_t(fp[f])].type == "empty") continue;
      const FaceRef F = face(c, f);
      BoundaryPiece<2> q{{-1, -1}, fp[f]};
      std::size_t m = 0;
      for (std::size_t i = 0; i < F.k; ++i)
        if (node_of[std::size_t(F[i])] >= 0){
          if (m == 2) partrac::fail(c.dir, ": side face ", f, " has more than two points on the front plane");
          q.key[m++] = F[i];
        }
      if (m != 2) partrac::fail(c.dir, ": side face ", f, " has ", m, " points on the front plane, not two");
      std::sort(q.key.begin(), q.key.end());
      emit(std::size_t(c.owner[f]), q);
    }
  });
  boundary_facets<3>(pieces, c, master, s);
}

}  // namespace

std::vector<std::int32_t> face_bases(const CaseData& c, std::size_t* failures, const std::vector<char>* only){
  const std::size_t nf = c.nfaces(), ni = std::size_t(c.n_internal);
  std::vector<std::int32_t> base(nf, -1);
  std::size_t failed = 0;
  const std::vector<std::int32_t> fp = face_patches(c);
  // The centre across a cyclic face, in this side's frame
  const auto across = [&](const std::size_t f){
    const auto& p = c.patches[std::size_t(fp[f])];
    const auto& q = c.patches[std::size_t(p.neighbour)];
    P3 x = centre(c, std::size_t(c.owner[f - std::size_t(p.start) + std::size_t(q.start)]));
    for (std::size_t d = 0; d < 3; ++d) x[d] -= p.separation[d];
    return x;
  };
#pragma omp parallel for schedule(static) reduction(+:failed)
  for (std::size_t f = 0; f < nf; ++f){
    if (only && !(*only)[f]) continue;
    const FaceRef F = face(c, f);
    const P3 own = centre(c, std::size_t(c.owner[f]));
    bool shared = false;
    P3 nbr{};
    if (f < ni){
      shared = true;
      nbr = centre(c, std::size_t(c.neighbour[f]));
    }
    else if (c.patches[std::size_t(fp[f])].type == "cyclic"){
      if (!c.patches[std::size_t(fp[f])].owner) continue;
      nbr = across(f);
      shared = true;
    }
    double best_q = -1e300;
    std::size_t best = 0;
    for (std::size_t b = 0; b < F.k; ++b){
      double q = min_quality(c, own, F, true, b);
      if (shared) q = std::min(q, min_quality(c, nbr, F, false, b));
      if (q > min_tet_quality){
        base[f] = std::int32_t(b);
        break;
      }
      if (q > best_q){
        best_q = q;
        best = b;
      }
    }
    // No good base: the best, where OpenFOAM takes point 0, whose fan need not conform
    if (base[f] < 0){
      base[f] = std::int32_t(best);
      ++failed;
    }
  }
  // A cyclic's other side takes the owner side's base, reindexed
  for (const auto& p : c.patches){
    if (p.type != "cyclic" || p.owner) continue;
    const auto& q = c.patches[std::size_t(p.neighbour)];
    for (std::size_t i = 0; i < std::size_t(p.size); ++i){
      const std::int32_t b = base[std::size_t(q.start) + i];
      if (b < 0) continue;
      const std::size_t f = std::size_t(p.start) + i;
      base[f] = b < 1 ? b : std::int32_t(face(c, f).k) - b;
    }
  }
  if (failures) *failures = failed;
  return base;
}

std::vector<std::int32_t> cyclic_masters(const CaseData& c){
  const std::size_t np = c.npoints();
  std::vector<std::int32_t> up(np);
  std::iota(up.begin(), up.end(), 0);
  for (const auto& a : c.patches){
    if (a.type != "cyclic" || !a.owner) continue;
    const auto& b = c.patches[std::size_t(a.neighbour)];
    for (std::size_t i = 0; i < std::size_t(a.size); ++i){
      const FaceRef fa = face(c, std::size_t(a.start) + i), fb = face(c, std::size_t(b.start) + i);
      for (std::size_t j = 0; j < fa.k; ++j){
        const std::int32_t r1 = root(up, fa[j]), r2 = root(up, fb[(fa.k - j) % fa.k]);
        if (r1 != r2) up[std::max(r1, r2)] = std::min(r1, r2);
      }
    }
  }
  // The root is the smallest of its set: every union links the larger root below the smaller
  std::vector<std::int32_t> m(np);
  for (std::size_t p = 0; p < np; ++p) m[p] = root(up, std::int32_t(p));
  return m;
}

SplitData split(const CaseData& c, const int tets_per_hex){
  if (tets_per_hex != 6 && tets_per_hex != 12)
    partrac::fail("split is 6 or 12 tets a hex, not ", tets_per_hex);
  SplitData s;
  if (c.empty_axis >= 0) split_2d(c, tets_per_hex, s);
  else split_3d(c, tets_per_hex, s);
  return s;
}

}  // namespace openfoam_split
