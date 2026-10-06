// OpenFOAM meshes as openfoam_load::CaseData without OpenFOAM, for the unit
// tests of the split and the node values: hex lattices, single cells and the
// checked-in fixtures' polyMesh files (ascii or binary, plain or gzipped),
// with the geometry by OpenFOAM's formulas (derive_geometry). A test's
// reader, not the loader's: the loader reads through OpenFOAM.
#pragma once

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

#include "openfoam_load.hpp"

namespace foam_arrays {

using openfoam_load::CaseData;
using openfoam_load::Patch;

using P3 = std::array<double, 3>;

inline P3 sub(const P3& a, const P3& b){ return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline P3 cross(const P3& a, const P3& b){
  return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
}
inline double dot(const P3& a, const P3& b){ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

// The geometry OpenFOAM derives from the points and faces, by its formulas
// (face::areaAndCentre, primitiveMesh::makeCellCentresAndVols): the cell
// centres and volumes, the boundary faces' centres and areas, each cyclic
// pair's separation
inline void derive_geometry(CaseData& c){
  const std::size_t nf = c.nfaces();
  const auto point = [&](const std::int32_t p){
    const std::size_t i = 3*std::size_t(p);
    return P3{c.points[i], c.points[i + 1], c.points[i + 2]};
  };
  std::vector<P3> fc(nf), fa(nf);
  for (std::size_t f = 0; f < nf; ++f){
    const std::int32_t* F = c.face_points.data() + c.face_start[f];
    const std::size_t k = std::size_t(c.face_start[f + 1] - c.face_start[f]);
    if (k == 3){
      const P3 a = point(F[0]), b = point(F[1]), d = point(F[2]);
      const P3 n = cross(sub(b, a), sub(d, a));
      for (std::size_t i = 0; i < 3; ++i){
        fa[f][i] = 0.5*n[i];
        fc[f][i] = (a[i] + b[i] + d[i])/3.;
      }
      continue;
    }
    P3 pav{0., 0., 0.};
    for (std::size_t j = 0; j < k; ++j){
      const P3 x = point(F[j]);
      for (std::size_t i = 0; i < 3; ++i) pav[i] += x[i];
    }
    for (std::size_t i = 0; i < 3; ++i) pav[i] /= double(k);
    std::vector<P3> a(k);
    P3 sa{0., 0., 0.};
    for (std::size_t j = 0; j < k; ++j){
      const P3 x = point(F[j]), y = point(F[(j + 1) % k]);
      a[j] = cross(sub(y, x), sub(pav, x));
      for (std::size_t i = 0; i < 3; ++i) sa[i] += a[j][i];
    }
    const double mag = std::sqrt(dot(sa, sa));
    P3 hat{0., 0., 0.};
    if (mag > 0) for (std::size_t i = 0; i < 3; ++i) hat[i] = sa[i]/mag;
    P3 num{0., 0., 0.};
    double san = 0.;
    for (std::size_t j = 0; j < k; ++j){
      const P3 x = point(F[j]), y = point(F[(j + 1) % k]);
      const double an = dot(a[j], hat);
      san += an;
      for (std::size_t i = 0; i < 3; ++i) num[i] += an*(x[i] + y[i] + pav[i]);
    }
    for (std::size_t i = 0; i < 3; ++i){
      fc[f][i] = san > 1e-300 ? num[i]/(3*san) : pav[i];
      fa[f][i] = 0.5*sa[i];
    }
  }
  const std::size_t nc = std::size_t(c.ncells), ni = std::size_t(c.n_internal);
  std::vector<P3> est(nc, P3{0., 0., 0.});
  std::vector<int> cnt(nc, 0);
  for (int side = 0; side < 2; ++side)
    for (std::size_t f = 0; f < (side ? ni : nf); ++f){
      const std::size_t cell = std::size_t(side ? c.neighbour[f] : c.owner[f]);
      for (std::size_t i = 0; i < 3; ++i) est[cell][i] += fc[f][i];
      ++cnt[cell];
    }
  for (std::size_t i = 0; i < nc; ++i)
    for (std::size_t d = 0; d < 3; ++d) est[i][d] /= double(cnt[i]);
  std::vector<P3> ctr(nc, P3{0., 0., 0.});
  std::vector<double> vol(nc, 0.);
  for (int side = 0; side < 2; ++side)
    for (std::size_t f = 0; f < (side ? ni : nf); ++f){
      const std::size_t cell = std::size_t(side ? c.neighbour[f] : c.owner[f]);
      const double pyr3 = (side ? -1. : 1.)*dot(fa[f], sub(fc[f], est[cell]));
      for (std::size_t i = 0; i < 3; ++i) ctr[cell][i] += pyr3*(0.75*fc[f][i] + 0.25*est[cell][i]);
      vol[cell] += pyr3;
    }
  c.cell_centres.resize(3*nc);
  c.cell_volumes.resize(nc);
  for (std::size_t i = 0; i < nc; ++i){
    const bool ok = std::abs(vol[i]) > 1e-300;
    for (std::size_t d = 0; d < 3; ++d) c.cell_centres[3*i + d] = ok ? ctr[i][d]/vol[i] : est[i][d];
    c.cell_volumes[i] = vol[i]/3.;
  }
  const std::size_t nb = nf - ni;
  c.bface_centres.resize(3*nb);
  c.bface_areas.resize(3*nb);
  for (std::size_t b = 0; b < nb; ++b)
    for (std::size_t d = 0; d < 3; ++d){
      c.bface_centres[3*b + d] = fc[ni + b][d];
      c.bface_areas[3*b + d] = fa[ni + b][d];
    }
  for (auto& p : c.patches){
    if (p.type != "cyclic") continue;
    const auto& q = c.patches[std::size_t(p.neighbour)];
    P3 sep{0., 0., 0.};
    for (std::int64_t i = 0; i < p.size; ++i)
      for (std::size_t d = 0; d < 3; ++d) sep[d] += fc[std::size_t(q.start + i)][d] - fc[std::size_t(p.start + i)][d];
    for (std::size_t d = 0; d < 3; ++d) p.separation[d] = p.size ? sep[d]/double(p.size) : 0.;
  }
}

// A file's bytes, gunzipped where it is name.gz
inline std::string slurp(const std::string& path){
  gzFile f = gzopen(path.c_str(), "rb");
  if (!f) f = gzopen((path + ".gz").c_str(), "rb");
  if (!f) throw std::runtime_error("no file " + path);
  std::string out;
  char buf[65536];
  int n;
  while ((n = gzread(f, buf, sizeof buf)) > 0) out.append(buf, std::size_t(n));
  gzclose(f);
  return out;
}

// A FoamFile: its header's format and class, and the position after the header
struct FoamText {
  std::string s, format, cls;
  std::size_t pos = 0;
  explicit FoamText(const std::string& path) : s(slurp(path)){
    const std::size_t h = s.find("FoamFile");
    const std::size_t e = s.find('}', h);
    const std::string head = s.substr(h, e - h);
    format = word_after(head, "format");
    cls = word_after(head, "class");
    pos = e + 1;
  }
  static std::string word_after(const std::string& t, const char* key){
    std::size_t p = t.find(key);
    if (p == std::string::npos) return "";
    p += std::strlen(key);
    while (std::isspace(static_cast<unsigned char>(t[p]))) ++p;
    std::size_t q = p;
    while (q < t.size() && t[q] != ';' && !std::isspace(static_cast<unsigned char>(t[q]))) ++q;
    return t.substr(p, q - p);
  }
  void skip(){
    for (;;){
      while (pos < s.size() && std::isspace(static_cast<unsigned char>(s[pos]))) ++pos;
      if (s.compare(pos, 2, "//") == 0) pos = s.find('\n', pos);
      else if (s.compare(pos, 2, "/*") == 0) pos = s.find("*/", pos) + 2;
      else return;
    }
  }
  long long integer(){
    skip();
    std::size_t q = pos;
    while (q < s.size() && (std::isdigit(static_cast<unsigned char>(s[q])) || s[q] == '-')) ++q;
    const long long v = std::stoll(s.substr(pos, q - pos));
    pos = q;
    return v;
  }
  double number(){
    skip();
    std::size_t q = pos;
    while (q < s.size() && !std::isspace(static_cast<unsigned char>(s[q])) && s[q] != ')' && s[q] != '(') ++q;
    const double v = std::stod(s.substr(pos, q - pos));
    pos = q;
    return v;
  }
  void expect(const char c){
    skip();
    if (s[pos] != c) throw std::runtime_error(std::string("expected ") + c);
    ++pos;
  }
  bool binary() const { return format == "binary"; }
  // n raw values of type T between parentheses
  template<typename T>
  std::vector<T> raw(const std::size_t n){
    expect('(');
    std::vector<T> v(n);
    std::memcpy(v.data(), s.data() + pos, n*sizeof(T));
    pos += n*sizeof(T);
    expect(')');
    return v;
  }
  std::vector<std::int64_t> labels(){
    const std::size_t n = std::size_t(integer());
    if (binary()){
      const auto v = raw<std::int32_t>(n);
      return std::vector<std::int64_t>(v.begin(), v.end());
    }
    std::vector<std::int64_t> v(n);
    expect('(');
    for (auto& x : v) x = integer();
    expect(')');
    return v;
  }
};

inline std::vector<double> read_points(const std::string& path){
  FoamText t(path);
  const std::size_t n = std::size_t(t.integer());
  if (t.binary()) return t.raw<double>(3*n);
  std::vector<double> x(3*n);
  t.expect('(');
  for (std::size_t i = 0; i < n; ++i){
    t.expect('(');
    for (int d = 0; d < 3; ++d) x[3*i + std::size_t(d)] = t.number();
    t.expect(')');
  }
  return x;
}

inline void read_faces(const std::string& path, CaseData& c){
  FoamText t(path);
  c.face_start.assign(1, 0);
  c.face_points.clear();
  if (t.cls == "faceCompactList"){
    const auto start = t.labels();
    const auto pts = t.labels();
    c.face_start.assign(start.begin(), start.end());
    c.face_points.assign(pts.begin(), pts.end());
    return;
  }
  const std::size_t n = std::size_t(t.integer());
  t.expect('(');
  for (std::size_t f = 0; f < n; ++f){
    const auto k = t.integer();
    t.expect('(');
    for (long long j = 0; j < k; ++j) c.face_points.push_back(std::int32_t(t.integer()));
    t.expect(')');
    c.face_start.push_back(std::int64_t(c.face_points.size()));
  }
}

// constant/polyMesh/boundary: the patches and their cyclic partners
inline std::vector<Patch> read_boundary(const std::string& path){
  FoamText t(path);
  const std::size_t n = std::size_t(t.integer());
  t.expect('(');
  std::vector<Patch> out(n);
  std::vector<std::string> nbr(n);
  for (std::size_t i = 0; i < n; ++i){
    t.skip();
    std::size_t q = t.pos;
    while (!std::isspace(static_cast<unsigned char>(t.s[q])) && t.s[q] != '{') ++q;
    out[i].name = t.s.substr(t.pos, q - t.pos);
    t.pos = q;
    t.expect('{');
    const std::size_t e = t.s.find('}', t.pos);
    const std::string body = t.s.substr(t.pos, e - t.pos);
    out[i].type = FoamText::word_after(body, "type");
    out[i].size = std::stoll(FoamText::word_after(body, "nFaces"));
    out[i].start = std::stoll(FoamText::word_after(body, "startFace"));
    nbr[i] = FoamText::word_after(body, "neighbourPatch");
    t.pos = e + 1;
  }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      if (out[i].type == "cyclic" && out[j].name == nbr[i]){
        out[i].neighbour = int(j);
        out[i].owner = i < j;
      }
  return out;
}

// The normal axis of the empty patches, -1 without any
inline int empty_axis(const CaseData& c){
  for (const auto& p : c.patches){
    if (p.type != "empty" || p.size == 0) continue;
    const double* a = &c.bface_areas[3*std::size_t(p.start - c.n_internal)];
    int ax = 0;
    for (int d = 1; d < 3; ++d) if (std::abs(a[d]) > std::abs(a[ax])) ax = d;
    return ax;
  }
  return -1;
}

// A case's polyMesh with the geometry derived
inline CaseData read_case(const std::string& dir){
  CaseData c;
  c.dir = dir;
  const std::string poly = dir + "/constant/polyMesh/";
  c.points = read_points(poly + "points");
  read_faces(poly + "faces", c);
  const auto own = FoamText(poly + "owner").labels();
  const auto nbr = FoamText(poly + "neighbour").labels();
  c.owner.assign(own.begin(), own.end());
  c.neighbour.assign(nbr.begin(), nbr.end());
  c.n_internal = std::int64_t(c.neighbour.size());
  c.ncells = 1 + std::max(*std::max_element(c.owner.begin(), c.owner.end()),
                          c.neighbour.empty() ? 0 : *std::max_element(c.neighbour.begin(), c.neighbour.end()));
  c.patches = read_boundary(poly + "boundary");
  derive_geometry(c);
  c.empty_axis = empty_axis(c);
  return c;
}

// A hex lattice on [0,1]^3, n cells an axis: patches xmin, xmax, ymin, ymax,
// zmin, zmax, a cyclic pair along each axis in `cyclic`, walls elsewhere.
// Internal faces point to the higher cell; each face starts at a corner
// drawn from rot_seed (a cyclic pair's two sides alike), the points are
// renumbered by a shuffle from perm_seed and moved by up to jitter of the
// spacing, interior points only; seeds of 0 leave either alone.
inline CaseData lattice(const std::array<int, 3> n, const std::array<bool, 3> cyclic = {false, false, false},
                        const unsigned perm_seed = 0, const unsigned rot_seed = 0, const double jitter = 0.){
  CaseData c;
  const int nx = n[0], ny = n[1], nz = n[2];
  const auto pid = [&](int i, int j, int k){ return i + (nx + 1)*(j + (ny + 1)*k); };
  const auto cid = [&](int i, int j, int k){ return i + nx*(j + ny*k); };
  const int np = (nx + 1)*(ny + 1)*(nz + 1);
  c.points.resize(3*std::size_t(np));
  std::mt19937 jit(12345u);
  std::uniform_real_distribution<double> u(-1., 1.);
  for (int k = 0; k <= nz; ++k)
    for (int j = 0; j <= ny; ++j)
      for (int i = 0; i <= nx; ++i){
        const int ijk[3] = {i, j, k};
        const bool inner = i > 0 && i < nx && j > 0 && j < ny && k > 0 && k < nz;
        for (int d = 0; d < 3; ++d)
          c.points[3*std::size_t(pid(i, j, k)) + std::size_t(d)] =
            double(ijk[d])/double(n[std::size_t(d)]) + (inner && jitter > 0 ? jitter*u(jit)/double(n[std::size_t(d)]) : 0.);
      }
  std::mt19937 rot(rot_seed);
  // The face normal to axis a at plane s, corner (u, v), normal +a
  const auto quad = [&](int a, int s, int uu, int vv){
    const int o0 = (a + 1) % 3, o1 = (a + 2) % 3;
    std::array<std::int32_t, 4> q;
    const int du[4] = {0, 1, 1, 0}, dv[4] = {0, 0, 1, 1};
    for (int m = 0; m < 4; ++m){
      int ijk[3];
      ijk[a] = s;
      ijk[o0] = uu + du[m];
      ijk[o1] = vv + dv[m];
      q[std::size_t(m)] = pid(ijk[0], ijk[1], ijk[2]);
    }
    return q;
  };
  const auto rotate = [&](std::array<std::int32_t, 4> q, int r){
    std::rotate(q.begin(), q.begin() + r, q.end());
    return q;
  };
  // Internal faces, then each patch
  std::vector<std::array<std::int32_t, 4>> faces;
  for (int a = 0; a < 3; ++a){
    const int o0 = (a + 1) % 3, o1 = (a + 2) % 3;
    for (int s = 1; s < n[std::size_t(a)]; ++s)
      for (int vv = 0; vv < n[std::size_t(o1)]; ++vv)
        for (int uu = 0; uu < n[std::size_t(o0)]; ++uu){
          faces.push_back(rotate(quad(a, s, uu, vv), rot_seed ? int(rot() % 4) : 0));
          int lo[3], hi[3];
          lo[a] = s - 1; hi[a] = s;
          lo[o0] = hi[o0] = uu;
          lo[o1] = hi[o1] = vv;
          c.owner.push_back(cid(lo[0], lo[1], lo[2]));
          c.neighbour.push_back(cid(hi[0], hi[1], hi[2]));
        }
  }
  c.n_internal = std::int64_t(faces.size());
  const char* names[6] = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
  std::vector<int> rots;
  for (int a = 0; a < 3; ++a){
    const int o0 = (a + 1) % 3, o1 = (a + 2) % 3;
    for (int side = 0; side < 2; ++side){
      Patch p;
      p.name = names[2*a + side];
      p.type = cyclic[std::size_t(a)] ? "cyclic" : "wall";
      p.start = std::int64_t(faces.size());
      int m = 0;
      for (int vv = 0; vv < n[std::size_t(o1)]; ++vv)
        for (int uu = 0; uu < n[std::size_t(o0)]; ++uu, ++m){
          auto q = quad(a, side ? n[std::size_t(a)] : 0, uu, vv);
          // the low side points out along -a: the same first point, the rest reversed
          if (!side) q = {q[0], q[3], q[2], q[1]};
          int r = rot_seed ? int(rot() % 4) : 0;
          if (cyclic[std::size_t(a)]){
            if (!side) rots.push_back(r);
            // face i there is face i here, its points reversed from the same first
            else r = (4 - rots[std::size_t(m)]) % 4;
          }
          faces.push_back(rotate(q, r));
          int cell[3];
          cell[a] = side ? n[std::size_t(a)] - 1 : 0;
          cell[o0] = uu;
          cell[o1] = vv;
          c.owner.push_back(cid(cell[0], cell[1], cell[2]));
        }
      p.size = std::int64_t(faces.size()) - p.start;
      if (cyclic[std::size_t(a)]){
        p.neighbour = int(c.patches.size()) + (side ? -1 : 1);
        p.owner = side == 0;
      }
      c.patches.push_back(p);
    }
    rots.clear();
  }
  // The points shuffled
  std::vector<std::int32_t> perm(static_cast<std::size_t>(np));
  for (int i = 0; i < np; ++i) perm[std::size_t(i)] = i;
  if (perm_seed){
    std::mt19937 g(perm_seed);
    std::shuffle(perm.begin(), perm.end(), g);
  }
  std::vector<double> x(c.points.size());
  for (int i = 0; i < np; ++i)
    for (int d = 0; d < 3; ++d) x[3*std::size_t(perm[std::size_t(i)]) + std::size_t(d)] = c.points[3*std::size_t(i) + std::size_t(d)];
  c.points.swap(x);
  c.face_start.assign(1, 0);
  for (const auto& q : faces){
    for (const auto p : q) c.face_points.push_back(perm[std::size_t(p)]);
    c.face_start.push_back(std::int64_t(c.face_points.size()));
  }
  c.ncells = std::int64_t(nx)*ny*nz;
  derive_geometry(c);
  return c;
}

// One cell from its points and faces (each pointing out), a wall all round
inline CaseData single_cell(const std::vector<std::array<double, 3>>& points,
                            const std::vector<std::vector<std::int32_t>>& faces){
  CaseData c;
  for (const auto& p : points) c.points.insert(c.points.end(), p.begin(), p.end());
  c.face_start.assign(1, 0);
  for (const auto& f : faces){
    c.face_points.insert(c.face_points.end(), f.begin(), f.end());
    c.face_start.push_back(std::int64_t(c.face_points.size()));
    c.owner.push_back(0);
  }
  c.ncells = 1;
  Patch p;
  p.name = "walls";
  p.type = "wall";
  p.size = std::int64_t(faces.size());
  c.patches.push_back(p);
  derive_geometry(c);
  return c;
}

// A hex lattice on the coordinates ax (each axis ascending): patches xmin,
// xmax, ymin, ymax, zmin, zmax of the types given (patch by default; a
// cyclic side pairs with its opposite), points x-fastest, internal faces to
// the higher cell; jitter moves the interior points by up to that fraction
// of the smallest spacing (seeded), an empty pair makes it a 2D case
inline CaseData box(const std::array<std::vector<double>, 3>& ax, const std::map<std::string, std::string>& types = {},
                    const double jitter = 0., const unsigned seed = 1){
  CaseData c;
  const int n[3] = {int(ax[0].size()) - 1, int(ax[1].size()) - 1, int(ax[2].size()) - 1};
  const auto pid = [&](int i, int j, int k){ return i + (n[0] + 1)*(j + (n[1] + 1)*k); };
  const auto cid = [&](int i, int j, int k){ return i + n[0]*(j + n[1]*k); };
  double h = 1e300;
  for (int d = 0; d < 3; ++d)
    for (int i = 0; i < n[d]; ++i) h = std::min(h, ax[std::size_t(d)][std::size_t(i + 1)] - ax[std::size_t(d)][std::size_t(i)]);
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(-1., 1.);
  c.points.resize(3*std::size_t((n[0] + 1)*(n[1] + 1)*(n[2] + 1)));
  for (int k = 0; k <= n[2]; ++k)
    for (int j = 0; j <= n[1]; ++j)
      for (int i = 0; i <= n[0]; ++i){
        const int ijk[3] = {i, j, k};
        const bool inner = i > 0 && i < n[0] && j > 0 && j < n[1] && k > 0 && k < n[2];
        for (int d = 0; d < 3; ++d)
          c.points[3*std::size_t(pid(i, j, k)) + std::size_t(d)] =
            ax[std::size_t(d)][std::size_t(ijk[d])] + (inner && jitter > 0 ? jitter*h*u(rng) : 0.);
      }
  const auto quad = [&](int a, int s, int uu, int vv){
    const int o0 = (a + 1) % 3, o1 = (a + 2) % 3;
    std::array<std::int32_t, 4> q;
    const int du[4] = {0, 1, 1, 0}, dv[4] = {0, 0, 1, 1};
    for (int m = 0; m < 4; ++m){
      int ijk[3];
      ijk[a] = s;
      ijk[o0] = uu + du[m];
      ijk[o1] = vv + dv[m];
      q[std::size_t(m)] = pid(ijk[0], ijk[1], ijk[2]);
    }
    return q;
  };
  std::vector<std::array<std::int32_t, 4>> faces;
  for (int a = 0; a < 3; ++a){
    const int o0 = (a + 1) % 3, o1 = (a + 2) % 3;
    for (int s = 1; s < n[a]; ++s)
      for (int vv = 0; vv < n[o1]; ++vv)
        for (int uu = 0; uu < n[o0]; ++uu){
          faces.push_back(quad(a, s, uu, vv));
          int lo[3], hi[3];
          lo[a] = s - 1; hi[a] = s;
          lo[o0] = hi[o0] = uu;
          lo[o1] = hi[o1] = vv;
          c.owner.push_back(cid(lo[0], lo[1], lo[2]));
          c.neighbour.push_back(cid(hi[0], hi[1], hi[2]));
        }
  }
  c.n_internal = std::int64_t(faces.size());
  const char* names[6] = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
  for (int a = 0; a < 3; ++a){
    const int o0 = (a + 1) % 3, o1 = (a + 2) % 3;
    for (int side = 0; side < 2; ++side){
      Patch p;
      p.name = names[2*a + side];
      const auto t = types.find(p.name);
      p.type = t == types.end() ? "patch" : t->second;
      p.start = std::int64_t(faces.size());
      for (int vv = 0; vv < n[o1]; ++vv)
        for (int uu = 0; uu < n[o0]; ++uu){
          auto q = quad(a, side ? n[a] : 0, uu, vv);
          if (!side) q = {q[0], q[3], q[2], q[1]};
          faces.push_back(q);
          int cell[3];
          cell[a] = side ? n[a] - 1 : 0;
          cell[o0] = uu;
          cell[o1] = vv;
          c.owner.push_back(cid(cell[0], cell[1], cell[2]));
        }
      p.size = std::int64_t(faces.size()) - p.start;
      if (p.type == "cyclic"){
        p.neighbour = int(c.patches.size()) + (side ? -1 : 1);
        p.owner = side == 0;
      }
      c.patches.push_back(p);
    }
  }
  c.face_start.assign(1, 0);
  for (const auto& q : faces){
    c.face_points.insert(c.face_points.end(), q.begin(), q.end());
    c.face_start.push_back(std::int64_t(c.face_points.size()));
  }
  c.ncells = std::int64_t(n[0])*n[1]*n[2];
  derive_geometry(c);
  c.empty_axis = empty_axis(c);
  return c;
}

// n + 1 coordinates on [0, 1], blockMesh's simpleGrading: the last cell over the first is ratio
inline std::vector<double> graded(const int n, const double ratio){
  std::vector<double> x(std::size_t(n) + 1, 0.);
  const double r = n > 1 ? std::pow(ratio, 1./double(n - 1)) : 1.;
  double sum = 0., w = 1.;
  for (int i = 0; i < n; ++i, w *= r) sum += w;
  w = 1.;
  for (int i = 0; i < n; ++i, w *= r) x[std::size_t(i) + 1] = x[std::size_t(i)] + w/sum;
  x[std::size_t(n)] = 1.;
  return x;
}

// n + 1 coordinates from a to b, equally spaced
inline std::vector<double> linspace(const double a, const double b, const int n){
  std::vector<double> x(std::size_t(n) + 1);
  for (int i = 0; i <= n; ++i) x[std::size_t(i)] = a + (b - a)*double(i)/double(n);
  x[std::size_t(n)] = b;
  return x;
}

}  // namespace foam_arrays
