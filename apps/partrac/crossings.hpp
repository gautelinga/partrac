#ifndef __CROSSINGS_HPP
#define __CROSSINGS_HPP

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include <omp.h>

#include "Error.hpp"
#include "Initializer.hpp"
#include "ParticleSet.hpp"

// A plane normal: x, y, z or nx,ny,nz, normalized; zero for none
inline Vector3d parse_normal(const std::string& value){
  Vector3d n = unit_vector(value);
  char rest;
  if (value == "none" || n.norm() > 0.)
    return n;
  if (std::sscanf(value.c_str(), "%lf,%lf,%lf%c", &n[0], &n[1], &n[2], &rest) != 3)
    partrac::fail("plane normal ", value, ": not x, y, z or nx,ny,nz");
  if (!(n.norm() > 0.))
    partrac::fail("plane normal ", value, ": zero");
  return n / n.norm();
}

// First crossings of the planes x.n = k dn: when a node's running max of x.n passes
// a plane, a record per edge of the node, six float64 after a text line: c, the time
// within the step, ln(dl/dl0) at its end, dl0/2, dl/2, k; without edges (particles),
// a record per node, three float64: id, the time, k. A new node starts from its
// neighbours' mean, so material on an edge split across a plane can be lost or
// counted twice. "<file>.filter": t, edges before, after. A restart appends, the
// running max afresh.
class Crossings {
public:
  Crossings(const ParticleSet& ps, const bool has_edges, const Vector3d& normal, const double spacing,
            const std::string& path, const bool append)
    : edges_(has_edges), n(normal), dn(spacing), out(path, append ? std::ios::binary | std::ios::app : std::ios::binary),
      flog(path + ".filter", append ? std::ios::app : std::ios::out) {
    if (!out || !flog)
      partrac::fail("cannot open ", path, " or ", path, ".filter");
    if (!append)
      out << "partrac crossings: float64 " << (edges_ ? "c t lnrho weight length k" : "id t k")
          << "; n = (" << n[0] << ", " << n[1] << ", " << n[2] << "), dn = " << dn << "\n";
    for (Uint i = 0; i < ps.N(); ++i){
      grow(ps.id(i));
      s_[ps.id(i)] = smax_[ps.id(i)] = ps.x(i).dot(n);
    }
  }

  // After the step from t to t + h
  void record(const ParticleSet& ps, const EdgesType& edges, const Node2EdgesType& node2edges,
              const double t, const double h){
    const Uint N = ps.N();
    if (edges_ && node2edges.size() != N)
      partrac::fail("crossings: node2edges has ", node2edges.size(), " rows for ", N, " nodes");
    Uint max_id = 0;
    std::vector<std::vector<Uint>> fresh(omp_get_max_threads());
    #pragma omp parallel for reduction(max:max_id)
    for (Uint i = 0; i < N; ++i){
      max_id = std::max(max_id, ps.id(i));
      if (ps.id(i) >= s_.size() || std::isnan(smax_[ps.id(i)]))
        fresh[omp_get_thread_num()].push_back(i);
    }
    grow(max_id);
    // New nodes: their seen neighbours' mean before the step
    for (const auto& f : fresh)
      for (const Uint i : f){
        double s = 0., m = 0.;
        int count = 0;
        if (edges_)
          for (const Uint e : node2edges[i]){
            const Uint j = ps.id(edges[e].first[0] == i ? edges[e].first[1] : edges[e].first[0]);
            if (!std::isnan(smax_[j])){ s += s_[j]; m += smax_[j]; ++count; }
          }
        s_[ps.id(i)] = count > 0 ? s / count : ps.x(i).dot(n);
        smax_[ps.id(i)] = count > 0 ? m / count : s_[ps.id(i)];
      }
    std::vector<std::vector<double>> found(omp_get_max_threads());
    #pragma omp parallel for
    for (Uint i = 0; i < N; ++i){
      const Uint id = ps.id(i);
      const double s = ps.x(i).dot(n), s0 = s_[id], m0 = smax_[id];
      s_[id] = s;
      if (s <= m0)
        continue;
      auto& f = found[omp_get_thread_num()];
      for (long k = static_cast<long>(std::floor(m0 / dn)) + 1; k * dn <= s; ++k){
        const double tk = t + h * std::max((k * dn - s0) / (s - s0), 0.);
        if (!edges_)
          f.insert(f.end(), {double(id), tk, double(k)});
        else
          for (const Uint e : node2edges[i]){
            const double dl = ps.dist(edges[e].first[0], edges[e].first[1]), dl0 = edges[e].second;
            f.insert(f.end(), {ps.c(i), tk, std::log(dl / dl0), dl0 / 2, dl / 2, double(k)});
          }
      }
      smax_[id] = s;
    }
    for (const auto& f : found)
      out.write(reinterpret_cast<const char*>(f.data()), f.size() * sizeof(f[0]));
  }

  void filtered(const double t, const Uint before, const Uint after){
    flog << t << " " << before << " " << after << "\n";
  }

private:
  const bool edges_;
  const Vector3d n;
  const double dn;
  std::ofstream out, flog;
  std::vector<double> s_, smax_;   // by node id, NaN: not seen; 16 bytes an id ever made

  void grow(const Uint id){
    if (id >= s_.size()){
      s_.resize(std::max<std::size_t>(id + 1, 2 * s_.size()), std::numeric_limits<double>::quiet_NaN());
      smax_.resize(s_.size(), std::numeric_limits<double>::quiet_NaN());
    }
  }
};

#endif
