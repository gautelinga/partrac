#include "mesh_tables.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <utility>
#include <omp.h>

#include "Error.hpp"

namespace mesh_tables {

namespace {

// The thread's slice of a sorted key array, moved out to group boundaries so no
// group is split between threads
void thread_range(const std::vector<std::uint64_t>& key, const int t, const int nt,
                  std::size_t& lo, std::size_t& hi){
  const std::size_t n = key.size();
  lo = n*std::size_t(t)/std::size_t(nt);
  hi = n*std::size_t(t+1)/std::size_t(nt);
  while (lo > 0 && lo < n && key[lo] == key[lo-1]) ++lo;
  while (hi > 0 && hi < n && key[hi] == key[hi-1]) ++hi;
  if (hi < lo) hi = lo;
}

// Entries fitting a 32-bit payload
void check_entry_count(const std::size_t m, const char* what){
  if (m > std::numeric_limits<std::uint32_t>::max()){
    partrac::fail(what, ": ", m, " entries, more than a 32-bit index can hold");
  }
}

// An exclusive prefix sum of a counting sort's histogram, a[0] left at zero:
// each thread's block total first, then its own block from that offset
void block_scan(std::vector<std::uint32_t>& a){
  const std::size_t n = a.size() - 1;
  const int nt = omp_get_max_threads();
  std::vector<std::uint32_t> base(std::size_t(nt) + 1, 0);
#pragma omp parallel num_threads(nt)
  {
    const int t = omp_get_thread_num();
    const std::size_t lo = 1 + (n*std::size_t(t))/std::size_t(nt);
    const std::size_t hi = 1 + (n*std::size_t(t+1))/std::size_t(nt);
    std::uint32_t run = 0;
    for (std::size_t i = lo; i < hi; ++i) run += a[i];
    base[std::size_t(t) + 1] = run;
#pragma omp barrier
#pragma omp single
    for (int u = 0; u < nt; ++u) base[std::size_t(u) + 1] += base[std::size_t(u)];
    std::uint32_t acc = base[std::size_t(t)];
    for (std::size_t i = lo; i < hi; ++i){ acc += a[i]; a[i] = acc; }
  }
}

// The key of a facet: its two lowest vertices, the third left as the tail
template<int NV>
void facet_key(const std::uint32_t* row, const int k, std::uint64_t& key, std::uint32_t& tail){
  std::array<std::uint32_t, NV-1> w;
  int m = 0;
  for (int j = 0; j < NV; ++j)
    if (j != k) w[m++] = row[j];
  std::sort(w.begin(), w.end());
  key = (std::uint64_t(w[0]) << 32) | w[1];
  tail = NV > 3 ? w[NV-2] : 0;
}

}  // namespace

void sort_by_key(std::vector<std::uint64_t>& key, std::vector<std::uint32_t>& payload){
  const std::size_t n = key.size();
  if (payload.size() != n){
    partrac::fail("internal: sort_by_key has ", n, " keys and ", payload.size(), " payloads");
  }
  if (n < 2) return;
  const int nt = omp_get_max_threads();
  constexpr int n_digits = 8;
  constexpr int n_bins = 256;
  // All digit histograms in one read pass
  std::vector<std::size_t> hist(n_digits*n_bins, 0);
  {
    std::vector<std::size_t> local(std::size_t(nt)*n_digits*n_bins, 0);
#pragma omp parallel
    {
      std::size_t* h = local.data() + std::size_t(omp_get_thread_num())*n_digits*n_bins;
#pragma omp for schedule(static)
      for (std::size_t i = 0; i < n; ++i){
        const std::uint64_t k = key[i];
        for (int d = 0; d < n_digits; ++d) ++h[d*n_bins + ((k >> (8*d)) & 0xff)];
      }
    }
    for (int t = 0; t < nt; ++t)
      for (std::size_t j = 0; j < std::size_t(n_digits)*n_bins; ++j)
        hist[j] += local[std::size_t(t)*n_digits*n_bins + j];
  }
  std::vector<std::uint64_t> key2(n);
  std::vector<std::uint32_t> pay2(n);
  std::vector<std::size_t> off(std::size_t(nt)*n_bins);
  std::uint64_t* src_k = key.data();
  std::uint32_t* src_p = payload.data();
  std::uint64_t* dst_k = key2.data();
  std::uint32_t* dst_p = pay2.data();
  for (int d = 0; d < n_digits; ++d){
    // A digit every key shares moves nothing
    bool single = false;
    for (int b = 0; b < n_bins; ++b) if (hist[d*n_bins + b] == n){ single = true; break; }
    if (single) continue;
#pragma omp parallel
    {
      const int t = omp_get_thread_num();
      const std::size_t lo = n*std::size_t(t)/std::size_t(nt);
      const std::size_t hi = n*std::size_t(t+1)/std::size_t(nt);
      std::size_t* c = off.data() + std::size_t(t)*n_bins;
      for (int b = 0; b < n_bins; ++b) c[b] = 0;
      for (std::size_t i = lo; i < hi; ++i) ++c[(src_k[i] >> (8*d)) & 0xff];
#pragma omp barrier
#pragma omp single
      {
        // Bin major, thread minor: the pass is stable
        std::size_t run = 0;
        for (int b = 0; b < n_bins; ++b)
          for (int u = 0; u < nt; ++u){
            const std::size_t cnt = off[std::size_t(u)*n_bins + b];
            off[std::size_t(u)*n_bins + b] = run;
            run += cnt;
          }
      }
      for (std::size_t i = lo; i < hi; ++i){
        const std::size_t o = c[(src_k[i] >> (8*d)) & 0xff]++;
        dst_k[o] = src_k[i];
        dst_p[o] = src_p[i];
      }
    }
    std::swap(src_k, dst_k);
    std::swap(src_p, dst_p);
  }
  if (src_k != key.data()){
    key.swap(key2);
    payload.swap(pay2);
  }
}

void sort_by_key_indexed(std::vector<std::uint64_t>& key, std::vector<std::uint32_t>& payload){
  const std::size_t n = key.size();
  check_entry_count(n, "internal: sort_by_key_indexed");
  payload.resize(n);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < n; ++i) payload[i] = std::uint32_t(i);
  sort_by_key(key, payload);
}

template<int NV>
std::size_t build_edge_table(const std::vector<std::uint32_t>& topo, const std::size_t ncells,
                             std::vector<std::uint32_t>& edges){
  constexpr auto loc = local_edges<NV>();
  constexpr int ne = n_local_edges<NV>;
  if (topo.size() < ncells*NV){
    partrac::fail("internal: topology of ", topo.size(), " entries for ", ncells, " cells");
  }
  const std::size_t m = ncells*std::size_t(ne);
  check_entry_count(m, "edge table");
  edges.assign(m, 0);
  if (ncells == 0) return 0;
  std::vector<std::uint64_t> key(m);
  std::vector<std::uint32_t> pay;
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < ncells; ++i){
    const std::uint32_t* row = topo.data() + i*NV;
    for (int e = 0; e < ne; ++e){
      const std::uint32_t a = row[loc[e][0]], b = row[loc[e][1]];
      const std::uint32_t l = std::min(a, b), h = std::max(a, b);
      key[i*ne + e] = (std::uint64_t(l) << 32) | h;
    }
  }
  sort_by_key_indexed(key, pay);
  // Equal keys are one edge; number the groups, then write each entry's id
  const int nt = omp_get_max_threads();
  std::vector<std::size_t> counts(std::size_t(nt) + 1, 0);
#pragma omp parallel
  {
    const int t = omp_get_thread_num();
    std::size_t lo, hi;
    thread_range(key, t, nt, lo, hi);
    std::size_t c = 0;
    for (std::size_t i = lo; i < hi; ++i)
      if (i == 0 || key[i] != key[i-1]) ++c;
    counts[std::size_t(t) + 1] = c;
  }
  for (int t = 0; t < nt; ++t) counts[std::size_t(t) + 1] += counts[std::size_t(t)];
  check_entry_count(counts[std::size_t(nt)], "edge numbering");
#pragma omp parallel
  {
    const int t = omp_get_thread_num();
    std::size_t lo, hi;
    thread_range(key, t, nt, lo, hi);
    std::uint32_t id = std::uint32_t(counts[std::size_t(t)]);
    for (std::size_t i = lo; i < hi; ++i){
      if (i > lo && key[i] != key[i-1]) ++id;
      edges[pay[i]] = id;
    }
  }
  return counts[std::size_t(nt)];
}

template<int NV>
void scatter_dofs_to_nodes(const std::vector<std::uint32_t>& topo,
                           const std::vector<std::uint32_t>& edges,
                           const std::size_t ncells, const std::size_t nverts,
                           const std::size_t nedges, const std::size_t ncomp,
                           const std::vector<std::uint32_t>& cell_dofs,
                           const std::vector<double>& vec,
                           std::vector<double>& values,
                           DofNodes& map){
  constexpr int ne = n_local_edges<NV>;
  const bool quadratic = !edges.empty();
  const std::size_t n_nodes = NV + (quadratic ? std::size_t(ne) : 0);
  const std::size_t n_total = nverts + (quadratic ? nedges : 0);
  if (ncomp == 0){
    partrac::fail("internal: scatter with no components");
  }
  if (cell_dofs.size() != ncells*n_nodes*ncomp){
    partrac::fail("dof table has ", cell_dofs.size(), " entries, not ", ncells*n_nodes*ncomp);
  }
  if (quadratic && edges.size() != ncells*std::size_t(ne)){
    partrac::fail("internal: edge table has ", edges.size(), " entries, not ", ncells*std::size_t(ne));
  }
  const std::size_t m = ncells*n_nodes;
  check_entry_count(m, "dof scatter");
  check_entry_count(n_total*ncomp, "the node numbering");
  values.assign(n_total*ncomp, 0.);
  map.start.assign(vec.size() + 1, 0);
  map.slot.clear();
  if (ncells == 0) return;
  // A dof outside the vector would be read in the grouped pass below
  bool bad_dof = false;
#pragma omp parallel for schedule(static) reduction(||: bad_dof)
  for (std::size_t i = 0; i < cell_dofs.size(); ++i)
    if (cell_dofs[i] >= vec.size()) bad_dof = true;
  if (bad_dof){
    const std::size_t i = std::find_if(cell_dofs.begin(), cell_dofs.end(),
                                       [&](const std::uint32_t d){ return d >= vec.size(); })
                          - cell_dofs.begin();
    partrac::fail("dof ", cell_dofs[i], " of cell ", i/(n_nodes*ncomp),
                  " is outside a vector of ", vec.size());
  }
  // The entries grouped by node, by a counting sort over the node ids
  const auto node_of = [&](const std::size_t i, const std::size_t j){
    return std::size_t(j) < std::size_t(NV)
      ? std::size_t(topo[i*NV + j])
      : nverts + std::size_t(edges[i*std::size_t(ne) + (j - NV)]);
  };
  std::vector<std::uint32_t> start(n_total + 1, 0);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < ncells; ++i)
    for (std::size_t j = 0; j < n_nodes; ++j){
      std::uint32_t* const slot = &start[node_of(i, j) + 1];
#pragma omp atomic update
      ++(*slot);
    }
  block_scan(start);
  std::vector<std::uint32_t> cursor(start.begin(), start.end() - 1);
  std::vector<std::uint32_t> entry(m);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < ncells; ++i)
    for (std::size_t j = 0; j < n_nodes; ++j){
      std::uint32_t* const slot = &cursor[node_of(i, j)];
      std::uint32_t at;
#pragma omp atomic capture
      at = (*slot)++;
      entry[at] = std::uint32_t(i*n_nodes + j);
    }
  // Each node's entries are now contiguous: the lowest cell writes, the rest
  // are read back and compared. A throw here would end the run, so the first
  // disagreement per thread is kept and reported after the region.
  // The stored dof of each node slot is kept for the mapping below; a node no
  // cell refers to keeps no_dof and feeds nothing.
  constexpr std::uint32_t no_dof = std::numeric_limits<std::uint32_t>::max();
  std::vector<std::uint32_t> dof_of(n_total*ncomp, no_dof);
  const int nt = omp_get_max_threads();
  std::vector<std::size_t> bad_node(std::size_t(nt), 0), bad_lo(std::size_t(nt), 0),
    bad_hi(std::size_t(nt), 0), bad_comp(std::size_t(nt), 0);
  std::vector<char> bad(std::size_t(nt), 0);
#pragma omp parallel
  {
    const int t = omp_get_thread_num();
#pragma omp for schedule(static)
    for (std::size_t node = 0; node < n_total; ++node){
      const std::size_t lo = start[node], hi = start[node + 1];
      if (lo == hi) continue;
      // The owner is the lowest entry of the group, whatever order it was placed in
      std::uint32_t e0 = entry[lo];
      for (std::size_t q = lo + 1; q < hi; ++q) e0 = std::min(e0, entry[q]);
      const std::size_t c0 = std::size_t(e0)/n_nodes, j0 = std::size_t(e0) % n_nodes;
      for (std::size_t c = 0; c < ncomp; ++c){
        const std::uint32_t d = cell_dofs[(c0*ncomp + c)*n_nodes + j0];
        values[node*ncomp + c] = vec[d];
        dof_of[node*ncomp + c] = d;
      }
      for (std::size_t q = lo; q < hi; ++q){
        const std::size_t e = entry[q];
        if (e == e0) continue;
        const std::size_t ci = e/n_nodes, ji = e % n_nodes;
        for (std::size_t c = 0; c < ncomp; ++c){
          const std::uint32_t d = cell_dofs[(ci*ncomp + c)*n_nodes + ji];
          if (vec[d] != values[node*ncomp + c] && !bad[std::size_t(t)]){
            bad[std::size_t(t)] = 1;
            bad_node[std::size_t(t)] = node;
            bad_lo[std::size_t(t)] = c0;
            bad_hi[std::size_t(t)] = ci;
            bad_comp[std::size_t(t)] = c;
          }
        }
      }
    }
  }
  for (int t = 0; t < nt; ++t)
    if (bad[std::size_t(t)]){
      const std::size_t node = bad_node[std::size_t(t)], c = bad_comp[std::size_t(t)];
      partrac::fail("node ", node, ", component ", c, ": cell ", bad_lo[std::size_t(t)],
                    " gives ", values[node*ncomp + c], " and cell ", bad_hi[std::size_t(t)],
                    " disagrees");
    }
  // The dof -> node slots mapping, by the same counting sort the groups came
  // from: a periodic-reduced space gives one stored dof a node and its images,
  // and a later stamp has to reach all of them.
  const std::size_t nslots = n_total*ncomp;
#pragma omp parallel for schedule(static)
  for (std::size_t sl = 0; sl < nslots; ++sl){
    if (dof_of[sl] == no_dof) continue;
    std::uint32_t* const at = &map.start[std::size_t(dof_of[sl]) + 1];
#pragma omp atomic update
    ++(*at);
  }
  block_scan(map.start);
  map.slot.resize(map.start.back());
  std::vector<std::uint32_t> place(map.start.begin(), map.start.end() - 1);
#pragma omp parallel for schedule(static)
  for (std::size_t sl = 0; sl < nslots; ++sl){
    if (dof_of[sl] == no_dof) continue;
    std::uint32_t* const at = &place[std::size_t(dof_of[sl])];
    std::uint32_t q;
#pragma omp atomic capture
    q = (*at)++;
    map.slot[q] = std::uint32_t(sl);
  }
  // Ascending within a group, so the mapping does not depend on the placement order
#pragma omp parallel for schedule(static)
  for (std::size_t d = 0; d < vec.size(); ++d)
    std::sort(map.slot.begin() + std::ptrdiff_t(map.start[d]),
              map.slot.begin() + std::ptrdiff_t(map.start[d + 1]));
}

template<int NV>
void build_facet_table(const std::vector<std::uint32_t>& topo, const std::size_t ncells,
                       const std::vector<double>& coords, const Uint gdim,
                       std::vector<std::int32_t>& across,
                       std::vector<std::pair<Vector3d, std::size_t>>& exterior){
  if (ncells > std::size_t(std::numeric_limits<std::int32_t>::max())){
    partrac::fail("mesh has ", ncells, " cells, more than a cell id can hold");
  }
  if (topo.size() < ncells*NV){
    partrac::fail("internal: topology of ", topo.size(), " entries for ", ncells, " cells");
  }
  const std::size_t m = ncells*NV;
  check_entry_count(m, "facet table");
  across.assign(m, facet_wall);
  exterior.clear();
  if (ncells == 0) return;
  std::vector<std::uint64_t> key(m);
  std::vector<std::uint32_t> tail(m);
  std::vector<std::uint32_t> pay;
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < ncells; ++i){
    const std::uint32_t* row = topo.data() + i*NV;
    for (int k = 0; k < NV; ++k)
      facet_key<NV>(row, k, key[i*NV + k], tail[i*NV + k]);
  }
  sort_by_key_indexed(key, pay);
  // A group shares the two lowest vertices; within it the third orders the
  // facets, and equal facets are neighbours. More than two is not a manifold.
  const int nt = omp_get_max_threads();
  std::vector<std::vector<std::size_t>> ext(static_cast<std::size_t>(nt));
  std::vector<std::size_t> bad_slot(std::size_t(nt), 0);
  std::vector<char> bad(std::size_t(nt), 0);
#pragma omp parallel
  {
    const int t = omp_get_thread_num();
    std::size_t lo, hi;
    thread_range(key, t, nt, lo, hi);
    std::vector<std::pair<std::uint32_t, std::uint32_t>> grp;   // (third vertex, slot)
    std::size_t i = lo;
    while (i < hi){
      std::size_t g = i + 1;
      while (g < hi && key[g] == key[i]) ++g;
      grp.clear();
      for (std::size_t q = i; q < g; ++q) grp.push_back({tail[pay[q]], pay[q]});
      std::sort(grp.begin(), grp.end());
      std::size_t a = 0;
      while (a < grp.size()){
        std::size_t b = a + 1;
        while (b < grp.size() && grp[b].first == grp[a].first) ++b;
        if (b - a == 1){
          ext[std::size_t(t)].push_back(grp[a].second);
        }
        else if (b - a == 2){
          across[grp[a].second] = std::int32_t(grp[a+1].second/NV);
          across[grp[a+1].second] = std::int32_t(grp[a].second/NV);
        }
        else if (!bad[std::size_t(t)]){
          bad[std::size_t(t)] = 1;
          bad_slot[std::size_t(t)] = grp[a].second;
        }
        a = b;
      }
      i = g;
    }
  }
  for (int t = 0; t < nt; ++t)
    if (bad[std::size_t(t)]){
      partrac::fail("the facet facing vertex ", bad_slot[std::size_t(t)] % NV, " of cell ",
                    bad_slot[std::size_t(t)]/NV, " is shared by more than two cells");
    }
  // Exterior facets in cell-then-facet order, with their midpoints
  std::vector<std::size_t> slots;
  for (int t = 0; t < nt; ++t)
    slots.insert(slots.end(), ext[std::size_t(t)].begin(), ext[std::size_t(t)].end());
  std::sort(slots.begin(), slots.end());
  exterior.resize(slots.size());
  const std::size_t nfv = NV - 1;
#pragma omp parallel for schedule(static)
  for (std::size_t s = 0; s < slots.size(); ++s){
    const std::size_t slot = slots[s];
    const std::uint32_t* row = topo.data() + (slot/NV)*NV;
    const std::size_t k = slot % NV;
    std::array<std::uint32_t, NV-1> w;
    int mm = 0;
    for (int j = 0; j < NV; ++j)
      if (std::size_t(j) != k) w[mm++] = row[j];
    std::sort(w.begin(), w.end());
    Vector3d p = Vector3d::Zero();
    for (std::size_t j = 0; j < nfv; ++j)
      for (Uint d = 0; d < gdim; ++d) p[d] += coords[std::size_t(w[j])*gdim + d];
    exterior[s] = {Vector3d(p/double(nfv)), slot};
  }
}

void match_periodic_facets(std::vector<std::int32_t>& across,
                           const std::vector<std::pair<Vector3d, std::size_t>>& exterior,
                           const std::vector<bool>& periodic,
                           const Vector3d& x_min, const Vector3d& x_max,
                           const Uint dim, const std::size_t nv, const double tol){
  // Each facet on the low face of a periodic axis, shifted onto the high face
  std::vector<std::vector<std::pair<Vector3d, std::size_t>>> lo(dim), hi(dim);
  for (const auto& f : exterior){
    const Vector3d& pt = f.first;
    for (Uint a = 0; a < dim; ++a){
      if (!periodic[a]) continue;
      if (pt[a] < x_min[a] + tol){
        Vector3d q = pt;
        q[a] += x_max[a] - x_min[a];
        lo[a].push_back({q, f.second});
        break;
      }
      if (pt[a] > x_max[a] - tol){
        hi[a].push_back({pt, f.second});
        break;
      }
    }
  }
  // Match within a sorted window along the next axis
  std::size_t unmatched = 0;
  for (Uint a = 0; a < dim; ++a){
    const Uint b = (a + 1) % dim;
    auto& h = hi[a];
    std::sort(h.begin(), h.end(), [b](const auto& p, const auto& q){ return p.first[b] < q.first[b]; });
    std::vector<bool> h_matched(h.size(), false);
    for (const auto& l : lo[a]){
      const double w = 2*tol + 4*std::numeric_limits<double>::epsilon()*(std::abs(l.first[b]) + 1.);
      auto it = std::lower_bound(h.begin(), h.end(), l.first[b] - w,
                                 [b](const auto& p, const double v){ return p.first[b] < v; });
      bool matched = false;
      for ( ; it != h.end() && it->first[b] <= l.first[b] + w; ++it){
        if ((l.first - it->first).norm() < tol){
          across[l.second] = facet_periodic(std::int32_t(it->second / nv));
          across[it->second] = facet_periodic(std::int32_t(l.second / nv));
          h_matched[it - h.begin()] = true;
          matched = true;
          break;
        }
      }
      if (!matched) ++unmatched;
    }
    unmatched += std::count(h_matched.begin(), h_matched.end(), false);
  }
  if (unmatched)
    std::cout << unmatched << " periodic facets have no partner and reflect as walls" << std::endl;
}

std::vector<std::uint32_t> match_periodic_vertices(const std::vector<double>& coords,
                                                   const std::size_t nverts, const Uint gdim,
                                                   const std::vector<bool>& periodic,
                                                   const Vector3d& x_min, const Vector3d& x_max,
                                                   const double tol){
  std::vector<std::uint32_t> master(nverts);
  for (std::size_t v = 0; v < nverts; ++v) master[v] = std::uint32_t(v);
  bool any = false;
  for (Uint a = 0; a < gdim; ++a) any = any || periodic[a];
  if (!any || nverts == 0) return master;
  const auto point = [&](const std::size_t v){
    Vector3d p = Vector3d::Zero();
    for (Uint d = 0; d < gdim; ++d) p[d] = coords[v*gdim + d];
    return p;
  };
  // One link per vertex, across the first periodic axis whose high face it is
  // on; a corner vertex reaches its master through the chain of links
  std::vector<std::uint32_t> link(nverts);
  for (std::size_t v = 0; v < nverts; ++v) link[v] = std::uint32_t(v);
  std::size_t unmatched = 0;
  for (Uint a = 0; a < gdim; ++a){
    if (!periodic[a]) continue;
    const Uint b = (a + 1) % gdim;
    std::vector<std::pair<Vector3d, std::size_t>> lo, hi;
    for (std::size_t v = 0; v < nverts; ++v){
      const Vector3d p = point(v);
      if (p[a] < x_min[a] + tol)
        lo.push_back({p, v});
      else if (p[a] > x_max[a] - tol && link[v] == v){
        Vector3d q = p;
        q[a] -= x_max[a] - x_min[a];
        hi.push_back({q, v});
      }
    }
    std::sort(lo.begin(), lo.end(), [b](const auto& p, const auto& q){ return p.first[b] < q.first[b]; });
    for (const auto& h : hi){
      const double w = 2*tol + 4*std::numeric_limits<double>::epsilon()*(std::abs(h.first[b]) + 1.);
      auto it = std::lower_bound(lo.begin(), lo.end(), h.first[b] - w,
                                 [b](const auto& p, const double v){ return p.first[b] < v; });
      bool matched = false;
      for ( ; it != lo.end() && it->first[b] <= h.first[b] + w; ++it){
        if ((h.first - it->first).norm() < tol){
          link[h.second] = std::uint32_t(it->second);
          matched = true;
          break;
        }
      }
      if (!matched) ++unmatched;
    }
  }
  if (unmatched)
    std::cout << unmatched << " vertices on a periodic face have no image and keep their own value"
              << std::endl;
  for (std::size_t v = 0; v < nverts; ++v){
    std::uint32_t m = std::uint32_t(v);
    for (std::size_t step = 0; step <= gdim && link[m] != m; ++step) m = link[m];
    master[v] = m;
  }
  return master;
}

template<int NV>
void build_facet_neighbours(const std::vector<std::uint32_t>& topo, const std::size_t ncells,
                            const std::vector<double>& coords, const Uint gdim,
                            const std::vector<bool>& periodic,
                            const Vector3d& x_min, const Vector3d& x_max,
                            const double tol,
                            std::vector<std::int32_t>& across){
  std::vector<std::pair<Vector3d, std::size_t>> exterior;
  build_facet_table<NV>(topo, ncells, coords, gdim, across, exterior);
  match_periodic_facets(across, exterior, periodic, x_min, x_max, gdim, NV, tol);
}

template std::size_t build_edge_table<3>(const std::vector<std::uint32_t>&, const std::size_t,
                                         std::vector<std::uint32_t>&);
template std::size_t build_edge_table<4>(const std::vector<std::uint32_t>&, const std::size_t,
                                         std::vector<std::uint32_t>&);
template void scatter_dofs_to_nodes<3>(const std::vector<std::uint32_t>&, const std::vector<std::uint32_t>&,
                                       const std::size_t, const std::size_t, const std::size_t,
                                       const std::size_t, const std::vector<std::uint32_t>&,
                                       const std::vector<double>&, std::vector<double>&,
                                       DofNodes&);
template void scatter_dofs_to_nodes<4>(const std::vector<std::uint32_t>&, const std::vector<std::uint32_t>&,
                                       const std::size_t, const std::size_t, const std::size_t,
                                       const std::size_t, const std::vector<std::uint32_t>&,
                                       const std::vector<double>&, std::vector<double>&,
                                       DofNodes&);
template void build_facet_table<3>(const std::vector<std::uint32_t>&, const std::size_t,
                                   const std::vector<double>&, const Uint,
                                   std::vector<std::int32_t>&,
                                   std::vector<std::pair<Vector3d, std::size_t>>&);
template void build_facet_table<4>(const std::vector<std::uint32_t>&, const std::size_t,
                                   const std::vector<double>&, const Uint,
                                   std::vector<std::int32_t>&,
                                   std::vector<std::pair<Vector3d, std::size_t>>&);
template void build_facet_neighbours<3>(const std::vector<std::uint32_t>&, const std::size_t,
                                        const std::vector<double>&, const Uint,
                                        const std::vector<bool>&, const Vector3d&, const Vector3d&,
                                        const double, std::vector<std::int32_t>&);
template void build_facet_neighbours<4>(const std::vector<std::uint32_t>&, const std::size_t,
                                        const std::vector<double>&, const Uint,
                                        const std::vector<bool>&, const Vector3d&, const Vector3d&,
                                        const double, std::vector<std::int32_t>&);

}  // namespace mesh_tables
