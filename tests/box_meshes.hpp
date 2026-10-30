// Structured simplex meshes of the unit square and cube, gdim doubles a point
// as the loaders store them; point (i, j[, k]) of the (n + 1)-point lattice is
// numbered (i m + j) m + k, m = n + 1
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

// n x n x n cubes, each cut into the six Kuhn tets (the tetrahedra 0 -> e_a ->
// e_a+e_b -> 1 over the six axis orders), which agree on every shared face, so
// the mesh is conforming and its outer faces pair up across the box
inline void kuhn_tets(const std::size_t n, std::vector<double>& coords, std::vector<std::uint32_t>& topo){
  const std::size_t m = n + 1;
  const auto vid = [m](const std::size_t i, const std::size_t j, const std::size_t k){
    return std::uint32_t((i*m + j)*m + k);
  };
  coords.assign(3*m*m*m, 0.);
  topo.clear();
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < m; ++j)
      for (std::size_t k = 0; k < m; ++k){
        double* c = coords.data() + 3*std::size_t(vid(i, j, k));
        c[0] = double(i)/double(n);
        c[1] = double(j)/double(n);
        c[2] = double(k)/double(n);
      }
  const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t k = 0; k < n; ++k)
        for (const auto& p : perm){
          std::array<std::size_t, 3> c = {i, j, k};
          topo.push_back(vid(c[0], c[1], c[2]));
          for (const int a : p){
            ++c[std::size_t(a)];
            topo.push_back(vid(c[0], c[1], c[2]));
          }
        }
}

// n x n squares, each cut along its (0,0)-(1,1) diagonal
inline void diagonal_triangles(const std::size_t n, std::vector<double>& coords, std::vector<std::uint32_t>& topo){
  const std::size_t m = n + 1;
  const auto vid = [m](const std::size_t i, const std::size_t j){ return std::uint32_t(i*m + j); };
  coords.assign(2*m*m, 0.);
  topo.clear();
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < m; ++j){
      double* c = coords.data() + 2*std::size_t(vid(i, j));
      c[0] = double(i)/double(n);
      c[1] = double(j)/double(n);
    }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j){
      const std::uint32_t v00 = vid(i, j), v10 = vid(i+1, j), v11 = vid(i+1, j+1), v01 = vid(i, j+1);
      topo.insert(topo.end(), {v00, v10, v11});
      topo.insert(topo.end(), {v00, v11, v01});
    }
}
