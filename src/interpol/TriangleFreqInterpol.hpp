#ifndef __TRIANGLEFREQINTERPOL_HPP
#define __TRIANGLEFREQINTERPOL_HPP

// A time series given as frequency components: one steady Taylor-Hood pair per
// component, each a dolfin HDF5 checkpoint on the same triangle mesh, summed
// with a cosine of the base frequency. The components are read the way
// SimplexInterpol reads a stamp; an evaluation needs every component at one
// cell, so each component's coefficients are held cell by cell.

#include <memory>
#include <string>
#include <vector>

#include "MeshCore.hpp"
#include "cell_tree.hpp"
#include "FreqStamps.hpp"
#include "strings.hpp"
#include "Triangle.hpp"
#include <omp.h>

class TriangleFreqInterpol final
  : public MeshCore<Triangle>
{
public:
  TriangleFreqInterpol(const std::string& infilename);
  void update(const double t);
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  double get_t_min() { return dolfin_params.get<double>("t_min"); };
  double get_t_max() { return dolfin_params.get<double>("t_max"); };
protected:
  // Out of the step loops: the cell tree, from nothing known
  bool locate_tree(const Vector3d& xx, CellPos& pos) override;

  FreqStamps fs; // frequencies holder

  // Per component, the cells' coefficients in a row each: the velocity's D
  // components blocked as the basis reads them
  std::vector<std::vector<double>> u_coefficients_;
  std::vector<std::vector<double>> p_coefficients_;

  double omega0 = 0.;

  // The arrays the tree addresses; the cells hold the geometry
  std::vector<std::uint32_t> topo_;
  std::vector<double> coords_;
  std::size_t ncells_ = 0, nverts_ = 0;
  std::unique_ptr<partrac::CellTree> tree_;
};

#endif
