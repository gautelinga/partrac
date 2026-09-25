#ifndef __OPENFOAM_PHASE_HPP
#define __OPENFOAM_PHASE_HPP

// A phase field on the split: the volume where its P1 field exceeds 1/2,
// per OpenFOAM cell, against the cell's volume times its volume fraction. A
// fanned cell's centre is read by no other cell, so its value is chosen to
// make the two agree; a cell without one keeps its error. Needs no OpenFOAM.

#include <cstdint>
#include <vector>

namespace openfoam_phase {

// The fraction of a simplex (nv 3 or 4 nodes) where the linear field of its
// node values v exceeds 1/2
double above_half(const double* v, int nv);

// The simplices of each OpenFOAM cell, and its centre node where fanned
struct Cells {
  Cells() = default;
  // cell_of: the OpenFOAM cell of each simplex; centre: a cell's centre node, -1 for none
  Cells(std::size_t ncells, const std::vector<std::int32_t>& cell_of, const std::vector<std::int32_t>& centre);
  std::vector<std::uint32_t> start;    // ncells + 1 into simplex
  std::vector<std::uint32_t> simplex;
  std::vector<std::int32_t> centre;
  std::size_t ncells() const { return centre.size(); }
};

// What a correction found and did; volumes in the split's measure (areas in 2D)
struct Report {
  double target = 0.;             // the sum of V_c clamp(alpha_c, 0, 1)
  double before = 0.;             // the enclosed volume with each centre its cell's value
  double after = 0.;              // with the centres corrected
  double worst_before = 0., worst_after = 0.;   // the largest cell error, of its volume
  std::int64_t worst_cell = -1;   // the OpenFOAM cell of worst_after
  // Fanned cells: within purity of 0 or 1; brought to the target; left at 0
  // or 1 short of it; a target of 0 or V_c no centre reaches
  std::size_t pure = 0, corrected = 0, clipped = 0, unreachable = 0;
  std::size_t unfanned = 0;       // cells without a centre, off by more than the tolerance
};

// A volume fraction this close to 0 or 1 is pure
constexpr double purity = 1e-6;

// Each fanned cell's centre value in values (by node) moved within [0, 1]
// until the volume where the P1 field exceeds 1/2 in the cell is
// V_c clamp(alpha_c, 0, 1), else to the end nearer it; a pure cell, or one
// whose target no centre reaches, keeps its centre. The simplices are topo's
// (nv a simplex) on coords (nv - 1 a node)
Report conserve(const Cells& cells, const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                int nv, const std::vector<double>& alpha, std::vector<double>& values);

// Each fanned cell's centre in grad (d a node) the mean of its other nodes'
void centre_means(const Cells& cells, const std::vector<std::uint32_t>& topo, int nv, int d,
                  std::vector<double>& grad);

}  // namespace openfoam_phase

#endif
