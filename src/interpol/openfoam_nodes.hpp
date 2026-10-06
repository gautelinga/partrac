#ifndef __OPENFOAM_NODES_HPP
#define __OPENFOAM_NODES_HPP

// Node values of an OpenFOAM field on its split, as sparse operators over
// [cell values; boundary face values] fixed at load. InverseDistance is
// volPointInterpolation's, so cellPoint's; LeastSquares a weighted linear
// fit over a point's cells (the next ring where they are nearly coplanar,
// mirrored at symmetry planes, across cyclics), or the value of a face whose
// condition fixes it: exact for linear fields. Needs no OpenFOAM.

#include <array>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "openfoam_load.hpp"
#include "openfoam_split.hpp"

namespace openfoam_nodes {

// Whether a patch's values are data for the reconstruction: its condition
// fixes the value and is none of those that carry the cell's value
bool is_data(const openfoam_load::FieldPatch& p);

// What a boundary patch is to the tracers, by the velocity's condition
enum class PatchClass : std::int8_t { wall, moving_wall, cyclic, other };

// Per patch: a no-slip wall where the condition is data and zero on every
// face, a moving wall where it is data, nonzero and along every face, cyclic,
// else other (inlets, outlets, symmetry, empty); zero to 1e-12 of the field's
// largest value
std::vector<PatchClass> classify_patches(const openfoam_load::CaseData& c,
                                         const openfoam_load::FieldData& u);
const char* class_name(PatchClass k);

// The cells every mesh point of which lies on a no-slip wall
std::size_t walled_cells(const openfoam_load::CaseData& c, const std::vector<PatchClass>& k);

// A case whose every cell is walled, as one cell thick between two walls:
// refused under split=6, where every node is then fixed at zero; a warning
// into log under split=12, whose cell centres alone carry the flow. Some
// cells walled: a warning with their count under either split
void check_walled(const openfoam_load::CaseData& c, const std::vector<PatchClass>& k, int tets_per_hex,
                  std::ostream& log);

// What the least-squares operator shares between fields of one mesh; the
// case must outlive it
struct Geometry {
  explicit Geometry(const openfoam_load::CaseData& c);
  const openfoam_load::CaseData* c = nullptr;
  std::vector<std::int32_t> master;           // a point's lowest cyclic image
  std::vector<double> shift;                  // 3 a point: x + shift = x[master]
  // CSR: a point's cells, its boundary faces (global ids), a master's images
  std::vector<std::size_t> pc_start, pc, pb_start, pb, img_start, img;
  // CSR: a cell's face neighbours, across cyclic faces too, with the shift
  // that carries the neighbour's centre to this side
  std::vector<std::size_t> link_start, link;
  std::vector<double> link_shift;
  std::vector<std::int32_t> bface_patch;
  std::vector<double> bface_normal;           // unit, 3 a boundary face
  std::vector<int> axes;                      // the fit's axes: all three, or the in-plane two
  double size = 0.;                           // the bounding box's diagonal
};

// Per mesh point, what the least-squares build did
struct LeastSquaresReport {
  std::vector<std::int8_t> ring;              // 0 fixed, 1 or 2 fitted, -1 deficient at the last, -2 no node
  std::vector<double> cond;                   // of the unit-diagonal normal matrix
  std::vector<std::int32_t> zero_wins;        // masters where a zero patch won
};

class LeastSquares {
public:
  static constexpr double default_rank_tol = 1e-3;
  LeastSquares() = default;
  // From the first stamp's field: its conditions, and its fixed patches'
  // values to tell uniform from varying and whether zero wins at a node
  LeastSquares(const Geometry& g, const openfoam_split::SplitData& s, const openfoam_load::FieldData& f,
               LeastSquaresReport* report = nullptr, double rank_tol = default_rank_tol);
  // The split's node values, the components comps of each, node after node;
  // a field whose patches are data where the first stamp's were not, or the
  // other way, is refused
  void apply(const openfoam_load::FieldData& f, const std::vector<int>& comps,
             std::vector<double>& out) const;
  std::size_t nnodes() const { return kind_.size(); }
  std::size_t nrows() const { return start_.empty() ? 0 : start_.size() - 1; }
  std::size_t nonzeros() const { return col_.size(); }
  std::size_t bytes() const;
  // What the build found, over the rows
  std::size_t fixed = 0, second_ring = 0, deficient = 0, zero_won = 0, mirrored = 0;
private:
  std::string field_;
  int ncomp_ = 1;
  std::size_t ncells_ = 0, n_internal_ = 0;
  // Per patch at load: whether data, its condition, its faces
  std::vector<char> data_;
  std::vector<std::string> name_, condition_;
  std::vector<std::size_t> patch_start_, patch_size_;
  // Node n: a centre (its cell), or a point whose row is ref_[n]
  std::vector<std::int8_t> kind_;
  std::vector<std::uint32_t> ref_;
  // Rows over the data [cell values; boundary face values], face f at ncells + f - n_internal
  std::vector<std::size_t> start_;
  std::vector<std::uint32_t> col_;
  std::vector<double> coef_;
  // A mirrored entry's reflection of a vector, into mirror_table_; empty without any
  std::vector<std::uint32_t> mirror_;
  std::vector<std::array<double, 9>> mirror_table_;
};

class InverseDistance {
public:
  InverseDistance() = default;
  InverseDistance(const Geometry& g, const openfoam_split::SplitData& s);
  // The split's node values, the components comps of each, node after node
  void apply(const openfoam_load::FieldData& f, const std::vector<int>& comps,
             std::vector<double>& out) const;
  std::size_t nnodes() const { return kind_.size(); }
  std::size_t bytes() const;
private:
  // Node n: a centre (its cell), or a point whose row of weights is row_of_[n]
  std::vector<std::int8_t> kind_;
  std::vector<std::size_t> ref_;
  // Rows over the data [cell values; boundary face values]: sources, weights, sum
  std::vector<std::size_t> start_;
  std::vector<std::size_t> src_;
  std::vector<double> w_;
  std::vector<double> den_;
  std::vector<std::int32_t> constraint_;     // per row, into t_, -1 for none
  std::vector<std::array<double, 9>> t_;
  // What a boundary face's value is: its patch's, or its owner's
  std::size_t ncells_ = 0, n_internal_ = 0;
  std::vector<std::int32_t> bface_patch_, bface_owner_;
  std::vector<std::array<double, 3>> bface_normal_;
  std::vector<char> symmetry_;               // per patch
  std::vector<std::size_t> patch_start_;
};

// The gradient of a scalar field at each point: G of the weighted linear fit
// LeastSquares makes at a point it does not fix, over the same cells, rings,
// mirrors and cyclic images, and never a face value, wherever the point is
class Gradient {
public:
  Gradient() = default;
  Gradient(const Geometry& g, const openfoam_split::SplitData& s,
           double rank_tol = LeastSquares::default_rank_tol);
  // The fit's axes (Geometry::axes) a node, node after node; a centre's zero
  void apply(const openfoam_load::FieldData& f, std::vector<double>& out) const;
  std::size_t nnodes() const { return ref_.size(); }
  std::size_t nrows() const { return start_.empty() ? 0 : start_.size() - 1; }
  std::size_t nonzeros() const { return col_.size(); }
  std::size_t bytes() const;
  // What the build found, over the rows
  std::size_t second_ring = 0, deficient = 0, mirrored = 0;
private:
  std::size_t d_ = 3, ncells_ = 0;
  std::vector<std::int32_t> ref_;             // a point's row, -1 for a centre
  std::vector<std::size_t> start_;
  std::vector<std::uint32_t> col_;            // cells
  std::vector<double> coef_;                  // d_ an entry
};

}  // namespace openfoam_nodes

#endif
