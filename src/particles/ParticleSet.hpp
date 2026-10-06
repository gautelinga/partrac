#ifndef __PARTICLESET_HPP
#define __PARTICLESET_HPP
#include <memory>
#include "geometry.hpp"
#include "typedefs.hpp"
#include <string>
#include "Interpol.hpp"
#include "TransportElement.hpp"
#include <random>

namespace H5 { class H5File; }

// Datasets a dump writes besides points and id; those of a field not recorded are skipped
struct OutputFields {
  bool u = false, c = false, p = false, rho = false;
  bool H = false, n = false, tau = false, t_loc = false;
  bool phi = false, J = false, S = false, cell_type = false;
};

// The cells of a mesh interpolator, 0 for the rest (stepping.cpp)
Uint cell_count(Interpol& intp);

class ParticleSet {
public:
    ParticleSet(std::shared_ptr<Interpol> intp, const Uint Nrw_max);
    void add(const std::vector<Vector3d> &pos_init, const Uint irw0);
    void add(const std::vector<Vector3d> &pos_init, const std::vector<Uint> &which, const Uint irw0);
    bool insert_node_between(const Uint, const Uint, const bool check_if_inside=true);
    double dist(const Uint inode, const Uint jnode) const { Vector3d dx = x_rw[inode]-x_rw[jnode]; return dx.norm(); };
    // Keep the given slots, in order, from slot `from` on
    void compact(const std::vector<Uint>& kept, const Uint from);
    double triangle_area(const Uint iface, const FacesType& faces, const EdgesType& edges) const;
    double triangle_area(const Uint iedge, const Uint jedge, const EdgesType& edges) const { return cross_product(iedge, jedge, edges).norm()/2; };
    Vector3d cross_product(const Uint iedge, const Uint jedge, const EdgesType& edges) const;
    void replace_nodes(Vector3d& x, const Uint inode, const Uint jnode);
    void collapse_nodes(const Uint inode, const Uint jnode, Node2EdgesType& node2edges);
    Vector3d x(const Uint i) const { return x_rw[i]; };
    void set_x(const Uint i, const Vector3d& pos) { x_rw[i] = pos; };
    // Moved: its cell unknown
    void move(const Uint i, const Vector3d& pos) { x_rw[i] = pos; cell_id_rw[i] = -1; };
    double t_loc(const Uint i) const { return t_loc_rw[i]; };
    void set_t_loc(const Uint i, const double val) { t_loc_rw[i] = val; };
    Vector3d u(const Uint i) const { return u_rw[i]; };
    void set_c(const Uint i, const double val) { c_rw[i] = val; };
    Vector3d facet_normal(const Uint iface,
                          const FacesType &faces,
                          const EdgesType &edges);
    void set_normals(const InteriorAnglesType& interior_angles, const std::vector<Vector3d> &face_normals);
    void compute_curvature(const EdgesType &edges, const Node2EdgesType &node2edges, std::vector<double>, std::vector<double>);
    void compute_strip_curvature(const EdgesType &edges, const Node2EdgesType &node2edges);
    bool has_space() const { return Nrw < Nrw_max; };
    bool has_space(const Uint n) const { return Nrw + n <= Nrw_max; };
    Uint N() const { return Nrw; };
    void set_N(Uint n) { Nrw=n; };
    // Old text checkpoints
    void load_scalar(const std::string filename, const std::string fieldname);
    void load_positions(const std::string filename);
    void dump_hdf5(H5::H5File& h5f, const std::string& groupname, const OutputFields& output_fields) const;
    // Fields at the particles (typed interpolator)
    template<typename Interp>
    void update_fields(Interp&, const double, const OutputFields&);
    int get_cell_id(const Uint irw) const { return cell_id_rw[irw]; };
    void set_cell_id(const Uint irw, const int cell_id) { cell_id_rw[irw] = cell_id; };
    // Carried and recorded fields, sized on request
    void carry(const TransportElement e);
    TransportElement carries() const { return element; };
    void record_J() { J_rw.resize(Nrw_max); has_J = true; };
    void record_phi() { phi_rw.resize(Nrw_max); has_phi = true; };
    void record_cell_type() { cell_type_rw.resize(Nrw_max); has_cell_type = true; };
    // Walker generation (dumped as w)
    void record_generation();
    bool records_generation() const { return has_generation; };
    double generation(const Uint i) const { return generation_rw[i]; };
    void set_generation(const Uint i, const double g) { generation_rw[i] = g; };
    // Dump a field under another name
    void dump_as(const std::string& field, const std::string& name);
    Vector3d rhohat(const Uint i) const { return rhohat_rw[i]; };
    void set_rhohat(const Uint i, const Vector3d& r) { rhohat_rw[i] = r; };
    double w(const Uint i) const { return w_rw[i]; };
    void set_w(const Uint i, const double v) { w_rw[i] = v; };
    double S(const Uint i) const { return S_rw[i]; };
    void set_S(const Uint i, const double v) { S_rw[i] = v; };
    // The deformation's stretching rates along its settled frame, diag(Q^T J Q)
    Vector3d S3(const Uint i) const { return S3_rw[i]; };
    // Deformation gradient F = Q diag(exp(logstretch)) U, U unit upper triangular:
    // the stretches as logs, the frame orthonormal
    Matrix3d F(const Uint i) const;
    void set_F(const Uint i, const Matrix3d& F) { factor_F(i, F); };
    Matrix3d frame(const Uint i) const { return Q_rw[i]; };
    // The frame after a step, M = (propagator) Q: kept while well conditioned,
    // else factored into F's factors
    void advance_frame(const Uint i, const Matrix3d& M);
    // The factors with the frame's growth folded in, as dumped
    void settled(const Uint i, Matrix3d& Q, Vector3d& s, Vector3d& u) const;
    Vector3d logstretch(const Uint i) const { Matrix3d Q; Vector3d s, u; settled(i, Q, s, u); return s; };
    double phi(const Uint i) const { return phi_rw[i]; };
    // Particle id
    Uint id(const Uint i) const { return id_rw[i]; };
    // Random unit rho-hat in the named directions
    void spin_rhohat(const std::string& dirs, std::mt19937& gen);
    // Sort by cell; returns old -> new, empty if unsorted
    std::vector<Uint> sort_by_cell();
    // Shuffle slots; returns old -> new
    std::vector<Uint> shuffle(std::mt19937& gen);
    void load_vector(const std::string filename, const std::string fieldname);
    void load_tensor(const std::string filename, const std::string fieldname);
    void load_ids(const std::string filename);
    // Checkpoint: positions, ids and next_id, the carried fields and cell ids; t_loc if asked
    void write_checkpoint(H5::H5File& h5f, const bool with_t_loc) const;
    void read_checkpoint(const H5::H5File& h5f, const bool with_t_loc);
  private:
    Uint Nrw = 0;
    Uint Nrw_max;
    std::shared_ptr<Interpol> intp;
    // Vector fields
    std::vector<Vector3d> x_rw;
    std::vector<Vector3d> u_rw;
    std::vector<Vector3d> n_rw;  // the sheet normal
    // Scalar fields
    std::vector<double> c_rw;
    std::vector<double> H_rw;
    std::vector<double> rho_rw;
    std::vector<double> p_rw;
    std::vector<double> t_loc_rw;  // eigentime
    std::vector<int> cell_id_rw; // for speed (if applicable)
    std::vector<Uint> id_rw;
    Uint next_id = 0;
    // Carried fields, sized by carry()
    TransportElement element = TransportElement::Point;
    std::vector<Vector3d> rhohat_rw;   // material line element, unit
    std::vector<double> w_rw;          // its log stretching
    std::vector<double> S_rw;          // its stretching rate rhohat^T J rhohat
    std::vector<Matrix3d> Q_rw;        // deformation gradient: its frame, orthonormal once settled,
    std::vector<Vector3d> logstretch_rw;   // the log of its triangular factor's diagonal
    std::vector<Vector3d> U_rw;        // and the unit triangular factor's entries 01, 02, 12
    std::vector<Vector3d> S3_rw;       // its stretching rates along the settled frame
    // F = G, factored into the fields of slot i
    void factor_F(const Uint i, const Matrix3d& G);
    // Recorded fields, sized by record_*()
    bool has_J = false, has_phi = false, has_cell_type = false, has_generation = false;
    std::vector<Matrix3d> J_rw;
    std::vector<double> phi_rw;
    std::vector<int> cell_type_rw;     // cell marker (XDMF)
    std::vector<double> generation_rw; // walker splits
    std::string rhohat_name = "rhohat";
    std::string t_loc_name = "t_loc";
    // Apply fn to every particle array
    template<typename Fn> void for_each_array(Fn&& fn);
    void init_carried(const Uint irw);
    void init_carried_fields(const Uint irw);
    // Reorder all arrays (slot k gets old slot order[k]); returns old -> new
    std::vector<Uint> reorder(const std::vector<Uint>& order);
    void interpolate_carried(const Uint k, const Uint inode, const Uint jnode);
};

// Per face: dumps, statistics, tau, remeshing
inline Vector3d ParticleSet::cross_product(const Uint iedge, const Uint jedge, const EdgesType& edges) const {
  Vector3d a = x_rw[edges[iedge].first[0]]-x_rw[edges[iedge].first[1]];
  Vector3d b = x_rw[edges[jedge].first[0]]-x_rw[edges[jedge].first[1]];
  return a.cross(b);
}

inline double ParticleSet::triangle_area(const Uint iface, const FacesType& faces, const EdgesType& edges) const {
  Uint iedge = faces[iface].first[0];
  Uint jedge = faces[iface].first[1];
  return triangle_area(iedge, jedge, edges);
}

inline Matrix3d ParticleSet::F(const Uint i) const {
  Matrix3d Uf = Matrix3d::Identity();
  Uf(0, 1) = U_rw[i][0]; Uf(0, 2) = U_rw[i][1]; Uf(1, 2) = U_rw[i][2];
  return Q_rw[i] * logstretch_rw[i].array().exp().matrix().asDiagonal() * Uf;
}

// F = M D U = Q' R' D U for any frame M, and R' D U = D' U' with
// D' = diag(r'_ii d_i), U'_ij = sum_k r'_ik d_k U_kj / (r'_ii d_i): the old
// stretches enter only as exp(s_k - s_i), k > i
inline void fold_frame(const Matrix3d& M, const Vector3d& s, const Vector3d& u,
                       Matrix3d& Q, Vector3d& s_out, Vector3d& u_out){
  Vector3d q0 = M.col(0), q1 = M.col(1), q2 = M.col(2);
  const double r00 = q0.norm();
  q0 /= r00;
  const double r01 = q0.dot(q1);
  q1 -= r01*q0;
  const double r11 = q1.norm();
  q1 /= r11;
  const double r02 = q0.dot(q2);
  q2 -= r02*q0;
  const double r12 = q1.dot(q2);
  q2 -= r12*q1;
  const double r22 = q2.norm();
  q2 /= r22;
  // Zero stays zero: no 0*inf
  auto term = [](const double a, const double ds){ return a == 0. ? 0. : a*exp(ds); };
  const double t01 = term(r01/r00, s[1] - s[0]), t02 = term(r02/r00, s[2] - s[0]), t12 = term(r12/r11, s[2] - s[1]);
  Q << q0, q1, q2;
  u_out = {u[0] + t01, u[1] + t01*u[2] + t02, u[2] + t12};
  s_out = {s[0] + log(r00), s[1] + log(r11), s[2] + log(r22)};
}

inline void ParticleSet::settled(const Uint i, Matrix3d& Q, Vector3d& s, Vector3d& u) const {
  fold_frame(Q_rw[i], logstretch_rw[i], U_rw[i], Q, s, u);
}

// Folded only when a column's squared length leaves [1/16, 16] or two columns
// come within 60 degrees: a well-conditioned frame carries F as it is
inline void ParticleSet::advance_frame(const Uint i, const Matrix3d& M){
  const Vector3d n = M.colwise().squaredNorm();
  const double c01 = M.col(0).dot(M.col(1)), c02 = M.col(0).dot(M.col(2)), c12 = M.col(1).dot(M.col(2));
  if (n.minCoeff() > 1./16 && n.maxCoeff() < 16.
      && 4*c01*c01 < n[0]*n[1] && 4*c02*c02 < n[0]*n[2] && 4*c12*c12 < n[1]*n[2]){
    Q_rw[i] = M;
    return;
  }
  fold_frame(M, logstretch_rw[i], U_rw[i], Q_rw[i], logstretch_rw[i], U_rw[i]);
}

#endif
