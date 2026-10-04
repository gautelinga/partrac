#ifndef __PARTICLESET_HPP
#define __PARTICLESET_HPP
#include <cstdint>
#include <memory>
#include <utility>
#include "Error.hpp"
#include "geometry.hpp"
#include "strings.hpp"
#include "typedefs.hpp"
#include <string>
#include "Interpol.hpp"
//#include "Integrator.hpp"  // remove!
#include "io.hpp"
#include "TransportElement.hpp"
#include <random>
#include <numeric>
#include <algorithm>


// The cells of a mesh interpolator, 0 for the rest (stepping.cpp)
Uint cell_count(Interpol& intp);

class ParticleSet {
public:
    //ParticleSet(std::shared_ptr<Interpol> intp, std::shared_ptr<Integrator> integrator, const Uint Nrw_max);
    ParticleSet(std::shared_ptr<Interpol> intp, const Uint Nrw_max);
    //ParticleSet(const Uint Nrw_max);
    void add(const std::vector<Vector3d> &pos_init, const Uint irw0);
    void add(const std::vector<Vector3d> &pos_init, const std::vector<Uint> &which, const Uint irw0);
    //template<typename T>
    bool insert_node_between(const Uint, const Uint, const bool check_if_inside);
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
    void dump_hdf5(H5::H5File& h5f, const std::string& groupname, std::map<std::string, bool> &output_fields) const;
    //bool integrate(const double t, const double dt);
    // void attach_integrator(std::shared_ptr<Integrator> integrator) { this->integrator = integrator; };
    //Uint get_accepted() { return integrator->get_accepted(); };
    //Uint get_declined() { return integrator->get_declined(); };
    // Fields at the particles (typed interpolator)
    template<typename Interp>
    void update_fields(Interp&, const double, std::map<std::string, bool>&);
    //void reduce(ParticleSet& psb, std::map<std::string, bool> &output_fields);
    std::shared_ptr<Interpol>& interpolator() { return intp; };
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
    //bool do_output_all = false;
    Uint Nrw = 0;
    Uint Nrw_max;
    std::shared_ptr<Interpol> intp;
    //std::shared_ptr<Integrator> integrator;
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
    // 
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

/*
ParticleSet::ParticleSet (std::shared_ptr<Interpol> intp, std::shared_ptr<Integrator> integrator, const Uint Nrw_max) : ParticleSet(intp, Nrw_max) {
    this->integrator = integrator;
}*/

inline ParticleSet::ParticleSet(std::shared_ptr<Interpol> intp, const Uint Nrw_max) {
    this->intp = intp;
    this->Nrw_max = Nrw_max;
    // Vector fields
    this->x_rw.resize(Nrw_max);
    this->u_rw.resize(Nrw_max);
    this->n_rw.resize(Nrw_max);
    // Scalar fields
    this->c_rw.resize(Nrw_max);
    this->H_rw.resize(Nrw_max);
    this->rho_rw.resize(Nrw_max);
    this->p_rw.resize(Nrw_max);
    this->t_loc_rw.resize(Nrw_max);  // eigentime
    //
    this->cell_id_rw.resize(Nrw_max);
    this->id_rw.resize(Nrw_max);
}

inline void ParticleSet::add(const std::vector<Vector3d> &pos_init, const Uint irw0) {
  for (Uint irw=irw0; irw < irw0+pos_init.size(); ++irw){
    // Assign initial position
    x_rw[irw] = pos_init[irw-irw0]; // could be done more efficiently

    // Does not work for injection:
    c_rw[irw] = double(irw)/(Nrw + pos_init.size() - 1);

    t_loc_rw[irw] = 0.;  // anything else?

    cell_id_rw[irw] = -1;
    init_carried(irw);
  }
  Nrw += pos_init.size();
}

// Add only the listed entries of pos_init
inline void ParticleSet::add(const std::vector<Vector3d> &pos_init,
                             const std::vector<Uint> &which, const Uint irw0) {
  for (Uint k=0; k < which.size(); ++k){
    const Uint irw = irw0 + k;
    x_rw[irw] = pos_init[which[k]];
    c_rw[irw] = double(irw)/(Nrw + which.size() - 1);
    t_loc_rw[irw] = 0.;
    cell_id_rw[irw] = -1;
    init_carried(irw);
  }
  Nrw += which.size();
}

//template<typename T>
inline bool ParticleSet::insert_node_between(const Uint inode, const Uint jnode, const bool check_if_inside=true){
  Vector3d x_rw_new = 0.5*(x_rw[inode]+x_rw[jnode]);
  int cell_id = cell_id_rw[inode];
  
  if (check_if_inside){
    double t0 = 0.; // not needed?
    
    bool inside = intp->locate(x_rw_new, t0, cell_id);
    if (!inside){
      std::cout << "Insertion failed! Need something more refined here." << std::endl;
      //return false;
      //exit(1);

      Vector3d dx_rw_new = x_rw[inode]-x_rw[jnode];
      double dx0 = dx_rw_new.norm();

      // tangent vector
      Vector3d tau0 = dx_rw_new / dx0;

      //Vector3d n0 = u_rw[inode]+u_rw[jnode];
      Vector3d n0 = intp->get_boundary_normal(x_rw[inode], cell_id_rw[inode]) + intp->get_boundary_normal(x_rw[jnode], cell_id_rw[jnode]);

      // Check that normal is valid.
      if (n0.norm() < 1e-2){
        return false;
      }

      n0 -= n0.dot(tau0) * tau0;
      n0 /= -n0.norm();
      
      double ddx = 1e-2 * dx0;
      double dx1 = 0;

      for (Uint iddx=1; iddx < 1000; ++iddx){
        dx1 = iddx * ddx;
        inside = intp->locate(x_rw_new + dx1 * n0, t0, cell_id);
        if (inside){
          break;
        }
      }

      if (!inside){
        return false;
      }


      x_rw_new += dx1*n0;
      // inside = intp->locate(x_rw_new, t0, cell_id);

    }

  }

  x_rw[Nrw] = x_rw_new;
  // The cell located, else unknown
  cell_id_rw[Nrw] = check_if_inside ? cell_id : -1;

  c_rw[Nrw] = 0.5*(c_rw[inode]+c_rw[jnode]);
  t_loc_rw[Nrw] = 0.5*(t_loc_rw[inode]+t_loc_rw[jnode]);

  H_rw[Nrw] = 0.5*(H_rw[inode]+H_rw[jnode]);
  n_rw[Nrw] = 0.5*(n_rw[inode]+n_rw[jnode]);
  n_rw[Nrw] /= n_rw[Nrw].norm();

  u_rw[Nrw] = {0., 0., 0.}; // 0.5*(u_rw[inode]+u_rw[jnode]);
  // Interpolate carried fields
  interpolate_carried(Nrw, inode, jnode);
  id_rw[Nrw] = next_id++;

  ++Nrw;
  return true;
}

inline void ParticleSet::compact(const std::vector<Uint>& kept, const Uint from){
  for_each_array([&](auto& a){
    for (Uint i = from; i < kept.size(); ++i)
      a[i] = a[kept[i]];
  });
}

inline double ParticleSet::triangle_area(const Uint iface,
                                         const FacesType& faces, const EdgesType& edges) const {
  Uint iedge = faces[iface].first[0];
  Uint jedge = faces[iface].first[1];
  return triangle_area(iedge, jedge, edges);
}

inline Vector3d ParticleSet::cross_product(const Uint iedge, const Uint jedge, const EdgesType& edges) const {
  Vector3d a = x_rw[edges[iedge].first[0]]-x_rw[edges[iedge].first[1]];
  Vector3d b = x_rw[edges[jedge].first[0]]-x_rw[edges[jedge].first[1]];
  return a.cross(b);
}

inline void ParticleSet::replace_nodes(Vector3d& x, const Uint inode, const Uint jnode){
  //intp->probe(x);  //

  Uint irws[2] = {inode, jnode};
  for (Uint i=0; i<2; ++i){
    Uint irw = irws[i];
    move(irw, x);
    c_rw[irw] = 0.5*(c_rw[inode]+c_rw[jnode]);
    t_loc_rw[irw] = 0.5*(t_loc_rw[inode]+t_loc_rw[jnode]);

    H_rw[irw] = 0.5*(H_rw[inode]+H_rw[jnode]);
    n_rw[irw] = 0.5*(n_rw[inode]+n_rw[jnode]);
    n_rw[irw] /= n_rw[irw].norm();
  }
  // Carried fields
  interpolate_carried(inode, inode, jnode);
  interpolate_carried(jnode, inode, jnode);
}

inline void ParticleSet::collapse_nodes(const Uint inode, const Uint jnode, Node2EdgesType& node2edges){
    bool inode_is_border = node2edges[inode].size() > 1;
    bool jnode_is_border = node2edges[jnode].size() > 1;

    Uint new_inode = std::min(inode, jnode);

    // can be made simpler!!

    Vector3d pos_new;
    if ((inode_is_border && jnode_is_border) || (!inode_is_border && !jnode_is_border)){
        pos_new = 0.5*(x_rw[inode] + x_rw[jnode]);
    }
    else if (inode_is_border){
        pos_new = x_rw[inode];
    }
    else {
        pos_new = x_rw[jnode];
    }
    move(inode, pos_new);
    move(jnode, pos_new);

    c_rw[new_inode] = 0.5*(c_rw[inode]+c_rw[jnode]);
    t_loc_rw[new_inode] = 0.5*(t_loc_rw[inode]+t_loc_rw[jnode]);

    H_rw[new_inode] = 0.5*(H_rw[inode]+H_rw[jnode]);
    n_rw[new_inode] = 0.5*(n_rw[inode]+n_rw[jnode]);
    n_rw[new_inode] /= n_rw[new_inode].norm();
    interpolate_carried(new_inode, inode, jnode);
}

inline Vector3d ParticleSet::facet_normal(const Uint iface,
                                          const FacesType &faces,
                                          const EdgesType &edges){
  Uint iedge = faces[iface].first[0];
  Uint jedge = faces[iface].first[1];
  Uint i00 = edges[iedge].first[0];
  Uint i01 = edges[iedge].first[1];
  Uint i10 = edges[jedge].first[0];
  Uint i11 = edges[jedge].first[1];

  Vector3d a = x_rw[i01]-x_rw[i00];
  Vector3d b = x_rw[i11]-x_rw[i10];
  Vector3d crossprod = a.cross(b);
  return crossprod/crossprod.norm();
}

inline void ParticleSet::set_normals(const InteriorAnglesType& interior_angles, const std::vector<Vector3d> &face_normals){
  for (Uint irw=0; irw<Nrw; ++irw){
    n_rw[irw] = {0., 0., 0.};
  }
  for (Uint iface=0; iface < interior_angles.size(); ++iface){
    std::map<Uint, double> angles = interior_angles[iface];
    for (std::map<Uint, double>::const_iterator angit=angles.begin();
         angit != angles.end(); ++angit){
      n_rw[angit->first] += angit->second * face_normals[iface];
    }
  }
  for (Uint irw=0; irw<Nrw; ++irw){
    n_rw[irw] /= n_rw[irw].norm();
  }
}

inline void ParticleSet::compute_curvature(const EdgesType &edges, const Node2EdgesType &node2edges,
                                           std::vector<double> edge_w, std::vector<double> mixed_areas){
  for (Uint inode=0; inode<Nrw; ++inode){
    Vector3d lapl_v(0., 0., 0.);
    for (auto edgeit=node2edges[inode].begin();
         edgeit != node2edges[inode].end(); ++edgeit){
      Uint jnode = get_other(edges[*edgeit].first[0],
                             edges[*edgeit].first[1], inode);
      Vector3d dv = x_rw[jnode]-x_rw[inode];
      lapl_v += edge_w[*edgeit]*dv;
    }
    lapl_v /= 2*mixed_areas[inode];
    H_rw[inode] = 0.5*lapl_v.dot(n_rw[inode]);
  }
}

inline void ParticleSet::compute_strip_curvature(const EdgesType &edges,
                                                 const Node2EdgesType &node2edges){
  for (Uint inode=0; inode<Nrw; ++inode){
    H_rw[inode] = 0.;
    n_rw[inode] = {1., 0., 0.};
    if (node2edges[inode].size() == 2){
      auto edgeit = node2edges[inode].begin();
      Uint jnode = get_other(edges[*edgeit].first[0],
                             edges[*edgeit].first[1], inode);
      ++edgeit;
      Uint knode = get_other(edges[*edgeit].first[0],
                             edges[*edgeit].first[1], inode);

      double R = circumcenter(x_rw[inode], x_rw[jnode], x_rw[knode]);

      H_rw[inode] = 1./R;

      //Vector3d n_loc = dij/dij.norm() + dik/dik.norm();
      //if (n_loc.norm() > 0)
      //  n_rw[inode] = n_loc/n_loc.norm();
    }
  }
}

inline void ParticleSet::load_scalar(const std::string filename, const std::string fieldname){
  if (fieldname == "c"){
    load_scalar_field(filename, c_rw, N());
  }
  else if (fieldname == "t_loc"){
    load_scalar_field(filename, t_loc_rw, N());
  }
  else if (fieldname == "w"){
    load_scalar_field(filename, w_rw, N());
  }
  else if (fieldname == "S"){
    load_scalar_field(filename, S_rw, N());
  }
  else if (fieldname == "generation"){
    load_scalar_field(filename, generation_rw, N());
  }
  else {
    // Unknown field
    partrac::fail("ParticleSet::load_scalar: no field '", fieldname, "'");
  }
}

inline void ParticleSet::load_positions(const std::string filename){
    assert(N() == 0);
    std::vector<Vector3d> pos;
    load_vector_field(filename, pos);
    if (pos.size() > Nrw_max)
      partrac::fail(filename, ": ", pos.size(), " particles, more than Nrw_max = ", Nrw_max);
    add(pos, 0);
}


inline void ParticleSet::dump_hdf5(H5::H5File& h5f, const std::string& groupname, std::map<std::string, bool> &output_fields) const {
    vector2hdf5(h5f, groupname + "/points", x_rw, N());
    if (output_fields["u"])
        vector2hdf5(h5f, groupname + "/u", u_rw, N());
    if (output_fields["rho"])
        scalar2hdf5(h5f, groupname + "/rho", rho_rw, N());
    if (output_fields["p"])
        scalar2hdf5(h5f, groupname + "/p", p_rw, N());
    if (output_fields["c"])
        scalar2hdf5(h5f, groupname + "/c", c_rw, N());
    if (output_fields["H"])
        scalar2hdf5(h5f, groupname + "/H", H_rw, N());
    if (output_fields["n"])
        vector2hdf5(h5f, groupname + "/n", n_rw, N());
    //if (faces.size() == 0)
    //  scalar2hdf5(h5f, groupname + "/e", ps.e_rw, ps.Nrw);
    if (output_fields["t_loc"])
        scalar2hdf5(h5f, groupname + "/" + t_loc_name, t_loc_rw, N());
    if (element == TransportElement::Vector){
        vector2hdf5(h5f, groupname + "/" + rhohat_name, rhohat_rw, N());
        scalar2hdf5(h5f, groupname + "/w", w_rw, N());
        scalar2hdf5(h5f, groupname + "/S", S_rw, N());
    }
    if (element == TransportElement::Tensor){
        std::vector<Matrix3d> Q(N());
        std::vector<Vector3d> ls(N()), U(N());
        for (Uint i = 0; i < N(); ++i) settled(i, Q[i], ls[i], U[i]);
        tensor2hdf5(h5f, groupname + "/Q", Q, N());
        vector2hdf5(h5f, groupname + "/logstretch", ls, N());
        vector2hdf5(h5f, groupname + "/U", U, N());
        if (output_fields["S"])
            vector2hdf5(h5f, groupname + "/S", S3_rw, N());
    }
    if (has_J && output_fields["J"])
        tensor2hdf5(h5f, groupname + "/J", J_rw, N());
    if (has_phi && output_fields["phi"])
        scalar2hdf5(h5f, groupname + "/phi", phi_rw, N());
    if (has_cell_type && output_fields["cell_type"])
        int2hdf5(h5f, groupname + "/cell_type", cell_type_rw, N());
    if (has_generation)
        scalar2hdf5(h5f, groupname + "/w", generation_rw, N());
    ulong2hdf5(h5f, groupname + "/id", id_rw, N());
}

inline void ParticleSet::dump_as(const std::string& field, const std::string& name){
  if (field == "rhohat") rhohat_name = name;
  else if (field == "t_loc") t_loc_name = name;
  else { partrac::fail("ParticleSet::dump_as: no field '", field, "'"); }
}

inline void ParticleSet::record_generation(){
  // Both dumped as w
  if (element == TransportElement::Vector){
    partrac::fail("ParticleSet: a line element's w and a walker's generation cannot share a dump");
  }
  generation_rw.resize(Nrw_max);
  has_generation = true;
}

inline void ParticleSet::carry(const TransportElement e){
  if (e == TransportElement::Vector && has_generation){
    partrac::fail("ParticleSet: a line element's w and a walker's generation cannot share a dump");
  }
  element = e;
  if (e == TransportElement::Vector){
    rhohat_rw.resize(Nrw_max);
    w_rw.resize(Nrw_max);
    S_rw.resize(Nrw_max);
  }
  if (e == TransportElement::Tensor){
    Q_rw.resize(Nrw_max);
    logstretch_rw.resize(Nrw_max);
    U_rw.resize(Nrw_max);
    S3_rw.resize(Nrw_max);
  }
  for (Uint i = 0; i < Nrw; ++i) init_carried_fields(i);
}

// Apply fn to every particle array; unused arrays are skipped
template<typename Fn>
inline void ParticleSet::for_each_array(Fn&& fn){
  fn(x_rw); fn(u_rw); fn(n_rw);
  fn(c_rw); fn(H_rw); fn(rho_rw); fn(p_rw); fn(t_loc_rw);
  fn(cell_id_rw); fn(id_rw);
  if (!rhohat_rw.empty()){ fn(rhohat_rw); fn(w_rw); fn(S_rw); }
  if (!Q_rw.empty()){ fn(Q_rw); fn(logstretch_rw); fn(U_rw); fn(S3_rw); }
  if (has_J) fn(J_rw);
  if (has_phi) fn(phi_rw);
  if (has_cell_type) fn(cell_type_rw);
  if (has_generation) fn(generation_rw);
}

inline void ParticleSet::init_carried(const Uint irw){
  id_rw[irw] = next_id++;
  init_carried_fields(irw);
}

inline void ParticleSet::init_carried_fields(const Uint irw){
  if (element == TransportElement::Vector){
    rhohat_rw[irw] = {0., 0., 0.};
    w_rw[irw] = 0.;
    S_rw[irw] = 0.;
  }
  if (element == TransportElement::Tensor){
    Q_rw[irw] = Matrix3d::Identity();
    logstretch_rw[irw] = Vector3d::Zero();
    U_rw[irw] = Vector3d::Zero();
    S3_rw[irw] = Vector3d::Zero();
  }
  if (has_J) J_rw[irw] = Matrix3d::Zero();
  if (has_phi) phi_rw[irw] = 0.;
  if (has_cell_type) cell_type_rw[irw] = 0;
  if (has_generation) generation_rw[irw] = 0.;
}

// Interpolate carried fields
inline void ParticleSet::interpolate_carried(const Uint k, const Uint inode, const Uint jnode){
  if (has_J) J_rw[k] = 0.5*(J_rw[inode] + J_rw[jnode]);
  if (has_phi) phi_rw[k] = 0.5*(phi_rw[inode] + phi_rw[jnode]);
  if (has_cell_type) cell_type_rw[k] = cell_type_rw[inode];   // label
  if (has_generation) generation_rw[k] = std::max(generation_rw[inode], generation_rw[jnode]);   // count
}

inline void ParticleSet::spin_rhohat(const std::string& dirs, std::mt19937& gen){
  std::normal_distribution<double> rnd_normal(0.0, 1.0);
  const bool dx = contains(dirs, "x"), dy = contains(dirs, "y"), dz = contains(dirs, "z");
  for (Uint irw = 0; irw < Nrw; ++irw){
    Vector3d r = {0., 0., 0.};
    if (dx) r[0] = rnd_normal(gen);
    if (dy) r[1] = rnd_normal(gen);
    if (dz) r[2] = rnd_normal(gen);
    rhohat_rw[irw] = r/r.norm();
  }
}

inline std::vector<Uint> ParticleSet::sort_by_cell(){
  // No cells
  if (Nrw == 0 || cell_id_rw[0] < 0)
    return {};
  std::vector<std::pair<int, Uint>> keys(Nrw);
  #pragma omp parallel for
  for (Uint i = 0; i < Nrw; ++i)
    keys[i] = {cell_id_rw[i], i};
  std::sort(keys.begin(), keys.end());
  std::vector<Uint> order(Nrw);
  for (Uint k = 0; k < Nrw; ++k)
    order[k] = keys[k].second;
  return reorder(order);
}

inline std::vector<Uint> ParticleSet::shuffle(std::mt19937& gen){
  std::vector<Uint> order(Nrw);
  std::iota(order.begin(), order.end(), 0);
  std::shuffle(order.begin(), order.end(), gen);
  return reorder(order);
}

inline std::vector<Uint> ParticleSet::reorder(const std::vector<Uint>& order){
  std::vector<Uint> old2new(Nrw);
  for (Uint k = 0; k < Nrw; ++k) old2new[order[k]] = k;
  for_each_array([&](auto& a){
    using T = typename std::remove_reference_t<decltype(a)>::value_type;
    std::unique_ptr<T[]> tmp(new T[Nrw]);
    #pragma omp parallel for
    for (Uint k = 0; k < Nrw; ++k) tmp[k] = a[order[k]];
    #pragma omp parallel for
    for (Uint k = 0; k < Nrw; ++k) a[k] = tmp[k];
  });
  return old2new;
}

inline void ParticleSet::load_vector(const std::string filename, const std::string fieldname){
  if (fieldname == "rhohat") load_vector_field(filename, rhohat_rw, N());
  else if (fieldname == "logstretch") load_vector_field(filename, logstretch_rw, N());
  else if (fieldname == "U") load_vector_field(filename, U_rw, N());
  else { partrac::fail("ParticleSet::load_vector: no field '", fieldname, "'"); }
}
// F whole, as a checkpoint from before the factors has it: factored on load
inline void ParticleSet::load_tensor(const std::string filename, const std::string fieldname){
  if (fieldname == "Q") load_tensor_field(filename, Q_rw, N());
  else if (fieldname == "F"){
    std::vector<Matrix3d> F_whole(N());
    load_tensor_field(filename, F_whole, N());
    for (Uint i = 0; i < N(); ++i) factor_F(i, F_whole[i]);
  }
  else { partrac::fail("ParticleSet::load_tensor: no field '", fieldname, "'"); }
}
inline Matrix3d ParticleSet::F(const Uint i) const {
  Matrix3d Uf = Matrix3d::Identity();
  Uf(0, 1) = U_rw[i][0]; Uf(0, 2) = U_rw[i][1]; Uf(1, 2) = U_rw[i][2];
  return Q_rw[i] * logstretch_rw[i].array().exp().matrix().asDiagonal() * Uf;
}

// Gram-Schmidt: G = Q R, R's diagonal positive
inline void ParticleSet::factor_F(const Uint i, const Matrix3d& G){
  Vector3d q0 = G.col(0), q1 = G.col(1), q2 = G.col(2);
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
  Q_rw[i] << q0, q1, q2;
  logstretch_rw[i] = {log(r00), log(r11), log(r22)};
  U_rw[i] = {r01/r00, r02/r00, r12/r11};
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

// Missing ids: identity
inline void ParticleSet::load_ids(const std::string filename){
  std::vector<Uint> ids;
  std::ifstream infile(filename);
  if (infile) load_list(filename, ids);
  if (ids.size() != N()) { ids.resize(N()); std::iota(ids.begin(), ids.end(), 0); }
  for (Uint i = 0; i < N(); ++i) id_rw[i] = ids[i];
  next_id = ids.empty() ? 0 : *std::max_element(ids.begin(), ids.end()) + 1;
}

// The frame and factors as carried
inline void ParticleSet::write_checkpoint(H5::H5File& h5f, const bool with_t_loc) const {
  vector2hdf5(h5f, "points", x_rw, N());
  ulong2hdf5(h5f, "id", id_rw, N());
  scalar2hdf5(h5f, "c", c_rw, N());
  if (with_t_loc)
    scalar2hdf5(h5f, "t_loc", t_loc_rw, N());
  if (element == TransportElement::Vector){
    vector2hdf5(h5f, "rhohat", rhohat_rw, N());
    scalar2hdf5(h5f, "w", w_rw, N());
    scalar2hdf5(h5f, "S", S_rw, N());
  }
  if (element == TransportElement::Tensor){
    tensor2hdf5(h5f, "Q", Q_rw, N());
    vector2hdf5(h5f, "logstretch", logstretch_rw, N());
    vector2hdf5(h5f, "U", U_rw, N());
  }
  if (has_generation)
    scalar2hdf5(h5f, "generation", generation_rw, N());
  int2hdf5(h5f, "cell_id", cell_id_rw, N());
  const std::uint64_t next = next_id;
  h5f.createAttribute("next_id", H5::PredType::NATIVE_UINT64, H5::DataSpace(H5S_SCALAR))
     .write(H5::PredType::NATIVE_UINT64, &next);
}

inline void ParticleSet::read_checkpoint(const H5::H5File& h5f, const bool with_t_loc){
  assert(N() == 0);
  const Uint n = hdf5_rows(h5f, "points");
  if (n > Nrw_max)
    partrac::fail(h5f.getFileName(), ": ", n, " particles, more than Nrw_max = ", Nrw_max);
  std::vector<Vector3d> pos;
  hdf52vectors(h5f, "points", n, pos);
  add(pos, 0);
  // Scalar fields
  std::vector<double> a;
  auto scalars = [&](const std::string& name, std::vector<double>& v){
    hdf52doubles(h5f, name, n, 1, a);
    std::copy(a.begin(), a.end(), v.begin());
  };
  scalars("c", c_rw);
  if (with_t_loc)
    scalars("t_loc", t_loc_rw);
  if (element == TransportElement::Vector){
    hdf52vectors(h5f, "rhohat", n, rhohat_rw);
    scalars("w", w_rw);
    scalars("S", S_rw);
  }
  if (element == TransportElement::Tensor){
    hdf52tensors(h5f, "Q", n, Q_rw);
    hdf52vectors(h5f, "logstretch", n, logstretch_rw);
    hdf52vectors(h5f, "U", n, U_rw);
  }
  if (has_generation)
    scalars("generation", generation_rw);
  std::vector<Uint> ids;
  hdf52ulongs(h5f, "id", n, 1, ids);
  std::copy(ids.begin(), ids.end(), id_rw.begin());
  next_id = ids.empty() ? 0 : *std::max_element(ids.begin(), ids.end()) + 1;
  // Ids of removed nodes stay unused; none saved in older checkpoints
  if (H5Aexists(h5f.getId(), "next_id") > 0){
    std::uint64_t next = 0;
    try { h5f.openAttribute("next_id").read(H5::PredType::NATIVE_UINT64, &next); }
    catch (const H5::Exception&){ partrac::fail(h5f.getFileName(), ": cannot read 'next_id'"); }
    next_id = std::max<Uint>(next_id, next);
  }
  // Cells of this mesh only; none in older checkpoints
  if (h5f.nameExists("cell_id")){
    std::vector<int> cell_ids;
    hdf52ints(h5f, "cell_id", n, 1, cell_ids);
    const Uint nc = intp ? cell_count(*intp) : 0;
    for (Uint i = 0; i < n; ++i)
      cell_id_rw[i] = cell_ids[i] >= 0 && Uint(cell_ids[i]) < nc ? cell_ids[i] : -1;
  }
}

// bool ParticleSet::integrate(const double t, const double dt){
//     Vector3d dx_rw = {0., 0., 0.};
//     bool all_nodes_are_fine = true;
//     for (Uint irw=0; irw < N(); ++irw){
//         //dx_rw = u_rw[irw]*dt;
//         dx_rw = integrator->integrate(x_rw[irw], t, dt);

//         //intp->probe(x_rw[irw] + dx_rw, t + dt);
//         //assert (intp->inside_domain());
//         x_rw[irw] += dx_rw;

//         if (integrator->is_stuck()){
//           all_nodes_are_fine = false;

//         }
//         /*
//         u_rw[irw] = intp->get_u();
//         */

//         /*if ((it+1) % int_dump_intv == 0){
//             rho_rw[irw] = intp->get_rho();
//             p_rw[irw] = intp->get_p();
//         }*/

//         /*
//         // Second-order terms
//         if (int_order >= 2){
//             a_rw[irw] = intp->get_Ju() + intp->get_a();
//         }
//         //return true;
//         */
//     }
//     return all_nodes_are_fine;
// }

// void ParticleSet::reduce(ParticleSet& psb, std::map<std::string, bool> &output_fields){
//   std::vector<double> xvals_(3*psb.N());
//   for (Uint i=0; i<psb.N(); ++i){
//     for (int d=0; d<3; ++d){
//       xvals_[3*i+d] = psb.x_rw[i][d];
//     }
//   }
//   auto all_xvals_ = gather_vector<double>(m_mpi, xvals_, MPI_DOUBLE);
//   Nrw = all_xvals_.size()/3;
//   for (Uint i=0; i<Nrw; ++i){
//     for (int d=0; d<3; ++d){
//       x_rw[i][d] = all_xvals_[3*i+d];
//     }
//   }
// }

#endif
