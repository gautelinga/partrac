#ifndef __PARTICLESET_HPP
#define __PARTICLESET_HPP
#include "typedefs.hpp"
#include <string>
#include "Interpol.hpp"
#include "Integrator.hpp"
#include "io.hpp"

class ParticleSet {
public:
    ParticleSet(Interpol* intp, const Uint Nrw_max);
    void add(const std::vector<Vector3d> &pos_init, const Uint irw0);
    bool insert_node_between(const Uint, const Uint);
    double dist(const Uint inode, const Uint jnode) const { Vector3d dx = x_rw[inode]-x_rw[jnode]; return dx.norm(); };
    void copy_node(const Uint, const Uint);
    double triangle_area(const Uint iface, const FacesType& faces, const EdgesType& edges) const;
    double triangle_area(const Uint iedge, const Uint jedge, const EdgesType& edges) const { return cross_product(iedge, jedge, edges).norm()/2; };
    std::vector<double> triangle_areas(const std::vector<Uint> &kfaces, const FacesType &faces, const EdgesType &edges);
    Vector3d cross_product(const Uint iedge, const Uint jedge, const EdgesType& edges) const;
    void replace_nodes(Vector3d& x, const Uint inode, const Uint jnode);
    void collapse_nodes(const Uint inode, const Uint jnode, Node2EdgesType& node2edges);
    Vector3d x(const Uint i) const { return x_rw[i]; };
    Vector3d u(const Uint i) const { return u_rw[i]; };
    Vector3d facet_normal(const Uint iface,
                          const FacesType &faces,
                          const EdgesType &edges);
    void set_normals(const InteriorAnglesType& interior_angles, const std::vector<Vector3d> &face_normals);
    void compute_curvature(const EdgesType &edges, const Node2EdgesType &node2edges, std::vector<double>, std::vector<double>);
    void compute_strip_curvature(const EdgesType &edges, const Node2EdgesType &node2edges);
    bool has_space() const { return Nrw < Nrw_max; };
    bool has_space(const Uint n) const { return Nrw + n < Nrw_max; };
    Uint N() const { return Nrw; };
    void set_N(Uint n) { Nrw=n; };
    void load_scalar(const std::string filename, const std::string fieldname);
    void dump_scalar(const std::string filename, const std::string fieldname) const;
    void load_positions(const std::string filename);
    void dump_positions(const std::string filename) const;
    void dump_hdf5(H5FilePtr h5f, const std::string groupname, std::map<std::string, bool> &output_fields) const;
    bool integrate(const double t, const double dt);
    void attach_integrator(Integrator* integrator) { this->integrator = integrator; };
    Uint get_accepted() { return integrator->get_accepted(); };
    Uint get_declined() { return integrator->get_declined(); };
    void update_fields(const double t, std::map<std::string, bool> &output_fields);
  private:
    bool do_output_all = false;
    Uint Nrw = 0;
    Uint Nrw_max;
    int int_order;
    Interpol* intp;
    Integrator* integrator;
    // Vector fields
    std::vector<Vector3d> x_rw;
    std::vector<Vector3d> u_rw;
    std::vector<Vector3d> a_rw;  // for second order integration
    std::vector<Vector3d> n_rw;  // useless?
    // Scalar fields
    std::vector<double> c_rw;
    std::vector<double> e_rw;
    std::vector<double> H_rw;
    std::vector<double> rho_rw;
    std::vector<double> p_rw;
    std::vector<double> tau_rw;  // eigentime
};

ParticleSet::ParticleSet (Interpol* intp, const Uint Nrw_max) {
    this->intp = intp;
    this->Nrw_max = Nrw_max;
    this->int_order = intp->get_int_order();
    // Vector fields
    this->x_rw.resize(Nrw_max);
    this->u_rw.resize(Nrw_max);
    this->a_rw.resize(Nrw_max);  // for second order integration
    this->n_rw.resize(Nrw_max);  // useless?
    // Scalar fields
    this->c_rw.resize(Nrw_max);
    this->e_rw.resize(Nrw_max);
    this->H_rw.resize(Nrw_max);
    this->rho_rw.resize(Nrw_max);
    this->p_rw.resize(Nrw_max);
    this->tau_rw.resize(Nrw_max);  // eigentime
}

void ParticleSet::add(const std::vector<Vector3d> &pos_init, const Uint irw0) {
  for (Uint irw=irw0; irw < irw0+pos_init.size(); ++irw){
    // Assign initial position
    x_rw[irw] = pos_init[irw-irw0]; // could be done more efficiently

    // Does not work for injection:
    c_rw[irw] = double(irw)/(Nrw + pos_init.size() - 1);

    tau_rw[irw] = 0.;  // anything else?

    /*
    intp->probe(x_rw[irw]);
    u_rw[irw] = intp->get_u();

    rho_rw[irw] = intp->get_rho();
    p_rw[irw] = intp->get_p();

    // Second-order terms
    if (int_order >= 2){
      a_rw[irw] = intp->get_Ju() + intp->get_a();
    }
    */
  }
  Nrw += pos_init.size();
}

bool ParticleSet::insert_node_between(const Uint inode, const Uint jnode){
  Vector3d x_rw_new = 0.5*(x_rw[inode]+x_rw[jnode]);
  int refinement_insertion_levels = 10;
  intp->probe(x_rw_new);
  if (!intp->inside_domain()){
    std::cout << "Insertion failed! Need something more refined here." << std::endl;
    exit(0);
    Vector3d dx_rw_new = x_rw[inode]-x_rw[jnode];
    double dx0 = dx_rw_new.norm();
    Vector3d n0 = u_rw[inode]+u_rw[jnode];
    n0 /= -n0.norm();
    intp->probe(x_rw_new + dx0*n0);

    if (!intp->inside_domain()){
      std::cout << "Insertion failed! Information:" << std::endl;
      std::cout << n0 << std::endl;
      std::cout << dx0 << std::endl;
      std::cout << x_rw_new << std::endl;
      std::cout << x_rw_new + dx0*n0 << std::endl;
      exit(0);

      return false;
    }
    double ddx = dx0/2;
    double dx1 = dx0;
    for (int i=2; i<(2+refinement_insertion_levels); ++i){
      intp->probe(x_rw_new + ddx*n0);
      if (intp->inside_domain()){
        ddx -= dx0/pow(2, i);
        dx1 = ddx;
      }
      else {
        ddx += dx0/pow(2, i);
      }
    }
    x_rw_new += dx1*n0;
    intp->probe(x_rw_new);
  }

  x_rw[Nrw] = x_rw_new;

  c_rw[Nrw] = 0.5*(c_rw[inode]+c_rw[jnode]);
  tau_rw[Nrw] = 0.5*(tau_rw[inode]+tau_rw[jnode]);

  H_rw[Nrw] = 0.5*(H_rw[inode]+H_rw[jnode]);
  n_rw[Nrw] = 0.5*(n_rw[inode]+n_rw[jnode]);
  n_rw[Nrw] /= n_rw[Nrw].norm();

  /*
  // u_rw[Nrw] = intp->get_u();  // For some reason this goes wrong?
  u_rw[Nrw] = 0.5*(u_rw[inode]+u_rw[jnode]);
  if (do_output_all){
    // rho_rw[Nrw] = intp->get_rho();
    // p_rw[Nrw] = intp->get_p();
    rho_rw[Nrw] = 0.5*(rho_rw[inode]+rho_rw[jnode]);
    p_rw[Nrw] = 0.5*(p_rw[inode]+p_rw[jnode]);
  }
  // Second-order terms
  if (int_order >= 2){
    // a_rw[Nrw] = intp->get_Ju() + intp->get_a();
    a_rw[Nrw] = 0.5*(a_rw[inode] + a_rw[jnode]);
  }*/
  ++Nrw;
  return true;
}

void ParticleSet::copy_node(const Uint i, const Uint j){
  x_rw[i] = x_rw[j];
  c_rw[i] = c_rw[j];
  tau_rw[i] = tau_rw[j];  // if it is used?

  /*
  u_rw[i] = u_rw[j];

  if (do_output_all){
    rho_rw[i] = rho_rw[j];
    p_rw[i] = p_rw[j];
  }

  H_rw[i] = H_rw[j];
  n_rw[i] = n_rw[j];

  if (int_order >= 2){
    a_rw[i] = a_rw[j];
  }
  */
}

double ParticleSet::triangle_area(const Uint iface,
                                  const FacesType& faces, const EdgesType& edges) const {
  Uint iedge = faces[iface].first[0];
  Uint jedge = faces[iface].first[1];
  return triangle_area(iedge, jedge, edges);
}

Vector3d ParticleSet::cross_product(const Uint iedge, const Uint jedge, const EdgesType& edges) const {
  Vector3d a = x_rw[edges[iedge].first[0]]-x_rw[edges[iedge].first[1]];
  Vector3d b = x_rw[edges[jedge].first[0]]-x_rw[edges[jedge].first[1]];
  return a.cross(b);
}

std::vector<double> ParticleSet::triangle_areas(const std::vector<Uint> &kfaces,
                                                const FacesType &faces, const EdgesType &edges){
  std::vector<double> a;
  for (std::vector<Uint>::const_iterator faceit=kfaces.begin();
       faceit != kfaces.end(); ++faceit){
    a.push_back(this->triangle_area(*faceit, faces, edges));
  }
  return a;
}

void ParticleSet::replace_nodes(Vector3d& x, const Uint inode, const Uint jnode){
  intp->probe(x);  //

  Uint irws[2] = {inode, jnode};
  for (Uint i=0; i<2; ++i){
    Uint irw = irws[i];
    x_rw[irw] = x;
    c_rw[irw] = 0.5*(c_rw[inode]+c_rw[jnode]);
    tau_rw[irw] = 0.5*(tau_rw[inode]+tau_rw[jnode]);

    H_rw[irw] = 0.5*(H_rw[inode]+H_rw[jnode]);
    n_rw[irw] = 0.5*(n_rw[inode]+n_rw[jnode]);
    n_rw[irw] /= n_rw[irw].norm();

    /*
    u_rw[irw] = intp->get_u();

    if (do_output_all){
      rho_rw[irw] = intp->get_rho();
      p_rw[irw] = intp->get_p();
    }
    // Second-order terms
    if (int_order >= 2){
      a_rw[irw] = intp->get_Ju() + intp->get_a();
    }

    */
  }
}

void ParticleSet::collapse_nodes(const Uint inode, const Uint jnode, Node2EdgesType& node2edges){
    bool inode_is_border = node2edges[inode].size() > 1;
    bool jnode_is_border = node2edges[jnode].size() > 1;

    Uint new_inode = std::min(inode, jnode);
    Uint old_inode = std::max(inode, jnode);

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
    x_rw[inode] = pos_new;
    x_rw[jnode] = pos_new;

    c_rw[new_inode] = 0.5*(c_rw[inode]+c_rw[jnode]);
    tau_rw[new_inode] = 0.5*(tau_rw[inode]+tau_rw[jnode]);

    H_rw[new_inode] = 0.5*(H_rw[inode]+H_rw[jnode]);
    n_rw[new_inode] = 0.5*(n_rw[inode]+n_rw[jnode]);
    n_rw[new_inode] /= n_rw[new_inode].norm();
    
    /*
    intp->probe(x_rw[new_inode]);
    u_rw[new_inode] = intp->get_u();
    if (do_output_all){
        rho_rw[new_inode] = intp->get_rho();
        p_rw[new_inode] = intp->get_p();
    }

    if (int_order >= 2){
        a_rw[new_inode] = intp->get_Ju() + intp->get_a();
    }*/
}

Vector3d ParticleSet::facet_normal(const Uint iface,
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

void ParticleSet::set_normals(const InteriorAnglesType& interior_angles, const std::vector<Vector3d> &face_normals){
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

void ParticleSet::compute_curvature(const EdgesType &edges, const Node2EdgesType &node2edges,
                                    std::vector<double> edge_w, std::vector<double> mixed_areas){
  for (Uint inode=0; inode<Nrw; ++inode){
    Vector3d lapl_v(0., 0., 0.);
    for (EdgesListType::const_iterator edgeit=node2edges[inode].begin();
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

void ParticleSet::compute_strip_curvature(const EdgesType &edges,
                                          const Node2EdgesType &node2edges){
  for (Uint inode=0; inode<Nrw; ++inode){
    H_rw[inode] = 0.;
    n_rw[inode] = {1., 0., 0.};
    if (node2edges[inode].size() == 2){
      EdgesListType::const_iterator edgeit = node2edges[inode].begin();
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

void ParticleSet::load_scalar(const std::string filename, const std::string fieldname){
  if (fieldname == "c"){
    load_scalar_field(filename, c_rw, N());
  }
  else if (fieldname == "tau"){
    load_scalar_field(filename, tau_rw, N());
  }
}

void ParticleSet::dump_scalar(const std::string filename, const std::string fieldname) const {
  if (fieldname == "c"){
    dump_scalar_field(filename, c_rw, N());
  }
  else if (fieldname == "tau"){
    dump_scalar_field(filename, tau_rw, N());
  }
}

void ParticleSet::load_positions(const std::string filename){
    assert(N() == 0);
    std::vector<Vector3d> pos;
    load_vector_field(filename, pos);
    add(pos, 0);
}

void ParticleSet::dump_positions(const std::string filename) const {
    dump_vector_field(filename, x_rw, N());
}

void ParticleSet::update_fields(const double t, std::map<std::string, bool> &output_fields){
  for (Uint irw=0; irw < N(); ++irw){
    intp->probe(x_rw[irw], t);
    if (output_fields["u"])
      u_rw[irw] = intp->get_u();
    if (output_fields["rho"])
      rho_rw[irw] = intp->get_rho();
    if (output_fields["p"])
      p_rw[irw] = intp->get_p();
  }
}

void ParticleSet::dump_hdf5(H5FilePtr h5f, const std::string groupname, std::map<std::string, bool> &output_fields) const {
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
    if (output_fields["tau"])
        scalar2hdf5(h5f, groupname + "/tau", tau_rw, N());
}

bool ParticleSet::integrate(const double t, const double dt){
    Vector3d dx_rw = {0., 0., 0.};
    for (Uint irw=0; irw < N(); ++irw){
        //dx_rw = u_rw[irw]*dt;
        dx_rw = integrator->integrate(x_rw[irw], t, dt);

        intp->probe(x_rw[irw] + dx_rw, t + dt);
        assert (intp->inside_domain());
        x_rw[irw] += dx_rw;
        /*
        u_rw[irw] = intp->get_u();
        */

        /*if ((it+1) % int_dump_intv == 0){
            rho_rw[irw] = intp->get_rho();
            p_rw[irw] = intp->get_p();
        }*/

        /*
        // Second-order terms
        if (int_order >= 2){
            a_rw[irw] = intp->get_Ju() + intp->get_a();
        }
        //return true;
        */
    }
    return true;
}

#endif