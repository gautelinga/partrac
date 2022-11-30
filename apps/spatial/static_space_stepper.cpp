#include <iostream>
#include <vector>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <set>
#include <iterator>
#include "H5Cpp.h"
//#include "hdf5.h"
#include <ctime>

#include "io.hpp"
#include "utils.hpp"
#include "Parameters.hpp"

#include "ParticleSet.hpp"
#include "Topology.hpp"
#include "Integrator.hpp"

//#include "Integrator.hpp"
//#include "ExplicitIntegrator.hpp"
//#include "RKIntegrator.hpp"
#include "helpers.hpp"
#include "MPIwrap.hpp"

class Integrator_Spatial : public Integrator {
public:
  Integrator_Spatial(const int int_order, const double u_min, const double dl_max, const double T);
  ~Integrator_Spatial() {};
  template<typename InterpolType, typename T>
  std::set<Uint> step(InterpolType&, T&, double t, double s);
  std::set<Uint> step(ParticleSet&, double t, double s) { std::set<Uint> dummy; return dummy; };
protected:
  double   m_u_min;
  double   m_dl_max;
  int      m_int_order;
  double   m_T;
};


class Integrator_Directional : public Integrator_Spatial {
public:
  Integrator_Directional(const Vector3d& direction, const int int_order, const double un_min, const double dl_max, const double T);
  ~Integrator_Directional() {};
  template<typename InterpolType, typename T>
  std::set<Uint> step(InterpolType&, T&, double t, double s);
protected:
  Vector3d m_direction;
};

Integrator_Spatial::Integrator_Spatial(const int int_order, const double u_min, const double dl_max, const double T)
  : Integrator(), m_int_order(int_order), m_u_min(u_min), m_dl_max(dl_max), m_T(T) {
    std::cout << "Choosing a spatial integrator of order " << int_order << "." << std::endl;
}

template<typename InterpolType, typename T>
std::set<Uint> Integrator_Spatial::step(InterpolType& intp, T& ps, const double t, const double ds) {
    std::set<Uint> outside_nodes;
    bool is_inside;
    double uabs_est, dt;
    Vector3d dx;

    for (Uint i=0; i < ps.N(); ++i){
        Vector3d x = ps.x(i);
        int cell_id = ps.get_cell_id(i);

        intp.probe(x, t, cell_id);

        Vector3d u_1 = intp.get_u();

        uabs_est = u_1.norm();

        is_inside = false;
        if (uabs_est > m_u_min && ps.t_loc(i) < m_T){
            dt = ds / uabs_est;

            dx = u_1 * dt;

            // Second-order terms
            if (m_int_order >= 2){
                dx += 0.5 * (intp.get_a() + intp.get_Ju()) * dt * dt;
            }

            if (dx.norm() < m_dl_max){
                intp.probe(x + dx, t, cell_id);  // Frozen time, otherwise: intp.probe(x+dx, t+dt);
                is_inside = intp.inside_domain();
            }
            else {
                std::cout << "Step too long (dl=" << dx.norm() << "), consider doing something smart!" << std::endl;
            }
        }
        // count things
        if (is_inside){
            ++n_accepted;
            ps.set_x(i, x + dx);
            ps.set_t_loc(i, ps.t_loc(i) + dt);
            ps.set_cell_id(i, cell_id);
        }
        else {
            outside_nodes.insert(i);
            ++n_declined;
        }
    }
    return outside_nodes;
}

Integrator_Directional::Integrator_Directional(const Vector3d& direction, const int int_order, const double u_min, const double dl_max, const double T)
  : Integrator_Spatial(int_order, u_min, dl_max, T), m_direction(direction) {
    std::cout << "Choosing a directional integrator." << std::endl;
}

template<typename InterpolType, typename T>
std::set<Uint> Integrator_Directional::step(InterpolType& intp, T& ps, const double t, const double s) {
    std::set<Uint> outside_nodes;
    bool is_inside;
    double s_prev, un_est, dt;
    Vector3d dx;

    for (Uint i=0; i < ps.N(); ++i){
        Vector3d x = ps.x(i);
        int cell_id = ps.get_cell_id(i);

        intp.probe(x, t, cell_id);

        Vector3d u_1 = intp.get_u();

        un_est = u_1.dot(m_direction);

        is_inside = false;
        if (un_est > m_u_min && ps.t_loc(i) < m_T){
            s_prev = x.dot(m_direction);
            dt = (s - s_prev) / un_est;
            dx = u_1 * dt;

            // Second-order terms
            if (m_int_order >= 2){
                dx += 0.5 * (intp.get_a() + intp.get_Ju()) * dt * dt;
            }

            if (dx.norm() < m_dl_max){
                intp.probe(x + dx, t, cell_id);  // Frozen time, otherwise: intp.probe(x+dx, t+dt);
                is_inside = intp.inside_domain();
            }
            else {
                std::cout << "Step too long (dl=" << dx.norm() << "), consider doing something smart!" << std::endl;
            }
        }
        // count things
        if (is_inside){
            ++n_accepted;
            ps.set_x(i, x + dx);
            ps.set_t_loc(i, ps.t_loc(i) + dt);
            ps.set_cell_id(i, cell_id);
        }
        else {
            outside_nodes.insert(i);
            ++n_declined;
        }
    }
    return outside_nodes;
}

int main(int argc, char* argv[])
{
  MPIwrap mpi(argc, argv);

  if (mpi.rank() == 0)
    std::cout << "Initialized spatial stepper with " << mpi.size() << " processes." << std::endl;
  mpi.barrier();
  std::cout << "This is process " << mpi.rank() << " out of " << mpi.size() << "." << std::endl;
  mpi.barrier();

  // Input parameters
  if (argc < 2 && mpi.rank() == 0) {
    std::cout << "Specify an input file." << std::endl;
    return 0;
  }
  Parameters prm(argc, argv);
  if (prm.restart_folder != ""){
    prm.parse_file(prm.restart_folder + "/Checkpoints/params.dat");
    prm.parse_cmd(argc, argv);
  }

  std::string infilename = std::string(argv[1]);

  std::cout << "Setting interpolator..." << std::endl;

  std::shared_ptr<Interpol> intp;
  set_interpolate_mode(intp, prm.mode, infilename);
  
  intp->set_U0(prm.U0);
  intp->set_int_order(prm.int_order);

  bool refine = prm.refine;
  bool coarsen = prm.coarsen;

  std::cout << "Creating folders..." << std::endl;

  std::string folder = intp->get_folder();
  std::string rwfolder = folder + "/StaticSpaceStepper/"; 
  if (mpi.rank() == 0)
    create_folder(rwfolder);
  std::string newfolder;
  if (prm.restart_folder != ""){
    newfolder = prm.folder;
  }
  else {
    newfolder = get_newfoldername(rwfolder, prm);
    mpi.barrier();
    if (mpi.rank() == 0)
      create_folder(newfolder);
    mpi.barrier();
  }
  newfolder = newfolder + "" + std::to_string(mpi.rank()) + "/";
  std::string posfolder = newfolder + "Positions/";
  std::string checkpointsfolder = newfolder + "Checkpoints/";
  //std::string histfolder = newfolder + "Histograms/";
  //if (mpi.rank() == 0){ // Might change in the future!
  {
    create_folder(newfolder);
    create_folder(posfolder);
    create_folder(checkpointsfolder);
    //create_folder(histfolder);
  }
  prm.folder = newfolder;

  if (mpi.rank() == 0)
    prm.print();

  std::mt19937 gen;
  if (prm.random) {
    std::random_device rd;
    gen.seed(rd());
  }
  else {
    std::seed_seq rd{prm.seed + mpi.rank()};
    gen.seed(rd);
  }

  // TODO: These should not be stored in particle tracker parameters.
  prm.Lx = intp->get_Lx();
  prm.Ly = intp->get_Ly();
  prm.Lz = intp->get_Lz();

  double t0 = std::max(intp->get_t_min(), prm.t0);
  prm.t0 = t0;

  // Higher-order time integration?
  if (prm.int_order > 2){
    if (mpi.rank() == 0)
      std::cout << "No support for such high temporal integration order." << std::endl;
    exit(0);
  }

  intp->update(t0);

  //Vector3d direction = {1., 0., 0.};

  //std::shared_ptr<Integrator> integrator;
  //integrator = std::make_shared<DirectionalIntegrator>(direction, prm.int_order);
  Integrator_Spatial integrator(prm.int_order, prm.u_eps, prm.dx_max, prm.T);

  ParticleSet ps(intp, prm.Nrw_max, mpi);
  Topology mesh(ps, prm, mpi);

  if (prm.restart_folder != ""){
    mesh.load_checkpoint(prm.restart_folder + "/Checkpoints", prm);
  }
  else {
    std::shared_ptr<Initializer> init_state;
    set_initial_state(init_state, intp, mpi, prm, gen);
    mesh.load_initial_state(init_state);
  }

  mesh.compute_maps();

  // Initial refinement
  if (refine && !prm.inject && mesh.dim() > 0){
    std::cout << "Initial refinement" << std::endl;
    Uint n_add = mesh.refine();

    std::cout << "Initial coarsening" << std::endl;
    Uint n_rem = mesh.coarsen();

    if (prm.verbose && mpi.rank() == 0)
      std::cout << "Added " << n_add << " edges and removed " << n_rem << " edges." << std::endl;
  }

  mesh.compute_interior();

  int it = 0;

  double xn = prm.x0;
  double dxn = prm.dxn;

  //if (mpi.rank() == 0)
  prm.dump(newfolder, xn);

  // Should not be taken from parameters
  //Uint n_accepted = prm.n_accepted;
  //Uint n_declined = prm.n_declined;

  std::string h5fname = newfolder + "/data_from_t" + std::to_string(xn) + ".h5";
  H5::H5File h5f(h5fname.c_str(), H5F_ACC_TRUNC);
  //h5f->openFile(h5fname.c_str(), H5F_ACC_TRUNC);
  //H5wrap h5file(mpi);
  //h5file.open(h5fname, "w");

  Uint int_stat_intv = int(prm.stat_intv/dxn);
  Uint int_dump_intv = int(prm.dump_intv/dxn);
  Uint int_checkpoint_intv = int(prm.checkpoint_intv/dxn);
  Uint int_chunk_intv = int_dump_intv*prm.dump_chunk_size;
  Uint int_refine_intv = int(prm.refine_intv/dxn);
  Uint int_coarsen_intv = int(prm.coarsen_intv/dxn);

  std::map<std::string, bool> output_fields;
  output_fields["u"] = !prm.minimal_output;
  output_fields["c"] = true;
  output_fields["p"] = !prm.minimal_output && prm.output_all_props;
  output_fields["rho"] = !prm.minimal_output && prm.output_all_props;        
  output_fields["H"] = !prm.minimal_output && mesh.dim() > 0;
  output_fields["n"] = !prm.minimal_output && mesh.dim() > 1;
  output_fields["t_loc"] = true;
  output_fields["tau"] = true;

  //std::string write_mode = prm.write_mode;

  std::ofstream statfile;
  //if (mpi.rank() == 0){
  {
    statfile.open(newfolder + "/tdata_from_t" + std::to_string(xn) + ".dat");
    write_stats_header(mpi, statfile, mesh.dim());
  }
  std::ofstream declinedfile(newfolder + "/declinedpos_from_t" + std::to_string(xn) + ".dat");

  // Simulation start
  std::clock_t clock_0 = std::clock();
  while (xn <= prm.Ln){
    // Statistics
    if (it % int_stat_intv == 0){
      std::cout << "Position = " << xn << std::endl;
      mesh.write_statistics(statfile, xn, prm.ds_max, integrator);
    }
    // Checkpoint
    if (it % int_checkpoint_intv == 0){
      //std::cout << "Writing checkpoint..." << std::endl;
      mesh.write_checkpoint(checkpointsfolder, xn, prm);
      //std::cout << "Done." << std::endl;
    }
    // Curvature computation
    if ((refine && it % int_refine_intv == 0) || (coarsen && it % int_coarsen_intv == 0) || it % int_dump_intv == 0){
      mesh.compute_interior();
    }

    // Refinement
    if (refine && it % int_refine_intv == 0 && it > 0){
      Uint n_add = mesh.refine();
      if (prm.verbose)
        std::cout << "Added " << n_add << " edges." << std::endl;
    }
    // Coarsening
    if (coarsen && it % int_coarsen_intv == 0){
      Uint n_rem = mesh.coarsen();
      if (prm.verbose)
        std::cout << "Removed " << n_rem << " edges." << std::endl;
    }

    // Dump detailed data
    if (it % int_dump_intv == 0){
      std::cout << "Dumping..." << std::endl;
      ps.update_fields(t0, output_fields);

      std::string groupname = std::to_string(xn);
      // Clear file if it exists, otherwise create
      if (int_chunk_intv > 0 && it % int_chunk_intv == 0 && it > 0){
        h5fname = newfolder + "/data_from_t" + std::to_string(xn) + ".h5";
        h5f.openFile(h5fname.c_str(), H5F_ACC_TRUNC);
      }
      else {
        h5f.openFile(h5fname.c_str(), H5F_ACC_RDWR);
      }
      h5f.createGroup(groupname + "/");
      mesh.dump_hdf5(h5f, groupname, output_fields);
      h5f.close();
    }

    xn += dxn;

    auto nodes_to_remove = integrator.step(*intp, ps, t0, dxn);

    if (nodes_to_remove.size() > 0){
      std::vector<bool> node_isactive(ps.N(), true);
      for (auto sit = nodes_to_remove.begin();
            sit != nodes_to_remove.end(); ++sit){
        node_isactive[*sit] = false;
      }
      mesh.remove_nodes_safe(node_isactive);
    }

    it += 1;
  }
  std::clock_t clock_1 = std::clock();
  double duration = (clock_1-clock_0) / (double) CLOCKS_PER_SEC;
  std::cout << "Total simulation time: " << duration << " seconds" << std::endl;

  mesh.write_checkpoint(checkpointsfolder, xn, prm);
  statfile.close();
  declinedfile.close();

  return 0;
}
