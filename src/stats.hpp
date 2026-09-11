#ifndef __STATS_HPP
#define __STATS_HPP

#include "utils.hpp"
#include "stats_columns.hpp"

inline std::vector<StatsColumn> stats_columns(
                 const double t,
                 const ParticleSet& ps,
                 const FacesType &faces,
                 const EdgesType &edges,
                 const double ds_max,
                 //const bool do_dump_hist,
                 //const std::string histfolder,
                 const unsigned long int n_accepted,
                 const unsigned long int n_declined,
                 const Uint mesh_dim)
{
  std::vector<StatsColumn> cols;
  Vector3d x_mean = {0., 0., 0.};
  Vector3d dx2_mean = {0., 0., 0.};
  Vector3d u_mean = {0., 0., 0.};
  Uint Nrw = ps.N();
  for (Uint irw=0; irw < Nrw; ++irw){
    // Sample mean
    x_mean += ps.x(irw); // /Nrw;
    u_mean += ps.u(irw); // /Nrw;
  }
  // the mean of nothing, and the spread of a single point, are both
  // undefined; report zero rather than a NaN that spreads downstream
  if (Nrw > 0){
    x_mean /= Nrw;
    u_mean /= Nrw;
  }

  if (Nrw > 1){
    for (Uint irw=0; irw < Nrw; ++irw){
      // Sample variance
      Vector3d dx = ps.x(irw)-x_mean;
      dx2_mean += dx.cwiseProduct(dx)/(Nrw-1);
    }
  }

  cols.push_back({"t", t});
  cols.push_back({"x_mean", x_mean[0]});
  cols.push_back({"dx2_mean", dx2_mean[0]});
  cols.push_back({"y_mean", x_mean[1]});
  cols.push_back({"dy2_mean", dx2_mean[1]});
  cols.push_back({"z_mean", x_mean[2]});
  cols.push_back({"dz2_mean", dx2_mean[2]});
  cols.push_back({"ux_mean", u_mean[0]});
  cols.push_back({"uy_mean", u_mean[1]});
  cols.push_back({"uz_mean", u_mean[2]});
  cols.push_back({"Nrw", double(Nrw), true});
  cols.push_back({"n_accepted", double(n_accepted), true});
  cols.push_back({"n_declined", double(n_declined), true});

  double s = 0.;
  double s0 = 0.;
  double A = 0.;
  double A0 = 0.;

  // Which columns this run writes is fixed by the dimension it settles into,
  // not by what the mesh happens to be right now: the header is written once,
  // before the loop, and an injecting run is still its inlet at that point.
  // A sheet that has not been swept yet reports zero area rather than the
  // length of the curve about to sweep it.
  const bool do_strip = mesh_dim == 1;
  const bool do_sheet = mesh_dim > 1;

  if (do_strip || do_sheet){
    Uint n_too_long = 0;
    for (EdgesType::const_iterator edgeit = edges.begin();
         edgeit != edges.end(); ++edgeit)
      if (ps.dist(edgeit->first[0], edgeit->first[1]) > ds_max)
        ++n_too_long;
    cols.push_back({"n_too_long", double(n_too_long), true});
  }

  if (do_strip){
    // Strip method
    double logelong_wmean = 0.;
    double logelong_w0mean = 0.;

    std::vector<std::array<double, 3>> logelong_vec;
    for (EdgesType::const_iterator edgeit = edges.begin();
         edgeit != edges.end(); ++edgeit){
      int inode = edgeit->first[0];
      int jnode = edgeit->first[1];
      double ds0 = edgeit->second;
      double ds = ps.dist(inode, jnode);
      double logelong = log(ds/ds0);
      logelong_wmean += logelong*ds;
      logelong_w0mean += logelong*ds0;
      logelong_vec.push_back({logelong, ds, ds0});
      s += ds;
      s0 += ds0;
    }
    logelong_wmean = s > 0. ? logelong_wmean/s : 0.;
    logelong_w0mean = s0 > 0. ? logelong_w0mean/s0 : 0.;

    double logelong_wvar = 0.;
    double logelong_w0var = 0.;
    for (std::vector<std::array<double, 3>>::const_iterator lit = logelong_vec.begin();
         lit != logelong_vec.end(); ++lit){
      logelong_wvar += pow((*lit)[0]-logelong_wmean, 2)*(*lit)[1];
      logelong_w0var += pow((*lit)[0]-logelong_w0mean, 2)*(*lit)[2];
    }
    logelong_wvar = s > 0. ? logelong_wvar/s : 0.;
    logelong_w0var = s0 > 0. ? logelong_w0var/s0 : 0.;
    cols.push_back({"s", s});
    cols.push_back({"s0", s0});
    cols.push_back({"logelong_wmean", logelong_wmean});
    cols.push_back({"logelong_wvar", logelong_wvar});
    cols.push_back({"logelong_w0mean", logelong_w0mean});
    cols.push_back({"logelong_w0var", logelong_w0var});
  }
  else if (do_sheet){
    // Sheet method
    double logelong_wmean = 0.;
    double logelong_w0mean = 0.;
    std::vector<std::array<double, 3>> logelong_vec;
    for (FacesType::const_iterator faceit = faces.begin();
         faceit != faces.end(); ++faceit){
      Uint iedge = faceit->first[0];
      Uint jedge = faceit->first[1];
      // Uint kedge = faceit->first[2];
      double dA0 = faceit->second;
      if (!(dA0 > 0.))
        continue;                 // a flat sweep, waiting to be culled
      double dA = ps.triangle_area(iedge, jedge, edges);
      double logelong = log(dA/dA0);
      logelong_wmean += logelong*dA;
      logelong_w0mean += logelong*dA0;
      logelong_vec.push_back({logelong, dA, dA0});
      A += dA;
      A0 += dA0;
    }
    logelong_wmean = A > 0. ? logelong_wmean/A : 0.;
    logelong_w0mean = A0 > 0. ? logelong_w0mean/A0 : 0.;
    double logelong_wvar = 0.;
    double logelong_w0var = 0.;
    for (std::vector<std::array<double, 3>>::const_iterator lit = logelong_vec.begin();
         lit != logelong_vec.end(); ++lit){
      logelong_wvar += pow((*lit)[0]-logelong_wmean, 2)*(*lit)[1];
      logelong_w0var += pow((*lit)[0]-logelong_w0mean, 2)*(*lit)[2];
    }
    logelong_wvar = A > 0. ? logelong_wvar/A : 0.;
    logelong_w0var = A0 > 0. ? logelong_w0var/A0 : 0.;
    cols.push_back({"A", A});
    cols.push_back({"A0", A0});
    cols.push_back({"logelong_wmean", logelong_wmean});
    cols.push_back({"logelong_wvar", logelong_wvar});
    cols.push_back({"logelong_w0mean", logelong_w0mean});
    cols.push_back({"logelong_w0var", logelong_w0var});
  }

  return cols;
}

#endif
