#ifndef __STATS_HPP
#define __STATS_HPP

#include "utils.hpp"
#include "stats_columns.hpp"

inline std::vector<StatsColumn> mesh_stats_columns(
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
  // Summed in whatever order the threads finish, so the last digits move with it
  {
    double xm0=0., xm1=0., xm2=0., um0=0., um1=0., um2=0.;
    #pragma omp parallel for reduction(+:xm0,xm1,xm2,um0,um1,um2)
    for (Uint irw=0; irw < Nrw; ++irw){
      // Sample mean
      const Vector3d xi = ps.x(irw), ui = ps.u(irw);
      xm0 += xi[0]; xm1 += xi[1]; xm2 += xi[2];
      um0 += ui[0]; um1 += ui[1]; um2 += ui[2];
    }
    x_mean = {xm0, xm1, xm2};
    u_mean = {um0, um1, um2};
  }
  // the mean of nothing, and the spread of a single point, are both
  // undefined; report zero rather than a NaN that spreads downstream
  if (Nrw > 0){
    x_mean /= Nrw;
    u_mean /= Nrw;
  }

  if (Nrw > 1){
    double v0=0., v1=0., v2=0.;
    #pragma omp parallel for reduction(+:v0,v1,v2)
    for (Uint irw=0; irw < Nrw; ++irw){
      // Sample variance
      Vector3d dx = ps.x(irw)-x_mean;
      dx = dx.cwiseProduct(dx)/(Nrw-1);
      v0 += dx[0]; v1 += dx[1]; v2 += dx[2];
    }
    dx2_mean = {v0, v1, v2};
  }

  cols = {{"t", t},
          {"x_mean", x_mean[0]}, {"dx2_mean", dx2_mean[0]},
          {"y_mean", x_mean[1]}, {"dy2_mean", dx2_mean[1]},
          {"z_mean", x_mean[2]}, {"dz2_mean", dx2_mean[2]},
          {"ux_mean", u_mean[0]},
          {"uy_mean", u_mean[1]},
          {"uz_mean", u_mean[2]},
          {"Nrw", double(Nrw), true},
          {"n_accepted", double(n_accepted), true},
          {"n_declined", double(n_declined), true}};

  double s = 0.;
  double s0 = 0.;
  double A = 0.;
  double A0 = 0.;

  // Fixed by the dimension the run settles into, before the first injection
  const bool do_strip = mesh_dim == 1;
  const bool do_sheet = mesh_dim > 1;

  if (do_sheet){
    // the sheet's own loop is over faces, so its edges are measured here
    Uint n_too_long = 0;
    #pragma omp parallel for reduction(+:n_too_long)
    for (Uint i = 0; i < edges.size(); ++i)
      if (ps.dist(edges[i].first[0], edges[i].first[1]) > ds_max)
        ++n_too_long;
    cols.push_back({"n_too_long", double(n_too_long), true});
  }

  if (do_strip){
    // Strip method
    double logelong_wmean = 0.;
    double logelong_w0mean = 0.;

    // the strip measures every edge anyway, so the count rides along
    Uint n_too_long = 0;
    // Written by index, so the pass can be split
    std::vector<std::array<double, 3>> logelong_vec(edges.size());
    #pragma omp parallel for reduction(+:logelong_wmean,logelong_w0mean,s,s0,n_too_long)
    for (Uint i = 0; i < edges.size(); ++i){
      Uint inode = edges[i].first[0];
      Uint jnode = edges[i].first[1];
      double ds0 = edges[i].second;
      double ds = ps.dist(inode, jnode);
      if (ds > ds_max)
        ++n_too_long;
      double logelong = log(ds/ds0);
      logelong_wmean += logelong*ds;
      logelong_w0mean += logelong*ds0;
      logelong_vec[i] = {logelong, ds, ds0};
      s += ds;
      s0 += ds0;
    }
    cols.push_back({"n_too_long", double(n_too_long), true});
    logelong_wmean = s > 0. ? logelong_wmean/s : 0.;
    logelong_w0mean = s0 > 0. ? logelong_w0mean/s0 : 0.;

    double logelong_wvar = 0.;
    double logelong_w0var = 0.;
    #pragma omp parallel for reduction(+:logelong_wvar,logelong_w0var)
    for (Uint i = 0; i < logelong_vec.size(); ++i){
      logelong_wvar += pow(logelong_vec[i][0]-logelong_wmean, 2)*logelong_vec[i][1];
      logelong_w0var += pow(logelong_vec[i][0]-logelong_w0mean, 2)*logelong_vec[i][2];
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
    // Written by index, so the pass can be split; a skipped face keeps its zeros
    std::vector<std::array<double, 3>> logelong_vec(faces.size(), {0., 0., 0.});
    #pragma omp parallel for reduction(+:logelong_wmean,logelong_w0mean,A,A0)
    for (Uint i = 0; i < faces.size(); ++i){
      Uint iedge = faces[i].first[0];
      Uint jedge = faces[i].first[1];
      // Uint kedge = faces[i].first[2];
      double dA0 = faces[i].second;
      if (!(dA0 > 0.))
        continue;                 // a flat sweep, waiting to be culled
      double dA = ps.triangle_area(iedge, jedge, edges);
      double logelong = log(dA/dA0);
      logelong_wmean += logelong*dA;
      logelong_w0mean += logelong*dA0;
      logelong_vec[i] = {logelong, dA, dA0};
      A += dA;
      A0 += dA0;
    }
    logelong_wmean = A > 0. ? logelong_wmean/A : 0.;
    logelong_w0mean = A0 > 0. ? logelong_w0mean/A0 : 0.;
    double logelong_wvar = 0.;
    double logelong_w0var = 0.;
    #pragma omp parallel for reduction(+:logelong_wvar,logelong_w0var)
    for (Uint i = 0; i < logelong_vec.size(); ++i){
      logelong_wvar += pow(logelong_vec[i][0]-logelong_wmean, 2)*logelong_vec[i][1];
      logelong_w0var += pow(logelong_vec[i][0]-logelong_w0mean, 2)*logelong_vec[i][2];
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
