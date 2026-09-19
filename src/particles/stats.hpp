#ifndef __STATS_HPP
#define __STATS_HPP

#include <cmath>
#include <vector>
#include "typedefs.hpp"
#include "ParticleSet.hpp"
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
  // Thread order changes the last digits
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
  // Zero rather than NaN for empty or single-point sets
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

  // Fixed by the settled dimension
  const bool do_strip = mesh_dim == 1;
  const bool do_sheet = mesh_dim > 1;

  if (do_sheet){
    // Count too long edges
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

    Uint n_too_long = 0;
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
    // A skipped face keeps its zeros
    std::vector<std::array<double, 3>> logelong_vec(faces.size(), {0., 0., 0.});
    #pragma omp parallel for reduction(+:logelong_wmean,logelong_w0mean,A,A0)
    for (Uint i = 0; i < faces.size(); ++i){
      Uint iedge = faces[i].first[0];
      Uint jedge = faces[i].first[1];
      // Uint kedge = faces[i].first[2];
      double dA0 = faces[i].second;
      if (!(dA0 > 0.))
        continue;                 // degenerate, culled later
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

// Line-element moments, split by the sign of phi
inline std::vector<StatsColumn> vector_stats_columns(
                 const double t,
                 const ParticleSet& ps,
                 const bool has_phi,
                 const unsigned long int n_declined)
{
  const Uint Nrw = ps.N();
  double xm0=0., xm1=0., xm2=0., um0=0., um1=0., um2=0., wm=0., Sm=0., phim=0.;
  double u1m0=0., u1m1=0., u1m2=0., w1m=0., S1m=0.;
  double u2m0=0., u2m1=0., u2m2=0., w2m=0., S2m=0.;
  Uint Nrw1 = 0;
  #pragma omp parallel for reduction(+:xm0,xm1,xm2,um0,um1,um2,wm,Sm,phim,u1m0,u1m1,u1m2,w1m,S1m,u2m0,u2m1,u2m2,w2m,S2m,Nrw1)
  for (Uint i = 0; i < Nrw; ++i){
    const Vector3d xi = ps.x(i), ui = ps.u(i);
    const double wi = ps.w(i), Si = ps.S(i), phii = has_phi ? ps.phi(i) : 1.;
    xm0 += xi[0]; xm1 += xi[1]; xm2 += xi[2];
    um0 += ui[0]; um1 += ui[1]; um2 += ui[2];
    wm += wi; Sm += Si; phim += phii;
    if (phii > 0){ ++Nrw1; u1m0 += ui[0]; u1m1 += ui[1]; u1m2 += ui[2]; w1m += wi; S1m += Si; }
    else         {         u2m0 += ui[0]; u2m1 += ui[1]; u2m2 += ui[2]; w2m += wi; S2m += Si; }
  }
  const Uint Nrw2 = Nrw - Nrw1;
  auto mean = [](double& a, const Uint n){ a = n > 0 ? a/n : 0.; };
  for (double* a : {&xm0, &xm1, &xm2, &um0, &um1, &um2, &wm, &Sm, &phim}) mean(*a, Nrw);
  for (double* a : {&u1m0, &u1m1, &u1m2, &w1m, &S1m}) mean(*a, Nrw1);
  for (double* a : {&u2m0, &u2m1, &u2m2, &w2m, &S2m}) mean(*a, Nrw2);

  double xv0=0., xv1=0., xv2=0., uv0=0., uv1=0., uv2=0., wv=0.;
  double u1v0=0., u1v1=0., u1v2=0., w1v=0., u2v0=0., u2v1=0., u2v2=0., w2v=0.;
  #pragma omp parallel for reduction(+:xv0,xv1,xv2,uv0,uv1,uv2,wv,u1v0,u1v1,u1v2,w1v,u2v0,u2v1,u2v2,w2v)
  for (Uint i = 0; i < Nrw; ++i){
    const Vector3d xi = ps.x(i), ui = ps.u(i);
    const double wi = ps.w(i), phii = has_phi ? ps.phi(i) : 1.;
    xv0 += pow(xi[0]-xm0, 2); xv1 += pow(xi[1]-xm1, 2); xv2 += pow(xi[2]-xm2, 2);
    uv0 += pow(ui[0]-um0, 2); uv1 += pow(ui[1]-um1, 2); uv2 += pow(ui[2]-um2, 2);
    wv += pow(wi-wm, 2);
    if (phii > 0){ u1v0 += pow(ui[0]-u1m0, 2); u1v1 += pow(ui[1]-u1m1, 2); u1v2 += pow(ui[2]-u1m2, 2); w1v += pow(wi-w1m, 2); }
    else         { u2v0 += pow(ui[0]-u2m0, 2); u2v1 += pow(ui[1]-u2m1, 2); u2v2 += pow(ui[2]-u2m2, 2); w2v += pow(wi-w2m, 2); }
  }
  // Unbiased; zero below two particles
  auto var = [](double& a, const Uint n){ a = n > 1 ? a/(n-1) : 0.; };
  for (double* a : {&xv0, &xv1, &xv2, &uv0, &uv1, &uv2, &wv}) var(*a, Nrw);
  for (double* a : {&u1v0, &u1v1, &u1v2, &w1v}) var(*a, Nrw1);
  for (double* a : {&u2v0, &u2v1, &u2v2, &w2v}) var(*a, Nrw2);

  return {{"t", t},
          {"x_mean", xm0}, {"y_mean", xm1}, {"z_mean", xm2},
          {"x_var", xv0}, {"y_var", xv1}, {"z_var", xv2},
          {"ux_mean", um0}, {"uy_mean", um1}, {"uz_mean", um2},
          {"ux_var", uv0}, {"uy_var", uv1}, {"uz_var", uv2},
          {"w_mean", wm}, {"w_var", wv}, {"S_mean", Sm}, {"phi_mean", phim},
          {"Nrw", double(Nrw), true}, {"n_declined", double(n_declined), true},
          {"Nrw1", double(Nrw1), true},
          {"u1x_mean", u1m0}, {"u1y_mean", u1m1}, {"u1z_mean", u1m2},
          {"u1x_var", u1v0}, {"u1y_var", u1v1}, {"u1z_var", u1v2},
          {"w1_mean", w1m}, {"w1_var", w1v}, {"S1_mean", S1m},
          {"Nrw2", double(Nrw2), true},
          {"u2x_mean", u2m0}, {"u2y_mean", u2m1}, {"u2z_mean", u2m2},
          {"u2x_var", u2v0}, {"u2y_var", u2v1}, {"u2z_var", u2v2},
          {"w2_mean", w2m}, {"w2_var", w2v}, {"S2_mean", S2m}};
}

// Cloud position and velocity moments
inline std::vector<StatsColumn> cloud_stats_columns(
                 const double t,
                 const ParticleSet& ps,
                 const unsigned long int n_declined)
{
  const Uint Nrw = ps.N();
  double xm0=0., xm1=0., xm2=0., um0=0., um1=0., um2=0.;
  #pragma omp parallel for reduction(+:xm0,xm1,xm2,um0,um1,um2)
  for (Uint i = 0; i < Nrw; ++i){
    const Vector3d xi = ps.x(i), ui = ps.u(i);
    xm0 += xi[0]; xm1 += xi[1]; xm2 += xi[2];
    um0 += ui[0]; um1 += ui[1]; um2 += ui[2];
  }
  // Zero, not NaN, for N < 2
  if (Nrw > 0){
    for (double* a : {&xm0, &xm1, &xm2, &um0, &um1, &um2}) *a /= Nrw;
  }
  double dx0=0., dx1=0., dx2=0.;
  if (Nrw > 1){
    #pragma omp parallel for reduction(+:dx0,dx1,dx2)
    for (Uint i = 0; i < Nrw; ++i){
      const Vector3d xi = ps.x(i);
      dx0 += pow(xi[0]-xm0, 2); dx1 += pow(xi[1]-xm1, 2); dx2 += pow(xi[2]-xm2, 2);
    }
    for (double* a : {&dx0, &dx1, &dx2}) *a /= (Nrw-1);
  }
  return {{"t", t},
          {"x_mean", xm0}, {"y_mean", xm1}, {"z_mean", xm2},
          {"dx2_mean", dx0}, {"dy2_mean", dx1}, {"dz2_mean", dx2},
          {"ux_mean", um0}, {"uy_mean", um1}, {"uz_mean", um2},
          {"Nrw", double(Nrw), true}, {"n_declined", double(n_declined), true}};
}

// Cloud statistics plus edge elongation with doublings
inline std::vector<StatsColumn> pair_stats_columns(
                 const double t,
                 const ParticleSet& ps,
                 const EdgesType& edges,
                 const std::vector<Uint>& doublings,
                 const unsigned long int n_declined)
{
  std::vector<StatsColumn> cols = cloud_stats_columns(t, ps, n_declined);
  const Uint nsum = edges.size();
  std::vector<double> logelong(nsum);
  double elong_mean = 0., elong2_mean = 0., logelong_mean = 0.;
  #pragma omp parallel for reduction(+:elong_mean,elong2_mean,logelong_mean)
  for (Uint i = 0; i < nsum; ++i){
    const double r = ps.dist(edges[i].first[0], edges[i].first[1])/edges[i].second;
    const double elong = r * exp2(doublings[i]);
    logelong[i] = log(r) + doublings[i] * log(2);
    elong_mean += elong;
    elong2_mean += elong*elong;
    logelong_mean += logelong[i];
  }
  double logelong_var = 0.;
  if (nsum > 0){
    elong_mean /= nsum; elong2_mean /= nsum; logelong_mean /= nsum;
    #pragma omp parallel for reduction(+:logelong_var)
    for (Uint i = 0; i < nsum; ++i){
      const double d = logelong[i] - logelong_mean;
      logelong_var += d*d;
    }
    logelong_var /= nsum;
  }
  cols.push_back({"n_edges", double(nsum), true});
  cols.push_back({"elong_mean", elong_mean});
  cols.push_back({"elong2_mean", elong2_mean});
  cols.push_back({"logelong_mean", logelong_mean});
  cols.push_back({"logelong_var", logelong_var});
  return cols;
}

#endif
