#ifndef __EXP_STATISTICS_HPP
#define __EXP_STATISTICS_HPP

//#include "utils.hpp"
#include "typedefs.hpp"

template<typename T>
void write_stats( MPIwrap& mpi
                , std::ofstream &statfile
                , const Real t
                , T& ps
                , const unsigned long int n_declined
                )
{
  Vector x_mean = {0., 0., 0.};
  Vector dx2_mean = {0., 0., 0.};
  Vector u_mean = {0., 0., 0.};

  Vector x_mean_glob = {0., 0., 0.};
  Vector dx2_mean_glob = {0., 0., 0.};
  Vector u_mean_glob = {0., 0., 0.};

  Uint Nrw = ps.particles().size();
  Uint Nrw_glob = mpi.allsum<Uint>(Nrw, MPI_UNSIGNED_LONG);

  Real logelong_mean = 0.;
  Real logelong_wmean = 0.;
  Real logelong_w0mean = 0.;
  Real nsum = 0.;
  Real wsum = 0.;
  Real w0sum = 0.;
  Real logelong_mean_glob = 0.;
  Real logelong_wmean_glob = 0.;
  Real logelong_w0mean_glob = 0.;
  Real nsum_glob = 0.;
  Real wsum_glob = 0.;
  Real w0sum_glob = 0.;

  for ( auto & particle : ps.particles() )
  {
    // Sample mean
    x_mean += particle.x(); // /Nrw;
    u_mean += particle.u(); // /Nrw;
  }
  for (int d=0; d<3; ++d)
  {
    auto xtmp = mpi.allsum<Real>(x_mean[d], MPI_Real);
    auto utmp = mpi.allsum<Real>(u_mean[d], MPI_Real);
    x_mean_glob[d] = xtmp / Nrw_glob;
    u_mean_glob[d] = utmp / Nrw_glob;
  }
  x_mean /= Nrw;
  u_mean /= Nrw;

  for ( auto & particle : ps.particles() )
  {
    // Sample variance
    Vector dx = particle.x()-x_mean;
    dx2_mean += dx.cwiseProduct(dx)/(Nrw-1);

    Vector dx_glob = particle.x() - x_mean_glob;
    dx2_mean_glob += dx_glob.cwiseProduct(dx_glob);
  }
  for (int d=0; d<3; ++d)
  {
    auto dx2tmp = mpi.allsum<Real>(dx2_mean_glob[d], MPI_Real);
    dx2_mean_glob[d] = dx2tmp / (Nrw_glob - 1);
  }
  
  if (ps.dim() > 0){
    switch (ps.dim())
    {
      case 1:
      {
        // Strip method
        for ( auto & edge : ps.edges() ){
          Real ds0 = edge.l0();
          Real ds = edge.length(ps);
          Real logelong = log(ds/ds0);
          logelong_mean += logelong;
          logelong_wmean += logelong * ds;
          logelong_w0mean += logelong * ds0;
          nsum += 1;
          wsum += ds;
          w0sum += ds0;
        }
        break;
      }
      case 2:
      {
        // Sheet method
        for ( auto & face : ps.faces() )
        {
          Real dA0 = face.A0();
          Real dA = face.area(ps);
          Real logelong = log(dA/dA0);
          logelong_mean += logelong;
          logelong_wmean += logelong * dA;
          logelong_w0mean += logelong * dA0;
          nsum += 1;
          wsum += dA;
          w0sum += dA0;
        }
        break;
      }
      default:
      {
        std::cout << "ERROR: Unknown topology!" << std::endl;
        exit(0);
      }
    }
    nsum_glob = mpi.sum<Real>(nsum, MPI_Real);
    wsum_glob = mpi.sum<Real>(wsum, MPI_Real);
    w0sum_glob = mpi.sum<Real>(w0sum, MPI_Real);
    logelong_mean_glob = mpi.sum<Real>(logelong_mean, MPI_Real);
    logelong_wmean_glob = mpi.sum<Real>(logelong_wmean, MPI_Real);
    logelong_w0mean_glob = mpi.sum<Real>(logelong_w0mean, MPI_Real);

    logelong_mean /= nsum;
    logelong_wmean /= wsum;
    logelong_w0mean /= w0sum;
    if (mpi.rank() == 0){
      logelong_mean_glob /= nsum_glob;
      logelong_wmean_glob /= wsum_glob;
      logelong_w0mean_glob /= w0sum_glob;
    }
  }
  statfile << t                       << "\t"           //  1
           << x_mean[0]               << "\t"           //  2
           << dx2_mean[0]             << "\t"           //  3
           << x_mean[1]               << "\t"           //  4
           << dx2_mean[1]             << "\t"           //  5
           << x_mean[2]               << "\t"           //  6
           << dx2_mean[2]             << "\t"           //  7
           << u_mean[0]               << "\t"           //  8
           << u_mean[1]               << "\t"           //  9
           << u_mean[2]               << "\t"           // 10
           << Nrw                     << "\t"           // 11
           << n_declined              << "\t";          // 12

  if (ps.dim() > 0){ 
    statfile << nsum                    << "\t"           // 13
             << wsum                    << "\t"           // 14
             << w0sum                   << "\t"           // 15
             << logelong_mean           << "\t"           // 16
             << logelong_wmean          << "\t"           // 17
             << logelong_w0mean         << "\t";          // 18
  }
  if (mpi.rank() == 0){
    statfile << x_mean_glob[0]        << "\t"           // 19
             << dx2_mean_glob[0]      << "\t"           // 20
             << x_mean_glob[1]        << "\t"           // 21
             << dx2_mean_glob[1]      << "\t"           // 22
             << x_mean_glob[2]        << "\t"           // 23
             << dx2_mean_glob[2]      << "\t";          // 24
    if (ps.dim() > 0){
      statfile << nsum_glob             << "\t"           // 25
               << wsum_glob             << "\t"           // 26
               << w0sum_glob            << "\t"           // 27
               << logelong_mean_glob    << "\t"           // 28
               << logelong_wmean_glob   << "\t"           // 29
               << logelong_w0mean_glob  << "\t";          // 30
    }
  }
  statfile << std::endl;
}

void write_stats_header(MPIwrap& mpi, std::ofstream &statfile, Uint mesh_dim){
  std::string wsumstr = "";
  if (mesh_dim == 1){
    wsumstr = "s";
  }
  else if (mesh_dim == 2){
    wsumstr = "A";
  }

  statfile << "# t" << "\t"                   //  1
           << "x_mean" << "\t"                //  2
           << "dx2_mean" << "\t"              //  3
           << "y_mean" << "\t"                //  4
           << "dy2_mean" << "\t"              //  5
           << "z_mean" << "\t"                //  6
           << "dz2_mean" << "\t"              //  7
           << "ux_mean" << "\t"               //  8
           << "uy_mean" << "\t"               //  9
           << "uz_mean" << "\t"               // 10
           << "Nrw" << "\t"                   // 11
           << "n_declined" << "\t";           // 12
  if (mesh_dim > 0){
    statfile << "n_edges \t"                                          // 13
             << wsumstr << "\t"                                    // 14
             << wsumstr << "0" << "\t"                             // 15
             << "logelong_mean" << "\t"                            // 16
             << "logelong_" << wsumstr << "mean" << "\t"           // 17
             << "logelong_" << wsumstr << "0mean" << "\t";         // 18
  }
  if (mpi.rank() == 0){
    statfile << "x_mean" << "\t"                                   //  19
             << "dx2_mean" << "\t"                                 //  20
             << "y_mean" << "\t"                                   //  21
             << "dy2_mean" << "\t"                                 //  22
             << "z_mean" << "\t"                                   //  23
             << "dz2_mean" << "\t";                                //  24
    if (mesh_dim > 0){
      statfile << "n_edges_glob \t"                                        // 25
               << wsumstr << "_glob" << "\t"                       // 26
               << wsumstr << "0_glob" << "\t"                      // 27
               << "logelong_mean" << "\t"                          // 28
               << "logelong_" << wsumstr << "mean_glob" << "\t"    // 29
               << "logelong_" << wsumstr << "0mean_glob" << "\t";  // 30
    }
  }
  statfile << std::endl;
}

#endif
