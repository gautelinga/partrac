#ifndef __STATS_HPP
#define __STATS_HPP

#include "utils.hpp"

void write_stats(MPIwrap& mpi,
                 std::ofstream &statfile,
                 const double t,
                 const ParticleSet& ps,
                 const FacesType &faces,
                 const EdgesType &edges,
                 const double ds_max,
                 //const bool do_dump_hist,
                 //const std::string histfolder,
                 const unsigned long int n_accepted,
                 const unsigned long int n_declined)
{
  Vector3d x_mean = {0., 0., 0.};
  Vector3d dx2_mean = {0., 0., 0.};
  Vector3d u_mean = {0., 0., 0.};
  Vector3d x_mean_glob = {0., 0., 0.};
  Vector3d dx2_mean_glob = {0., 0., 0.};
  Vector3d u_mean_glob = {0., 0., 0.};
  Uint Nrw = ps.N();
  Uint Nrw_glob = mpi.allsum<Uint>(Nrw, MPI_UNSIGNED_LONG);
  for (Uint irw=0; irw < Nrw; ++irw){
    // Sample mean
    x_mean += ps.x(irw); // /Nrw;
    u_mean += ps.u(irw); // /Nrw;
  }
  for (int d=0; d<3; ++d){
    auto xtmp = mpi.allsum<double>(x_mean[d], MPI_DOUBLE);
    auto utmp = mpi.allsum<double>(u_mean[d], MPI_DOUBLE);
    x_mean_glob[d] = xtmp / Nrw_glob;
    u_mean_glob[d] = utmp / Nrw_glob;
  }
  x_mean /= Nrw;
  u_mean /= Nrw;

  for (Uint irw=0; irw < Nrw; ++irw){
    // Sample variance
    Vector3d dx = ps.x(irw)-x_mean;
    dx2_mean += dx.cwiseProduct(dx)/(Nrw-1);

    Vector3d dx_glob = ps.x(irw) - x_mean_glob;
    dx2_mean_glob += dx_glob.cwiseProduct(dx_glob);
  }
  for (int d=0; d<3; ++d){
    auto dx2tmp = mpi.allsum<double>(dx2_mean_glob[d], MPI_DOUBLE);
    dx2_mean_glob[d] = dx2tmp / (Nrw_glob - 1);
  }

  statfile << t << "\t"                   //  1
           << x_mean[0] << "\t"           //  2
           << dx2_mean[0] << "\t"         //  3
           << x_mean[1] << "\t"           //  4
           << dx2_mean[1] << "\t"         //  5
           << x_mean[2] << "\t"           //  6
           << dx2_mean[2] << "\t"         //  7
           << u_mean[0] << "\t"           //  8
           << u_mean[1] << "\t"           //  9
           << u_mean[2] << "\t"           // 10
           << Nrw << "\t"                 // 11
           << n_accepted << "\t"          // 12
           << n_declined << "\t";         // 13

  double s = 0.;
  double s0 = 0.;
  double A = 0.;
  double A0 = 0.;

  bool do_strip = edges.size() > 0 && faces.size() == 0;
  bool do_sheet = faces.size() > 0;

  if (do_strip){
    // Strip method
    Uint n_too_long = 0;
    double logelong_wmean = 0.;
    double logelong_w0mean = 0.;

    std::vector<std::array<double, 3>> logelong_vec;
    for (EdgesType::const_iterator edgeit = edges.begin();
         edgeit != edges.end(); ++edgeit){
      int inode = edgeit->first[0];
      int jnode = edgeit->first[1];
      double ds0 = edgeit->second;
      double ds = ps.dist(inode, jnode);
      if (ds > ds_max){
        ++n_too_long;
      }
      double logelong = log(ds/ds0);
      logelong_wmean += logelong*ds;
      logelong_w0mean += logelong*ds0;
      logelong_vec.push_back({logelong, ds, ds0});
      s += ds;
      s0 += ds0;
    }
    logelong_wmean /= s;
    logelong_w0mean /= s0;
    statfile << n_too_long << "\t";          // 14

    if (faces.size() == 0){
      double logelong_wvar = 0.;
      double logelong_w0var = 0.;
      for (std::vector<std::array<double, 3>>::const_iterator lit = logelong_vec.begin();
           lit != logelong_vec.end(); ++lit){
        logelong_wvar += pow((*lit)[0]-logelong_wmean, 2)*(*lit)[1];
        logelong_w0var += pow((*lit)[0]-logelong_w0mean, 2)*(*lit)[2];
      }
      logelong_wvar /= s;
      logelong_w0var /= s0;
      /*double s_ = 0.;
      double s0_ = 0.;
      if (do_dump_hist){
        std::sort(logelong_vec.begin(), logelong_vec.end());
        std::ofstream logelongfile(histfolder + "/logelong_t" + std::to_string(t) + ".hist");
        for (std::vector<std::array<double, 3>>::const_iterator lit = logelong_vec.begin();
             lit != logelong_vec.end(); ++lit){
          s_ += (*lit)[1];
          s0_ += (*lit)[2];
          logelongfile << (*lit)[0] << " "
                       << (*lit)[1] << " " << s_ << " " << s_/s << " "
                       << (*lit)[2] << " " << s0_ << " " << s0_/s0 << " "
                       << std::endl;
        }
        logelongfile.close();
      }*/
      statfile << s << "\t"                   // 15
               << s0 << "\t"                  // 16
               << logelong_wmean << "\t"       // 17
               << logelong_wvar << "\t"        // 18
               << logelong_w0mean << "\t"      // 19
               << logelong_w0var << "\t";      // 20
    }
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
      double dA = ps.triangle_area(iedge, jedge, edges);
      double logelong = log(dA/dA0);
      logelong_wmean += logelong*dA;
      logelong_w0mean += logelong*dA0;
      logelong_vec.push_back({logelong, dA, dA0});
      A += dA;
      A0 += dA0;
    }
    logelong_wmean /= A;
    logelong_w0mean /= A0;
    double logelong_wvar = 0.;
    double logelong_w0var = 0.;
    for (std::vector<std::array<double, 3>>::const_iterator lit = logelong_vec.begin();
         lit != logelong_vec.end(); ++lit){
      logelong_wvar += pow((*lit)[0]-logelong_wmean, 2)*(*lit)[1];
      logelong_w0var += pow((*lit)[0]-logelong_w0mean, 2)*(*lit)[1];
    }
    logelong_wvar /= A;
    logelong_w0var /= A0;
    /*double A_ = 0.;
    double A0_ = 0.;
    if (do_dump_hist){
      std::sort(logelong_vec.begin(), logelong_vec.end());
      std::ofstream logelongfile(histfolder + "/logelong_t" + std::to_string(t) + ".hist");
      for (std::vector<std::array<double, 3>>::const_iterator lit=logelong_vec.begin();
           lit != logelong_vec.end(); ++lit){
        A_ += (*lit)[1];
        A0_ += (*lit)[2];
        logelongfile << (*lit)[0] << " "
                         << (*lit)[1] << " " << A_ << " " << A_/A << " "
                         << (*lit)[2] << " " << A0_ << " " << A0_/A0 << " "
                     << std::endl;
      }
      logelongfile.close();
    }*/
    statfile << A << "\t"                   // 15
             << A0 << "\t"                  // 16
             << logelong_wmean << "\t"  // 17
             << logelong_wvar << "\t"   // 18
             << logelong_w0mean << "\t" // 19
             << logelong_w0var << "\t"; // 20
  }
  double s_glob = mpi.sum<double>(s, MPI_DOUBLE);
  double s0_glob = mpi.sum<double>(s0, MPI_DOUBLE);
  double A_glob = mpi.sum<double>(A, MPI_DOUBLE);
  double A0_glob = mpi.sum<double>(A0, MPI_DOUBLE);
  if (mpi.rank() == 0){
    statfile << x_mean_glob[0] << "\t"           //  21
             << dx2_mean_glob[0] << "\t"         //  22
             << x_mean_glob[1] << "\t"           //  23
             << dx2_mean_glob[1] << "\t"         //  24
             << x_mean_glob[2] << "\t"           //  25
             << dx2_mean_glob[2] << "\t";        //  26
    if (do_strip){
      statfile << s_glob << "\t"           // 27
               << s0_glob << "\t";         // 28
    }
    else if (do_sheet){
      statfile << A_glob << "\t"           // 27
               << A0_glob << "\t";         // 28
    }
  }
 
  statfile << std::endl;
}

void write_stats_header(MPIwrap& mpi, std::ofstream &statfile, Uint mesh_dim){
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
           << "n_accepted" << "\t"            // 12
           << "n_declined" << "\t";           // 13
  if (mesh_dim > 0){
    statfile << "n_too_long" << "\t";         // 14
  }
  if (mesh_dim > 1){
    statfile << "s" << "\t"                   // 15
             << "s0" << "\t"                  // 16
             << "logelong_wmean" << "\t"      // 17
             << "logelong_wvar" << "\t"       // 18
             << "logelong_w0mean" << "\t"     // 19
             << "logelong_w0var" << "\t";     // 20
  }
  else {
    statfile << "A" << "\t"                   // 15
             << "A0" << "\t"                  // 16
             << "logelong_wmean" << "\t"      // 17
             << "logelong_wvar" << "\t"       // 18
             << "logelong_w0mean" << "\t"     // 19
             << "logelong_w0var" << "\t";     // 20
  }
  if (mpi.rank() == 0){
    statfile << "x_mean" << "\t"                //  21
             << "dx2_mean" << "\t"              //  22
             << "y_mean" << "\t"                //  23
             << "dy2_mean" << "\t"              //  24
             << "z_mean" << "\t"                //  25
             << "dz2_mean" << "\t";             //  26
    if (mesh_dim == 1){
      statfile << "s_glob" << "\t"           // 27
               << "s0_glob" << "\t";         // 28
    }
    else if (mesh_dim == 2){
      statfile << "A_glob" << "\t"           // 27
               << "A0_glob" << "\t";         // 28
    }
  }
  statfile << std::endl;
}

#endif
