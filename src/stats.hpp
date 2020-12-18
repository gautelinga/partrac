#include "utils.hpp"

#ifndef __STATS_HPP
#define __STATS_HPP

//using namespace std;

void write_stats(std::ofstream &statfile,
                 const double t,
                 Vector3d* x_rw,
                 Vector3d* u_rw,
                 const Uint Nrw,
                 FacesType &faces,
                 EdgesType &edges,
                 const double ds_max,
                 const bool do_dump_hist,
                 const std::string histfolder,
                 const unsigned long int n_accepted,
                 const unsigned long int n_declined
                 ){
  double x_mean = 0.;
  double y_mean = 0.;
  double z_mean = 0.;
  double dx2_mean = 0.;
  double dy2_mean = 0.;
  double dz2_mean = 0.;
  double ux_mean = 0.;
  double uy_mean = 0.;
  double uz_mean = 0.;
  for (Uint irw=0; irw < Nrw; ++irw){
    // Sample mean
    x_mean += x_rw[irw][0]/Nrw;
    y_mean += x_rw[irw][1]/Nrw;
    z_mean += x_rw[irw][2]/Nrw;
    ux_mean += u_rw[irw][0]/Nrw;
    uy_mean += u_rw[irw][1]/Nrw;
    uz_mean += u_rw[irw][2]/Nrw;
  }
  for (Uint irw=0; irw < Nrw; ++irw){
    // Sample variance
    dx2_mean += pow(x_rw[irw][0]-x_mean, 2)/(Nrw-1);
    dy2_mean += pow(x_rw[irw][1]-y_mean, 2)/(Nrw-1);
    dz2_mean += pow(x_rw[irw][2]-z_mean, 2)/(Nrw-1);
  }
  statfile << t << "\t"                   //  1
           << x_mean << "\t"              //  2
           << dx2_mean << "\t"            //  3
           << y_mean << "\t"              //  4
           << dy2_mean << "\t"            //  5
           << z_mean << "\t"              //  6
           << dz2_mean << "\t"            //  7
           << ux_mean << "\t"             //  8
           << uy_mean << "\t"             //  9
           << uz_mean << "\t"             // 10
           << Nrw << "\t"                 // 11
           << n_accepted << "\t"          // 12
           << n_declined << "\t";          // 13

  if (edges.size() > 0){
    // Strip method
    double s = 0.;
    Uint n_too_long = 0;
    double logelong_wmean = 0.;
    double logelong_w0mean = 0.;
    double s0 = 0.;
    std::vector<std::array<double, 3>> logelong_vec;
    for (EdgesType::const_iterator edgeit = edges.begin();
         edgeit != edges.end(); ++edgeit){
      int inode = edgeit->first[0];
      int jnode = edgeit->first[1];
      double ds0 = edgeit->second;
      double ds = dist(inode, jnode, x_rw);
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
      double s_ = 0.;
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
      }
      statfile << s << "\t"                   // 15
               << s0 << "\t"                  // 16
               << logelong_wmean << "\t"       // 17
               << logelong_wvar << "\t"        // 18
               << logelong_w0mean << "\t"      // 19
               << logelong_w0var << "\t";      // 20
    }
  }
  if (faces.size() > 0){
    // Sheet method
    double A = 0.;
    double A0 = 0.;
    double logelong_wmean = 0.;
    double logelong_w0mean = 0.;
    std::vector<std::array<double, 3>> logelong_vec;
    for (FacesType::const_iterator faceit = faces.begin();
         faceit != faces.end(); ++faceit){
      Uint iedge = faceit->first[0];
      Uint jedge = faceit->first[1];
      // Uint kedge = faceit->first[2];
      double dA0 = faceit->second;
      double dA = area(iedge, jedge, x_rw, edges);
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
    double A_ = 0.;
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
    }
    statfile << A << "\t"                   // 15
             << A0 << "\t"                  // 16
             << logelong_wmean << "\t"  // 17
             << logelong_wvar << "\t"   // 18
             << logelong_w0mean << "\t" // 19
             << logelong_w0var << "\t"; // 20
  }
  statfile << std::endl;
}

void write_stats_header(std::ofstream &statfile,
                        const FacesType &faces,
                        const EdgesType &edges){
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
  if (edges.size() > 0){
    statfile << "n_too_long" << "\t";         // 14
  }
  if (faces.size() == 0){
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
             << "logelong_wmean" << "\t"  // 17
             << "logelong_wvar" << "\t"   // 18
             << "logelong_w0mean" << "\t" // 19
             << "logelong_w0var" << "\t"; // 20
  }
  statfile << std::endl;
}

#endif
