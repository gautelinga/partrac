#ifndef __EXP_STATISTICS_HPP
#define __EXP_STATISTICS_HPP

//#include "utils.hpp"
#include "typedefs.hpp"

template<typename T>
void write_stats( std::ofstream &statfile
                , const Real t
                , T& ps
                , const unsigned long int n_declined
                )
{
  Vector x_mean = {0., 0., 0.};
  Vector dx2_mean = {0., 0., 0.};
  Vector u_mean = {0., 0., 0.};


  Uint Nrw = ps.particles().size();

  Real logelong_mean = 0.;
  Real logelong_var = 0.;
  Real elong_mean = 0.;
  Real elong2_mean = 0.;
  //Real logelong_wmean = 0.;
  //Real logelong_w0mean = 0.;
  //Real nsum = 0.;
  Uint nsum = 0;

  for ( auto & particle : ps.particles() )
  {
    // Sample mean
    x_mean += particle.x(); // /Nrw;
    u_mean += particle.u(); // /Nrw;
  }
  x_mean /= Nrw;
  u_mean /= Nrw;

  for ( auto & particle : ps.particles() )
  {
    // Sample variance
    Vector dx = particle.x()-x_mean;
    dx2_mean += dx.cwiseProduct(dx)/(Nrw-1);

  }
  
  if (ps.dim() > 0){
    switch (ps.dim())
    {
      case 1:
      {
        // Strip method
        for ( auto & edge : ps.edges() ){
          //Real ds0 = edge.l0();
          //Real ds = edge.length(ps);
          //Real logelong = log(ds/ds0);
          Real elong = edge.elong(ps);
          Real logelong = edge.logelong(ps);
          elong_mean += elong;
          elong2_mean += elong*elong;
          logelong_mean += logelong;
          //logelong_wmean += logelong * ds;
          //logelong_w0mean += logelong * ds0;
          //nsum += 1;
          //wsum += ds;
          //w0sum += ds0;
        }
        nsum = ps.edges().size();
        logelong_mean /= nsum;
        for ( auto & edge : ps.edges() )
        {
          Real Dlogelong = edge.logelong(ps) - logelong_mean;
          logelong_var += Dlogelong*Dlogelong;
        }
        logelong_var /= nsum;
        break;
      }
      case 2:
      {
        // Sheet method
        for ( auto & face : ps.faces() )
        {
          //Real dA0 = face.A0();
          //Real dA = face.area(ps);
          //Real logelong = log(dA/dA0);
          Real elong = face.elong(ps);
          Real logelong = face.elong(ps);
          elong_mean += elong;
          elong2_mean += elong*elong;
          logelong_mean += logelong;
          //logelong_wmean += logelong * dA;
          //logelong_w0mean += logelong * dA0;
          //nsum += 1;
          //wsum += dA;
          //w0sum += dA0;
        }
        nsum = ps.faces().size();
        logelong_mean /= nsum;
        for ( auto & face : ps.faces() )
        {
          Real Dlogelong = face.elong(ps)-logelong_mean;
          logelong_var += Dlogelong*Dlogelong;
        }
        logelong_var /= nsum;
        break;
      }
      default:
      {
        std::cout << "ERROR: Unknown topology!" << std::endl;
        exit(1);
      }
    }

    elong_mean /= nsum;
    elong2_mean /= nsum;

  }
  statfile << t                       << "\t"           //  1
           << x_mean[0]               << "\t"           //  2
           << x_mean[1]               << "\t"           //  3
           << x_mean[2]               << "\t"           //  4
           << dx2_mean[0]             << "\t"           //  5
           << dx2_mean[1]             << "\t"           //  6
           << dx2_mean[2]             << "\t"           //  7
           << u_mean[0]               << "\t"           //  8
           << u_mean[1]               << "\t"           //  9
           << u_mean[2]               << "\t"           // 10
           << Nrw                     << "\t"           // 11
           << n_declined              << "\t";          // 12

  if (ps.dim() > 0){ 
    statfile << nsum                    << "\t"           // 13
             //<< wsum                    << "\t"           // 14
             //<< w0sum                   << "\t"           // 15
             << elong_mean              << "\t"            // 14
             << elong2_mean             << "\t"            // 15
             << logelong_mean           << "\t"            // 16
             << logelong_var            << "\t";           // 17
             //<< logelong_wmean          << "\t"           // 17
             //<< logelong_w0mean         << "\t";          // 18
  }
  statfile << std::endl;
}

void write_stats_header(std::ofstream &statfile, Uint mesh_dim){
  std::string wsumstr = "";
  if (mesh_dim == 1){
    wsumstr = "s";
  }
  else if (mesh_dim == 2){
    wsumstr = "A";
  }

  statfile << "# t" << "\t"                   //  1
           << "x_mean" << "\t"                //  2
           << "y_mean" << "\t"                //  3
           << "z_mean" << "\t"                //  4
           << "dx2_mean" << "\t"              //  5
           << "dy2_mean" << "\t"              //  6
           << "dz2_mean" << "\t"              //  7
           << "ux_mean" << "\t"               //  8
           << "uy_mean" << "\t"               //  9
           << "uz_mean" << "\t"               // 10
           << "Nrw" << "\t"                   // 11
           << "n_declined" << "\t";           // 12
  if (mesh_dim > 0){
    statfile << "n_edges \t"                                          // 13
             << "elong_mean" << "\t"                                  // 14
             << "elong2_mean" << "\t"                                 // 15
             //<< wsumstr << "\t"                                    // 14
             //<< wsumstr << "0" << "\t"                             // 15
             << "logelong_mean" << "\t";                            // 16
            // << "logelong_" << wsumstr << "mean" << "\t"           // 17
            // << "logelong_" << wsumstr << "0mean" << "\t";         // 18
  }
  statfile << std::endl;
}

#endif
