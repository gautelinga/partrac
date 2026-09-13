#ifndef __EXP_STATISTICS_HPP
#define __EXP_STATISTICS_HPP

//#include "utils.hpp"
#include "typedefs.hpp"
#include "stats_columns.hpp"

// One list, read once for the header and once per row. ps.dim() decides the
// columns in both, so the header cannot promise what the rows do not fill.
template<typename T>
std::vector<StatsColumn> particle_stats_columns( const Real t
                                                , T& ps
                                                , const unsigned long int n_declined
                                                )
{
  std::vector<StatsColumn> cols;
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
  // the mean of nothing, and the spread of a single point, are both
  // undefined; report zero rather than a NaN that spreads downstream
  if (Nrw > 0){
    x_mean /= Nrw;
    u_mean /= Nrw;
  }

  if (Nrw > 1){
    for ( auto & particle : ps.particles() )
    {
      // Sample variance
      Vector dx = particle.x()-x_mean;
      dx2_mean += dx.cwiseProduct(dx)/(Nrw-1);
    }
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
  cols.push_back({"t", t});
  cols.push_back({"x_mean", x_mean[0]});
  cols.push_back({"y_mean", x_mean[1]});
  cols.push_back({"z_mean", x_mean[2]});
  cols.push_back({"dx2_mean", dx2_mean[0]});
  cols.push_back({"dy2_mean", dx2_mean[1]});
  cols.push_back({"dz2_mean", dx2_mean[2]});
  cols.push_back({"ux_mean", u_mean[0]});
  cols.push_back({"uy_mean", u_mean[1]});
  cols.push_back({"uz_mean", u_mean[2]});
  cols.push_back({"Nrw", double(Nrw), true});
  cols.push_back({"n_declined", double(n_declined), true});

  if (ps.dim() > 0){
    cols.push_back({"n_edges", double(nsum), true});
    cols.push_back({"elong_mean", elong_mean});
    cols.push_back({"elong2_mean", elong2_mean});
    cols.push_back({"logelong_mean", logelong_mean});
    cols.push_back({"logelong_var", logelong_var});
  }
  return cols;
}

#endif
