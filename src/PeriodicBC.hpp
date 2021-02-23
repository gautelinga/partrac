#ifndef PeriodicBC_HPP
#define PeriodicBC_HPP

#include <dolfin.h>


class PeriodicBC : public dolfin::SubDomain {
public:
  PeriodicBC(const std::vector<bool> &periodic,
      const Vector3d &x_min, const Vector3d &x_max, const Uint dim) {
    periodic_x = periodic[0] && dim > 0;
    periodic_y = periodic[1] && dim > 1;
    periodic_z = periodic[2] && dim > 2;
    this->x_min = x_min;
    this->x_max = x_max;
  };
  bool inside(const dolfin::Array<double>& x, bool on_boundary) const {
    return (on_boundary &&
            ((periodic_x && x[0] < x_min[0] + DOLFIN_EPS_LARGE) ||
             (periodic_y && x[1] < x_min[1] + DOLFIN_EPS_LARGE) ||
             (periodic_z && x[2] < x_min[2] + DOLFIN_EPS_LARGE)) &&
            !(
              (periodic_x && periodic_y &&
               x[0] < x_min[0] + DOLFIN_EPS_LARGE &&
               x[1] > x_max[1] - DOLFIN_EPS_LARGE) ||
              (periodic_x && periodic_y &&
               x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
               x[1] < x_min[1] + DOLFIN_EPS_LARGE) ||
              (periodic_x && periodic_z &&
               x[0] < x_min[0] + DOLFIN_EPS_LARGE &&
               x[2] > x_max[2] - DOLFIN_EPS_LARGE) ||
              (periodic_x && periodic_z &&
               x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
               x[2] < x_min[2] + DOLFIN_EPS_LARGE) ||
              (periodic_y && periodic_z &&
               x[1] < x_min[1] + DOLFIN_EPS_LARGE &&
               x[2] > x_max[2] - DOLFIN_EPS_LARGE) ||
              (periodic_y && periodic_z &&
               x[1] > x_max[1] - DOLFIN_EPS_LARGE &&
               x[2] < x_min[2] + DOLFIN_EPS_LARGE)
              )
            );
  }
  void map(const dolfin::Array<double>& x, dolfin::Array<double>& y) const {
    if (periodic_x && periodic_y && periodic_z &&
        x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
        x[1] > x_max[1] - DOLFIN_EPS_LARGE &&
        x[2] > x_max[2] - DOLFIN_EPS_LARGE){
      y[0] = x[0] - (x_max[0]-x_min[0]);
      y[1] = x[1] - (x_max[1]-x_min[1]);
      y[2] = x[2] - (x_max[2]-x_min[2]);
    }
    else if (periodic_x && periodic_y &&
             x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
             x[1] > x_max[1] - DOLFIN_EPS_LARGE){
      y[0] = x[0] - (x_max[0]-x_min[0]);
      y[1] = x[1] - (x_max[1]-x_min[1]);
      y[2] = x[2];
    }
    else if (periodic_x && periodic_z &&
             x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
             x[2] > x_max[2] - DOLFIN_EPS_LARGE){
      y[0] = x[0] - (x_max[0]-x_min[0]);
      y[1] = x[1];
      y[2] = x[2] - (x_max[2]-x_min[2]);
    }
    else if (periodic_y && periodic_z &&
             x[1] > x_max[1] - DOLFIN_EPS_LARGE &&
             x[2] > x_max[2] - DOLFIN_EPS_LARGE){
      y[0] = x[0];
      y[1] = x[1] - (x_max[1]-x_min[1]);
      y[2] = x[2] - (x_max[2]-x_min[2]);
    }
    else if (periodic_x && x[0] > x_max[0] - DOLFIN_EPS_LARGE){
      y[0] = x[0] - (x_max[0]-x_min[0]);
      y[1] = x[1];
      y[2] = x[2];
    }
    else if (periodic_y && x[1] > x_max[1] - DOLFIN_EPS_LARGE){
      y[0] = x[0];
      y[1] = x[1] - (x_max[1]-x_min[1]);
      y[2] = x[2];
    }
    else if (periodic_z && x[2] > x_max[2] - DOLFIN_EPS_LARGE){
      y[0] = x[0];
      y[1] = x[1];
      y[2] = x[2] - (x_max[2]-x_min[2]);
    }
    else {
      y[0] = x[0]-1000;
      y[1] = x[1]-1000;
      y[2] = x[2]-1000;
    }
  }
private:
  double periodic_x;
  double periodic_y;
  double periodic_z;
  Vector3d x_min;
  Vector3d x_max;
};

#endif
