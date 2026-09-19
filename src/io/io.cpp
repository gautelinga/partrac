#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <limits>
#include "H5Cpp.h"
#include "io.hpp"

// using namespace std;

// recently moved here
void load_vector_field(const std::string& input_file,
                       std::vector<Vector3d> &pos_init){
  std::ifstream infile(input_file);
  double x, y, z;
  while (infile >> x >> y >> z){
    pos_init.push_back({x, y, z});
  }
  infile.close();
}

// Full precision for checkpoints
const int checkpoint_precision = std::numeric_limits<double>::max_digits10;

void dump_vector_field(const std::string& output_file,
                       const std::vector<Vector3d>& x_rw,
                       const Uint Nrw){
  std::ofstream outfile(output_file);
  for (Uint irw=0; irw<Nrw; ++irw){
    outfile << std::setprecision(checkpoint_precision)
            << x_rw[irw][0] << " "
            << x_rw[irw][1] << " "
            << x_rw[irw][2] << "\n";
  }
  outfile.close();
}

void dump_vector_field(const std::string& output_file,
                       const std::vector<Vector3d> &pos){
  std::ofstream outfile(output_file);
  for (std::vector<Vector3d>::const_iterator posit=pos.begin();
       posit != pos.end(); ++posit){
    outfile << std::setprecision(checkpoint_precision)
            << (*posit)[0] << " "
            << (*posit)[1] << " "
            << (*posit)[2] << "\n";
  }
  outfile.close();
}

// tau and rho_prev optional (older checkpoints)
void load_faces(const std::string& input_file,
                FacesType& faces){
  std::ifstream infile(input_file);
  std::string line;
  while (std::getline(infile, line)){
    std::istringstream ss(line);
    Uint first, second, third;
    double dA0, tau = 0., rho_prev = 1.;
    if (!(ss >> first >> second >> third >> dA0)) continue;
    ss >> tau >> rho_prev;
    faces.push_back({{first, second, third}, dA0, tau, rho_prev});
  }
  infile.close();
}

void load_edges(const std::string& input_file,
                EdgesType &edges){
  std::ifstream infile(input_file);
  std::string line;
  while (std::getline(infile, line)){
    std::istringstream ss(line);
    Uint first, second;
    double ds0, tau = 0., rho_prev = 1.;
    if (!(ss >> first >> second >> ds0)) continue;
    ss >> tau >> rho_prev;
    edges.push_back({{first, second}, ds0, tau, rho_prev});
  }
  infile.close();
}

void load_list(const std::string& input_file,
               std::vector<Uint> &li){
  std::ifstream infile(input_file);
  Uint a;
  while (infile >> a){
    li.push_back(a);
  }
  infile.close();
}

void dump_list(const std::string& output_file,
               const std::vector<Uint> &li){
  dump_list(output_file, li, li.size());
}

void dump_list(const std::string& output_file,
               const std::vector<Uint> &li, const Uint n){
  std::ofstream outfile(output_file);
  for (Uint k = 0; k < n; ++k){
    outfile << li[k] << "\n";
  }
  outfile.close();
}

void dump_faces(const std::string& output_file,
                const FacesType &faces){
  std::ofstream outfile(output_file);
  for (FacesType::const_iterator faceit = faces.begin();
       faceit != faces.end(); ++faceit){
    outfile << faceit->first[0] << " " << faceit->first[1] << " " << faceit->first[2]
            << " " << std::setprecision(checkpoint_precision) << faceit->second
            << " " << faceit->tau << " " << faceit->rho_prev << "\n";
  }
  outfile.close();
}

void dump_edges(const std::string& output_file,
                const EdgesType &edges){
  std::ofstream outfile(output_file);
  for (auto edgeit = edges.begin();
       edgeit != edges.end(); ++edgeit){
    outfile << edgeit->first[0] << " " << edgeit->first[1] << " "
            << std::setprecision(checkpoint_precision) << edgeit->second
            << " " << edgeit->tau << " " << edgeit->rho_prev << "\n";
  }
  outfile.close();
}

void load_scalar_field(const std::string& input_file,
                       std::vector<double>& c_rw, const Uint Nrw){
  // TODO: to hdf5
  std::ifstream infile(input_file);
  for (Uint irw=0; irw < Nrw; ++irw){
    infile >> c_rw[irw];
  }
  infile.close();
}

void dump_scalar_field(const std::string& output_file,
                       const std::vector<double>& c_rw, const Uint Nrw){
                   // TODO: to hdf5
  std::ofstream outfile(output_file);
  for (Uint irw=0; irw < Nrw; ++irw){
    outfile << std::setprecision(checkpoint_precision) << c_rw[irw] << "\n";
  }
  outfile.close();
}

void tensor2hdf5(H5::H5File& h5f, const std::string& dsetname,
                 const std::vector<double>& axx_rw, const std::vector<double>& axy_rw, const std::vector<double>& axz_rw,
                 const std::vector<double>& ayx_rw, const std::vector<double>& ayy_rw, const std::vector<double>& ayz_rw,
                 const std::vector<double>& azx_rw, const std::vector<double>& azy_rw, const std::vector<double>& azz_rw,
                 const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 3*3;
  H5::DataSpace dspace(2, dims);
  std::vector<double> data(Nrw*3*3);
  for (Uint irw=0; irw < Nrw; ++irw){
    data[irw*3*3+0] = axx_rw[irw];
    data[irw*3*3+1] = axy_rw[irw];
    data[irw*3*3+2] = axz_rw[irw];
    data[irw*3*3+3] = ayx_rw[irw];
    data[irw*3*3+4] = ayy_rw[irw];
    data[irw*3*3+5] = ayz_rw[irw];
    data[irw*3*3+6] = azx_rw[irw];
    data[irw*3*3+7] = azy_rw[irw];
    data[irw*3*3+8] = azz_rw[irw];
  }
  H5::DataSet dset = h5f.createDataSet(dsetname,
                                    H5::PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}

void vector2hdf5(H5::H5File& h5f, const std::string& dsetname,
                 const std::vector<double>& ax_rw, const std::vector<double>& ay_rw, const std::vector<double>& az_rw,
                 const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 3;
  H5::DataSpace dspace(2, dims);
  std::vector<double> data(Nrw*3);
  for (Uint irw=0; irw < Nrw; ++irw){
    data[irw*3+0] = ax_rw[irw];
    data[irw*3+1] = ay_rw[irw];
    data[irw*3+2] = az_rw[irw];
  }
  H5::DataSet dset = h5f.createDataSet(dsetname,
                                    H5::PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}

// Write storage directly
void vector2hdf5(H5::H5File& h5f, const std::string& dsetname,
                 const std::vector<Vector3d>& a_rw, const Uint Nrw){
  static_assert(sizeof(Vector3d) == 3*sizeof(double), "Vector3d is not three bare doubles");
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 3;
  H5::DataSpace dspace(2, dims);
  H5::DataSet dset = h5f.createDataSet(dsetname,
                                    H5::PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(Nrw > 0 ? a_rw[0].data() : nullptr, H5::PredType::NATIVE_DOUBLE);
}


// Row-major (Eigen is column-major)
void tensor2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<Matrix3d>& M_rw,
                 const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 9;
  H5::DataSpace dspace(2, dims);
  std::vector<double> data(Nrw*9);
  #pragma omp parallel for
  for (Uint irw=0; irw < Nrw; ++irw){
    for (Uint i=0; i<3; ++i)
      for (Uint j=0; j<3; ++j)
        data[irw*9 + 3*i + j] = M_rw[irw](i, j);
  }
  H5::DataSet dset = h5f.createDataSet(dsetname, H5::PredType::NATIVE_DOUBLE, dspace);
  dset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}

void ulong2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<Uint>& a, const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 1;
  H5::DataSpace dspace(2, dims);
  H5::DataSet dset = h5f.createDataSet(dsetname, H5::PredType::NATIVE_ULONG, dspace);
  dset.write(a.data(), H5::PredType::NATIVE_ULONG);
}

void int2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<int>& a, const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 1;
  H5::DataSpace dspace(2, dims);
  H5::DataSet dset = h5f.createDataSet(dsetname, H5::PredType::NATIVE_INT, dspace);
  dset.write(a.data(), H5::PredType::NATIVE_INT);
}

void dump_tensor_field(const std::string& output_file, const std::vector<Matrix3d>& M_rw, const Uint Nrw){
  std::ofstream outfile(output_file);
  outfile << std::setprecision(checkpoint_precision);
  for (Uint irw=0; irw < Nrw; ++irw){
    for (Uint i=0; i<3; ++i)
      for (Uint j=0; j<3; ++j)
        outfile << M_rw[irw](i, j) << (i == 2 && j == 2 ? "\n" : " ");
  }
  outfile.close();
}

void load_tensor_field(const std::string& input_file, std::vector<Matrix3d>& M_rw, const Uint Nrw){
  std::ifstream infile(input_file);
  for (Uint irw=0; irw < Nrw; ++irw)
    for (Uint i=0; i<3; ++i)
      for (Uint j=0; j<3; ++j)
        infile >> M_rw[irw](i, j);
  infile.close();
}

// Into a pre-sized array
void load_vector_field(const std::string& input_file, std::vector<Vector3d>& a_rw, const Uint Nrw){
  std::ifstream infile(input_file);
  for (Uint irw=0; irw < Nrw; ++irw)
    infile >> a_rw[irw][0] >> a_rw[irw][1] >> a_rw[irw][2];
  infile.close();
}

void scalar2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<double>& c_rw,
                 const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 1;
  H5::DataSpace dspace(2, dims);
  H5::DataSet dset = h5f.createDataSet(dsetname,
                                    H5::PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(c_rw.data(), H5::PredType::NATIVE_DOUBLE);
}

