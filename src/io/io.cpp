#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <limits>
#include "H5Cpp.h"
#include "Error.hpp"
#include "h5direct.hpp"
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

void load_scalar_field(const std::string& input_file,
                       std::vector<double>& c_rw, const Uint Nrw){
  // TODO: to hdf5
  std::ifstream infile(input_file);
  for (Uint irw=0; irw < Nrw; ++irw){
    infile >> c_rw[irw];
  }
  infile.close();
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

void hdf52vectors(const H5::H5File& h5f, const std::string& dsetname, const Uint rows,
                  std::vector<Vector3d>& a_rw){
  std::vector<double> data;
  hdf52doubles(h5f, dsetname, rows, 3, data);
  if (a_rw.size() < rows) a_rw.resize(rows);
  for (Uint irw=0; irw < rows; ++irw)
    a_rw[irw] = {data[3*irw], data[3*irw+1], data[3*irw+2]};
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

// Row-major
void hdf52tensors(const H5::H5File& h5f, const std::string& dsetname, const Uint rows,
                  std::vector<Matrix3d>& M_rw){
  std::vector<double> data;
  hdf52doubles(h5f, dsetname, rows, 9, data);
  if (M_rw.size() < rows) M_rw.resize(rows);
  for (Uint irw=0; irw < rows; ++irw)
    for (Uint i=0; i<3; ++i)
      for (Uint j=0; j<3; ++j)
        M_rw[irw](i, j) = data[irw*9 + 3*i + j];
}

// Uint, as the readers take it
void ulong2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<Uint>& a,
                const Uint rows, const Uint cols){
  const hsize_t dims[2] = {rows, cols};
  const hid_t type = partrac::h5_native_type<Uint>();
  const partrac::H5Id space(H5Screate_simple(2, dims, nullptr), H5Sclose);
  const partrac::H5Id dset(H5Dcreate2(h5f.getId(), dsetname.c_str(), type, space,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
  if (!dset.valid() || H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, a.data()) < 0)
    throw H5::DataSetIException("ulong2hdf5", "cannot write '" + dsetname + "'");
}

void ulong2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<Uint>& a, const Uint Nrw){
  ulong2hdf5(h5f, dsetname, a, Nrw, 1);
}

void int2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<int>& a, const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 1;
  H5::DataSpace dspace(2, dims);
  H5::DataSet dset = h5f.createDataSet(dsetname, H5::PredType::NATIVE_INT, dspace);
  dset.write(a.data(), H5::PredType::NATIVE_INT);
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


// Dataset shape; missing is an error
static partrac::H5DatasetInfo hdf5_info(const H5::H5File& h5f, const std::string& dsetname){
  if (H5Lexists(h5f.getId(), dsetname.c_str(), H5P_DEFAULT) <= 0)
    partrac::fail(h5f.getFileName(), ": no dataset '", dsetname, "'");
  const partrac::H5DatasetInfo info = partrac::h5_dataset_info(h5f.getId(), dsetname);
  if (info.rank() < 1 || info.rank() > 2)
    partrac::fail(h5f.getFileName(), ": '", dsetname, "' has rank ", info.rank(), ", expected 1 or 2");
  return info;
}

Uint hdf5_rows(const H5::H5File& h5f, const std::string& dsetname){
  return hdf5_info(h5f, dsetname).rows();
}

static void hdf5_check(const H5::H5File& h5f, const std::string& dsetname, const Uint rows, const Uint cols){
  const partrac::H5DatasetInfo info = hdf5_info(h5f, dsetname);
  if (info.rows() != rows || info.cols() != cols)
    partrac::fail(h5f.getFileName(), ": '", dsetname, "' is ", info.rows(), " x ", info.cols(),
                  ", expected ", rows, " x ", cols);
}

void hdf52doubles(const H5::H5File& h5f, const std::string& dsetname, const Uint rows, const Uint cols,
                  std::vector<double>& a){
  hdf5_check(h5f, dsetname, rows, cols);
  partrac::h5_read(h5f.getId(), dsetname, a);
}

void hdf52ulongs(const H5::H5File& h5f, const std::string& dsetname, const Uint rows, const Uint cols,
                 std::vector<Uint>& a){
  hdf5_check(h5f, dsetname, rows, cols);
  partrac::h5_read(h5f.getId(), dsetname, a);
}

void hdf52ints(const H5::H5File& h5f, const std::string& dsetname, const Uint rows, const Uint cols,
               std::vector<int>& a){
  hdf5_check(h5f, dsetname, rows, cols);
  partrac::h5_read(h5f.getId(), dsetname, a);
}
