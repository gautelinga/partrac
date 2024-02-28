#ifndef __EXP_IO_HPP
#define __EXP_IO_HPP

#include "typedefs.hpp"

void tensor_to_h5(H5::H5File& h5f, const std::string dsetname,
                 const std::vector<Real>& axx_rw, const std::vector<Real>& axy_rw, const std::vector<Real>& axz_rw,
                 const std::vector<Real>& ayx_rw, const std::vector<Real>& ayy_rw, const std::vector<Real>& ayz_rw,
                 const std::vector<Real>& azx_rw, const std::vector<Real>& azy_rw, const std::vector<Real>& azz_rw,
                 const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 3*3;
  H5::DataSpace dspace(2, dims);
  std::vector<Real> data(Nrw*3*3);
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

void vector_to_h5(H5::H5File& h5f, const std::string dsetname,
                 const std::vector<Real>& ax_rw, const std::vector<Real>& ay_rw, const std::vector<Real>& az_rw,
                 const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 3;
  H5::DataSpace dspace(2, dims);
  std::vector<Real> data(Nrw*3);
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

void vector_to_h5(H5::H5File& h5f, const std::string dsetname,
                 const std::vector<Vector>& a_rw, const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 3;
  H5::DataSpace dspace(2, dims);
  std::vector<Real> data(Nrw*3);
  for (Uint irw=0; irw < Nrw; ++irw){
    for (Uint d=0; d<3; ++d){
      data[irw*3+d] = a_rw[irw][d];
    }
  }
  H5::DataSet dset = h5f.createDataSet(dsetname,
                                    H5::PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}

void vector_to_h5(H5::H5File& h5f, const std::string dsetname,
                 const std::vector<Real>& a, const Uint dim = 3){
  hsize_t dims[2];
  dims[0] = a.size()/dim;
  dims[1] = dim;
  H5::DataSpace dspace(2, dims);
  H5::DataSet dset = h5f.createDataSet(dsetname, H5::PredType::NATIVE_DOUBLE, dspace);
  dset.write(a.data(), H5::PredType::NATIVE_DOUBLE);
}

void vector_to_h5(H5::H5File& h5f, const std::string dsetname,
                 const std::vector<int>& a, const Uint dim = 3){
  hsize_t dims[2];
  dims[0] = a.size()/dim;
  dims[1] = dim;
  H5::DataSpace dspace(2, dims);
  H5::DataSet dset = h5f.createDataSet(dsetname, H5::PredType::NATIVE_INT, dspace);
  dset.write(a.data(), H5::PredType::NATIVE_INT);
}

void scalar_to_h5(H5::H5File& h5f
                , const std::string dsetname
                , const std::vector<Real>& c_rw
                , const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 1;
  H5::DataSpace dspace(2, dims);
  std::vector<Real> data(Nrw);
  for (Uint irw=0; irw < Nrw; ++irw){
    data[irw] = c_rw[irw];
  }
  H5::DataSet dset = h5f.createDataSet(dsetname,
                                    H5::PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}

void scalar_to_h5(H5::H5File& h5f, const std::string dsetname, const std::vector<Real>& c){
  vector_to_h5(h5f, dsetname, c, 1);
}

void scalar_to_h5(H5::H5File& h5f, const std::string dsetname, const std::vector<int>& c){
  vector_to_h5(h5f, dsetname, c, 1);
}

#endif