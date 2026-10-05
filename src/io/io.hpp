#ifndef __IO_HPP
#define __IO_HPP

#include <iostream>
#include <vector>
#include <map>
#include <fstream>

#include "H5Cpp.h"
#include "typedefs.hpp"

typedef std::shared_ptr<H5::H5File> H5FilePtr;

// Recently moved here

void load_vector_field(const std::string& input_file, std::vector<Vector3d> &pos_init);
void load_faces(const std::string& input_file,
                FacesType& faces);
void load_edges(const std::string& input_file,
                EdgesType &edges);
void load_list(const std::string& input_file,
               std::vector<Uint> &li);

void load_scalar_field(const std::string& input_file,
                       std::vector<double>& c_rw, const Uint Nrw);
void tensor2hdf5(H5::H5File& h5f, const std::string& dsetname,
                 const std::vector<double>& axx_rw, const std::vector<double>& axy_rw, const std::vector<double>& axz_rw,
                 const std::vector<double>& ayx_rw, const std::vector<double>& ayy_rw, const std::vector<double>& ayz_rw,
                 const std::vector<double>& azx_rw, const std::vector<double>& azy_rw, const std::vector<double>& azz_rw,
                 const Uint Nrw);
void vector2hdf5(H5::H5File& h5f, const std::string& dsetname,
                 const std::vector<double>& ax_rw, const std::vector<double>& ay_rw, const std::vector<double>& az_rw,
                 const Uint Nrw);
void vector2hdf5(H5::H5File& h5f, const std::string& dsetname,
                 const std::vector<Vector3d>& a_rw, const Uint Nrw);
void tensor2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<Matrix3d>& M_rw,
                 const Uint Nrw);
void ulong2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<Uint>& a,
                const Uint rows, const Uint cols);
void ulong2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<Uint>& a, const Uint Nrw);
void int2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<int>& a, const Uint Nrw);
void load_tensor_field(const std::string& input_file, std::vector<Matrix3d>& M_rw, const Uint Nrw);
void load_vector_field(const std::string& input_file, std::vector<Vector3d>& a_rw, const Uint Nrw);
void scalar2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<double>& c_rw,
                 const Uint Nrw);
Uint hdf5_rows(const H5::H5File& h5f, const std::string& dsetname);
// Read back as rows x cols
void hdf52doubles(const H5::H5File& h5f, const std::string& dsetname, const Uint rows, const Uint cols,
                  std::vector<double>& a);
void hdf52ulongs(const H5::H5File& h5f, const std::string& dsetname, const Uint rows, const Uint cols,
                 std::vector<Uint>& a);
void hdf52ints(const H5::H5File& h5f, const std::string& dsetname, const Uint rows, const Uint cols,
               std::vector<int>& a);
// Vectors and tensors, one a row
void hdf52vectors(const H5::H5File& h5f, const std::string& dsetname, const Uint rows,
                  std::vector<Vector3d>& a_rw);
void hdf52tensors(const H5::H5File& h5f, const std::string& dsetname, const Uint rows,
                  std::vector<Matrix3d>& M_rw);



#endif
