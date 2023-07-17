#ifndef __IO_HPP
#define __IO_HPP

#include <iostream>
#include <vector>
#include <map>
#include <filesystem>
#include <fstream>

#include "H5Cpp.h"
#include "typedefs.hpp"

//using namespace std;
//using namespace H5;

std::map<std::string, std::string> load_settings(const std::string& casefile);
std::vector<std::string> get_files(const std::string& s);
std::vector<std::vector<std::string>> load_grid(const std::string& infile);
std::vector<std::vector<std::string>> load_fields(const std::string& infile);
void print_grid(bool** grid, const int nx, const int ny);
void copy_arr(double*** f,
	      double*** f_prev,
	      bool** grid,
	      const int nx,
	      const int ny,
	      const int nc);
void copy_arr(double** a_from,
	      double** a_to,
	      bool** grid,
	      const int nx,
	      const int ny);
void dump2file(const std::string& filename,
	       double** rho,
	       double** m_x,
	       double** m_y,
	       bool** grid,
	       const int nx,
	       const int ny);
std::string create_folder(const std::string& folder);
void verify_file_exists(const std::string& infilename);
void print_param(const std::string& key, const double val);
void print_param(const std::string& key, const int val);
void print_param(const std::string& key, const std::string& val);
void write_param(std::ofstream &ofile, const std::string& key, double val);
void write_param(std::ofstream &ofile, const std::string& key, int val);
void write_param(std::ofstream &ofile, const std::string& key, long int val);
void write_param(std::ofstream &ofile, const std::string& key, const std::string& val);
/*void load_field(H5::H5File &h5file,
                //double*** u,
                std::vector<double>&,
                const std::string field
                //, const int nx, const int ny, const int nz
                );
void load_int_field(H5::H5File &h5file,
                    //int*** u,
                    std::vector<int>&,
                    const std::string field
                    //, const int nx, const int ny, const int nz
                    );
void load_h5(const std::string h5filename,
             //double*** ux,
             //double*** uy,
             //double*** uz,
             //double*** rho,
             //double*** p,
             //const int nx,
             //const int ny,
             //const int nz,
             std::vector<double>& ux,
             std::vector<double>& uy,
             std::vector<double>& uz,
             std::vector<double>& rho,
             std::vector<double>& p,
             const bool verbose, const bool, const bool, const bool);*/

// Recently moved here

void load_vector_field(const std::string& input_file, std::vector<Vector3d> &pos_init);
void dump_vector_field(const std::string& output_file,
                       const std::vector<Vector3d>& x_rw,
                       const Uint Nrw);
void dump_vector_field(const std::string& output_file,
                       const std::vector<Vector3d> &pos);
void load_faces(const std::string& input_file,
                FacesType& faces);
void load_edges(const std::string& input_file,
                EdgesType &edges);
void load_list(const std::string& input_file,
               std::list<Uint> &li);
void dump_list(const std::string& output_file,
               const std::list<Uint> &li);
void dump_faces(const std::string& output_file,
                const FacesType &faces);
void dump_edges(const std::string& output_file,
                const EdgesType &edges);

void load_scalar_field(const std::string& input_file,
                       std::vector<double>& c_rw, const Uint Nrw);
void dump_scalar_field(const std::string& output_file,
                       const std::vector<double>& c_rw, const Uint Nrw);
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
void scalar2hdf5(H5::H5File& h5f, const std::string& dsetname, const std::vector<double>& c_rw,
                 const Uint Nrw);
void print_mesh(const FacesType& faces, const EdgesType& edges,
                const Edge2FacesType& edge2faces,
                const std::vector<Vector3d>& x_rw, const Uint Nrw);
void dump_mesh(const FacesType& faces, const EdgesType& edges,
               const Edge2FacesType& edge2faces,
               const std::vector<Vector3d>& x_rw, const Uint Nrw);
void posdata2txt(std::ofstream &pos_out,
                 std::vector<Vector3d>& x_rw,
                 std::vector<Vector3d>& u_rw,
                 const Uint Nrw);



#endif
