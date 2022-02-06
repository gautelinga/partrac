#include <iostream>
#include <vector>
#include <map>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <filesystem>
#include "H5Cpp.h"
#include "io.hpp"

// using namespace std;

std::map<std::string, std::string> load_settings(std::string casefile){
  std::ifstream input(casefile);
  std::string line;
  std::string el;
  std::map<std::string, std::string> settings;

  while (getline(input, line)){
    size_t pos = line.find("=");

    // std::cout << line.substr(0, pos) << ": " << line.substr(pos+1) << std::endl;
    settings[line.substr(0, pos)] = line.substr(pos+1);
  }
  return settings;
}

std::vector<std::string> get_files(const std::string& s)
{
    std::vector<std::string> r;
    for(auto& p : std::filesystem::directory_iterator(s))
        if (!p.is_directory())
            r.push_back(p.path().string());
    return r;
}

std::vector<std::vector<std::string>> load_grid(std::string infile){
  std::ifstream input(infile);

  std::vector<char> x;
  std::string line;

  std::vector<std::vector<std::string>> sites;
  int nx = 0;
  while (getline(input, line)){
    // std::cout << line << std::endl;
    std::vector<std::string> sites_loc;
    boost::split(sites_loc, line, boost::is_any_of(" "));
    int sites_loc_size = sites_loc.size();
    nx = std::max(nx, sites_loc_size);
    sites.push_back(sites_loc);
  }
  int ny = sites.size();

  std::cout << "Size of grid: " << nx << " x " << ny << std::endl;
  for (int iy=0; iy < ny; ++iy){
    int nx_loc = sites[iy].size();
    for (int j=0; j < nx-nx_loc; ++j){
      sites[iy].push_back("s");
    }
  }
  return sites;
}

std::vector<std::vector<std::string>> load_fields(std::string infile){
  std::ifstream input(infile);

  std::vector<char> x;
  std::string line;

  std::vector<std::vector<std::string>> data;
  int nx = 0;
  while (getline(input, line)){
    // std::cout << line << std::endl;
    std::vector<std::string> data_loc;
    boost::split(data_loc, line, boost::is_any_of(" "));
    int data_loc_size = data_loc.size();
    nx = std::max(nx, data_loc_size);
    data.push_back(data_loc);
  }
  return data;
}

void print_grid(bool** grid, const int nx, const int ny){
  for (int iy=0; iy < ny; ++iy){
    for (int ix=0; ix < nx; ++ix){
      std::cout << grid[iy][ix];
    }
    std::cout << std::endl;
  }
}

void copy_arr(double*** f,
	      double*** f_prev,
	      bool** grid,
	      const int nx,
	      const int ny,
	      const int nc){
  for (int iy=0; iy < ny; ++iy){
    for (int ix=0; ix < nx; ++ix){
      if (grid[iy][ix]){
        for (int ic=0; ic < nc; ++ic){
          f_prev[iy][ix][ic] = f[iy][ix][ic];
        }
      }
    }
  }
}

void copy_arr(double** a_from,
	      double** a_to,
	      bool** grid,
	      const int nx,
	      const int ny){
  for (int iy=0; iy < ny; ++iy){
    for (int ix=0; ix < nx; ++ix){
      if (grid[iy][ix]){
        a_to[iy][ix] = a_from[iy][ix];
      }
    }
  }
}

void dump2file(std::string filename,
	       double** rho,
	       double** m_x,
	       double** m_y,
	       bool** grid,
	       const int nx,
	       const int ny){
  std::ofstream outfile(filename);
  outfile.precision(17);
  for (int iy=0; iy < ny; ++iy){
    for (int ix=0; ix < nx; ++ix){
      if (grid[iy][ix]){
        outfile << ix << " " << iy << " "
                << rho[iy][ix] << " " << m_x[iy][ix] << " " << m_y[iy][ix] << std::endl;
      }
    }
  }
  outfile.close();
}

std::string create_folder(const std::string folder){
  if (!std::filesystem::is_directory(folder)){
    std::filesystem::create_directory(folder);
  }
  return folder;
}

void verify_file_exists(const std::string infilename){
  if (!std::filesystem::exists(infilename)){
    std::cout << "No such file: " << infilename << std::endl;
    exit(0);
  }
}

// template this...
void print_param(const std::string key, const double val){
  std::cout << key << " = " << val << std::endl;
}

void print_param(const std::string key, const int val){
  std::cout << key << " = " << val << std::endl;
}

void print_param(const std::string key, const std::string val){
  std::cout << key << " = " << val << std::endl;
}

void write_param(std::ofstream &ofile, std::string key, double val){
  ofile << key << "=" << val << std::endl;
}

void write_param(std::ofstream &ofile, std::string key, int val){
  ofile << key << "=" << val << std::endl;
}

void write_param(std::ofstream &ofile, std::string key, long int val){
  ofile << key << "=" << val << std::endl;
}

void write_param(std::ofstream &ofile, std::string key, std::string val){
  ofile << key << "=" << val << std::endl;
}

// recently moved here
void load_vector_field(const std::string input_file,
                       std::vector<Vector3d> &pos_init){
  std::ifstream infile(input_file);
  double x, y, z;
  while (infile >> x >> y >> z){
    pos_init.push_back({x, y, z});
  }
  infile.close();
}

void dump_vector_field(const std::string output_file,
                       const std::vector<Vector3d>& x_rw,
                       const Uint Nrw){
  std::ofstream outfile(output_file);
  for (Uint irw=0; irw<Nrw; ++irw){
    outfile << std::setprecision(12)
            << x_rw[irw][0] << " "
            << x_rw[irw][1] << " "
            << x_rw[irw][2] << std::endl;
  }
  outfile.close();
}

void dump_vector_field(const std::string output_file,
                       const std::vector<Vector3d> &pos){
  std::ofstream outfile(output_file);
  for (std::vector<Vector3d>::const_iterator posit=pos.begin();
       posit != pos.end(); ++posit){
    outfile << std::setprecision(12)
            << (*posit)[0] << " "
            << (*posit)[1] << " "
            << (*posit)[2] << std::endl;
  }
  outfile.close();
}

void load_faces(const std::string input_file,
                FacesType& faces){
  std::ifstream infile(input_file);
  Uint first, second, third;
  double fourth;
  while (infile >> first >> second >> third >> fourth){
    faces.push_back({{first, second, third}, fourth});
  }
  infile.close();
}

void load_edges(const std::string input_file,
                EdgesType &edges){
  std::ifstream infile(input_file);
  Uint first, second;
  double third;
  while (infile >> first >> second >> third){
    edges.push_back({{first, second}, third});
  }
  infile.close();
}

void load_list(const std::string input_file,
               std::list<Uint> &li){
  std::ifstream infile(input_file);
  Uint a;
  while (infile >> a){
    li.push_back(a);
  }
  infile.close();
}

void dump_list(const std::string output_file,
               const std::list<Uint> &li){
  std::ofstream outfile(output_file);
  for (std::list<Uint>::const_iterator lit=li.begin();
       lit != li.end(); ++lit){
    outfile << *lit << std::endl;
  }
  outfile.close();
}

void dump_faces(const std::string output_file,
                const FacesType &faces){
  std::ofstream outfile(output_file);
  for (FacesType::const_iterator faceit = faces.begin();
       faceit != faces.end(); ++faceit){
    outfile << faceit->first[0] << " " << faceit->first[1] << " " << faceit->first[2]
            << " " << faceit->second << std::endl;
  }
  outfile.close();
}

void dump_edges(const std::string output_file,
                const EdgesType &edges){
  std::ofstream outfile(output_file);
  for (auto edgeit = edges.begin();
       edgeit != edges.end(); ++edgeit){
    outfile << edgeit->first[0] << " " << edgeit->first[1] << " " << edgeit->second << std::endl;
  }
  outfile.close();
}

void load_scalar_field(const std::string input_file,
                       std::vector<double>& c_rw, const Uint Nrw){
  // TODO: to hdf5
  std::ifstream infile(input_file);
  for (Uint irw=0; irw < Nrw; ++irw){
    infile >> c_rw[irw];
  }
  infile.close();
}

void dump_scalar_field(const std::string output_file,
                       const std::vector<double>& c_rw, const Uint Nrw){
                   // TODO: to hdf5
  std::ofstream outfile(output_file);
  for (Uint irw=0; irw < Nrw; ++irw){
    outfile << std::setprecision(12) << c_rw[irw] << std::endl;
  }
  outfile.close();
}

void tensor2hdf5(H5File& h5f, const std::string dsetname,
                 const std::vector<double>& axx_rw, const std::vector<double>& axy_rw, const std::vector<double>& axz_rw,
                 const std::vector<double>& ayx_rw, const std::vector<double>& ayy_rw, const std::vector<double>& ayz_rw,
                 const std::vector<double>& azx_rw, const std::vector<double>& azy_rw, const std::vector<double>& azz_rw,
                 const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 3*3;
  DataSpace dspace(2, dims);
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
  DataSet dset = h5f.createDataSet(dsetname,
                                    PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(data.data(), PredType::NATIVE_DOUBLE);
}

void vector2hdf5(H5File& h5f, const std::string dsetname,
                 const std::vector<double>& ax_rw, const std::vector<double>& ay_rw, const std::vector<double>& az_rw,
                 const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 3;
  DataSpace dspace(2, dims);
  std::vector<double> data(Nrw*3);
  for (Uint irw=0; irw < Nrw; ++irw){
    data[irw*3+0] = ax_rw[irw];
    data[irw*3+1] = ay_rw[irw];
    data[irw*3+2] = az_rw[irw];
  }
  DataSet dset = h5f.createDataSet(dsetname,
                                    PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(data.data(), PredType::NATIVE_DOUBLE);
}

void vector2hdf5(H5File& h5f, const std::string dsetname,
                 const std::vector<Vector3d>& a_rw, const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 3;
  DataSpace dspace(2, dims);
  std::vector<double> data(Nrw*3);
  for (Uint irw=0; irw < Nrw; ++irw){
    for (Uint d=0; d<3; ++d){
      data[irw*3+d] = a_rw[irw][d];
    }
  }
  DataSet dset = h5f.createDataSet(dsetname,
                                    PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(data.data(), PredType::NATIVE_DOUBLE);
}


void scalar2hdf5(H5File& h5f, const std::string dsetname, const std::vector<double>& c_rw,
                 const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 1;
  DataSpace dspace(2, dims);
  std::vector<double> data(Nrw);
  for (Uint irw=0; irw < Nrw; ++irw){
    data[irw] = c_rw[irw];
  }
  DataSet dset = h5f.createDataSet(dsetname,
                                    PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(data.data(), PredType::NATIVE_DOUBLE);
}

void print_mesh(const FacesType& faces, const EdgesType& edges,
                const Edge2FacesType& edge2faces,
                const std::vector<Vector3d>& x_rw, const Uint Nrw){
  std::cout << "==================" << std::endl;
  std::cout << "nodes" << std::endl;
  for (Uint irw=0; irw < Nrw; ++irw){
    std::cout << irw << " "
              << x_rw[irw][0] << " "
              << x_rw[irw][1] << " "
              << x_rw[irw][2] << std::endl;
  }

  std::cout << "edges" << std::endl;
  for (size_t ie=0; ie < edges.size(); ++ie){
    std::cout << ie << ": "
              << edges[ie].first[0] << " "
              << edges[ie].first[1] << " "
              << edges[ie].second << std::endl;
  }

  std::cout << "faces" << std::endl;
  for (size_t ifa=0; ifa < faces.size(); ++ifa){
    std::cout << ifa << ": "
              << faces[ifa].first[0] << " "
              << faces[ifa].first[1] << " "
              << faces[ifa].first[2] << " "
              << faces[ifa].second << std::endl;
  }

  std::cout << "edge2faces" << std::endl;
  for (size_t ie2f = 0; ie2f < edge2faces.size(); ++ie2f){
    std::cout << ie2f << ": ";
    for (FacesListType::const_iterator vit=edge2faces[ie2f].begin();
         vit != edge2faces[ie2f].end(); ++vit){
      std::cout << *vit << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "==================" << std::endl;
}

void dump_mesh(const FacesType& faces, const EdgesType& edges,
               const Edge2FacesType& edge2faces,
               const std::vector<Vector3d>& x_rw, const Uint Nrw){
  std::ofstream nodesf("mesh.node");
  for (Uint irw=0; irw < Nrw; ++irw){
    nodesf << x_rw[irw][0] << " "
           << x_rw[irw][1] << " "
           << x_rw[irw][2] << std::endl;
  }
  nodesf.close();

  std::ofstream edgesf("mesh.edge");
  for (size_t ie=0; ie < edges.size(); ++ie){
    edgesf << edges[ie].first[0] << " " << edges[ie].first[1] << " " << edges[ie].second << std::endl;
  }
  nodesf.close();

  std::ofstream facesf("mesh.face");
  for (size_t ifa=0; ifa < faces.size(); ++ifa){
    facesf << faces[ifa].first[0] << " " << faces[ifa].first[1] << " "
           << faces[ifa].first[2] << " " << faces[ifa].second << std::endl;
  }
  facesf.close();

  std::ofstream e2ff("mesh.e2f");
  for (size_t ie2f = 0; ie2f < edge2faces.size(); ++ie2f){
    for (FacesListType::const_iterator vit=edge2faces[ie2f].begin();
         vit != edge2faces[ie2f].end(); ++vit){
      e2ff << *vit << " ";
    }
    e2ff << std::endl;
  }
  e2ff.close();
}

void posdata2txt(std::ofstream &pos_out,
                 std::vector<Vector3d>& x_rw,
                 std::vector<Vector3d>& u_rw,
                 const Uint Nrw){
  for (Uint irw=0; irw < Nrw; ++irw){
    pos_out << irw << " "
            << x_rw[irw][0] << " "
            << x_rw[irw][1] << " "
            << x_rw[irw][2] << " "
            << u_rw[irw][0] << " "
            << u_rw[irw][1] << " "
            << u_rw[irw][2] << std::endl;
  }
}
