#include "xdmf_helpers.hpp"
#include <filesystem>

namespace pt = boost::property_tree;

void read_dataset_scalar(std::string& h5filename, std::string& field, std::vector<double>& data){
    H5::H5File h5file(h5filename, H5F_ACC_RDONLY);
    H5::DataSet dataset = h5file.openDataSet(field.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    // Move out, probably
    int rank = 2; // dataspace.getSimpleExtentNdims();
    std::vector<hsize_t> shape(rank);
    int ndims = dataspace.getSimpleExtentDims( shape.data(), NULL);
    shape[1] = 1;
    std::pair<std::int64_t, std::int64_t> range = dolfin::MPI::local_range(MPI_COMM_WORLD, shape[0]);

    // Hyperslab selection
    std::vector<hsize_t> offset(rank, 0);
    std::vector<hsize_t> count = shape;

    offset[0] = range.first;
    count[0] = range.second - range.first;

    // Allocate data for shape
    dataspace.selectHyperslab( H5S_SELECT_SET, count.data(), offset.data() );
    H5::DataSpace memspace( rank, count.data() );
    // memspace.selectHyperslab( H5S_SELECT_SET, count.data(), offset.data() );

    std::size_t data_size = 1;
    for (std::size_t i = 0; i < count.size(); ++i)
    {
      data_size *= count[i];
    }
    data.resize(data_size);

    dataset.read( data.data(), H5::PredType::NATIVE_DOUBLE, memspace, dataspace );
}

void reorder_indices(std::vector<double>& data_, const std::vector<double>& xdata, const std::vector<Uint>& j2i, const int dim){
  for ( Uint j=0; j < j2i.size(); ++j ){
    Uint i = j2i[j];
    for ( Uint k=0; k < dim; ++k)
      data_[dim*j+k] = xdata[dim*i+k];
  }
}

void read_dataset_vector(std::string& h5filename_u, std::string& field, std::vector<double>& data, const int dim){
    H5::H5File h5file_u(h5filename_u, H5F_ACC_RDONLY);
    H5::DataSet dataset = h5file_u.openDataSet(field.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    // Move out, probably
    int rank = 2; // dataspace.getSimpleExtentNdims();
    std::vector<hsize_t> shape(rank);
    int ndims = dataspace.getSimpleExtentDims( shape.data(), NULL);
    shape[1] = dim;
    std::pair<std::int64_t, std::int64_t> range = dolfin::MPI::local_range(MPI_COMM_WORLD, shape[0]);

    // Hyperslab selection
    std::vector<hsize_t> offset(rank, 0);
    std::vector<hsize_t> count = shape;

    offset[0] = range.first;
    count[0] = range.second - range.first;

    // Allocate data for shape
    dataspace.selectHyperslab( H5S_SELECT_SET, count.data(), offset.data() );
    H5::DataSpace memspace( rank, count.data() );
    // memspace.selectHyperslab( H5S_SELECT_SET, count.data(), offset.data() );

    std::size_t data_size = 1;
    for (std::size_t i = 0; i < count.size(); ++i)
    {
      data_size *= count[i];
    }
    data.resize(data_size);

    dataset.read( data.data(), H5::PredType::NATIVE_DOUBLE, memspace, dataspace );
}

std::vector<std::pair<double, std::vector<std::string>>> parse_xdmf(const std::string& xdmffilename, 
                                                                    std::string& h5filename,
                                                                    std::string& topology_path, 
                                                                    std::string& geometry_path)
{
  // Create empty property tree object
  pt::ptree tree;

  std::filesystem::path ppath = xdmffilename;
  std::string dirname = std::string(ppath.parent_path()) + "/";

  // Parse the XML into the property tree.
  pt::read_xml(xdmffilename, tree);

  topology_path = tree.get<std::string>("Xdmf.Domain.Grid.Grid.Topology.DataItem");
  geometry_path = tree.get<std::string>("Xdmf.Domain.Grid.Grid.Geometry.DataItem");

  auto topology_pos = topology_path.find(":");
  h5filename = topology_path.substr(0, topology_pos);
  topology_path.erase(0, topology_pos + 1);
  if (h5filename != geometry_path.substr(0, topology_pos)){
    std::cout << "XDMF error: Not matching filenames." << std::endl;
    exit(0);
  }
  geometry_path.erase(0, topology_pos + 1);

  h5filename = dirname + h5filename;

  std::vector<std::pair<double, std::vector<std::string>>> titems;
  for (auto & p : tree.get_child("Xdmf.Domain.Grid")) {
    //std :: cout << "[" << p.first << "]" << std :: endl;    
    if (p.first == "Grid"){
      double time;
      std::string filename;
      std::string location;

      for (auto & pp : p.second) {
        if ( pp.first == "Time" ){
          time = pp.second.get<double>("<xmlattr>.Value");
        }
        else if ( pp.first == "Attribute")
        {
          location = pp.second.get<std::string>("DataItem");
          auto pos = location.find(":");
          filename = location.substr(0, pos);
          location.erase(0, pos + 1);
          // Check filename == h5filename too
        }
      }
      //std::cout << " " << time << " " << location << std::endl;
      //titems.push_back({time, {filename, location}});
      std::vector<std::string> path = {dirname + filename, location};
      titems.push_back({time, path});
    }
  }
  return titems;
}

std::vector<std::pair<double, std::vector<std::string>>> parse_xdmf(const std::string& xdmffilename)
{
  std::string _dummy0, _dummy1, _dummy2;
  return parse_xdmf(xdmffilename, _dummy0, _dummy1, _dummy2);
}
