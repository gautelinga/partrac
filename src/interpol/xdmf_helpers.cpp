#include "Error.hpp"
#include "h5direct.hpp"
#include "xdmf_helpers.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <filesystem>
#include <iostream>

namespace pt = boost::property_tree;

void read_dataset_columns(const std::string& h5filename, const std::string& field,
                          std::vector<double>& data, const int ncols){
  const partrac::H5Id file = partrac::h5_open_read(h5filename);
  const partrac::H5DatasetInfo info = partrac::h5_dataset_info(file, field);
  if (info.rank() != 2){
    partrac::fail("XDMF: '", field, "' in ", h5filename, " has rank ", info.rank(), ", expected 2.");
  }
  partrac::h5_read(file, field, data, static_cast<std::size_t>(ncols));
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
    partrac::fail("XDMF: the velocity and pressure files do not match");
  }
  geometry_path.erase(0, topology_pos + 1);

  h5filename = dirname + h5filename;

  std::vector<std::pair<double, std::vector<std::string>>> titems;
  for (auto & p : tree.get_child("Xdmf.Domain.Grid")) {
    //std :: cout << "[" << p.first << "]" << std :: endl;    
    if (p.first == "Grid"){
      double time = 0.;
      bool has_time = false;
      std::string filename;
      std::string location;

      for (auto & pp : p.second) {
        if ( pp.first == "Time" ){
          time = pp.second.get<double>("<xmlattr>.Value");
          has_time = true;
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
      if (!has_time){
        partrac::fail("XDMF: a grid without a time");
      }
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
