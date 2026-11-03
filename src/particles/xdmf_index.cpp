#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <map>
#include <tuple>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include "Error.hpp"
#include "h5direct.hpp"
#include "io.hpp"
#include "xdmf_index.hpp"

namespace {

const std::string header =
  "<?xml version=\"1.0\"?>\n"
  "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
  "<Xdmf Version=\"3.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n"
  "  <Domain>\n"
  "    <Grid Name=\"Timeseries\" GridType=\"Collection\" CollectionType=\"Temporal\">";
const std::string footer = "\n    </Grid>\n  </Domain>\n</Xdmf>\n";
const std::string grid_open = "<Grid Name=\"mesh\" GridType=\"Uniform\">";

// Names in a group
std::vector<std::string> links(const hid_t loc, const std::string& group){
  std::vector<std::string> out;
  H5G_info_t gi;
  if (H5Gget_info_by_name(loc, group.c_str(), &gi, H5P_DEFAULT) < 0)
    return out;
  for (hsize_t k = 0; k < gi.nlinks; ++k){
    const ssize_t n = H5Lget_name_by_idx(loc, group.c_str(), H5_INDEX_NAME, H5_ITER_INC, k, nullptr, 0, H5P_DEFAULT);
    if (n <= 0) continue;
    std::string name(std::size_t(n), '\0');
    H5Lget_name_by_idx(loc, group.c_str(), H5_INDEX_NAME, H5_ITER_INC, k, &name[0], std::size_t(n) + 1, H5P_DEFAULT);
    out.push_back(name);
  }
  return out;
}

using Dsets = std::map<std::string, WrittenDataset>;

// Datasets of a group in a file, by name; per node where the rows are the
// points', else per face or edge where they match those
Dsets datasets(const hid_t file, const std::string& group){
  Dsets out;
  for (const auto& name : links(file, group)){
    const std::string p = group + "/" + name;
    {
      const partrac::H5Id obj(H5Oopen(file, p.c_str(), H5P_DEFAULT), H5Oclose);
      if (!obj.valid() || H5Iget_type(obj) != H5I_DATASET) continue;
    }
    const partrac::H5DatasetInfo d = partrac::h5_dataset_info(file, p);
    if (d.rank() < 1 || d.rank() > 2) continue;
    if (d.type_class != H5T_FLOAT && d.type_class != H5T_INTEGER) continue;
    const char kind = d.type_class == H5T_FLOAT ? 'f' : d.type_signed ? 'i' : 'u';
    out[name] = {name, Uint(d.rows()), Uint(d.cols()), kind, unsigned(d.type_size), 0};
  }
  const auto rows = [&](const char* name){ const auto it = out.find(name); return it == out.end() ? -1 : long(it->second.rows); };
  const long nodes = rows("points"), faces = rows("faces"), edges = rows("edges");
  for (auto& [name, d] : out)
    d.centre = long(d.rows) == nodes ? 'n' : long(d.rows) == faces ? 'f' : long(d.rows) == edges ? 'e' : 0;
  return out;
}

// Empty if the group is no dump
std::string grid_of(const Dsets& ds, const std::string& rel, const std::string& group){
  const auto points = ds.find("points");
  if (points == ds.end() || points->second.cols != 3)
    return "";
  const auto item = [&](const std::string& name, const WrittenDataset& d){
    const char* type = d.kind == 'f' ? "Float" : d.kind == 'i' ? "Int" : "UInt";
    return "<DataItem Dimensions=\"" + std::to_string(d.rows) + " " + std::to_string(d.cols)
      + "\" NumberType=\"" + type + "\" Precision=\"" + std::to_string(d.size) + "\" Format=\"HDF\">"
      + rel + ":/" + group + "/" + name + "</DataItem>";
  };
  const Uint nodes = points->second.rows;
  const auto faces = ds.find("faces"), edges = ds.find("edges");
  const auto cells = faces != ds.end() ? faces : edges;
  std::string s = "\n      " + grid_open;
  if (cells != ds.end()){
    const bool tri = cells == faces;
    s += "\n        <Topology NumberOfElements=\"" + std::to_string(cells->second.rows)
      + (tri ? "\" TopologyType=\"Triangle\" NodesPerElement=\"3\">"
             : "\" TopologyType=\"Polyline\" NodesPerElement=\"2\">")
      + "\n          " + item(cells->first, cells->second) + "\n        </Topology>";
  }
  else
    s += "\n        <Topology TopologyType=\"Polyvertex\" NodesPerElement=\"" + std::to_string(nodes)
      + "\">\n        </Topology>";
  s += "\n        <Geometry GeometryType=\"XYZ\">\n          " + item("points", points->second)
    + "\n        </Geometry>";
  s += "\n        <Time Value=\"" + group + "\" />";
  // Per node, or per face or edge
  for (const auto& [name, d] : ds){
    if (name == "points" || name == "faces" || name == "edges")
      continue;
    const bool cell = d.centre != 'n';
    if (cell ? (cells == ds.end() || d.centre != (cells == faces ? 'f' : 'e') || d.rows != cells->second.rows)
             : d.rows != nodes)
      continue;
    const char* kind = d.cols == 1 ? "Scalar" : d.cols == 3 ? "Vector" : d.cols == 9 ? "Tensor" : nullptr;
    if (!kind)
      continue;
    s += "\n        <Attribute Name=\"" + name + "\" AttributeType=\"" + kind + "\" Center=\""
      + (cell ? "Cell" : "Node") + "\">\n          " + item(name, d) + "\n        </Attribute>";
  }
  return s + "\n      </Grid>";
}

// HDF5's error printing off in its scope
class H5Quiet {
public:
  H5Quiet(){
    H5Eget_auto2(H5E_DEFAULT, &func_, &data_);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  }
  ~H5Quiet(){ H5Eset_auto2(H5E_DEFAULT, func_, data_); }
  H5Quiet(const H5Quiet&) = delete;
  H5Quiet& operator=(const H5Quiet&) = delete;
private:
  H5E_auto2_t func_ = nullptr;
  void* data_ = nullptr;
};

bool blank(const char c){ return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }

}

XdmfIndex::XdmfIndex(const std::string& folder)
  : folder(folder), path(folder + "/mesh.xdmf") {}

XdmfIndex::~XdmfIndex(){
  if (fd >= 0)
    close(fd);
}

void XdmfIndex::put(const std::string& s, const off_t at){
  std::size_t done = 0;
  while (done < s.size()){
    const ssize_t n = pwrite(fd, s.data() + done, s.size() - done, at + off_t(done));
    if (n < 0 && errno == EINTR) continue;
    if (n <= 0)
      partrac::fail("cannot write ", path, ": ", std::strerror(errno));
    done += std::size_t(n);
  }
}

void XdmfIndex::start(){
  fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (fd < 0)
    partrac::fail("cannot create ", path, ": ", std::strerror(errno));
  put(header + footer, 0);
  end = off_t(header.size());
}

// A reader during the write may see it incomplete
void XdmfIndex::append(const std::vector<WrittenDataset>& written, const std::string& file,
                       const std::string& group){
  if (fd < 0)
    return;
  Dsets ds;
  for (const auto& d : written)
    if (d.path.size() > group.size() + 1 && d.path.compare(0, group.size() + 1, group + "/") == 0)
      ds[d.path.substr(group.size() + 1)] = d;
  const std::string grid = grid_of(ds, file, group);
  // The restart's grid, dumped again
  if (at_tr >= 0 && group == tr){
    if (ftruncate(fd, at_tr) != 0)
      partrac::fail("cannot truncate ", path, ": ", std::strerror(errno));
    end = at_tr;
  }
  at_tr = -1;
  put(grid + footer, end);
  end += off_t(grid.size());
}

void XdmfIndex::resume(const std::string& tname, const std::string& created){
  tr = tname;
  const double t_restart = std::strtod(tname.c_str(), nullptr);
  fd = open(path.c_str(), O_RDWR | O_CREAT, 0644);
  if (fd < 0)
    partrac::fail("cannot open ", path, ": ", std::strerror(errno));
  if (!cut(t_restart, created))
    rebuild(t_restart, created);
}

// The grids from the end back to the last kept, the tail read in blocks; false
// where it is not as written
bool XdmfIndex::cut(const double t_restart, const std::string& created){
  const off_t size = lseek(fd, 0, SEEK_END);
  if (size < off_t(header.size() + footer.size()))
    return false;
  // The file from base on
  std::string s;
  off_t base = size;
  const auto more = [&]{
    if (base == 0)
      return false;
    const off_t n = std::min(base, std::max(off_t(1) << 16, off_t(s.size())));
    std::string b(std::size_t(n), '\0');
    for (off_t done = 0; done < n; ){
      const ssize_t r = pread(fd, &b[std::size_t(done)], std::size_t(n - done), base - n + done);
      if (r < 0 && errno == EINTR) continue;
      if (r <= 0)
        return false;
      done += r;
    }
    s.insert(0, b);
    base -= n;
    return true;
  };
  const auto at = [&](const off_t a){ return std::size_t(a - base); };
  // Back over the blanks before offset a
  const auto back = [&](const off_t a){
    while (true){
      std::size_t r = at(a);
      while (r > 0 && blank(s[r - 1])) --r;
      if (r > 0 || !more())
        return base + off_t(r);
    }
  };
  more();
  if (s.size() < footer.size() || s.compare(s.size() - footer.size(), footer.size(), footer) != 0)
    return false;
  off_t hi = size - off_t(footer.size());   // the end of what is left
  while (true){
    // The last grid that opens before hi
    off_t k = -1;
    while (true){
      if (hi - base >= off_t(grid_open.size())){
        const std::size_t q = s.rfind(grid_open, at(hi) - grid_open.size());
        if (q != std::string::npos){
          k = base + off_t(q);
          break;
        }
      }
      if (!more())
        break;
    }
    if (k < 0){
      // None: back to the header
      const std::string temporal = "CollectionType=\"Temporal\">";
      const std::size_t p = s.find(temporal);
      if (base != 0 || p == std::string::npos || off_t(p + temporal.size()) != hi)
        return false;
      break;
    }
    const std::size_t e = s.find("</Grid>", at(k)), tv = s.find("<Time Value=\"", at(k));
    const std::size_t hdf = s.find("Format=\"HDF\"", at(k));
    const std::size_t gt = s.find('>', hdf), colon = s.find(':', gt);
    if (e == std::string::npos || base + off_t(e + 7) != hi || tv > e || hdf > e || colon > e)
      return false;
    const double t = std::strtod(s.c_str() + tv + 13, nullptr);
    std::size_t f0 = gt + 1, f1 = colon;
    while (f0 < f1 && blank(s[f0])) ++f0;
    while (f1 > f0 && blank(s[f1 - 1])) --f1;
    if (t > t_restart || s.compare(f0, f1 - f0, created) == 0){
      hi = back(k);
      continue;
    }
    if (t == t_restart)
      at_tr = back(k);
    break;
  }
  if (ftruncate(fd, hi) != 0)
    partrac::fail("cannot truncate ", path, ": ", std::strerror(errno));
  put(footer, hi);
  end = hi;
  return true;
}

// For each time before the restart, or at it outside the file created: the
// group of the file started last but not after it, the history that holds it
void XdmfIndex::rebuild(const double t_restart, const std::string& created){
  namespace fs = std::filesystem;
  const H5Quiet quiet;
  std::map<double, std::tuple<double, std::string, std::string>> found;
  std::error_code ec;
  const std::string prefix = "data_from_t";
  for (const auto& entry : fs::directory_iterator(folder, ec)){
    const std::string name = entry.path().filename().string();
    if (name.rfind(prefix, 0) != 0 || name == created)
      continue;
    char* stop = nullptr;
    const double start = std::strtod(name.c_str() + prefix.size(), &stop);
    if (stop == name.c_str() + prefix.size() || std::string(stop) != ".h5")
      continue;
    const partrac::H5Id file(H5Fopen(entry.path().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    if (!file.valid())
      continue;
    for (const auto& group : links(file, "/")){
      const double t = std::strtod(group.c_str(), &stop);
      if (stop == group.c_str() || t > t_restart || start > t)
        continue;
      const auto it = found.find(t);
      if (it == found.end() || std::tie(std::get<0>(it->second), std::get<1>(it->second)) < std::tie(start, name))
        found[t] = {start, name, group};
    }
  }
  std::string s = header;
  for (const auto& [t, where] : found){
    const auto& [start, name, group] = where;
    const partrac::H5Id file(H5Fopen((folder + "/" + name).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    if (!file.valid())
      continue;
    const std::string grid = grid_of(datasets(file, group), name, group);
    if (t == t_restart && !grid.empty())
      at_tr = off_t(s.size());
    s += grid;
  }
  if (ftruncate(fd, 0) != 0)
    partrac::fail("cannot truncate ", path, ": ", std::strerror(errno));
  put(s + footer, 0);
  end = off_t(s.size());
}
