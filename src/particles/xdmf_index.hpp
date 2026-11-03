#ifndef __XDMF_INDEX_HPP
#define __XDMF_INDEX_HPP

#include <string>
#include <sys/types.h>
#include <vector>
#include "io.hpp"

// mesh.xdmf: the run's dumps as one temporal collection, a grid appended at each
class XdmfIndex {
public:
  explicit XdmfIndex(const std::string& folder);
  ~XdmfIndex();
  XdmfIndex(const XdmfIndex&) = delete;
  XdmfIndex& operator=(const XdmfIndex&) = delete;
  // New: the header and the closing tags
  void start();
  // Restart at the dump time tname: the grids before it, and the one at it unless
  // in the file just created, until a dump at it replaces it; from mesh.xdmf,
  // else the dump files
  void resume(const std::string& tname, const std::string& created);
  // The grid of a dump group from the datasets written, over the closing tags,
  // the grid's file closed; file relative to the folder
  void append(const std::vector<WrittenDataset>& written, const std::string& file, const std::string& group);
private:
  std::string folder, path;
  int fd = -1;
  off_t end = 0;     // where the closing tags begin
  std::string tr;    // the restart's group
  off_t at_tr = -1;  // where its grid begins, until the first dump
  bool cut(const double t_restart, const std::string& created);
  void rebuild(const double t_restart, const std::string& created);
  void put(const std::string& s, const off_t at);
};

#endif
