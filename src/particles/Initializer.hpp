#ifndef __DISTRIBUTE_HPP
#define __DISTRIBUTE_HPP

#include <algorithm>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include "typedefs.hpp"
#include "Interpol.hpp"
#include "strings.hpp"
#include "Params.hpp"

// true if init_mode starts with one of the given kinds
inline bool init_mode_is(const partrac::Params& prm, const std::vector<std::string>& kinds){
  const std::string kind = split_string(prm.get<std::string>("init_mode"), "_")[0];
  return std::find(kinds.begin(), kinds.end(), kind) != kinds.end();
}

// The initializers index the init_mode tokens without bounds checks
inline bool init_mode_shape_ok(const std::string& init_mode){
  if (init_mode.rfind("from", 0) == 0)
    return init_mode.rfind("from_file:", 0) == 0 && init_mode.size() > 10;
  const std::vector<std::string> key = split_string(init_mode, "_");
  if (!init_mode_dirs_ok(key)) return false;
  if (key[0] == "point") return key.size() == 1;
  if (key[0] == "pairs") return key.size() == 2 || key.size() == 3;
  if (key[0] == "randomgaussianstrip") return key.size() == 3;
  if (key[0] == "randomgaussiancircle") return key.size() == 2 || key.size() == 3;
  return key.size() == 2;
}

// Initializer parameters; La, Lb, ds_init and init_weight only for some init_modes
inline void add_initializer_params(partrac::Schema& s){
  s.require<std::string>("init_mode", "initial distribution");
  s.require<Uint>("Nrw", "number of particles");
  // Actual counts; Nrw is the request
  s.runtime<Uint>("Nrw_init", 0, "particles the initializer placed");
  s.runtime<Uint>("Nrw_current", 0, "particles in the set when this was written");
  s.require<Uint>("Nrw_max", "max number of particles");
  // ds bounds only for edges or remeshing
  auto has_edges_or_remeshes = [](const partrac::Params& p){
    const bool cloud = init_mode_is(p, {"point", "points", "randomgaussianstrip", "randomgaussiancircle"});
    const bool remeshes = (p.has("refine") && p.get<bool>("refine")) || (p.has("coarsen") && p.get<bool>("coarsen"));
    return !cloud || remeshes;
  };
  s.require_if<double>("ds_max", 0.0, has_edges_or_remeshes,
                       "the initial state has edges, or refine/coarsen is on", "max edge length");
  s.require_if<double>("ds_min", 0.0, has_edges_or_remeshes,
                       "the initial state has edges, or refine/coarsen is on", "min edge length");
  s.require_if<double>("La",
                       [](const partrac::Params& p){
                         return init_mode_is(p, {"strip", "sheet", "ellipsoid",
                                                 "randomgaussianstrip",
                                                 "randomgaussiancircle"});
                       },
                       "init_mode is a strip, sheet, ellipsoid or gaussian",
                       "principal extent");
  s.require_if<double>("Lb",
                       [](const partrac::Params& p){
                         return init_mode_is(p, {"sheet", "ellipsoid",
                                                 "randomgaussianstrip",
                                                 "randomgaussiancircle"});
                       },
                       "init_mode is a sheet, ellipsoid or gaussian",
                       "second extent");
  s.require_if<double>("ds_init",
                       [](const partrac::Params& p){
                         return init_mode_is(p, {"sheet", "pair", "pairs", "points"});
                       },
                       "init_mode is a sheet, pair or points distribution",
                       "initial edge length");
  // uniform steps across the domain
  s.check([](const partrac::Params& p){
            return !init_mode_is(p, {"uniform"}) || p.get<Uint>("Nrw") >= 2;
          },
          "init_mode uniform needs Nrw of 2 or more");
  s.require_if<std::string>("init_weight",
                            [](const partrac::Params& p){
                              return init_mode_is(p, {"points"});
                            },
                            "init_mode is a points distribution",
                            "sampling weight");
  s.opt<double>("t0", 0.0, "start time");
  s.opt<double>("x0", 0.0, "initial position");
  s.opt<double>("y0", 0.0, "initial position");
  s.opt<double>("z0", 0.0, "initial position");
  s.opt<bool>("inject", false, "inject new particles");
  s.opt<bool>("clear_initial_edges", false, "drop the initial edges");
  s.check([](const partrac::Params& p){
            return init_mode_shape_ok(p.get<std::string>("init_mode"));
          },
          "init_mode has the wrong number of directions: most modes take one,"
            " as in uniform_x; point takes none; pairs takes one or two;"
            " randomgaussianstrip takes two, as in randomgaussianstrip_x_y;"
            " randomgaussiancircle one or two, the second the spread's directions;"
            " from_file takes a path, as in from_file:positions.h5");
}

class Initializer {
public:
  Initializer(std::shared_ptr<Interpol> intp, partrac::Params& prm) : intp(intp), prm(prm) {
    x0 = {prm.get<double>("x0"), prm.get<double>("y0"), prm.get<double>("z0")};
    x_min = intp->get_x_min();
    x_max = intp->get_x_max();
    L = x_max - x_min;
    inject = prm.get<bool>("inject");
    clear_initial_edges = prm.get<bool>("clear_initial_edges");
  };
  ~Initializer() { nodes.clear(); edges.clear(); faces.clear(); };
  /*std::vector<Vector3d>::const_iterator node_begin() const { return nodes.begin(); };
  std::vector<Vector3d>::const_iterator node_end() const { return nodes.end(); };
  EdgesType::const_iterator edge_begin() const { return edges.begin(); };
  EdgesType::const_iterator edge_end() const { return edges.end(); };
  FacesType::const_iterator face_begin() const { return faces.begin(); };
  FacesType::const_iterator face_end() const { return faces.end(); };*/
  std::vector<Vector3d> nodes;
  EdgesType edges;
  FacesType faces;
  bool inject;
  bool clear_initial_edges;
protected:
  std::shared_ptr<Interpol> intp;
  partrac::Params& prm;
  Vector3d x0;
  Vector3d x_min;
  Vector3d x_max;
  Vector3d L;
};

// Initial state for init_mode
void set_initial_state(std::shared_ptr<Initializer>& init_state, std::shared_ptr<Interpol> intp, partrac::Params& prm, std::mt19937& gen);

// A gaussian strip (dim 2) or circle (dim 3), whatever init_mode's shape
std::shared_ptr<Initializer> make_gaussian_initializer(const Uint dim, const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, partrac::Params& prm, std::mt19937& gen);

#endif
