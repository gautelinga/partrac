#ifndef __INITIALIZER_HPP
#define __INITIALIZER_HPP

#include <algorithm>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include "typedefs.hpp"
#include "Interpol.hpp"
#include "strings.hpp"
#include "Params.hpp"

// Nodes, edges and faces to start from
struct InitialState {
  std::vector<Vector3d> nodes;
  EdgesType edges;
  FacesType faces;
};

// Initial states, one per init_mode; key is init_mode split at '_'

// Nrw copies of (x0, y0, z0)
InitialState init_point(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen);
// Nrw points across the domain, joined
InitialState init_uniform(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen);
// Nrw points on a line of length La, joined
InitialState init_strip(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen);
// An La by Lb rectangle, triangulated
InitialState init_sheet(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen);
// An ellipsoid, La along the normal and Lb in the plane, triangulated
InitialState init_ellipsoid(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen);
// Randomly oriented pairs ds_init apart
InitialState init_pairs(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen);
// Nrw random points weighted by init_weight; joined in order along one axis
InitialState init_points(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen);
// Nrw random points on a line of length La, spread by Lb
InitialState init_gaussian_strip(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen);
// Nrw random points on a disc of diameter La, spread by Lb
InitialState init_gaussian_circle(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen);
// The rows of an HDF5 file's nodes, joined
InitialState init_file(const std::string& path, std::shared_ptr<Interpol> intp, const partrac::Params& prm);

// An init_mode
struct InitMode {
  std::string name;
  std::vector<Uint> n_tokens;          // allowed token counts, the name included
  std::vector<std::string> needs;      // of La, Lb, ds_init and init_weight
  bool cloud;                          // no initial edges: ds_max and ds_min only for remeshing
  InitialState (*build)(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen);
};

// All init_modes but file:<path>
inline const std::vector<InitMode> init_modes = {
  // name                  tokens  needs                     cloud  build
  // point
  {"point",                {1},    {},                       true,  init_point},
  // uniform_x
  {"uniform",              {2},    {},                       false, init_uniform},
  // strip_x
  {"strip",                {2},    {"La"},                   false, init_strip},
  // sheet_xy
  {"sheet",                {2},    {"La", "Lb", "ds_init"},  false, init_sheet},
  // ellipsoid_xy
  {"ellipsoid",            {2},    {"La", "Lb"},             false, init_ellipsoid},
  // pair_xyz: one pair at (x0, y0, z0)
  {"pair",                 {2},    {"ds_init"},              false, init_pairs},
  // pairs_xyz, or pairs_xyz_xy with the centres drawn along x and y
  {"pairs",                {2, 3}, {"ds_init"},              false, init_pairs},
  // points_xy; along one axis, as points_x, also ds_init and edges
  {"points",               {2},    {"init_weight"},          true,  init_points},
  // randomgaussianstrip_x_y: along x, spread along y
  {"randomgaussianstrip",  {3},    {"La", "Lb"},             true,  init_gaussian_strip},
  // randomgaussiancircle_x or _x_yz: normal x, spread along y and z
  {"randomgaussiancircle", {2, 3}, {"La", "Lb"},             true,  init_gaussian_circle},
};

// The table entry named name, or nullptr
inline const InitMode* find_init_mode(const std::string& name){
  for (const auto& mode : init_modes)
    if (mode.name == name)
      return &mode;
  return nullptr;
}

// The table entry for init_mode, or nullptr
inline const InitMode* find_init_mode(const partrac::Params& prm){
  return find_init_mode(split_string(prm.get<std::string>("init_mode"), "_")[0]);
}

// points_x, points_y or points_z
inline bool points_along_one_axis(const std::vector<std::string>& key){
  return key.size() == 2 && key[0] == "points" && key[1].size() == 1;
}

// true if init_mode needs param; points along one axis also needs ds_init
inline bool init_mode_needs(const partrac::Params& prm, const std::string& param){
  const InitMode* mode = find_init_mode(prm);
  if (param == "ds_init" && points_along_one_axis(split_string(prm.get<std::string>("init_mode"), "_")))
    return true;
  return mode && std::find(mode->needs.begin(), mode->needs.end(), param) != mode->needs.end();
}

// true if the initial state has edges; points along one axis has
inline bool init_mode_has_edges(const partrac::Params& prm){
  const InitMode* mode = find_init_mode(prm);
  return !mode || !mode->cloud || points_along_one_axis(split_string(prm.get<std::string>("init_mode"), "_"));
}

// "a, b or c"
inline std::string join_words(const std::vector<std::string>& words, const std::string& last){
  std::string text;
  for (Uint i=0; i < words.size(); ++i)
    text += (i == 0 ? "" : i+1 == words.size() ? " " + last + " " : ", ") + words[i];
  return text;
}

// The init_modes that need param, as "strip, sheet or ellipsoid"
inline std::string init_modes_needing(const std::string& param){
  std::vector<std::string> names;
  for (const auto& mode : init_modes)
    if (std::find(mode.needs.begin(), mode.needs.end(), param) != mode.needs.end())
      names.push_back(mode.name);
  return join_words(names, "or");
}

// Directions each init_mode takes, as "point takes none; uniform and strip take one"
inline std::string init_mode_directions(){
  const std::vector<std::string> counts = {"none", "one", "two"};
  std::string text;
  std::vector<bool> listed(init_modes.size(), false);
  for (Uint i=0; i < init_modes.size(); ++i){
    if (listed[i])
      continue;
    std::vector<std::string> names;
    for (Uint j=i; j < init_modes.size(); ++j){
      if (init_modes[j].n_tokens == init_modes[i].n_tokens){
        names.push_back(init_modes[j].name);
        listed[j] = true;
      }
    }
    std::vector<std::string> n;
    for (const Uint n_tokens : init_modes[i].n_tokens)
      n.push_back(counts[n_tokens-1]);
    text += (text.empty() ? "" : "; ") + join_words(names, "and")
            + (names.size() == 1 ? " takes " : " take ") + join_words(n, "or");
  }
  return text;
}

// Positions read from a file: file:<path>, the path anything
inline bool init_mode_is_file(const std::string& init_mode){
  return init_mode.rfind("file:", 0) == 0;
}

// The initializers index the init_mode tokens without bounds checks; empty
// is refused by its own check
inline bool init_mode_shape_ok(const std::string& init_mode){
  if (init_mode.empty())
    return true;
  if (init_mode_is_file(init_mode))
    return init_mode.size() > 5;
  const std::vector<std::string> key = split_string(init_mode, "_");
  if (!init_mode_dirs_ok(key)) return false;
  const InitMode* mode = find_init_mode(key[0]);
  if (!mode) return key.size() == 2;
  return std::find(mode->n_tokens.begin(), mode->n_tokens.end(), key.size()) != mode->n_tokens.end();
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
    const bool remeshes = (p.has("refine") && p.get<bool>("refine")) || (p.has("coarsen") && p.get<bool>("coarsen"));
    return init_mode_has_edges(p) || remeshes;
  };
  s.require_if<double>("ds_max", 0.0, has_edges_or_remeshes,
                       "the initial state has edges, or refine/coarsen is on", "max edge length");
  s.require_if<double>("ds_min", 0.0, has_edges_or_remeshes,
                       "the initial state has edges, or refine/coarsen is on", "min edge length");
  s.require_if<double>("La",
                       [](const partrac::Params& p){ return init_mode_needs(p, "La"); },
                       "init_mode is " + init_modes_needing("La"),
                       "principal extent");
  s.require_if<double>("Lb",
                       [](const partrac::Params& p){ return init_mode_needs(p, "Lb"); },
                       "init_mode is " + init_modes_needing("Lb"),
                       "second extent");
  s.require_if<double>("ds_init",
                       [](const partrac::Params& p){ return init_mode_needs(p, "ds_init"); },
                       "init_mode is " + init_modes_needing("ds_init") + ", or points along one axis, as points_x",
                       "initial edge length");
  s.require_if<std::string>("init_weight",
                            [](const partrac::Params& p){ return init_mode_needs(p, "init_weight"); },
                            "init_mode is " + init_modes_needing("init_weight"),
                            "sampling weight");
  // uniform and strip step from end to end
  for (const std::string name : {"uniform", "strip"})
    s.check([name](const partrac::Params& p){
              const InitMode* mode = find_init_mode(p);
              return !mode || mode->name != name || p.get<Uint>("Nrw") >= 2;
            },
            "init_mode " + name + " needs Nrw of 2 or more");
  s.opt<double>("t0", 0.0, "start time");
  s.opt<double>("x0", 0.0, "initial position");
  s.opt<double>("y0", 0.0, "initial position");
  s.opt<double>("z0", 0.0, "initial position");
  s.opt<bool>("inject", false, "inject new particles");
  s.opt<bool>("clear_initial_edges", false, "drop the initial edges");
  s.check([](const partrac::Params& p){ return !p.get<std::string>("init_mode").empty(); },
          "init_mode is empty");
  s.check([](const partrac::Params& p){
            return init_mode_shape_ok(p.get<std::string>("init_mode"));
          },
          "init_mode has the wrong number of directions, each some of x, y and z as in uniform_x"
            " or sheet_xy: " + init_mode_directions() +
            "; file takes a path, as in file:positions.h5 (was from_file:positions.h5)");
}

// Initial state for init_mode
InitialState set_initial_state(std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen);

#endif
