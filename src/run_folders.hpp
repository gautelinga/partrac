#ifndef __RUN_FOLDERS_HPP
#define __RUN_FOLDERS_HPP

#include <string>
#include <sstream>
#include <iomanip>
#include "typedefs.hpp"
#include "Params.hpp"
#include "io.hpp"

// Output layout: <base>/<name>/<parameters>/0/{Positions,Checkpoints}/

inline std::string get_newfoldername(const std::string& rwfolder, const partrac::Params& prm){
  std::ostringstream ss_Dm, ss_dt, ss_Nrw, ss_seed;
  ss_Dm << std::scientific << std::setprecision(7) << prm.get<double>("Dm");
  ss_dt << std::scientific << std::setprecision(7) << prm.get<double>("dt");
  ss_Nrw << prm.get<Uint>("Nrw");
  ss_seed << prm.get<int>("seed");
  return rwfolder +
         "/Dm" + ss_Dm.str() +
         "_dt" + ss_dt.str() +
         "_Nrw" + ss_Nrw.str() +
         "_seed" + ss_seed.str() +
         prm.get<std::string>("tag") +
         "/";
}

// Layout options (bitmask)
enum RunFolderOpt : unsigned {
  DefaultLayout = 0,
  DryRun        = 1,   // work out the names, create nothing
  NoSubfolders  = 2,   // no Positions/Checkpoints
  NoRunIndex    = 4    // no trailing 0/
};

struct RunFolders {
  std::string run;
  std::string positions;
  std::string checkpoints;
};

inline RunFolders make_run_folders(const std::string& base,
                                   const std::string& name,
                                   partrac::Params& prm,
                                   const unsigned opts = DefaultLayout){
  const bool dry_run = opts & DryRun;
  const std::string appfolder = base + "/" + name + "/";
  if (!dry_run)
    create_folder(appfolder);

  std::string newfolder;
  if (prm.get<std::string>("restart_folder") != ""){
    newfolder = prm.get<std::string>("folder");
  }
  else {
    newfolder = get_newfoldername(appfolder, prm);
    if (!dry_run)
      create_folder(newfolder);
    if (!(opts & NoRunIndex))
      newfolder = newfolder + "0/";
  }

  RunFolders f{newfolder, newfolder + "Positions/", newfolder + "Checkpoints/"};
  if (!dry_run){
    create_folder(f.run);
    if (!(opts & NoSubfolders)){
      create_folder(f.positions);
      create_folder(f.checkpoints);
    }
  }
  prm.set<std::string>("folder", f.run);
  return f;
}

#endif
