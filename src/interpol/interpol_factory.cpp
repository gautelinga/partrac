#include <iostream>
#include "Error.hpp"
#include "interpol_factory.hpp"
#include "interpolators.hpp"

// Each interpolator built here is a row of the table in src/CMakeLists.txt
void set_interpolate_mode(std::shared_ptr<Interpol>& intp, const std::string& mode, const std::string& path){
  // A bare file name is in the working directory
  const std::string infilename = path.find('/') == std::string::npos ? "./" + path : path;
  if (mode == "analytic"){
    std::cout << "AnalyticInterpol initiated." << std::endl;
    intp = std::make_shared<AnalyticInterpol>(infilename);
  }
  else if (mode == "tet"){
    if (partrac::peek_bool(infilename, "divfree"))
      intp = std::make_shared<SplitTetInterpol>(infilename);
    else
      intp = std::make_shared<TetInterpol>(infilename);
  }
  else if (mode == "triangle"){
    if (partrac::peek_bool(infilename, "divfree"))
      intp = std::make_shared<SplitTriangleInterpol>(infilename);
    else
      intp = std::make_shared<TriangleInterpol>(infilename);
  }
  else if (mode == "trianglefreq"){
    intp = std::make_shared<TriangleFreqInterpol>(infilename);
  }
  else if (mode == "tetfreq"){
    intp = std::make_shared<TetFreqInterpol>(infilename);
  }
  else if (mode == "xdmftet"){
    intp = std::make_shared<XDMFTetInterpol>(infilename);
  }
  else if (mode == "xdmftriangle"){
    intp = std::make_shared<XDMFTriangleInterpol>(infilename);
  }
  else if (mode == "openfoam"){
    if (!openfoam_load::available())
      partrac::fail("mode ", mode, " needs a build with PARTRAC_ENABLE_OPENFOAM=ON.");
    // The cell type is the case's: triangles between an empty pair
    if (openfoam_load::has_empty_patches(infilename.substr(0, infilename.find_last_of('/'))))
      intp = std::make_shared<OpenFoamTriangleInterpol>(infilename);
    else
      intp = std::make_shared<OpenFoamTetInterpol>(infilename);
  }
  else if (mode == "fenics"){
#ifdef USE_DOLFIN
    // The cell type is the mesh's; a file without one fails in the tet loader's checks
    if (dolfin_mesh_dim(infilename) == 2)
      intp = std::make_shared<DolfTriangleInterpol>(infilename);
    else
      intp = std::make_shared<DolfTetInterpol>(infilename);
#else
    partrac::fail("mode ", mode, " needs a build with PARTRAC_ENABLE_DOLFIN=ON.");
#endif
  }
  else if (mode == "structured" || mode == "lbm" || mode == "felbm"){
    if (partrac::peek_file(infilename, "interpolation") == "constant")
      intp = std::make_shared<StructuredConstInterpol>(infilename);
    else
      intp = std::make_shared<StructuredInterpol>(infilename);
  }
  else {
    partrac::fail("mode not supported.");
  }
}
