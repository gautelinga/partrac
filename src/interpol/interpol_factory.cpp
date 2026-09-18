#include <iostream>
#include "interpol_factory.hpp"
#include "AnalyticInterpol.hpp"
#include "StructuredInterpol.hpp"
#ifdef USE_DOLFIN
#include "DolfTriangleInterpol.hpp"
#include "DolfTetInterpol.hpp"
#include "TetInterpol.hpp"
#include "TriangleInterpol.hpp"
#include "TriangleFreqInterpol.hpp"
#include "XDMFTriangleInterpol.hpp"
#include "XDMFTetInterpol.hpp"
#endif

// Each interpolator built here is also in with_concrete and in PARTRAC_STEP_INTERPOLATORS
void set_interpolate_mode(std::shared_ptr<Interpol>& intp, const std::string& mode, const std::string& infilename){
  if (mode == "analytic"){
    std::cout << "AnalyticInterpol initiated." << std::endl;
    intp = std::make_shared<AnalyticInterpol>(infilename);
  }
  else if (mode == "unstructured" || mode == "fenics" || mode == "xdmf" ||
           mode == "tet" || mode == "triangle" || mode == "trianglefreq" || mode == "xdmftriangle" || mode == "xdmftet"){
#ifdef USE_DOLFIN
    if (mode == "tet"){
      intp = std::make_shared<TetInterpol>(infilename);
    }
    else if (mode == "triangle"){
      intp = std::make_shared<TriangleInterpol>(infilename);
    }
    else if (mode == "trianglefreq"){
      intp = std::make_shared<TriangleFreqInterpol>(infilename);
    }
    else if (mode == "xdmftet"){
      intp = std::make_shared<XDMFTetInterpol>(infilename);
    }
    else if (mode == "xdmftriangle"){
      intp = std::make_shared<XDMFTriangleInterpol>(infilename);
    }
    else if (mode == "fenics"){
      // The cell type is the mesh's; a file without one fails in the tet loader's checks
      if (dolfin_mesh_dim(infilename) == 2)
        intp = std::make_shared<DolfTriangleInterpol>(infilename);
      else
        intp = std::make_shared<DolfTetInterpol>(infilename);
    }
    else if (mode == "xdmf"){
      std::cout << "XDMF format is not implemented yet." << std::endl;
      exit(1);
    }
    else {
      std::cout << "Mode should be 'fenics', 'tet' or 'triangle'." << std::endl;
      exit(1);
    }
#else
    std::cout << "You have to compile with PARTRAC_ENABLE_DOLFIN=ON." << std::endl;
    exit(1);
#endif
  }
  else if (mode == "structured" || mode == "lbm" || mode == "felbm"){
    if (partrac::peek_file(infilename, "interpolation") == "constant")
      intp = std::make_shared<StructuredConstInterpol>(infilename);
    else
      intp = std::make_shared<StructuredInterpol>(infilename);
  }
  else {
    std::cout << "Mode not supported." << std::endl;
    exit(1);
  }
}
