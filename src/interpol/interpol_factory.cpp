#include <iostream>
#include "Error.hpp"
#include "interpol_factory.hpp"
#include "AnalyticInterpol.hpp"
#include "StructuredInterpol.hpp"
#include "TetInterpol.hpp"
#include "TriangleInterpol.hpp"
#include "TriangleFreqInterpol.hpp"
#include "TetFreqInterpol.hpp"
#include "XDMFTriangleInterpol.hpp"
#include "XDMFTetInterpol.hpp"
#ifdef USE_DOLFIN
#include "DolfTriangleInterpol.hpp"
#include "DolfTetInterpol.hpp"
#endif

// Each interpolator built here is also in with_concrete and in PARTRAC_STEP_INTERPOLATORS
void set_interpolate_mode(std::shared_ptr<Interpol>& intp, const std::string& mode, const std::string& infilename){
  if (mode == "analytic"){
    std::cout << "AnalyticInterpol initiated." << std::endl;
    intp = std::make_shared<AnalyticInterpol>(infilename);
  }
  else if (mode == "tet"){
    intp = std::make_shared<TetInterpol>(infilename);
  }
  else if (mode == "triangle"){
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
  else if (mode == "unstructured" || mode == "fenics" || mode == "xdmf"){
    if (mode == "xdmf")
      partrac::fail("XDMF format is not implemented yet.");
    if (mode == "unstructured")
      partrac::fail("mode unstructured does not name a loader; a mesh is read by 'fenics', "
                    "'tet', 'triangle', 'trianglefreq', 'tetfreq', 'xdmftet' or 'xdmftriangle'.");
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
