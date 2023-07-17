#include <iostream>
#include <vector>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <set>
#include <iterator>
#include "H5Cpp.h"
//#include "hdf5.h"
#include <ctime>

#include "MPIwrap.hpp"
#include "StructuredInterpol.hpp"
#include "TriangleInterpol.hpp"


int main(int argc, char* argv[])
{
    MPIwrap mpi(argc, argv);

    if (mpi.rank() == 0)
    {
        std::cout << "======================================================================\n"
                  << "   Initialized mpitest with " << mpi.size() << " processes. \n"
                  << "======================================================================" << std::endl;
    }
    mpi.barrier();
    
    std::cout << "This is process " << mpi.rank() << " out of " << mpi.size() << "." << std::endl;
    mpi.barrier();

    // Input parameters
    if (argc < 2 && mpi.rank() == 0) {
        std::cout << "Please specify an input file." << std::endl;
        return 0;
    }

    std::string infilename = std::string(argv[1]);

    {
        TriangleInterpol intp(infilename);
    }
    mpi.finalize();

    return EXIT_SUCCESS;
}
