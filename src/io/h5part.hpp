#ifndef __H5PART_HPP
#define __H5PART_HPP

#include <string>
#include <vector>
#include <hdf5.h>
#include "Error.hpp"
#include "h5direct.hpp"

inline void write_h5part(const std::string& path,
                         const std::vector<std::string>& ptheader,
                         const std::vector<double>& ptdata_,
                         const double t0 = 0.0) {

    const hsize_t num_cols = ptheader.size();
    if (num_cols == 0)
        partrac::fail("write_h5part: no columns to write to ", path);
    const hsize_t num_rows = ptdata_.size() / num_cols;

    // Errors ours, not HDF5's stack
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    const partrac::H5Id file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);
    if (!file.valid())
        partrac::fail("cannot create ", path);
    const partrac::H5Id grp(H5Gcreate2(file, "Step#0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
    if (!grp.valid())
        partrac::fail(path, ": cannot create the group Step#0");

    std::vector<double> col_data(num_rows);

    for (hsize_t ic = 0; ic < num_cols; ++ic) {
        #pragma omp parallel for
        for (hsize_t ir = 0; ir < num_rows; ++ir)
            col_data[ir] = ptdata_[ir * num_cols + ic];
        const partrac::H5Id space(H5Screate_simple(1, &num_rows, nullptr), H5Sclose);
        const partrac::H5Id dset(H5Dcreate2(grp, ptheader[ic].c_str(), H5T_NATIVE_DOUBLE, space,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
        if (!dset.valid() || H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, col_data.data()) < 0)
            partrac::fail(path, ": cannot write the column ", ptheader[ic]);
    }

    const partrac::H5Id s(H5Screate(H5S_SCALAR), H5Sclose);
    const partrac::H5Id a(H5Acreate2(grp, "TimeValue", H5T_NATIVE_DOUBLE, s,
                                     H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (!a.valid() || H5Awrite(a, H5T_NATIVE_DOUBLE, &t0) < 0)
        partrac::fail(path, ": cannot write the attribute TimeValue");
    if (H5Fflush(file, H5F_SCOPE_LOCAL) < 0)
        partrac::fail("cannot write ", path);
}

#endif
