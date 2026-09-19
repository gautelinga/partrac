#ifndef __H5PART_HPP
#define __H5PART_HPP

#include <string>
#include <vector>
#include <hdf5.h>

inline void write_h5part(const std::string& path,
                         const std::vector<std::string>& ptheader,
                         const std::vector<double>& ptdata_,
                         const double t0 = 0.0) {

    const hsize_t num_cols = ptheader.size();
    const hsize_t num_rows = ptdata_.size() / num_cols;

    auto file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    auto grp  = H5Gcreate2(file, "Step#0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    std::vector<double> col_data(num_rows);

    for (hsize_t ic = 0; ic < num_cols; ++ic) {
        #pragma omp parallel for
        for (hsize_t ir = 0; ir < num_rows; ++ir)
            col_data[ir] = ptdata_[ir * num_cols + ic];
        auto space = H5Screate_simple(1, &num_rows, nullptr);
        auto dset  = H5Dcreate2(grp, ptheader[ic].c_str(), H5T_NATIVE_DOUBLE, space,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, col_data.data());
        H5Dclose(dset);
        H5Sclose(space);
    }

    auto s = H5Screate(H5S_SCALAR);
    auto a = H5Acreate2(grp, "TimeValue", H5T_NATIVE_DOUBLE, s,
                         H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(a, H5T_NATIVE_DOUBLE, &t0);
    H5Aclose(a);
    H5Sclose(s);

    H5Gclose(grp);
    H5Fclose(file);
}

#endif
