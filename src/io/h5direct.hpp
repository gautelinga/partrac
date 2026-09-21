#ifndef __H5DIRECT_HPP
#define __H5DIRECT_HPP

// Reading HDF5 datasets with the library on the metadata only: shape, type and
// byte offset through the C API, then the bytes of a contiguous unfiltered
// dataset with pread from OpenMP threads, a slab each, and one parallel pass
// into the caller's type. H5Dread reads anything else.

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

#include <hdf5.h>
#include <omp.h>

#include "Error.hpp"

namespace partrac {

// Closes its hid_t; a fail must not leak one
class H5Id {
public:
  H5Id(const hid_t id, herr_t (* const closer)(hid_t)) : id_(id), closer_(closer) {}
  ~H5Id(){ if (id_ >= 0) closer_(id_); }
  H5Id(const H5Id&) = delete;
  H5Id& operator=(const H5Id&) = delete;
  operator hid_t() const { return id_; }
  bool valid() const { return id_ >= 0; }
private:
  hid_t id_;
  herr_t (*closer_)(hid_t);
};

// The file an object belongs to, for messages
inline std::string h5_file_name(const hid_t loc_id){
  const ssize_t n = H5Fget_name(loc_id, nullptr, 0);
  if (n <= 0) return std::string();
  std::string s(static_cast<std::size_t>(n), '\0');
  H5Fget_name(loc_id, &s[0], static_cast<std::size_t>(n) + 1);
  return s;
}

inline H5Id h5_open_read(const std::string& path){
  const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) fail("h5: cannot open ", path);
  return H5Id(file, H5Fclose);
}

// What a dataset is on disk
struct H5DatasetInfo {
  std::vector<hsize_t> shape;
  H5T_class_t type_class = H5T_NO_CLASS;
  std::size_t type_size = 0;
  bool type_signed = false;
  bool type_native = false;   // the stored bytes are the machine's own
  bool contiguous = false;
  bool filtered = false;
  haddr_t offset = HADDR_UNDEF;

  std::size_t rank() const { return shape.size(); }
  std::size_t rows() const { return shape.empty() ? 0 : static_cast<std::size_t>(shape[0]); }
  std::size_t cols() const { return shape.size() < 2 ? 1 : static_cast<std::size_t>(shape[1]); }
  std::size_t count() const { return rows() * cols(); }
  // the bytes can be taken as they lie
  bool raw_readable() const { return contiguous && !filtered && type_native && offset != HADDR_UNDEF; }
};

// The machine's own type of a class, size and sign, or -1 for none
inline hid_t h5_native_of(const H5T_class_t cls, const std::size_t size, const bool sgn){
  if (cls == H5T_FLOAT){
    if (size == 4) return H5T_NATIVE_FLOAT;
    if (size == 8) return H5T_NATIVE_DOUBLE;
    return -1;
  }
  if (cls == H5T_INTEGER){
    switch (size){
      case 1: return sgn ? H5T_NATIVE_INT8 : H5T_NATIVE_UINT8;
      case 2: return sgn ? H5T_NATIVE_INT16 : H5T_NATIVE_UINT16;
      case 4: return sgn ? H5T_NATIVE_INT32 : H5T_NATIVE_UINT32;
      case 8: return sgn ? H5T_NATIVE_INT64 : H5T_NATIVE_UINT64;
      default: return -1;
    }
  }
  return -1;
}

// The type a std::vector<T> needs from HDF5
template<typename T>
inline hid_t h5_native_type(){
  static_assert(std::is_arithmetic<T>::value, "h5: only arithmetic element types");
  if constexpr (std::is_same_v<T, double>) return H5T_NATIVE_DOUBLE;
  else if constexpr (std::is_same_v<T, float>) return H5T_NATIVE_FLOAT;
  else {
    constexpr bool sgn = std::is_signed_v<T>;
    if constexpr (sizeof(T) == 1) return sgn ? H5T_NATIVE_INT8 : H5T_NATIVE_UINT8;
    else if constexpr (sizeof(T) == 2) return sgn ? H5T_NATIVE_INT16 : H5T_NATIVE_UINT16;
    else if constexpr (sizeof(T) == 4) return sgn ? H5T_NATIVE_INT32 : H5T_NATIVE_UINT32;
    else return sgn ? H5T_NATIVE_INT64 : H5T_NATIVE_UINT64;
  }
}

inline H5DatasetInfo h5_dataset_info(const hid_t loc_id, const std::string& path){
  const H5Id dset(H5Dopen2(loc_id, path.c_str(), H5P_DEFAULT), H5Dclose);
  if (!dset.valid()) fail("h5: no dataset '", path, "' in ", h5_file_name(loc_id));
  H5DatasetInfo info;
  {
    const H5Id space(H5Dget_space(dset), H5Sclose);
    const int ndims = H5Sget_simple_extent_ndims(space);
    if (ndims < 0) fail("h5: no shape for '", path, "' in ", h5_file_name(loc_id));
    info.shape.assign(static_cast<std::size_t>(ndims), 0);
    if (ndims > 0) H5Sget_simple_extent_dims(space, info.shape.data(), nullptr);
  }
  {
    const H5Id type(H5Dget_type(dset), H5Tclose);
    info.type_class = H5Tget_class(type);
    info.type_size = H5Tget_size(type);
    if (info.type_class == H5T_INTEGER) info.type_signed = H5Tget_sign(type) == H5T_SGN_2;
    const hid_t native = h5_native_of(info.type_class, info.type_size, info.type_signed);
    info.type_native = native >= 0 && H5Tequal(type, native) > 0;
  }
  {
    const H5Id dcpl(H5Dget_create_plist(dset), H5Pclose);
    info.contiguous = H5Pget_layout(dcpl) == H5D_CONTIGUOUS;
    info.filtered = H5Pget_nfilters(dcpl) > 0;
  }
  // defined for a contiguous unfiltered dataset whose space is allocated
  if (info.contiguous && !info.filtered) info.offset = H5Dget_offset(dset);
  return info;
}

// The descriptor pread reads from. The sec2 driver's own where the file has one
// -- no second open, and it is the file the id refers to whatever the path says
// now -- else a read-only open of the name the file reports, which is what a
// parallel HDF5 build opened with another driver leaves us.
class H5RawFile {
public:
  explicit H5RawFile(const hid_t loc_id) : file_(H5Iget_file_id(loc_id), H5Fclose){
    if (!file_.valid()) return;
    unsigned intent = 0;
    // the bytes must be in the file, not in the library's cache
    if (H5Fget_intent(file_, &intent) >= 0 && intent != H5F_ACC_RDONLY)
      H5Fflush(file_, H5F_SCOPE_LOCAL);
    const H5Id fapl(H5Fget_access_plist(file_), H5Pclose);
    if (fapl.valid() && H5Pget_driver(fapl) == H5FD_SEC2){
      int* handle = nullptr;
      if (H5Fget_vfd_handle(file_, H5P_DEFAULT, reinterpret_cast<void**>(&handle)) >= 0 && handle != nullptr)
        fd_ = *handle;
    }
    if (fd_ < 0){
      const std::string name = h5_file_name(file_);
      if (!name.empty()){
        fd_ = ::open(name.c_str(), O_RDONLY);
        own_ = fd_ >= 0;
      }
    }
  }
  ~H5RawFile(){ if (own_ && fd_ >= 0) ::close(fd_); }
  H5RawFile(const H5RawFile&) = delete;
  H5RawFile& operator=(const H5RawFile&) = delete;
  int fd() const { return fd_; }
private:
  H5Id file_;
  int fd_ = -1;
  bool own_ = false;
};

// One contiguous slab a thread, whole elements each
inline void h5_pread(const int fd, const haddr_t offset, char* const dst, const std::size_t nbytes,
                     const std::size_t esize, const std::string& what){
  if (nbytes == 0) return;
  int bad_errno = 0;
  int truncated = 0;
  #pragma omp parallel reduction(max:bad_errno) reduction(max:truncated)
  {
    const std::size_t nt = static_cast<std::size_t>(omp_get_num_threads());
    const std::size_t it = static_cast<std::size_t>(omp_get_thread_num());
    const std::size_t nelem = nbytes / esize;
    const std::size_t per = (nelem + nt - 1) / nt;
    const std::size_t begin = std::min(nbytes, per * esize * it);
    const std::size_t end = std::min(nbytes, begin + per * esize);
    std::size_t done = begin;
    while (done < end){
      const ssize_t n = ::pread(fd, dst + done, end - done, static_cast<off_t>(offset + done));
      if (n < 0){
        if (errno == EINTR) continue;
        bad_errno = errno;
        break;
      }
      if (n == 0){ truncated = 1; break; }
      done += static_cast<std::size_t>(n);
    }
  }
  if (bad_errno != 0) fail("h5: reading '", what, "' failed: ", std::strerror(bad_errno));
  if (truncated != 0) fail("h5: the file ends inside '", what, "'");
}

// v survives the conversion to T
template<typename T, typename S>
inline bool h5_in_range(const S v){
  static_assert(std::is_arithmetic<T>::value && std::is_arithmetic<S>::value, "h5: arithmetic only");
  if constexpr (std::is_floating_point_v<T>){
    return true;
  }
  else if constexpr (std::is_floating_point_v<S>){
    // the integer limits are exact in S, so max + 1 is the first value that does not fit
    return v >= static_cast<S>(std::numeric_limits<T>::lowest())
        && v < static_cast<S>(std::numeric_limits<T>::max()) + S(1);
  }
  else {
    if constexpr (std::is_signed_v<S>){
      if (v < S(0)){
        if constexpr (std::is_unsigned_v<T>) return false;
        else return static_cast<std::intmax_t>(v) >= static_cast<std::intmax_t>(std::numeric_limits<T>::lowest());
      }
    }
    return static_cast<std::uintmax_t>(v) <= static_cast<std::uintmax_t>(std::numeric_limits<T>::max());
  }
}

// The file's values into the caller's type, keeping the first cols of each row
template<typename T, typename S>
void h5_convert(const S* const src, std::vector<T>& out, const std::size_t rows,
                const std::size_t cols_file, const std::size_t cols, const std::string& path){
  constexpr std::size_t none = std::numeric_limits<std::size_t>::max();
  std::size_t bad = none;
  S bad_value = S(0);
  // a throw inside the region would terminate, so the first offender is carried out
  #pragma omp parallel for
  for (std::size_t i = 0; i < rows; ++i){
    for (std::size_t j = 0; j < cols; ++j){
      const S v = src[i * cols_file + j];
      if (h5_in_range<T>(v)){
        out[i * cols + j] = static_cast<T>(v);
      }
      else {
        #pragma omp critical
        {
          if (i * cols_file + j < bad){ bad = i * cols_file + j; bad_value = v; }
        }
      }
    }
  }
  if (bad != none)
    fail("h5: '", path, "' holds ", +bad_value, " at element ", bad, ", which the reader's type cannot hold");
}

// Raw bytes where the dataset allows it, H5Dread of the file's own type otherwise
template<typename S>
void h5_fill(const hid_t loc_id, const std::string& path, const H5DatasetInfo& info,
             S* const dst, const std::size_t n){
  if (info.raw_readable()){
    const H5RawFile raw(loc_id);
    if (raw.fd() >= 0){
      h5_pread(raw.fd(), info.offset, reinterpret_cast<char*>(dst), n * sizeof(S), sizeof(S), path);
      return;
    }
  }
  const H5Id dset(H5Dopen2(loc_id, path.c_str(), H5P_DEFAULT), H5Dclose);
  if (!dset.valid()) fail("h5: no dataset '", path, "' in ", h5_file_name(loc_id));
  if (H5Dread(dset, h5_native_type<S>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, dst) < 0)
    fail("h5: cannot read '", path, "' in ", h5_file_name(loc_id));
}

template<typename T, typename S>
void h5_read_typed(const hid_t loc_id, const std::string& path, const H5DatasetInfo& info,
                   const std::size_t rows, const std::size_t cols_file, const std::size_t cols,
                   std::vector<T>& out){
  const std::size_t n = rows * cols_file;
  if constexpr (std::is_same_v<T, S>){
    // no conversion and no dropped column: straight into the caller's vector
    if (cols == cols_file){
      h5_fill<S>(loc_id, path, info, out.data(), n);
      return;
    }
  }
  std::vector<S> buf(n);
  h5_fill<S>(loc_id, path, info, buf.data(), n);
  h5_convert<T, S>(buf.data(), out, rows, cols_file, cols, path);
}

// A type this reader does not know: HDF5 converts it
template<typename T>
void h5_read_converted(const hid_t loc_id, const std::string& path, const std::size_t rows,
                       const std::size_t cols_file, const std::size_t cols, std::vector<T>& out){
  const H5Id dset(H5Dopen2(loc_id, path.c_str(), H5P_DEFAULT), H5Dclose);
  if (!dset.valid()) fail("h5: no dataset '", path, "' in ", h5_file_name(loc_id));
  herr_t err = 0;
  if (cols == cols_file){
    err = H5Dread(dset, h5_native_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
  }
  else {
    const hsize_t count[2] = {rows, cols};
    const hsize_t offset[2] = {0, 0};
    const H5Id space(H5Dget_space(dset), H5Sclose);
    H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    const H5Id memspace(H5Screate_simple(2, count, nullptr), H5Sclose);
    err = H5Dread(dset, h5_native_type<T>(), memspace, space, H5P_DEFAULT, out.data());
  }
  if (err < 0) fail("h5: cannot read '", path, "' in ", h5_file_name(loc_id));
}

// A 1D or 2D dataset into out. cols, if given, keeps that many leading columns
// of each row and drops the rest, as a padded vector field has.
template<typename T>
void h5_read(const hid_t loc_id, const std::string& path, std::vector<T>& out,
             const std::size_t ncols = 0){
  const H5DatasetInfo info = h5_dataset_info(loc_id, path);
  if (info.rank() < 1 || info.rank() > 2)
    fail("h5: '", path, "' in ", h5_file_name(loc_id), " has rank ", info.rank(), ", expected 1 or 2");
  const std::size_t rows = info.rows();
  const std::size_t cols_file = info.cols();
  const std::size_t cols = ncols == 0 ? cols_file : ncols;
  if (cols > cols_file)
    fail("h5: '", path, "' in ", h5_file_name(loc_id), " has ", cols_file, " columns, ", cols, " asked for");
  out.resize(rows * cols);
  if (out.empty()) return;

  if (info.type_class == H5T_FLOAT && info.type_size == 8)
    h5_read_typed<T, double>(loc_id, path, info, rows, cols_file, cols, out);
  else if (info.type_class == H5T_FLOAT && info.type_size == 4)
    h5_read_typed<T, float>(loc_id, path, info, rows, cols_file, cols, out);
  else if (info.type_class == H5T_INTEGER && info.type_size == 8 && info.type_signed)
    h5_read_typed<T, std::int64_t>(loc_id, path, info, rows, cols_file, cols, out);
  else if (info.type_class == H5T_INTEGER && info.type_size == 8)
    h5_read_typed<T, std::uint64_t>(loc_id, path, info, rows, cols_file, cols, out);
  else if (info.type_class == H5T_INTEGER && info.type_size == 4 && info.type_signed)
    h5_read_typed<T, std::int32_t>(loc_id, path, info, rows, cols_file, cols, out);
  else if (info.type_class == H5T_INTEGER && info.type_size == 4)
    h5_read_typed<T, std::uint32_t>(loc_id, path, info, rows, cols_file, cols, out);
  else if (info.type_class == H5T_INTEGER && info.type_size == 2 && info.type_signed)
    h5_read_typed<T, std::int16_t>(loc_id, path, info, rows, cols_file, cols, out);
  else if (info.type_class == H5T_INTEGER && info.type_size == 2)
    h5_read_typed<T, std::uint16_t>(loc_id, path, info, rows, cols_file, cols, out);
  else if (info.type_class == H5T_INTEGER && info.type_size == 1 && info.type_signed)
    h5_read_typed<T, std::int8_t>(loc_id, path, info, rows, cols_file, cols, out);
  else if (info.type_class == H5T_INTEGER && info.type_size == 1)
    h5_read_typed<T, std::uint8_t>(loc_id, path, info, rows, cols_file, cols, out);
  else
    h5_read_converted<T>(loc_id, path, rows, cols_file, cols, out);
}

// A string attribute of a group or a dataset, fixed or variable length. An empty
// obj_path reads the attribute of loc_id itself.
inline std::string h5_read_string_attribute(const hid_t loc_id, const std::string& obj_path,
                                            const std::string& name){
  struct ObjId {
    hid_t id; bool own;
    ~ObjId(){ if (own && id >= 0) H5Oclose(id); }
  } obj{obj_path.empty() ? loc_id : H5Oopen(loc_id, obj_path.c_str(), H5P_DEFAULT), !obj_path.empty()};
  if (obj.id < 0) fail("h5: no object '", obj_path, "' in ", h5_file_name(loc_id));

  const H5Id attr(H5Aopen(obj.id, name.c_str(), H5P_DEFAULT), H5Aclose);
  if (!attr.valid()) fail("h5: no attribute '", name, "' on '", obj_path, "' in ", h5_file_name(loc_id));
  const H5Id type(H5Aget_type(attr), H5Tclose);
  if (H5Tget_class(type) != H5T_STRING) fail("h5: attribute '", name, "' is not a string");
  const H5Id space(H5Aget_space(attr), H5Sclose);
  if (H5Sget_simple_extent_npoints(space) != 1)
    fail("h5: attribute '", name, "' is not a single string");

  std::string value;
  if (H5Tis_variable_str(type) > 0){
    const H5Id mem(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(mem, H5T_VARIABLE);
    H5Tset_cset(mem, H5Tget_cset(type));
    char* p = nullptr;
    if (H5Aread(attr, mem, &p) < 0) fail("h5: cannot read attribute '", name, "'");
    if (p != nullptr) value.assign(p);
    H5Dvlen_reclaim(mem, space, H5P_DEFAULT, &p);
  }
  else {
    const std::size_t size = H5Tget_size(type);
    const H5Id mem(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(mem, size);
    H5Tset_strpad(mem, H5Tget_strpad(type));
    H5Tset_cset(mem, H5Tget_cset(type));
    std::vector<char> buf(size + 1, '\0');
    if (H5Aread(attr, mem, buf.data()) < 0) fail("h5: cannot read attribute '", name, "'");
    std::size_t len = std::strlen(buf.data());   // the buffer carries one byte more, always zero
    if (H5Tget_strpad(type) == H5T_STR_SPACEPAD)
      while (len > 0 && buf[len - 1] == ' ') --len;
    value.assign(buf.data(), len);
  }
  return value;
}

}  // namespace partrac

#endif
