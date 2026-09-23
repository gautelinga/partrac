// The raw HDF5 readers (io/h5direct.hpp). What they have to get right is that
// the bytes read with pread from several threads are the bytes HDF5 would have
// handed over, whatever the dataset's type; that a layout pread cannot follow
// (chunked, filtered) is read by the library instead; that a value the caller's
// type cannot hold stops the run rather than wrapping silently; and that a
// string attribute comes back as it was written, fixed or variable length.
#include <catch2/catch.hpp>

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include <hdf5.h>
#include <omp.h>

#include "Error.hpp"
#include "h5direct.hpp"
#include "case_dir.hpp"

namespace {

// A file of its own per case, removed when it goes out of scope
struct TempH5 {
  std::string name;
  explicit TempH5(const std::string& tag)
    : name(temp_path(tag + ".h5").string()) {}
  ~TempH5(){ std::remove(name.c_str()); }
};

// Several threads, so a read is split into slabs; the count is put back after
struct Threads {
  const int before;
  explicit Threads(const int n) : before(omp_get_max_threads()) { omp_set_num_threads(n); }
  ~Threads(){ omp_set_num_threads(before); }
};

// HDF5 prints its error stack to stderr, so a case that asks for something
// absent would otherwise bury the test output in it
struct QuietH5 {
  H5E_auto2_t printer = nullptr;
  void* data = nullptr;
  QuietH5(){
    H5Eget_auto2(H5E_DEFAULT, &printer, &data);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  }
  ~QuietH5(){ H5Eset_auto2(H5E_DEFAULT, printer, data); }
};

// gzip is optional in an HDF5 build, so the filtered case is only asked for
// where the filter is there; chunking alone already forces the fallback
bool have_deflate(){
  return H5Zfilter_avail(H5Z_FILTER_DEFLATE) > 0;
}

template<typename T>
void write_dataset(const hid_t file, const std::string& name, const hid_t type,
                   const std::vector<T>& v, const std::vector<hsize_t>& shape,
                   const std::vector<hsize_t>& chunk){
  const hid_t space = H5Screate_simple(static_cast<int>(shape.size()), shape.data(), nullptr);
  hid_t dcpl = H5P_DEFAULT;
  if (!chunk.empty()){
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, static_cast<int>(chunk.size()), chunk.data());
    if (have_deflate()) H5Pset_deflate(dcpl, 6);
  }
  const hid_t dset = H5Dcreate2(file, name.c_str(), type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
  H5Dclose(dset);
  if (dcpl != H5P_DEFAULT) H5Pclose(dcpl);
  H5Sclose(space);
}

void write_fixed_attribute(const hid_t loc, const std::string& name, const std::string& value,
                           const std::size_t width){
  const hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, width);
  H5Tset_strpad(type, H5T_STR_NULLPAD);
  const hid_t space = H5Screate(H5S_SCALAR);
  const hid_t attr = H5Acreate2(loc, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT);
  std::vector<char> buf(width, '\0');
  value.copy(buf.data(), std::min(width, value.size()));
  H5Awrite(attr, type, buf.data());
  H5Aclose(attr);
  H5Sclose(space);
  H5Tclose(type);
}

void write_variable_attribute(const hid_t loc, const std::string& name, const std::string& value){
  const hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, H5T_VARIABLE);
  const hid_t space = H5Screate(H5S_SCALAR);
  const hid_t attr = H5Acreate2(loc, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT);
  const char* p = value.c_str();
  H5Awrite(attr, type, &p);
  H5Aclose(attr);
  H5Sclose(space);
  H5Tclose(type);
}

// What HDF5 itself makes of a dataset, the reference every raw read is held to
template<typename T>
std::vector<T> read_with_hdf5(const std::string& file, const std::string& path, const hid_t type,
                              const std::size_t n){
  const hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  const hid_t d = H5Dopen2(f, path.c_str(), H5P_DEFAULT);
  std::vector<T> v(n);
  H5Dread(d, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
  H5Dclose(d);
  H5Fclose(f);
  return v;
}

const std::size_t n_ints = 1000;
const std::size_t n_rows = 200;
const std::size_t n_cols = 3;
const std::string fixed_text = "Lagrange triangle 2";
const std::string variable_text = "FiniteElement('Lagrange', tetrahedron, 2)";

// One file holding every layout and type the readers have to cope with
void write_fixture(const std::string& path){
  std::vector<std::int32_t> ints(n_ints);
  for (std::size_t i = 0; i < n_ints; ++i) ints[i] = static_cast<std::int32_t>(i) * 37 - 4000;
  std::vector<double> grid(n_rows * n_cols);
  for (std::size_t i = 0; i < grid.size(); ++i) grid[i] = 0.25 * static_cast<double>(i) - 1. / 3.;
  const std::vector<std::int32_t> wide = {0, 7, 100000, -5, 12};

  const hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  write_dataset(file, "ints", H5T_NATIVE_INT32, ints, {n_ints}, {});
  write_dataset(file, "grid", H5T_NATIVE_DOUBLE, grid, {n_rows, n_cols}, {});
  write_dataset(file, "chunked", H5T_NATIVE_DOUBLE, grid, {n_rows, n_cols}, {16, n_cols});
  write_dataset(file, "wide", H5T_NATIVE_INT32, wide, {wide.size()}, {});
  write_dataset(file, "wide_chunked", H5T_NATIVE_INT32, wide, {wide.size()}, {2});

  const hid_t group = H5Gcreate2(file, "meta", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  write_fixed_attribute(group, "element", fixed_text, fixed_text.size() + 4);
  write_variable_attribute(group, "signature", variable_text);
  H5Gclose(group);

  const hid_t dset = H5Dopen2(file, "grid", H5P_DEFAULT);
  write_variable_attribute(dset, "signature", variable_text);
  H5Dclose(dset);
  H5Fclose(file);
}

std::vector<double> reference_grid(){
  std::vector<double> grid(n_rows * n_cols);
  for (std::size_t i = 0; i < grid.size(); ++i) grid[i] = 0.25 * static_cast<double>(i) - 1. / 3.;
  return grid;
}

}  // namespace

TEST_CASE("h5direct reports how a dataset lies in the file", "[h5direct]") {
  TempH5 f("info");
  write_fixture(f.name);
  const partrac::H5Id file = partrac::h5_open_read(f.name);

  SECTION("a contiguous unfiltered dataset has a byte offset and can be read raw"){
    const partrac::H5DatasetInfo info = partrac::h5_dataset_info(file, "ints");
    REQUIRE(info.rank() == 1);
    REQUIRE(info.rows() == n_ints);
    REQUIRE(info.cols() == 1);
    REQUIRE(info.type_class == H5T_INTEGER);
    REQUIRE(info.type_size == 4);
    REQUIRE(info.type_signed);
    REQUIRE(info.contiguous);
    REQUIRE_FALSE(info.filtered);
    REQUIRE(info.offset != HADDR_UNDEF);
    REQUIRE(info.raw_readable());
  }
  SECTION("a 2D double dataset reports both extents"){
    const partrac::H5DatasetInfo info = partrac::h5_dataset_info(file, "grid");
    REQUIRE(info.rank() == 2);
    REQUIRE(info.rows() == n_rows);
    REQUIRE(info.cols() == n_cols);
    REQUIRE(info.type_class == H5T_FLOAT);
    REQUIRE(info.type_size == 8);
    REQUIRE(info.count() == n_rows * n_cols);
    REQUIRE(info.raw_readable());
  }
  SECTION("a chunked dataset has no offset to read from"){
    const partrac::H5DatasetInfo info = partrac::h5_dataset_info(file, "chunked");
    REQUIRE_FALSE(info.contiguous);
    REQUIRE(info.filtered == have_deflate());
    REQUIRE(info.offset == HADDR_UNDEF);
    REQUIRE_FALSE(info.raw_readable());
  }
  SECTION("a dataset that is not there is named in the failure"){
    const QuietH5 quiet;
    REQUIRE_THROWS_AS(partrac::h5_dataset_info(file, "absent"), partrac::Error);
  }
}

TEST_CASE("the raw read gives what H5Dread gives", "[h5direct]") {
  Threads threads(4);
  TempH5 f("raw");
  write_fixture(f.name);
  const partrac::H5Id file = partrac::h5_open_read(f.name);

  SECTION("a 1D dataset read as its own type"){
    const std::vector<std::int32_t> want =
      read_with_hdf5<std::int32_t>(f.name, "ints", H5T_NATIVE_INT32, n_ints);
    std::vector<std::int32_t> got;
    partrac::h5_read(file, "ints", got);
    REQUIRE(got == want);
  }
  SECTION("a 1D dataset converted to another type"){
    const std::vector<double> want = read_with_hdf5<double>(f.name, "ints", H5T_NATIVE_DOUBLE, n_ints);
    std::vector<double> got;
    partrac::h5_read(file, "ints", got);
    REQUIRE(got == want);
    std::vector<std::int64_t> wide;
    partrac::h5_read(file, "ints", wide);
    REQUIRE(wide.size() == n_ints);
    REQUIRE(static_cast<double>(wide[n_ints - 1]) == want[n_ints - 1]);
  }
  SECTION("a 2D dataset, every column"){
    const std::vector<double> want = reference_grid();
    std::vector<double> got;
    partrac::h5_read(file, "grid", got);
    REQUIRE(got == want);
  }
  SECTION("a 2D dataset, the leading columns of each row"){
    const std::vector<double> want = reference_grid();
    std::vector<double> got;
    partrac::h5_read(file, "grid", got, 2);
    REQUIRE(got.size() == 2 * n_rows);
    bool same = true;
    for (std::size_t i = 0; i < n_rows; ++i)
      for (std::size_t j = 0; j < 2; ++j)
        same = same && got[2 * i + j] == want[n_cols * i + j];
    REQUIRE(same);
  }
  SECTION("more columns than the file holds is a failure, not a read past the row"){
    std::vector<double> got;
    REQUIRE_THROWS_AS(partrac::h5_read(file, "grid", got, 4), partrac::Error);
  }
}

TEST_CASE("a chunked dataset falls back to the library and still reads", "[h5direct]") {
  Threads threads(4);
  TempH5 f("chunked");
  write_fixture(f.name);
  const partrac::H5Id file = partrac::h5_open_read(f.name);
  REQUIRE_FALSE(partrac::h5_dataset_info(file, "chunked").raw_readable());

  const std::vector<double> want = reference_grid();
  std::vector<double> got;
  partrac::h5_read(file, "chunked", got);
  REQUIRE(got == want);

  // the same dataset through the raw path, so the fallback is not a different answer
  std::vector<double> raw;
  partrac::h5_read(file, "grid", raw);
  REQUIRE(got == raw);

  std::vector<double> two;
  partrac::h5_read(file, "chunked", two, 2);
  REQUIRE(two.size() == 2 * n_rows);
  REQUIRE(two[0] == want[0]);
  REQUIRE(two[1] == want[1]);
  REQUIRE(two[2] == want[n_cols]);
}

TEST_CASE("a value the caller's type cannot hold stops the run", "[h5direct]") {
  Threads threads(4);
  TempH5 f("range");
  write_fixture(f.name);
  const partrac::H5Id file = partrac::h5_open_read(f.name);

  SECTION("an overflow on the raw path"){
    std::vector<std::uint16_t> got;
    REQUIRE_THROWS_AS(partrac::h5_read(file, "wide", got), partrac::Error);
  }
  SECTION("a negative value into an unsigned type"){
    std::vector<std::uint32_t> got;
    REQUIRE_THROWS_AS(partrac::h5_read(file, "wide", got), partrac::Error);
  }
  SECTION("the fallback path checks the same"){
    std::vector<std::uint32_t> got;
    REQUIRE_THROWS_AS(partrac::h5_read(file, "wide_chunked", got), partrac::Error);
  }
  SECTION("a type that holds every value does not"){
    std::vector<std::int64_t> got;
    REQUIRE_NOTHROW(partrac::h5_read(file, "wide", got));
    REQUIRE(got == std::vector<std::int64_t>({0, 7, 100000, -5, 12}));
  }
}

TEST_CASE("a string attribute comes back as it was written", "[h5direct]") {
  TempH5 f("attr");
  write_fixture(f.name);
  const partrac::H5Id file = partrac::h5_open_read(f.name);

  SECTION("fixed length, padded in the file"){
    REQUIRE(partrac::h5_read_string_attribute(file, "meta", "element") == fixed_text);
  }
  SECTION("variable length, on a group and on a dataset"){
    REQUIRE(partrac::h5_read_string_attribute(file, "meta", "signature") == variable_text);
    REQUIRE(partrac::h5_read_string_attribute(file, "grid", "signature") == variable_text);
  }
  SECTION("an attribute that is not there is a failure"){
    const QuietH5 quiet;
    REQUIRE_THROWS_AS(partrac::h5_read_string_attribute(file, "meta", "absent"), partrac::Error);
  }
}
