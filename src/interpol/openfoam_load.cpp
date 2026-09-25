// The OpenFOAM reader: the only translation unit that includes an OpenFOAM
// header, built as the plugin libpartrac_openfoam.so, whose C entry points
// fill the plain arrays of openfoam_load.hpp. The polyMesh is read, copied out
// and freed before the case is returned; a field is read from its file alone,
// each patch's condition asked whether it fixes the value on a one-cell probe
// mesh. OpenFOAM's errors are thrown, caught and returned as text.

#include "openfoam_load.hpp"

#include <algorithm>
#include <cmath>
#include <dlfcn.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>
#include <memory>
#include <numeric>
#include <sstream>

#include "Time.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "cyclicPolyPatch.H"
#include "emptyPolyPatch.H"
#include "wallPolyPatch.H"
#include "polyBoundaryMeshEntries.H"
#include "volFields.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "OStringStream.H"

namespace fs = std::filesystem;
using namespace openfoam_load;

namespace {

// A refusal of this reader, as a message
struct Refused {
  std::string what;
};

template<typename... Args>
[[noreturn]] void refuse(const Args&... args){
  std::ostringstream s;
  (s << ... << args);
  throw Refused{s.str()};
}

const char* const convert_cmd =
  "foamDictionary -entry writeFormat -set ascii system/controlDict && foamFormatConvert -case .";

// Patch types taken
const std::vector<std::string> accepted = {"patch", "wall", "mappedWall", "symmetry",
                                           "symmetryPlane", "cyclic", "empty"};

bool contains(const std::vector<std::string>& v, const std::string& s){
  return std::find(v.begin(), v.end(), s) != v.end();
}

std::string joined(const std::vector<std::string>& v){
  std::string out;
  for (const auto& s : v) out += (out.empty() ? "" : ", ") + s;
  return out;
}

// A file or its gzipped twin
bool exists_either(const fs::path& p){
  return fs::exists(p) || fs::exists(p.string() + ".gz");
}

bool is_number(const std::string& s){
  if (s.empty()) return false;
  char* end = nullptr;
  std::strtod(s.c_str(), &end);
  return end && *end == '\0';
}

// Quiet OpenFOAM's banners, and make its errors exceptions
void quiet(){
  Foam::messageStream::level = 0;
  Foam::FatalError.throwExceptions();
  Foam::FatalIOError.throwExceptions();
}

// The header of a file; binary with 64-bit labels is refused, since this
// OpenFOAM reads it as 32-bit
void check_header(const fs::path& p){
  // IFstream reads name.gz for name
  std::string name = p.string();
  if (name.size() > 3 && name.compare(name.size() - 3, 3, ".gz") == 0) name.resize(name.size() - 3);
  Foam::IFstream is(name);
  if (!is.good()) refuse(p.string(), ": cannot open");
  Foam::token t(is);
  while (t.good() && !(t.isWord() && t.wordToken() == "FoamFile")) t = Foam::token(is);
  if (!t.good()) return;
  const Foam::dictionary head(is);
  const Foam::word format = head.lookupOrDefault<Foam::word>("format", "ascii");
  const Foam::string arch = head.lookupOrDefault<Foam::string>("arch", "");
  if (format == "binary" && arch.find("label=64") != std::string::npos)
    refuse(p.string(), " is binary with 64-bit labels, which this OpenFOAM reads as 32-bit; "
           "convert the case to ascii first: ", convert_cmd);
}

// The case's root and name, as Time takes them
std::unique_ptr<Foam::Time> open_time(const std::string& dir){
  const Foam::fileName c = Foam::fileName(fs::absolute(dir).lexically_normal().string());
  return std::make_unique<Foam::Time>(Foam::Time::controlDictName, c.path(), c.name(), false);
}

// A case the reader takes: system/, a reconstructed mesh, no moving mesh
void check_case(const std::string& dir){
  const fs::path d(dir);
  if (!fs::is_directory(d)) refuse(dir, ": no such directory");
  if (!fs::is_directory(d/"system")) refuse(dir, ": no system/ directory; not an OpenFOAM case");
  const fs::path poly = d/"constant"/"polyMesh";
  if (!exists_either(poly/"faces")){
    for (const auto& e : fs::directory_iterator(d))
      if (e.is_directory() && e.path().filename().string().rfind("processor", 0) == 0)
        refuse(dir, " is decomposed and has no reconstructed mesh; run reconstructPar -case . "
               "(and reconstructParMesh if the mesh changed) first");
    refuse(dir, ": no constant/polyMesh");
  }
  for (const auto& e : fs::directory_iterator(d)){
    const std::string n = e.path().filename().string();
    if (e.is_directory() && is_number(n) && exists_either(e.path()/"polyMesh"/"points"))
      refuse(dir, ": ", n, "/polyMesh/points: a moving or changed mesh is not handled");
  }
  for (const char* f : {"points", "faces", "owner", "neighbour"}){
    const fs::path p = poly/f;
    check_header(fs::exists(p) ? p : fs::path(p.string() + ".gz"));
  }
}

// Every patch type is one the loader takes
void check_patch_types(const std::string& dir, const std::vector<Patch>& patches){
  for (const auto& p : patches)
    if (!contains(accepted, p.type))
      refuse(dir, ": patch ", p.name, " has type ", p.type, ", which is not handled (the accepted "
             "types: ", joined(accepted), ")");
}

// Two boundary faces on one set of points: a baffle
void check_baffles(const CaseData& c){
  const std::size_t n0 = std::size_t(c.n_internal), nb = c.nfaces() - n0;
  // Each boundary face's points sorted
  const std::int64_t off = c.face_start[n0];
  std::vector<std::int32_t> pts(c.face_points.begin() + off, c.face_points.end());
  const auto first = [&](const std::size_t b){ return pts.begin() + (c.face_start[n0 + b] - off); };
  for (std::size_t b = 0; b < nb; ++b) std::sort(first(b), first(b + 1));
  // The faces by their smallest point, then sorted by their points within each
  std::vector<std::size_t> start(c.npoints() + 1, 0);
  for (std::size_t b = 0; b < nb; ++b) ++start[std::size_t(*first(b)) + 1];
  std::partial_sum(start.begin(), start.end(), start.begin());
  std::vector<std::size_t> order(nb), at(start.begin(), start.end() - 1);
  for (std::size_t b = 0; b < nb; ++b) order[at[std::size_t(*first(b))]++] = b;
  const auto less = [&](const std::size_t a, const std::size_t b){
    return std::lexicographical_compare(first(a), first(a + 1), first(b), first(b + 1))
      || (std::equal(first(a), first(a + 1), first(b), first(b + 1)) && a < b);
  };
  for (std::size_t p = 0; p < c.npoints(); ++p){
    const auto o0 = order.begin() + std::ptrdiff_t(start[p]), o1 = order.begin() + std::ptrdiff_t(start[p + 1]);
    if (o1 - o0 < 2) continue;
    std::sort(o0, o1, less);
    for (auto o = o0 + 1; o < o1; ++o){
      if (!std::equal(first(o[-1]), first(o[-1] + 1), first(o[0]), first(o[0] + 1))) continue;
      const std::int64_t f = std::int64_t(n0 + o[-1]);
      std::string name;
      for (const auto& q : c.patches)
        if (f >= q.start && f < q.start + q.size) name = q.name;
      refuse(c.dir, ": two boundary faces on the points of face ", f, " (patch ", name,
             "): a baffle, which is not handled");
    }
  }
}

// Each cyclic pair a translation between planar patches, by its addressing:
// face i of one side is face i of the other, its points reversed
void check_cyclics(CaseData& c){
  const double tol = 1e-9*box_diagonal(c);
  const std::int64_t ni = c.n_internal;
  for (auto& a : c.patches){
    if (a.type != "cyclic") continue;
    const Patch& b = c.patches[std::size_t(a.neighbour)];
    if (a.size != b.size)
      refuse(c.dir, ": cyclic ", a.name, " has ", a.size, " faces, ", b.name, " ", b.size);
    std::array<double, 3> sep{0., 0., 0.}, n{0., 0., 0.};
    for (std::int64_t i = 0; i < a.size; ++i)
      for (int d = 0; d < 3; ++d){
        sep[d] += c.bface_centres[3*(b.start + i - ni) + d] - c.bface_centres[3*(a.start + i - ni) + d];
        n[d] += c.bface_areas[3*(a.start + i - ni) + d];
      }
    double nn = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    for (int d = 0; d < 3; ++d){
      if (a.size) sep[d] /= double(a.size);
      n[d] /= std::max(nn, 1e-300);
    }
    a.separation = sep;
    double worst = 0., plane_lo = 1e300, plane_hi = -1e300;
    for (std::int64_t i = 0; i < a.size; ++i){
      const std::int64_t fa = a.start + i, fb = b.start + i;
      const std::int64_t k = c.face_start[fa + 1] - c.face_start[fa];
      if (c.face_start[fb + 1] - c.face_start[fb] != k)
        refuse(c.dir, ": cyclic ", a.name, " faces differ in size from ", b.name, "'s");
      for (std::int64_t j = 0; j < k; ++j){
        const std::int32_t pa = c.face_points[std::size_t(c.face_start[fa] + j)];
        const std::int32_t pb = c.face_points[std::size_t(c.face_start[fb] + (k - j) % k)];
        double e = 0., h = 0.;
        for (int d = 0; d < 3; ++d){
          const double x = c.points[3*std::size_t(pa) + d];
          const double m = x + sep[d] - c.points[3*std::size_t(pb) + d];
          e += m*m;
          h += x*n[d];
        }
        worst = std::max(worst, std::sqrt(e));
        plane_lo = std::min(plane_lo, h);
        plane_hi = std::max(plane_hi, h);
      }
    }
    if (a.size && plane_hi - plane_lo > tol)
      refuse(c.dir, ": cyclic ", a.name, " is not planar; only a translation between planar "
             "patches is handled (checkMesh reports it; snap the sides and pair them again with "
             "createPatch)");
    if (worst > tol)
      refuse(c.dir, ": cyclic ", a.name, "/", b.name, " is not a translation: a point misses its "
             "image by ", worst, " (the tolerance is ", tol, ", 1e-9 of the mesh size)");
  }
}

// The normal axis of an empty pair one cell thick, -1 without empty faces
int empty_axis(const CaseData& c){
  std::vector<std::int64_t> f;
  for (const auto& p : c.patches)
    if (p.type == "empty")
      for (std::int64_t i = 0; i < p.size; ++i) f.push_back(p.start + i);
  if (f.empty()) return -1;
  const std::int64_t ni = c.n_internal;
  int axis = -1;
  for (const std::int64_t g : f){
    const double* A = &c.bface_areas[3*std::size_t(g - ni)];
    const double mag = std::sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
    int ax = 0;
    for (int d = 1; d < 3; ++d) if (std::abs(A[d]) > std::abs(A[ax])) ax = d;
    if (axis < 0) axis = ax;
    if (ax != axis || std::abs(std::abs(A[ax]) - mag) > 1e-9*mag)
      refuse(c.dir, ": the empty faces are not normal to one axis; only an axis-aligned "
             "empty pair is handled");
  }
  double lo = 1e300, hi = -1e300;
  for (const std::int64_t g : f){
    const double z = c.bface_centres[3*std::size_t(g - ni) + std::size_t(axis)];
    lo = std::min(lo, z);
    hi = std::max(hi, z);
  }
  const double tol = 1e-9*box_diagonal(c);
  std::vector<int> front(std::size_t(c.ncells), 0), back(std::size_t(c.ncells), 0);
  for (const std::int64_t g : f){
    const double z = c.bface_centres[3*std::size_t(g - ni) + std::size_t(axis)];
    const bool at_lo = std::abs(z - lo) <= tol, at_hi = std::abs(z - hi) <= tol;
    if (hi - lo <= tol || !(at_lo || at_hi))
      refuse(c.dir, ": the empty faces do not lie on two planes normal to ", "xyz"[axis]);
    ++(at_lo ? front : back)[std::size_t(c.owner[std::size_t(g)])];
  }
  for (std::int64_t i = 0; i < c.ncells; ++i)
    if (front[std::size_t(i)] != 1 || back[std::size_t(i)] != 1)
      refuse(c.dir, ": the case between the empty planes is not one cell thick; a 3D case "
             "with empty patches is not handled");
  return axis;
}

// The time directories that hold the field, by value; 0 only when alone
void read_times(const Foam::Time& runTime, const std::string& dir, const std::string& field,
                CaseData& c){
  const Foam::instantList all = runTime.times();
  std::vector<std::pair<double, std::string>> ts;
  for (const Foam::instant& t : all){
    const std::string n = t.name();
    if (n == "constant") continue;
    if (exists_either(fs::path(dir)/n/field)) ts.push_back({t.value(), n});
  }
  std::sort(ts.begin(), ts.end());
  std::vector<std::pair<double, std::string>> nz;
  for (const auto& t : ts) if (t.first != 0.) nz.push_back(t);
  if (!nz.empty()) ts.swap(nz);
  if (ts.empty()){
    for (const auto& e : fs::directory_iterator(dir))
      if (e.is_directory() && e.path().filename().string().rfind("processor", 0) == 0)
        refuse(dir, ": no time directory holds ", field, "; the times are under ",
               e.path().filename().string(), ": run reconstructPar -case . first");
    refuse(dir, ": no time directory holds ", field);
  }
  for (const auto& t : ts){
    c.time_values.push_back(t.first);
    c.times.push_back(t.second);
  }
}

// The patches of constant/polyMesh/boundary as written, (name, type, nFaces)
std::vector<Patch> boundary_entries(const Foam::Time& runTime){
  const Foam::polyBoundaryMeshEntries entries(
    Foam::IOobject("boundary", runTime.findInstance(Foam::polyMesh::meshSubDir, "boundary"),
                   Foam::polyMesh::meshSubDir, runTime, Foam::IOobject::MUST_READ,
                   Foam::IOobject::NO_WRITE, false));
  std::vector<Patch> out;
  for (const Foam::entry& e : entries){
    Patch p;
    p.name = e.keyword();
    p.type = e.dict().lookup<Foam::word>("type");
    p.size = e.dict().lookupOrDefault<Foam::label>("nFaces", 0);
    out.push_back(p);
  }
  return out;
}

void do_read_case(const std::string& dir, const std::string& velocity, CaseData& c){
  quiet();
  check_case(dir);
  c.dir = dir;
  const std::unique_ptr<Foam::Time> runTime = open_time(dir);
  // Refused by name before OpenFOAM builds a patch it would need more for
  check_patch_types(dir, boundary_entries(*runTime));
  {
    const Foam::polyMesh mesh(Foam::IOobject(Foam::polyMesh::defaultRegion, runTime->name(),
                                             *runTime, Foam::IOobject::MUST_READ));
    const Foam::pointField& X = mesh.points();
    c.points.resize(3*std::size_t(X.size()));
    for (Foam::label i = 0; i < X.size(); ++i)
      for (int d = 0; d < 3; ++d) c.points[3*std::size_t(i) + std::size_t(d)] = X[i][d];
    const Foam::faceList& F = mesh.faces();
    c.face_start.resize(std::size_t(F.size()) + 1);
    c.face_start[0] = 0;
    for (Foam::label f = 0; f < F.size(); ++f)
      c.face_start[std::size_t(f) + 1] = c.face_start[std::size_t(f)] + F[f].size();
    c.face_points.resize(std::size_t(c.face_start.back()));
    for (Foam::label f = 0; f < F.size(); ++f)
      std::copy(F[f].begin(), F[f].end(), c.face_points.begin() + c.face_start[std::size_t(f)]);
    c.owner.assign(mesh.faceOwner().begin(), mesh.faceOwner().end());
    c.neighbour.assign(mesh.faceNeighbour().begin(), mesh.faceNeighbour().end());
    c.n_internal = mesh.nInternalFaces();
    c.ncells = mesh.nCells();
    for (const Foam::polyPatch& pp : mesh.boundary()){
      Patch p;
      p.name = pp.name();
      p.type = pp.type();
      p.start = pp.start();
      p.size = pp.size();
      for (const Foam::word& g : pp.inGroups()) p.groups.push_back(g);
      if (Foam::isA<Foam::cyclicPolyPatch>(pp)){
        const auto& cp = Foam::refCast<const Foam::cyclicPolyPatch>(pp);
        p.neighbour = cp.nbrPatch().index();
        p.owner = cp.owner();
      }
      c.patches.push_back(std::move(p));
    }
    check_patch_types(dir, c.patches);
    const Foam::vectorField& C = mesh.cellCentres();
    const Foam::scalarField& V = mesh.cellVolumes();
    c.cell_centres.resize(3*std::size_t(C.size()));
    c.cell_volumes.assign(V.begin(), V.end());
    for (Foam::label i = 0; i < C.size(); ++i)
      for (int d = 0; d < 3; ++d) c.cell_centres[3*std::size_t(i) + std::size_t(d)] = C[i][d];
    const Foam::vectorField& fc = mesh.faceCentres();
    const Foam::vectorField& fa = mesh.faceAreas();
    const std::size_t nb = std::size_t(F.size() - c.n_internal);
    c.bface_centres.resize(3*nb);
    c.bface_areas.resize(3*nb);
    for (std::size_t i = 0; i < nb; ++i)
      for (int d = 0; d < 3; ++d){
        c.bface_centres[3*i + std::size_t(d)] = fc[Foam::label(i) + Foam::label(c.n_internal)][d];
        c.bface_areas[3*i + std::size_t(d)] = fa[Foam::label(i) + Foam::label(c.n_internal)][d];
      }
  }
  check_baffles(c);
  check_cyclics(c);
  c.empty_axis = empty_axis(c);
  read_times(*runTime, dir, velocity, c);
}

// A one-cell mesh whose patches each condition is built on to ask it
// whether it fixes the value: one of each patch type the conditions differ by
class Probe {
public:
  explicit Probe(const Foam::Time& runTime){
    Foam::pointField X(8);
    for (int i = 0; i < 8; ++i) X[i] = Foam::vector(i & 1, (i >> 1) & 1, (i >> 2) & 1);
    Foam::faceList F(6);
    F[0] = Foam::face(Foam::labelList({0, 2, 3, 1}));
    F[1] = Foam::face(Foam::labelList({4, 5, 7, 6}));
    F[2] = Foam::face(Foam::labelList({0, 1, 5, 4}));
    F[3] = Foam::face(Foam::labelList({2, 6, 7, 3}));
    F[4] = Foam::face(Foam::labelList({0, 4, 6, 2}));
    F[5] = Foam::face(Foam::labelList({1, 3, 7, 5}));
    mesh_ = std::make_unique<Foam::fvMesh>(
      Foam::IOobject("partracProbe", runTime.name(), runTime, Foam::IOobject::NO_READ,
                     Foam::IOobject::NO_WRITE, false),
      std::move(X), std::move(F), Foam::labelList(6, 0), Foam::labelList(), false);
    const Foam::polyBoundaryMesh& bm = static_cast<const Foam::polyMesh&>(*mesh_).boundary();
    Foam::List<Foam::polyPatch*> pp(3);
    pp[0] = new Foam::polyPatch("patch", 1, 0, 0, bm);
    pp[1] = new Foam::wallPolyPatch("wall", 1, 1, 1, bm);
    pp[2] = new Foam::polyPatch("rest", 4, 2, 2, bm);
    mesh_->addFvPatches(pp, false);
  }
  // The probe patch standing for a case patch type
  const Foam::fvPatch& patch(const std::string& type) const {
    return mesh_->boundary()[type == "wall" || type == "mappedWall" ? 1 : 0];
  }
  const Foam::fvMesh& mesh() const { return *mesh_; }
private:
  std::unique_ptr<Foam::fvMesh> mesh_;
};

template<typename Type>
struct ProbeField {
  explicit ProbeField(const Foam::fvMesh& m)
    : iF(Foam::IOobject("partracProbeField", m.time().name(), m, Foam::IOobject::NO_READ,
                        Foam::IOobject::NO_WRITE, false),
         m, Foam::dimensioned<Type>(Foam::dimless, Foam::Zero)) {}
  Foam::DimensionedField<Type, Foam::fvMesh> iF;
};

// The case's Time and probe, kept for the next field of the same case
struct Session {
  std::string dir;
  std::unique_ptr<Foam::Time> runTime;
  std::unique_ptr<Probe> probe;
  std::map<std::string, bool> fixes;   // (patch kind, condition) -> fixesValue()
};
Session& session(const std::string& dir){
  static Session s;
  if (s.dir != dir || !s.runTime){
    s.probe.reset();
    s.runTime = open_time(dir);
    s.probe = std::make_unique<Probe>(*s.runTime);
    s.fixes.clear();
    s.dir = dir;
  }
  return s;
}

// The entry of a patch in a boundaryField, as GeometricBoundaryField finds it:
// its name, then a group (the dictionary read backwards), then a pattern;
// nullptr for an empty patch without one
const Foam::entry* patch_entry(const Foam::dictionary& bf, const Patch& p){
  for (auto it = bf.begin(); it != bf.end(); ++it)
    if (it().isDict() && !it().keyword().isPattern() && it().keyword() == p.name) return &it();
  for (auto it = bf.rbegin(); it != bf.rend(); ++it)
    if (it().isDict() && !it().keyword().isPattern()
        && contains(p.groups, std::string(it().keyword())))
      return &it();
  if (p.type == "empty") return nullptr;
  return bf.lookupEntryPtr(p.name, false, true);
}

// Whether a condition fixes the value, asked of it on the probe: built by
// name, or where it has no such constructor from its dictionary with a value
// of the probe's size
template<typename Type>
bool probe_fixes(const Session& s, const Patch& p, const std::string& condition, const Foam::dictionary& pd,
                 const Foam::DimensionedField<Type, Foam::fvMesh>& iF){
  const Foam::fvPatch& patch = s.probe->patch(p.type);
  try {
    return Foam::fvPatchField<Type>::New(condition, patch, iF)().fixesValue();
  } catch (const Foam::error&) {}
  Foam::dictionary d(pd);
  Foam::OStringStream os;
  // Its other fields of the patch's size uniform too, as a contactAngle's gradient
  for (const Foam::entry& e : pd){
    if (!e.isStream() || e.keyword() == "value") continue;
    const Foam::ITstream& ts = e.stream();
    if (ts.size() < 2 || !ts[0].isWord() || ts[0].wordToken() != "nonuniform" || !ts[1].isCompound()) continue;
    const Foam::word t = ts[1].compoundToken().type();
    if (t == "List<scalar>") os << e.keyword() << " uniform 0;";
    else if (t == "List<vector>") os << e.keyword() << " uniform (0 0 0);";
  }
  os << "value uniform " << Foam::pTraits<Type>::zero << ";";
  Foam::IStringStream is(os.str());
  while (is.good()){
    Foam::autoPtr<Foam::entry> e(Foam::entry::New(is));
    if (!e.valid()) break;
    d.set(e.ptr());
  }
  return Foam::fvPatchField<Type>::New(patch, iF, d)().fixesValue();
}

// Whether OpenFOAM can build a condition here: its type in the table once the
// libraries its entry names are open
template<typename Type>
bool buildable(const std::string& condition, const Foam::dictionary& pd){
  Foam::libs.open(pd, "libs", Foam::fvPatchField<Type>::dictionaryConstructorTablePtr_);
  const auto* table = Foam::fvPatchField<Type>::dictionaryConstructorTablePtr_;
  return table && table->found(condition);
}

// Whether a file holds key, read in 1 MB blocks overlapping by the key's length less one
bool file_holds(const fs::path& file, const std::string& key){
  std::ifstream in(file, std::ios::binary);
  const std::size_t block = std::size_t(1) << 20, keep = key.size() - 1;
  std::string buf;
  std::vector<char> chunk(block);
  while (in.read(chunk.data(), std::streamsize(block)) || in.gcount() > 0){
    buf.append(chunk.data(), std::size_t(in.gcount()));
    if (buf.find(key) != std::string::npos) return true;
    if (buf.size() > keep) buf.erase(0, buf.size() - keep);
  }
  return false;
}

// The library beside OpenFOAM's own that holds a type's name as a string of
// its own, empty if none does
std::string library_of(const std::string& type){
  Dl_info info;
  if (!dladdr(reinterpret_cast<void*>(&Foam::fvPatchField<Foam::vector>::dictionaryConstructorTablePtr_), &info)
      || !info.dli_fname)
    return "";
  const std::string key = std::string(1, '\0') + type + std::string(1, '\0');
  std::vector<fs::path> libs;
  std::error_code ec;
  for (const auto& e : fs::directory_iterator(fs::path(info.dli_fname).parent_path(), ec))
    if (e.path().extension() == ".so") libs.push_back(e.path());
  std::sort(libs.begin(), libs.end());
  for (const fs::path& l : libs)
    if (file_holds(l, key)) return l.filename().string();
  return "";
}

template<typename Type>
void read_typed(const CaseData& c, Session& s, const Foam::dictionary& dict, const fs::path& file,
                FieldData& out){
  constexpr int ncomp = int(Foam::pTraits<Type>::nComponents);
  out.ncomp = ncomp;
  const Foam::Field<Type> in("internalField", dict, Foam::label(c.ncells));
  out.internal.resize(std::size_t(ncomp)*std::size_t(in.size()));
  for (Foam::label i = 0; i < in.size(); ++i)
    for (int d = 0; d < ncomp; ++d)
      out.internal[std::size_t(ncomp)*std::size_t(i) + std::size_t(d)] =
        Foam::component(in[i], Foam::direction(d));
  const Foam::dictionary& bf = dict.subDict("boundaryField");
  ProbeField<Type> pf(s.probe->mesh());
  for (const Patch& p : c.patches){
    FieldPatch fp;
    const Foam::entry* e = patch_entry(bf, p);
    if (!e && p.type != "empty")
      refuse(file.string(), ": boundaryField has no entry for patch ", p.name);
    if (!e){
      fp.condition = "empty";
      out.patches.push_back(std::move(fp));
      continue;
    }
    const Foam::dictionary& pd = e->dict();
    fp.condition = pd.lookup<Foam::word>("type");
    const bool constraint = p.type == "cyclic" || p.type == "empty" || p.type == "symmetry"
                         || p.type == "symmetryPlane";
    if (!constraint){
      const std::string key = std::string(p.type == "wall" || p.type == "mappedWall" ? "wall" : "patch")
                            + ":" + fp.condition;
      fp.unknown = !buildable<Type>(fp.condition, pd);
      if (fp.unknown && ncomp != 1){
        const std::string lib = library_of(fp.condition);
        refuse(file.string(), ": patch ", p.name, " is ", fp.condition, ", which OpenFOAM cannot build here: ",
               lib.empty() ? "no library beside OpenFOAM's names it; add the one that defines it"
                           : "its library is not loaded, probably " + lib + "; add libs (\"" + lib + "\");",
               " to system/controlDict");
      }
      if (fp.unknown) fp.fixes_value = false;
      else {
        auto it = s.fixes.find(key);
        if (it == s.fixes.end()) it = s.fixes.emplace(key, probe_fixes<Type>(s, p, fp.condition, pd, pf.iF)).first;
        fp.fixes_value = it->second;
      }
      if (pd.found("value")){
        const Foam::Field<Type> v("value", pd, Foam::label(p.size));
        fp.has_value = true;
        fp.values.resize(std::size_t(ncomp)*std::size_t(v.size()));
        for (Foam::label i = 0; i < v.size(); ++i)
          for (int d = 0; d < ncomp; ++d)
            fp.values[std::size_t(ncomp)*std::size_t(i) + std::size_t(d)] =
              Foam::component(v[i], Foam::direction(d));
      }
      else if (fp.fixes_value){
        // A condition that sets its value without writing it, as noSlip
        const Foam::tmp<Foam::fvPatchField<Type>> probe =
          Foam::fvPatchField<Type>::New(s.probe->patch(p.type), pf.iF, pd);
        fp.has_value = true;
        fp.values.resize(std::size_t(ncomp)*std::size_t(p.size));
        for (std::int64_t i = 0; i < p.size; ++i)
          for (int d = 0; d < ncomp; ++d)
            fp.values[std::size_t(ncomp)*std::size_t(i) + std::size_t(d)] =
              Foam::component(probe()[0], Foam::direction(d));
      }
    }
    out.patches.push_back(std::move(fp));
  }
}

void do_read_field(const CaseData& c, const std::string& time, const std::string& name,
                   FieldData& out){
  quiet();
  fs::path file = fs::path(c.dir)/time/name;
  if (!fs::exists(file)){
    if (!fs::exists(file.string() + ".gz")) refuse(file.string(), ": no such file");
    file = fs::path(file.string() + ".gz");
  }
  check_header(file);
  Session& s = session(c.dir);
  Foam::IFstream is((fs::path(c.dir)/time/name).string());
  Foam::IOobject io(name, time, *s.runTime, Foam::IOobject::NO_READ, Foam::IOobject::NO_WRITE, false);
  if (!io.readHeader(is)) refuse(file.string(), ": no FoamFile header");
  const Foam::dictionary dict(is);
  out.name = name;
  out.time = time;
  const std::string cls = io.headerClassName();
  if (cls == "volVectorField") read_typed<Foam::vector>(c, s, dict, file, out);
  else if (cls == "volScalarField") read_typed<Foam::scalar>(c, s, dict, file, out);
  else if (cls.rfind("surface", 0) == 0)
    refuse(file.string(), " is a ", cls, ", values on the faces, as OpenFOAM's phi, the volume flux through "
           "them; not a volScalarField or a volVectorField");
  else refuse(file.string(), " is a ", cls, ", not a volScalarField or a volVectorField");
}

bool do_has_empty(const std::string& dir){
  quiet();
  check_case(dir);
  const std::unique_ptr<Foam::Time> runTime = open_time(dir);
  for (const Patch& p : boundary_entries(*runTime))
    if (p.type == "empty" && p.size > 0) return true;
  return false;
}

// The host's streams as they were: a Time sets their precision to the case's
struct StreamState {
  std::ios_base::fmtflags out_flags = std::cout.flags(), err_flags = std::cerr.flags();
  std::streamsize out_prec = std::cout.precision(), err_prec = std::cerr.precision();
  ~StreamState(){
    std::cout.flags(out_flags);
    std::cout.precision(out_prec);
    std::cerr.flags(err_flags);
    std::cerr.precision(err_prec);
  }
};

// Runs body, any error into err: 0 done, 1 refused
template<typename Body>
int guarded(std::string* err, Body&& body){
  const StreamState streams;
  try {
    body();
    return 0;
  } catch (const Refused& r) {
    *err = r.what;
  } catch (const Foam::error& e) {
    *err = "OpenFOAM: " + std::string(e.message());
  } catch (const std::exception& e) {
    *err = e.what();
  }
  return 1;
}

}  // namespace

extern "C" {

int partrac_openfoam_read_case(const char* dir, const char* velocity, CaseData* out, std::string* err){
  return guarded(err, [&]{ do_read_case(dir, velocity, *out); });
}

int partrac_openfoam_read_field(const CaseData* c, const char* time, const char* name,
                                FieldData* out, std::string* err){
  return guarded(err, [&]{ do_read_field(*c, time, name, *out); });
}

int partrac_openfoam_has_empty(const char* dir, int* out, std::string* err){
  return guarded(err, [&]{ *out = do_has_empty(dir) ? 1 : 0; });
}

std::uint64_t partrac_openfoam_layout(){
  return layout_stamp;
}

}
