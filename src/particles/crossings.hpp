#ifndef __CROSSINGS_HPP
#define __CROSSINGS_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>
#include <omp.h>
#include "H5Cpp.h"

#include "Error.hpp"
#include "Initializer.hpp"
#include "ParticleSet.hpp"
#include "io.hpp"
#include "schedule.hpp"

// What a crossing record is of: a strip's edges, points, or points carrying a
// line element
enum class CrossingRecords { Edges, Points, Vectors };

// First crossings of the levels of X = n.(x - o): when a node's running max
// of X passes one. The levels: cross_intv a bare number dn, k dn for every
// integer k; a schedule (schedule.hpp, read as distances), its event values
// above 0, k = 1, 2, ... in turn, listed in <file>.levels as "k level". The
// origin o: a point, or start, each node's own position when first seen.
// A text line names the float64 columns, then the records:
//   Edges    c t lnrho weight length k, per edge of the node: c; t; ln(dl/dl0)
//            at the step's end; dl0/2; dl/2; k
//   Points   id t k x y z
//   Vectors  id t k x y z lnrho nx ny nz: lnrho is w, n is rhohat
// t, x (unwrapped, as dumped) and w are linear in the step at the crossing's
// fraction; rhohat too, renormalised. A new strip node starts from its seen
// neighbours' mean, so an edge split across a plane can be lost or counted
// twice. <file>.filter logs a strip's filtering: t, edges before, after.
// A checkpoint holds the files' lengths and each node's running max (a strip
// node's X too, a start its origin); a restart cuts the files to those lengths
// and appends, as the run never stopped would. A missing or empty file starts
// afresh, the running max from X; one of other records, or unknown to the
// checkpoint, is refused
class Crossings {
public:
  // intv and origin: cross_intv and cross_origin; checkpoint: the
  // checkpoint.h5 a restart resumes from, else empty
  Crossings(ParticleSet& ps, const CrossingRecords records, const Vector3d& normal, const std::string& intv,
            const std::string& origin, const std::string& path, const std::string& checkpoint)
    : kind(records), n(normal), path(path), header(header_line(records, normal, intv, origin)) {
    if (kind == CrossingRecords::Vectors && ps.carries() != TransportElement::Vector)
      partrac::fail("crossings: line element records need the particles to carry one");
    Schedule sc;
    if (const ScheduleProblem p = parse_schedule(intv, sc))
      partrac::fail(schedule_message("cross_intv", p));
    uniform_ = sc.plain();
    if (uniform_)
      dn = sc.segs[0].X;
    else
      levels_ = ScheduleLevels(sc);
    if (!try_parse_origin(origin, from_start_, o_))
      partrac::fail("cross_origin ", origin, ": not start or x,y,z");
    if (from_start_ && kind == CrossingRecords::Edges)
      partrac::fail("cross_origin=start follows points and line elements: a strip's new nodes have no "
                    "start of their own");
    max_ = ps.add_column();
    if (from_start_)
      for (Uint& c : o_cols_) c = ps.add_column();
    if (kind == CrossingRecords::Edges)
      xn_ = ps.add_column();
    else
      for (Uint& c : x0_) c = ps.add_column();
    if (kind == CrossingRecords::Vectors){
      w0_ = ps.add_column();
      for (Uint& c : r0_) c = ps.add_column();
    }
    const bool resumes = resume(ps, checkpoint);
    const auto mode = std::ios::binary | (resumes ? std::ios::app : std::ios::trunc);
    out.open(path, mode);
    if (!out)
      partrac::fail("crossings: cannot open ", path);
    if (kind == CrossingRecords::Edges){
      flog.open(path + ".filter", mode);
      if (!flog)
        partrac::fail("crossings: cannot open ", path, ".filter");
    }
    if (!uniform_){
      llog.open(path + ".levels", mode);
      if (!llog)
        partrac::fail("crossings: cannot open ", path, ".levels");
    }
    const State st = state(ps);
    if (!resumes){
      out << header << "\n";
      bytes_ = header.size() + 1;
      for (Uint i = 0; i < ps.N(); ++i){
        if (from_start_)
          for (int d = 0; d < 3; ++d) st.o[d][i] = ps.x(i)[d];
        st.max[i] = level(ps, st, i);
        if (st.xn) st.xn[i] = st.max[i];
      }
    }
    for (Uint i = 0; i < ps.N(); ++i)
      remember(ps, st, i);
    write_levels();
  }

  // Origins for the new nodes, where they are: before their first step
  // (cross_origin=start; their birth, as injected)
  void place(ParticleSet& ps){
    if (!from_start_)
      return;
    const State st = state(ps);
    for (Uint i = 0; i < ps.N(); ++i)
      if (std::isnan(st.o[0][i]))
        for (int d = 0; d < 3; ++d) st.o[d][i] = ps.x(i)[d];
  }

  static const char* columns(const CrossingRecords records){
    switch (records){
      case CrossingRecords::Edges: return "c t lnrho weight length k";
      case CrossingRecords::Points: return "id t k x y z";
      default: return "id t k x y z lnrho nx ny nz";
    }
  }

  // The file's text line
  static std::string header_line(const CrossingRecords records, const Vector3d& n, const std::string& intv,
                                 const std::string& origin){
    std::ostringstream os;
    os << "partrac crossings: float64 " << columns(records)
       << "; n = (" << n[0] << ", " << n[1] << ", " << n[2] << "), cross_intv = " << intv
       << ", cross_origin = " << origin;
    return os.str();
  }

  // After the step from t to t + h
  void record(ParticleSet& ps, const EdgesType& edges, const Node2EdgesType& node2edges,
              const double t, const double h){
    const Uint N = ps.N();
    const bool edges_ = kind == CrossingRecords::Edges;
    if (edges_ && node2edges.size() != N)
      partrac::fail("crossings: node2edges has ", node2edges.size(), " rows for ", N, " nodes");
    const State st = state(ps);
    double* const m = st.max;
    double* const xn = st.xn;
    // New nodes only where ids were made since
    if (ps.ids_made() != ids_seen_){
      start_new(ps, st, edges, node2edges);
      ids_seen_ = ps.ids_made();
    }
    // A schedule's levels made beyond the farthest node
    if (!uniform_){
      double far = -std::numeric_limits<double>::infinity();
      #pragma omp parallel for schedule(static) reduction(max:far) if(N >= parallel_from)
      for (Uint i = 0; i < N; ++i)
        far = std::max(far, level(ps, st, i));
      levels_.reach(far);
      write_levels();
    }
    const std::vector<double>& L = levels_.levels();
    found_.resize(std::max<std::size_t>(found_.size(), omp_get_max_threads()));
    for (auto& f : found_)
      f.clear();
    #pragma omp parallel for schedule(static) if(N >= parallel_from)
    for (Uint i = 0; i < N; ++i){
      const double s = level(ps, st, i), m0 = m[i];
      // Before the step, then after it
      Vector3d x0 = Vector3d::Zero(), r0 = Vector3d::Zero();
      double w0 = 0., s0;
      if (edges_){
        s0 = xn[i];
        xn[i] = s;
      }
      else {
        x0 = {st.x0[0][i], st.x0[1][i], st.x0[2][i]};
        s0 = (x0 - origin(st, i)).dot(n);
      }
      if (kind == CrossingRecords::Vectors){
        w0 = st.w0[i];
        r0 = {st.r0[0][i], st.r0[1][i], st.r0[2][i]};
      }
      remember(ps, st, i);
      if (s <= m0)
        continue;
      auto& f = found_[omp_get_thread_num()];
      // Level lk, of index k
      const auto cross = [&](const double lk, const double k){
        const double a = std::max((lk - s0) / (s - s0), 0.), tk = t + h * a;
        if (edges_){
          // A strip node's two edges in ascending order, as a restart rebuilds node2edges
          const AdjRowType& row = node2edges[i];
          const bool flip = row.size() == 2 && row[1] < row[0];
          for (Uint q = 0; q < row.size(); ++q){
            const Uint e = row[flip ? 1 - q : q];
            const double dl = ps.dist(edges[e].first[0], edges[e].first[1]), dl0 = edges[e].second;
            f.insert(f.end(), {ps.c(i), tk, std::log(dl / dl0), dl0 / 2, dl / 2, k});
          }
          return;
        }
        const Vector3d xk = x0 + a * (ps.x(i) - x0);
        f.insert(f.end(), {double(ps.id(i)), tk, k, xk[0], xk[1], xk[2]});
        if (kind == CrossingRecords::Vectors){
          Vector3d rk = r0 + a * (ps.rhohat(i) - r0);
          rk /= rk.norm();
          f.insert(f.end(), {w0 + a * (ps.w(i) - w0), rk[0], rk[1], rk[2]});
        }
      };
      if (uniform_)
        for (long k = static_cast<long>(std::floor(m0 / dn)) + 1; k * dn <= s; ++k)
          cross(k * dn, double(k));
      else
        for (std::size_t j = levels_.above(m0); j < L.size() && L[j] <= s; ++j)
          cross(L[j], double(j + 1));
      m[i] = s;
    }
    for (const auto& f : found_){
      out.write(reinterpret_cast<const char*>(f.data()), f.size() * sizeof(f[0]));
      bytes_ += f.size() * sizeof(f[0]);
    }
    logged_ = skip_log_ = false;
  }

  // A strip's filtering. A restart from a checkpoint written after its step's
  // filtering repeats that step: the filtering, done, is not logged again
  void filtered(const double t, const Uint before, const Uint after){
    logged_ = true;
    if (skip_log_)
      return;
    std::ostringstream os;
    os << t << " " << before << " " << after << "\n";
    const std::string line = os.str();
    flog.write(line.data(), line.size());
    filter_bytes_ += line.size();
  }

  // Into the checkpoint, the files flushed: the record file's, the filter
  // log's and the level list's lengths, whether this step's filtering is
  // logged, the levels listed; by slot the running max of X, a strip node's
  // X, a start
  void checkpoint(H5::H5File& h5f, const ParticleSet& ps){
    out.flush();
    for (std::ofstream* f : {&flog, &llog})
      if (f->is_open())
        f->flush();
    if (!out || (flog.is_open() && !flog) || (llog.is_open() && !llog))
      partrac::fail("crossings: cannot write ", path);
    const H5::StrType text(H5::PredType::C_S1, header.size());
    h5f.createAttribute("crossings", text, H5::DataSpace(H5S_SCALAR)).write(text, header);
    const std::uint64_t lengths[5] = {bytes_, filter_bytes_, logged_, levels_bytes_, levels_written_};
    const hsize_t five = 5;
    h5f.createAttribute("crossings_bytes", H5::PredType::NATIVE_UINT64, H5::DataSpace(1, &five))
       .write(H5::PredType::NATIVE_UINT64, lengths);
    scalar2hdf5(h5f, "crossings_max", ps.column(max_), ps.N());
    if (kind == CrossingRecords::Edges)
      scalar2hdf5(h5f, "crossings_xn", ps.column(xn_), ps.N());
    if (from_start_){
      std::vector<Vector3d> o(ps.N());
      for (Uint i = 0; i < ps.N(); ++i)
        o[i] = {ps.column(o_cols_[0])[i], ps.column(o_cols_[1])[i], ps.column(o_cols_[2])[i]};
      vector2hdf5(h5f, "crossings_start", o, ps.N());
    }
  }

private:
  const CrossingRecords kind;
  const Vector3d n;
  const std::string path, header;
  bool uniform_ = true;      // levels k dn, else a schedule's
  double dn = 0.;
  ScheduleLevels levels_;
  bool from_start_ = false;  // the origin each node's start, else o_
  Vector3d o_ = Vector3d::Zero();
  std::ofstream out, flog, llog;
  std::uint64_t bytes_ = 0, filter_bytes_ = 0, levels_bytes_ = 0;   // the files' lengths
  std::uint64_t levels_written_ = 0;   // levels in <file>.levels
  bool logged_ = false;     // a filtering since the last step
  bool skip_log_ = false;   // a restart's: its step's filtering logged
  // ParticleSet columns: the running max of x.n, NaN for a node not seen yet;
  // before the step, x.n (edges), else the position, and w and rhohat (vectors)
  Uint max_ = 0, xn_ = 0, w0_ = 0;
  std::array<Uint, 3> x0_{}, r0_{}, o_cols_{};   // and the start, NaN until placed

  // Nodes from which record() runs in parallel: below, waking the threads
  // costs more than the loop
  static constexpr Uint parallel_from = 4000;
  Uint ids_seen_ = Uint(-1);   // ps.ids_made() at the last record(); none yet
  std::vector<std::vector<double>> found_;   // records by thread
  std::vector<std::vector<Uint>> new_;       // new nodes by thread

  // The columns' rows; null where not kept
  struct State {
    double* max;
    double* xn;
    std::array<double*, 3> x0, r0, o;
    double* w0;
  };
  State state(ParticleSet& ps) const {
    const auto col = [&](const Uint k){ return ps.column(k).data(); };
    State st{col(max_), nullptr, {}, {}, {}, nullptr};
    if (from_start_)
      for (int d = 0; d < 3; ++d) st.o[d] = col(o_cols_[d]);
    if (kind == CrossingRecords::Edges)
      st.xn = col(xn_);
    else
      for (int d = 0; d < 3; ++d) st.x0[d] = col(x0_[d]);
    if (kind == CrossingRecords::Vectors){
      st.w0 = col(w0_);
      for (int d = 0; d < 3; ++d) st.r0[d] = col(r0_[d]);
    }
    return st;
  }

  // The origin of the node in slot i
  Vector3d origin(const State& st, const Uint i) const {
    return from_start_ ? Vector3d(st.o[0][i], st.o[1][i], st.o[2][i]) : o_;
  }
  // Its X now
  double level(const ParticleSet& ps, const State& st, const Uint i) const {
    return (ps.x(i) - origin(st, i)).dot(n);
  }

  // The levels made since into <file>.levels, "k level" a line
  void write_levels(){
    const std::vector<double>& L = levels_.levels();
    if (uniform_ || levels_written_ >= L.size())
      return;
    std::ostringstream os;
    os.precision(17);
    for (std::size_t j = levels_written_; j < L.size(); ++j)
      os << j + 1 << " " << L[j] << "\n";
    const std::string lines = os.str();
    llog.write(lines.data(), lines.size());
    levels_bytes_ += lines.size();
    levels_written_ = L.size();
  }

  // The state before the next step of the node in slot i
  void remember(const ParticleSet& ps, const State& st, const Uint i) const {
    if (kind == CrossingRecords::Edges)
      return;
    const Vector3d x = ps.x(i);
    for (int d = 0; d < 3; ++d) st.x0[d][i] = x[d];
    if (kind == CrossingRecords::Vectors){
      st.w0[i] = ps.w(i);
      const Vector3d r = ps.rhohat(i);
      for (int d = 0; d < 3; ++d) st.r0[d][i] = r[d];
    }
  }

  // New nodes, not seen yet (NaN max), in slot order: their seen neighbours'
  // mean before the step (edges), else their own position; from their start,
  // placed before the step, X = 0 there. A new neighbour earlier in that order
  // counts as seen, whatever the threads
  void start_new(const ParticleSet& ps, const State& st, const EdgesType& edges,
                 const Node2EdgesType& node2edges){
    const Uint N = ps.N();
    double* const m = st.max;
    new_.resize(std::max<std::size_t>(new_.size(), omp_get_max_threads()));
    for (auto& f : new_)
      f.clear();
    #pragma omp parallel for schedule(static) if(N >= parallel_from)
    for (Uint i = 0; i < N; ++i)
      if (std::isnan(m[i]))
        new_[omp_get_thread_num()].push_back(i);
    for (const auto& f : new_)
      for (const Uint i : f){
        double s = 0., ms = 0.;
        int count = 0;
        if (kind == CrossingRecords::Edges)
          for (const Uint e : node2edges[i]){
            const Uint j = edges[e].first[0] == i ? edges[e].first[1] : edges[e].first[0];
            if (!std::isnan(m[j])){ s += st.xn[j]; ms += m[j]; ++count; }
          }
        if (from_start_ && !std::isnan(st.o[0][i])){
          m[i] = 0.;
          for (int d = 0; d < 3; ++d) st.x0[d][i] = st.o[d][i];
          if (kind == CrossingRecords::Vectors){
            st.w0[i] = ps.w(i);
            for (int d = 0; d < 3; ++d) st.r0[d][i] = ps.rhohat(i)[d];
          }
          continue;
        }
        if (from_start_)
          for (int d = 0; d < 3; ++d) st.o[d][i] = ps.x(i)[d];
        const double s0 = count > 0 ? s / count : level(ps, st, i);
        m[i] = count > 0 ? ms / count : s0;
        if (kind == CrossingRecords::Edges)
          st.xn[i] = s0;
        remember(ps, st, i);
      }
  }

  // Whether a restart continues the files: then cut to the lengths its
  // checkpoint holds, the state restored. False on a new run, or if the file
  // is missing or empty; fails on other records, or records its checkpoint
  // does not know
  bool resume(ParticleSet& ps, const std::string& checkpoint){
    if (checkpoint.empty())
      return false;
    std::string line;
    {
      std::ifstream in(path, std::ios::binary);
      if (!in || !std::getline(in, line))
        return false;
    }
    if (line != header)
      partrac::fail("crossings: ", path, " holds other records (", line, "), not ", header);
    std::uint64_t lengths[5] = {0, 0, 0, 0, 0};
    std::vector<double> m, xn;
    std::vector<Vector3d> o;
    if (!read_state(checkpoint, ps.N(), lengths, m, xn, o))
      partrac::fail("crossings: the checkpoint holds no state of ", path,
                    "; move it away to start the records afresh");
    cut(path, lengths[0]);
    if (kind == CrossingRecords::Edges)
      cut(path + ".filter", lengths[1]);
    if (!uniform_)
      cut(path + ".levels", lengths[3]);
    bytes_ = lengths[0];
    filter_bytes_ = lengths[1];
    skip_log_ = lengths[2] != 0;
    levels_bytes_ = lengths[3];
    levels_written_ = lengths[4];
    // The levels listed, made again
    while (levels_.levels().size() < levels_written_)
      levels_.reach(levels_.levels().empty() ? 0. : levels_.levels().back());
    std::copy(m.begin(), m.end(), ps.column(max_).begin());
    if (kind == CrossingRecords::Edges)
      std::copy(xn.begin(), xn.end(), ps.column(xn_).begin());
    for (Uint i = 0; i < o.size(); ++i)
      for (int d = 0; d < 3; ++d) ps.column(o_cols_[d])[i] = o[i][d];
    return true;
  }

  // The state of these records in checkpoint, of N nodes; false if none
  bool read_state(const std::string& checkpoint, const Uint N, std::uint64_t* lengths,
                  std::vector<double>& m, std::vector<double>& xn, std::vector<Vector3d>& o) const {
    if (!std::ifstream(checkpoint))
      return false;
    try {
      H5::H5File h5f(checkpoint, H5F_ACC_RDONLY);
      if (H5Aexists(h5f.getId(), "crossings") <= 0)
        return false;
      H5::Attribute a = h5f.openAttribute("crossings");
      std::string held;
      a.read(a.getStrType(), held);
      if (held != header)
        return false;
      h5f.openAttribute("crossings_bytes").read(H5::PredType::NATIVE_UINT64, lengths);
      if (hdf5_rows(h5f, "crossings_max") != N)
        partrac::fail(checkpoint, ": crossings_max has ", hdf5_rows(h5f, "crossings_max"), " rows for ",
                      N, " nodes");
      hdf52doubles(h5f, "crossings_max", N, 1, m);
      if (kind == CrossingRecords::Edges)
        hdf52doubles(h5f, "crossings_xn", N, 1, xn);
      if (from_start_)
        hdf52vectors(h5f, "crossings_start", N, o);
    } catch (const H5::Exception&){
      partrac::fail("crossings: cannot read the state in ", checkpoint);
    }
    return true;
  }

  // The file at path cut to length; fails if shorter
  static void cut(const std::string& file, const std::uint64_t length){
    std::error_code ec;
    const std::uintmax_t size = std::filesystem::exists(file) ? std::filesystem::file_size(file, ec) : 0;
    if (ec || size < length)
      partrac::fail("crossings: ", file, " holds ", size, " bytes, fewer than the ", length,
                    " its checkpoint counted");
    std::filesystem::resize_file(file, length, ec);
    if (ec)
      partrac::fail("crossings: cannot cut ", file, " to ", length, " bytes: ", ec.message());
  }
};

#endif
