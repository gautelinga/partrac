#ifndef __STEPS_IMPL_HPP
#define __STEPS_IMPL_HPP

// Loop bodies of the steps; included only by the per-interpolator step sources

#include <algorithm>
#include <cmath>
#include <vector>
#include <omp.h>
#include "ParticleSet.hpp"
#include "RKIntegrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "SpatialIntegrator.hpp"
#include "tableau.hpp"
#include <type_traits>
#include <utility>

// A step reads velocity, acceleration and their gradients: an interpolator
// with an evaluate_motion leaves the rest out; called on the concrete type,
// so never through the vtable
template<typename Interp, typename = void>
struct has_evaluate_motion : std::false_type {};
template<typename Interp>
struct has_evaluate_motion<Interp, std::void_t<decltype(std::declval<Interp&>().evaluate_motion(
  std::declval<const Vector3d&>(), 0., std::declval<const CellPos&>(), std::declval<PointValues&>()))>> : std::true_type {};

template<typename Interp>
inline void evaluate_motion(Interp& intp, const Vector3d& x, const double t, const CellPos& pos, PointValues& ptvals){
  if constexpr (has_evaluate_motion<Interp>::value) intp.evaluate_motion(x, t, pos, ptvals);
  else intp.evaluate(x, t, pos, ptvals);
}

// A stage located from scratch: inside if a cell holds it
template<typename Interp>
struct LocatedEval {
  Interp& intp;
  PointValues& ptvals;
  CellPos& pos;
  bool operator()(const Vector3d& x, const double t){
    if (!intp.locate(x, t, pos))
      return false;
    evaluate_motion(intp, x, t, pos, ptvals);
    return true;
  }
  Vector3d u(){ return ptvals.get_u(); }
  Matrix3d J(){ return ptvals.get_J(); }
  bool end(const Vector3d& x, const double t){ return intp.locate(x, t, pos); }
};

// Pieces of a step that left the fluid: halves, down to this fraction
constexpr int max_pieces = 64;

// A piece: taken, failed (halved again), or stopped at its start, which halving keeps
enum class Piece { taken, failed, stopped };

// The step from t over dt in halves, a failed piece halved again down to
// dt/max_pieces; piece(t, h) takes one. Failed: a smallest piece failed
template<typename Fn>
inline Piece in_pieces(Fn&& piece, const double t, const double dt){
  constexpr int units = max_pieces;
  const double unit = dt/units;
  int u = 0, size = units/2;
  while (u < units){
    const Piece r = piece(t + u*unit, size*unit);
    if (r == Piece::taken){
      u += size;
      // The largest piece aligned here
      while (size < units/2 && u % (2*size) == 0)
        size *= 2;
    }
    else if (r == Piece::stopped || size == 1)
      return r;
    else
      size /= 2;
  }
  return Piece::taken;
}

// The displacement over h from the fields at its start
__attribute__((always_inline))
inline Vector3d drift(PointValues& pv, const int int_order, const double h){
  Vector3d dx = pv.get_u() * h;
  if (int_order >= 2)
    dx += 0.5 * (pv.get_a() + pv.get_Ju()) * h * h;
  return dx;
}

// The carried element after h from the gradients at its start: el from n, or
// the frame Fel from F
template<typename M>
__attribute__((always_inline))
inline M carried(PointValues& pv, const int int_order, const double h, const M& n){
  const Matrix3d J1 = pv.get_J();
  M e = n;
  e += J1 * e * h;
  if (int_order >= 2)
    e += 0.5*(J1*(J1*n) + pv.get_grada()*n) * h * h;
  return e;
}

template<TransportElement E>
__attribute__((always_inline))
inline void carry(PointValues& pv, const int int_order, const double h, [[maybe_unused]] const Vector3d& n,
                  [[maybe_unused]] const Matrix3d& F, [[maybe_unused]] Vector3d& el,
                  [[maybe_unused]] Matrix3d& Fel){
  if constexpr (E == TransportElement::Vector) el = carried(pv, int_order, h, n);
  if constexpr (E == TransportElement::Tensor) Fel = carried(pv, int_order, h, F);
}

// Particle i through the pieces of a step; pos: the cell of the piece's start
template<TransportElement E>
struct PieceState {
  Vector3d x;
  CellPos pos;
  Vector3d n = Vector3d::Zero();
  double w = 0.;
  Matrix3d F = Matrix3d::Identity();
  PieceState(ParticleSet& ps, const Uint i) : x(ps.x(i)) {
    pos.id = ps.get_cell_id(i);
    if constexpr (E == TransportElement::Vector){ n = ps.rhohat(i); w = ps.w(i); }
    if constexpr (E == TransportElement::Tensor) F = ps.frame(i);
  }
  // A piece over dx, its end in p
  void advance(const Vector3d& dx, const CellPos& p, const Vector3d& el, const Matrix3d& Fel){
    x += dx;
    pos = p;
    if constexpr (E == TransportElement::Vector){
      const double len = el.norm();
      n = el/len;
      w += log(len);
    }
    if constexpr (E == TransportElement::Tensor) F = Fel;
  }
  void store(ParticleSet& ps, const Uint i, const double dt_loc) const {
    ps.set_x(i, x);
    ps.set_t_loc(i, ps.t_loc(i) + dt_loc);
    ps.set_cell_id(i, pos.id);
    if constexpr (E == TransportElement::Vector){ ps.set_rhohat(i, n); ps.set_w(i, w); }
    if constexpr (E == TransportElement::Tensor) ps.advance_frame(i, F);
  }
};

// A step with a stage or its end outside, again in pieces with every stage
// inside; false: the particle leaves the fluid. Out of the loop: rare and bulky
template<TransportElement E, typename Interp>
__attribute__((noinline))
bool rk4_pieces(Interp& intp, ParticleSet& ps, const Uint i, const double t, const double dt){
  PieceState<E> s(ps, i);
  PointValues ptvals(intp.get_U0());
  // The first stage at the start, kept while its piece is halved
  bool have_stage1 = false, located = false;
  Vector3d u1;
  [[maybe_unused]] Matrix3d J1;
  const Piece r = in_pieces([&](const double ta, const double h){
    if (!have_stage1){
      // The first start located here; later ones are the last piece's end
      if (!located && !intp.locate(s.x, ta, s.pos))
        return Piece::stopped;
      evaluate_motion(intp, s.x, ta, s.pos, ptvals);
      u1 = ptvals.get_u();
      if constexpr (E != TransportElement::Point) J1 = ptvals.get_J();
      have_stage1 = located = true;
    }
    CellPos p = s.pos;
    LocatedEval<Interp> ev{intp, ptvals, p};
    Vector3d dx, el;
    Matrix3d dF;
    if (!rk_stages_from<RK4Tableau, E, true>(ev, u1, J1, s.x, s.n, s.F, ta, h, dx, el, dF))
      return Piece::failed;
    [[maybe_unused]] Matrix3d Fel;
    if constexpr (E == TransportElement::Tensor) Fel = s.F + dF;
    s.advance(dx, p, el, Fel);
    have_stage1 = false;
    return Piece::taken;
  }, t, dt);
  if (r == Piece::taken)
    s.store(ps, i, dt);
  return r == Piece::taken;
}

// An explicit step with its end outside, again in pieces of the same order,
// each with its share of the step's noise; start: the fields at x, dx_rw: the
// step. Declined: a start outside, or noisy on an interpolator that
// reflects. Out of the loop: rare and bulky
template<TransportElement E, typename Interp>
__attribute__((noinline))
bool explicit_pieces(Interp& intp, ParticleSet& ps, const Uint i, const double t, const double dt,
                     const int int_order, const bool noisy, PointValues& start,
                     const Vector3d& dx_rw){
  if (noisy && intp.can_reflect)
    return false;
  PieceState<E> s(ps, i);
  // A start outside: declined
  if (!intp.locate(s.x, t, s.pos))
    return false;
  // Noise: the step less its drift
  Vector3d noise = Vector3d::Zero();
  if (noisy)
    noise = dx_rw - drift(start, int_order, dt);
  PointValues ptvals(intp.get_U0());
  // The fields at the piece's start, kept while it is halved
  PointValues* at = &start;
  const Piece r = in_pieces([&](const double ta, const double h){
    if (!at){
      evaluate_motion(intp, s.x, ta, s.pos, ptvals);
      at = &ptvals;
    }
    Vector3d dx = drift(*at, int_order, h);
    if (noisy)
      dx += noise * (h/dt);
    CellPos p = s.pos;
    if (!intp.locate(s.x + dx, ta + h, p))
      return Piece::failed;
    Vector3d el = s.n;
    Matrix3d Fel = s.F;
    carry<E>(*at, int_order, h, s.n, s.F, el, Fel);
    s.advance(dx, p, el, Fel);
    at = nullptr;
    return Piece::taken;
  }, t, dt);
  if (r == Piece::taken)
    s.store(ps, i, dt);
  return r == Piece::taken;
}

// A spatial step from a start inside with its end outside, again in pieces of
// the path length, fields frozen at t; start: the fields at x. left: a
// smallest piece's end outside, not a start at rest. Out of the loop: rare
// and bulky
template<TransportElement E, typename Interp>
__attribute__((noinline))
bool spatial_pieces(Interp& intp, ParticleSet& ps, const Uint i, const double t, const double ds,
                    const int int_order, const double u_min, PointValues& start, bool& left){
  PieceState<E> s(ps, i);
  PointValues ptvals(intp.get_U0());
  // The fields at the piece's start, kept while it is halved
  PointValues* at = &start;
  double t_loc = 0.;
  const Piece r = in_pieces([&](const double, const double dsp){
    if (!at){
      evaluate_motion(intp, s.x, t, s.pos, ptvals);
      at = &ptvals;
    }
    const double uabs = at->get_u().norm();
    // At rest: no path
    if (!(uabs > u_min))
      return Piece::stopped;
    const double h = dsp / uabs;
    const Vector3d dx = drift(*at, int_order, h);
    CellPos p = s.pos;
    if (!intp.locate(s.x + dx, t, p))
      return Piece::failed;
    Vector3d el = s.n;
    Matrix3d Fel = s.F;
    carry<E>(*at, int_order, h, s.n, s.F, el, Fel);
    s.advance(dx, p, el, Fel);
    t_loc += h;
    at = nullptr;
    return Piece::taken;
  }, 0., ds);
  left = r == Piece::failed;
  if (r == Piece::taken)
    s.store(ps, i, t_loc);
  return r == Piece::taken;
}

template<TransportElement E, typename Interp>
PARTRAC_HOT_LOOP
std::vector<Uint> RK4Integrator::step(Interp& intp, ParticleSet& ps, const double t, const double dt) {
    std::vector<Uint> outside_nodes;
    // Parallel over particles
    #pragma omp parallel
    {
    std::vector<Uint> outside_nodes_loc;
    Uint n_accepted_loc = 0;
    Uint n_declined_loc = 0;

    #pragma omp for
    for (Uint i=0; i < ps.N(); ++i){
        Vector3d x = ps.x(i);
        CellPos pos;
        pos.id = ps.get_cell_id(i);
        PointValues ptvals(intp.get_U0());
        // Carried element
        [[maybe_unused]] const Vector3d n = E == TransportElement::Vector ? ps.rhohat(i) : Vector3d::Zero();
        [[maybe_unused]] const Matrix3d F = E == TransportElement::Tensor ? ps.frame(i) : Matrix3d::Zero();
        Vector3d dx;
        [[maybe_unused]] Vector3d el;
        [[maybe_unused]] Matrix3d dF;

        LocatedEval<Interp> ev{intp, ptvals, pos};
        // An outside stage or end: pieces
        if (rk_stages<RK4Tableau, E, false>(ev, x, n, F, t, dt, dx, el, dF)){
            ps.set_x(i, x + dx);
            ps.set_t_loc(i, ps.t_loc(i) + dt);
            if constexpr (E == TransportElement::Vector){
                const double len = el.norm();
                ps.set_rhohat(i, el/len);
                ps.set_w(i, ps.w(i) + log(len));
            }
            if constexpr (E == TransportElement::Tensor)
                ps.advance_frame(i, F + dF);
            ps.set_cell_id(i, pos.id);
            ++n_accepted_loc;
        }
        else if (rk4_pieces<E>(intp, ps, i, t, dt))
            ++n_accepted_loc;
        else {
            outside_nodes_loc.push_back(i);
            ++n_declined_loc;
        }
    }
    #pragma omp critical
    {
        outside_nodes.insert(outside_nodes.end(), outside_nodes_loc.begin(), outside_nodes_loc.end());
        n_accepted += n_accepted_loc;
        n_declined += n_declined_loc;
    }
    }
    // Sort: threads merge out of order
    std::sort(outside_nodes.begin(), outside_nodes.end());
    return outside_nodes;
}

template<TransportElement E, typename Interp>
PARTRAC_HOT_LOOP
std::vector<Uint> ExplicitIntegrator::step(Interp& intp, ParticleSet& ps, const double t, const double dt) {
    std::vector<Uint> outside_nodes;

    #pragma omp parallel
    {
        std::normal_distribution<double> _rnd_normal(0., 1.0);
        std::mt19937& gen = gens[omp_get_thread_num()];
        std::vector<Uint> outside_nodes_loc;
        Uint n_accepted_loc = 0;
        Uint n_declined_loc = 0;

        double sqrt2Dmdt = sqrt(2*Dm*dt);
        const bool reflects = Dm > 0.0 && intp.can_reflect;

        # pragma omp for
        for (Uint i=0; i < ps.N(); ++i){
            Vector3d x = ps.x(i);
            CellPos pos;
            pos.id = ps.get_cell_id(i);

            PointValues ptvals(intp.get_U0());

            // Carried element
            [[maybe_unused]] Vector3d n0, el;
            [[maybe_unused]] Matrix3d F0, Fel;
            if constexpr (E == TransportElement::Vector){ n0 = ps.rhohat(i); el = n0; }
            if constexpr (E == TransportElement::Tensor){ F0 = ps.frame(i); Fel = F0; }

            const bool x_inside = intp.locate(x, t, pos);
            Vector3d dx_rw = Vector3d::Zero();
            if (x_inside){
                // Outside: zero velocity
                evaluate_motion(intp, x, t, pos, ptvals);
                dx_rw = ptvals.get_u() * dt;
                if (int_order >= 2){
                    dx_rw += 0.5 * (ptvals.get_a() + ptvals.get_Ju()) * dt * dt;
                }
                // The gradient, for what carries it
                if constexpr (E != TransportElement::Point){
                    const Matrix3d J1 = ptvals.get_J();
                    if constexpr (E == TransportElement::Vector){
                        el += J1 * el * dt;
                        if (int_order >= 2)
                            el += 0.5*(J1*(J1*n0) + ptvals.get_grada()*n0) * dt * dt;
                    }
                    if constexpr (E == TransportElement::Tensor){
                        Fel += J1 * Fel * dt;
                        if (int_order >= 2)
                            Fel += 0.5*(J1*(J1*F0) + ptvals.get_grada()*F0) * dt * dt;
                    }
                }
            }
            if (Dm > 0.0){
                // TODO: Consider trying multiple times
                Vector3d eta = {_rnd_normal(gen),
                                _rnd_normal(gen),
                                _rnd_normal(gen)};
                dx_rw += sqrt2Dmdt * eta;
            }
            // Diffusive steps walk off the walls
            const bool is_inside = (reflects && x_inside)
                ? intp.reflect(x, dx_rw, pos)
                : intp.locate(x+dx_rw, t+dt, pos);
            if (is_inside){
                ps.set_x(i, x + dx_rw);
                ps.set_t_loc(i, ps.t_loc(i) + dt);
                ps.set_cell_id(i, pos.id);
                if constexpr (E == TransportElement::Vector){
                    const double len = el.norm();
                    ps.set_rhohat(i, el/len);
                    ps.set_w(i, ps.w(i) + log(len));
                }
                if constexpr (E == TransportElement::Tensor)
                    ps.advance_frame(i, Fel);
                ++n_accepted_loc;
            }
            // An end outside: pieces
            else if (explicit_pieces<E>(intp, ps, i, t, dt, int_order, Dm > 0.0, ptvals, dx_rw))
                ++n_accepted_loc;
            else {
                outside_nodes_loc.push_back(i);
                ++n_declined_loc;
            }
        }
        #pragma omp critical
        {
            outside_nodes.insert(outside_nodes.end(), outside_nodes_loc.begin(), outside_nodes_loc.end());
            n_accepted += n_accepted_loc;
            n_declined += n_declined_loc;
        }
    }
    // Sort: threads merge out of order
    std::sort(outside_nodes.begin(), outside_nodes.end());
    return outside_nodes;
}

template<TransportElement E, typename InterpolType>
PARTRAC_HOT_LOOP
std::vector<Uint> SpatialIntegrator::step(InterpolType& intp, ParticleSet& ps, const double t, const double ds) {
    std::vector<Uint> outside_nodes;
    // Nodes are independent; only the tally is shared
    #pragma omp parallel
    {
    std::vector<Uint> outside_nodes_loc;
    Uint n_accepted_loc = 0;
    Uint n_declined_loc = 0;
    Uint n_left_loc = 0;
    bool is_inside, tried;
    double uabs_est, dt;
    Vector3d dx;
    Vector3d el;
    Matrix3d F;

    #pragma omp for
    for (Uint i=0; i < ps.N(); ++i){
        // Done
        if (ps.t_loc(i) >= m_T){
            outside_nodes_loc.push_back(i);
            ++n_declined_loc;
            continue;
        }
        Vector3d x = ps.x(i);
        CellPos pos;
        pos.id = ps.get_cell_id(i);

        PointValues ptvals(intp.get_U0());
        // Outside: zero velocity, so declined below, as left
        bool left = !intp.locate(x, t, pos);
        if (!left)
            evaluate_motion(intp, x, t, pos, ptvals);

        Vector3d u_1 = ptvals.get_u();

        uabs_est = u_1.norm();

        is_inside = tried = false;
        if (uabs_est > m_u_min){
            dt = ds / uabs_est;

            dx = u_1 * dt;

            // Second-order terms
            if (m_int_order >= 2){
                dx += 0.5 * (ptvals.get_a() + ptvals.get_Ju()) * dt * dt;
            }
            [[maybe_unused]] Vector3d n;
            [[maybe_unused]] Matrix3d F0;
            if constexpr (E == TransportElement::Vector){ n = ps.rhohat(i); el = n; }
            if constexpr (E == TransportElement::Tensor){ F0 = ps.frame(i); F = F0; }
            carry<E>(ptvals, m_int_order, dt, n, F0, el, F);

            if (dx.norm() < m_dl_max){
                // Frozen time, otherwise: locate(x+dx, t+dt, pos)
                is_inside = intp.locate(x + dx, t, pos);
                tried = true;
            }
            else {
                #pragma omp critical
                std::cout << "Step too long (dl=" << dx.norm() << "), consider doing something smart!" << std::endl;
            }
        }
        // count things
        if (is_inside){
            ++n_accepted_loc;
            ps.set_x(i, x + dx);
            ps.set_t_loc(i, ps.t_loc(i) + dt);
            ps.set_cell_id(i, pos.id);
            if constexpr (E == TransportElement::Vector){
                const double len = el.norm();
                ps.set_rhohat(i, el/len);
                ps.set_w(i, ps.w(i) + log(len));
            }
            if constexpr (E == TransportElement::Tensor)
                ps.advance_frame(i, F);
        }
        // An end outside: pieces
        else if (tried && spatial_pieces<E>(intp, ps, i, t, ds, m_int_order, m_u_min, ptvals, left))
            ++n_accepted_loc;
        else {
            outside_nodes_loc.push_back(i);
            ++n_declined_loc;
            n_left_loc += left;
        }
    }
    #pragma omp critical
    {
        outside_nodes.insert(outside_nodes.end(), outside_nodes_loc.begin(), outside_nodes_loc.end());
        n_accepted += n_accepted_loc;
        n_declined += n_declined_loc;
        n_left += n_left_loc;
    }
    }
    // Sort: threads merge out of order
    std::sort(outside_nodes.begin(), outside_nodes.end());
    return outside_nodes;
}

template<typename Interp>
PARTRAC_HOT_LOOP
void ParticleSet::update_fields(Interp& intp, const double t, const OutputFields& output_fields){
  const bool do_rho = output_fields.rho;
  const bool do_p = output_fields.p;
  const bool do_J = has_J && output_fields.J;
  // phi always updated (vector statistics)
  const bool do_phi = has_phi;
  const bool do_cell_type = has_cell_type;
  const bool do_S = element == TransportElement::Vector;
  const bool do_S3 = element == TransportElement::Tensor;

  #pragma omp parallel for
  for (Uint irw=0; irw < N(); ++irw){
    CellPos pos;
    pos.id = get_cell_id(irw);
    PointValues ptvals(intp.get_U0());
    // Outside: zero fields
    if (intp.locate(x_rw[irw], t, pos))
      intp.evaluate(x_rw[irw], t, pos, ptvals);
    // Stretching rate
    if (do_S){
      const Matrix3d J = ptvals.get_J();
      S_rw[irw] = rhohat_rw[irw].transpose() * J * rhohat_rw[irw];
    }
    // Stretching rates along the frame
    if (do_S3){
      Matrix3d Q;
      Vector3d s, U;
      settled(irw, Q, s, U);
      S3_rw[irw] = (Q.transpose() * ptvals.get_J() * Q).diagonal();
    }
    // Always, for the statistics
    u_rw[irw] = ptvals.get_u();
    if (do_rho){
      rho_rw[irw] = ptvals.get_rho();
    }
    if (do_p){
      p_rw[irw] = ptvals.get_p();
    }
    if (do_J)
      J_rw[irw] = ptvals.get_J();
    if (do_phi)
      phi_rw[irw] = ptvals.get_phi();
    if (do_cell_type)
      cell_type_rw[irw] = ptvals.get_cell_type();
  }
}

#endif
