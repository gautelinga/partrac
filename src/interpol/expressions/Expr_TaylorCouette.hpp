#include "Expr.hpp"

#ifndef __EXPR_TAYLORCOUETTE_HPP
#define __EXPR_TAYLORCOUETTE_HPP

// Martinez-Ruiz et al., J. Fluid Mech. 837 (2018), eq. (3.2); lengths in inner radii, time in 1/Omega
class Expr_TaylorCouette final : public Expr {
public:
  Expr_TaylorCouette(const partrac::Params& expr_params) : Expr(expr_params) {
    R = expr_params.get<double>("R");
    H = expr_params.get<double>("H");
    K = expr_params.get<double>("K");
    a = expr_params.get<double>("a");
    c1 = expr_params.get<double>("c1");
    c2 = expr_params.get<double>("c2");
    c3 = expr_params.get<double>("c3");
    rho_inf = expr_params.get<double>("rho");
    p_inf = expr_params.get<double>("p_inf");
    x0 = {expr_params.get<double>("x0"),
          expr_params.get<double>("y0"),
          expr_params.get<double>("z0")};
    kr = M_PI/(R-1);
    kz = 2*M_PI/H;
  };
  // Keys this expression reads
  static void add_params(partrac::Schema& s) {
    s.require<double>("R", "outer radius, in inner radii");
    s.require<double>("H", "height");
    s.require<double>("K", "roll amplitude");
    s.require<double>("a", "second harmonic of the rolls");
    s.require<double>("c1", "decay of the azimuthal velocity from the inner cylinder");
    s.require<double>("c2", "end-plate decay");
    s.require<double>("c3", "azimuthal modulation by the rolls");
    s.require<double>("p_inf", "pressure");
    s.require<double>("rho", "density");
    s.require<double>("x0", "centre x");
    s.require<double>("y0", "centre y");
    s.require<double>("z0", "centre z");
  };
  // Distance to the nearest of cylinders and plates, positive in the fluid
  bool has_wall() const { return true; };
  double sdf(const Vector3d &x) const {
    const Vector3d r = x-x0;
    const double s = std::hypot(r[0], r[1]);
    return std::min({s - 1., R - s, H/2 - std::abs(r[2])});
  };
  Vector3d sdf_grad(const Vector3d &x) const {
    const Vector3d r = x-x0;
    const double s = std::hypot(r[0], r[1]);
    const double inner = s - 1., outer = R - s, plate = H/2 - std::abs(r[2]);
    if (plate <= inner && plate <= outer) return {0., 0., r[2] > 0. ? -1. : 1.};
    const Vector3d e = s > 0. ? Vector3d(r[0]/s, r[1]/s, 0.) : Vector3d::Zero();
    return inner <= outer ? e : Vector3d(-e);
  };
  bool inside(const Vector3d &x, const double t __attribute__((unused))) {
    Vector3d r = x-x0;
    double s = sqrt(r[0]*r[0]+r[1]*r[1]);
    return s >= 1. && s <= R && 2*std::abs(r[2]) <= H;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused)), PointValues& ptvals) {
    compute(x, ptvals.U, ptvals.gradU, ptvals.convGradU);
    ptvals.P = p_inf;
    ptvals.Rho = rho_inf;
  };
private:
  // No member writes: called inside omp for
  void compute(const Vector3d &x, Vector3d& U_, Matrix3d& gradU_, Matrix3d& convGradU_) const {
    Vector3d d = x-x0;
    double s = sqrt(d[0]*d[0]+d[1]*d[1]);
    double z = d[2];

    // Clamp to the walls
    double dr = std::max(s-1., 0.);
    double t1 = tanh(c2*std::max(H/2-z, 0.));
    double t2 = tanh(c2*std::max(H/2+z, 0.));
    double co = 1./t1 + 1./t2;                       // infinite at either plate
    double eE = dr > 0. ? exp(-c1*dr*co) : 1.;

    double q = kr*dr, p = kz*z;
    double cq = cos(q), sq = sin(q);
    double F = cos(p) + a*cos(2*p);                  // dG/dz / kz
    double G = sin(p) + a*sin(2*p)/2;
    double Fz = -kz*(sin(p) + 2*a*sin(2*p));
    double W = 1. - 4*z*z/(H*H);
    double A = K*kz, B = K*kr;

    double u_r = A*sq*F/s;
    double u_z = -B*cq*G/s;
    double u_t = eE + c3*sq*W;

    double dur_dr = A*F*(kr*cq/s - sq/(s*s));
    double dur_dz = A*sq*Fz/s;
    double duz_dr = B*G*(kr*sq/s + cq/(s*s));
    double duz_dz = -B*kz*cq*F/s;
    double dut_dr = c3*kr*cq*W;
    double dut_dz = -8*c3*sq*z/(H*H);
    // eE's derivatives, in s and z
    const bool graded = eE > 0. && std::isfinite(co);
    double Er = 0., Ez = 0., co_z = 0.;
    if (graded){
      co_z = c2*((1-t1*t1)/(t1*t1) - (1-t2*t2)/(t2*t2));
      Er = -c1*co*eE;
      Ez = -c1*dr*co_z*eE;
      dut_dr += Er;
      dut_dz += Ez;
    }

    double c = d[0]/s, sn = d[1]/s;
    U_ = {u_r*c - u_t*sn, u_r*sn + u_t*c, u_z};
    // gradU_(i, j) is dU_i/dx_j
    gradU_ <<
      c*c*dur_dr - c*sn*dut_dr + (sn*sn*u_r + sn*c*u_t)/s,
      sn*c*dur_dr - sn*sn*dut_dr - (sn*c*u_r + c*c*u_t)/s,
      c*dur_dz - sn*dut_dz,
      c*sn*dur_dr + c*c*dut_dr - (sn*c*u_r - sn*sn*u_t)/s,
      sn*sn*dur_dr + sn*c*dut_dr + (c*c*u_r - sn*c*u_t)/s,
      sn*dur_dz + c*dut_dz,
      c*duz_dr,
      sn*duz_dr,
      duz_dz;
    if (convective)
      convective_terms({d[0], d[1], d[2], dr, t1, t2, co, co_z, eE, Er, Ez, graded, sq, cq, sin(p), cos(p)}, convGradU_);
  };
  // What convective_terms reads of the point
  struct At {
    double r_0, r_1, r_2;                // from the centre
    double dr, t1, t2, co, co_z;         // eE's decay
    double eE, Er, Ez;
    bool graded;
    double sq, cq, sp, cp;               // sin and cos of kr*dr and kz*z
  };
  // (U.grad) gradU, out of the hot eval
  __attribute__((noinline))
  void convective_terms(const At& at, Matrix3d& out) const {
    double Err = 0., Erz = 0., Ezz = 0.;
    if (at.graded){
      const double t1 = at.t1, t2 = at.t2, co = at.co, co_z = at.co_z, dr = at.dr;
      const double co_zz = 2*c2*c2*((1-t1*t1)/(t1*t1*t1) + (1-t2*t2)/(t2*t2*t2));
      Err = c1*c1*co*co*at.eE;
      Erz = (-c1*co_z + c1*c1*dr*co*co_z)*at.eE;
      Ezz = (-c1*dr*co_zz + c1*c1*dr*dr*co_z*co_z)*at.eE;
    }
    // (U.grad) gradU, generated by Sympy (python/codegen/conv_gradu.py)
    const double w0 = ((at.r_0)*(at.r_0));
    const double w1 = ((at.r_1)*(at.r_1));
    const double w2 = w0 + w1;
    const double w3 = (1.0/(w2));
    const double w4 = K*at.sq;
    const double w5 = w3*w4;
    const double w6 = ((kz)*(kz));
    const double w7 = a*at.cp;
    const double w8 = 4*w7 + 1;
    const double w9 = w6*w8;
    const double w10 = at.sp*w9;
    const double w11 = w10*w5;
    const double w12 = (1.0/((H)*(H)));
    const double w13 = c3*w12;
    const double w14 = w3*(Erz - 8*at.cq*kr*at.r_2*w13);
    const double w15 = at.r_0*at.r_1;
    const double w16 = w14*w15;
    const double w17 = 8*at.sq;
    const double w18 = w13*w17;
    const double w19 = at.Ez - at.r_2*w18;
    const double w20 = (1.0/((w2)*sqrt(w2)));
    const double w21 = at.r_1*w20;
    const double w22 = at.r_0*w19*w21;
    const double w23 = (1.0/((w2)*(w2)));
    const double w24 = at.cq*w20;
    const double w25 = K*kr;
    const double w26 = at.sp*w25;
    const double w27 = w24*w26*w9;
    const double w28 = 2*K*at.sp*at.sq*w0*w23*w6*w8 - w0*w27 - w11 - w16 + w22;
    const double w29 = (1.0/sqrt(w2));
    const double w30 = w7 + 1;
    const double w31 = w26*w30;
    const double w32 = at.cq*w29*w31;
    const double w33 = a*(2*((at.cp)*(at.cp)) - 1);
    const double w34 = at.cp + w33;
    const double w35 = at.r_1*w34;
    const double w36 = kz*w5;
    const double w37 = -4*((at.r_2)*(at.r_2))*w12 + 1;
    const double w38 = c3*w37;
    const double w39 = at.sq*w38;
    const double w40 = at.eE + w39;
    const double w41 = at.r_0*w40;
    const double w42 = w29*w41 + w35*w36;
    const double w43 = kr*w38;
    const double w44 = at.Er + at.cq*w43;
    const double w45 = w3*w44;
    const double w46 = at.r_0*w45;
    const double w47 = w20*w41;
    const double w48 = at.Er*w20;
    const double w49 = w24*w43;
    const double w50 = ((kr)*(kr));
    const double w51 = w3*w39*w50;
    const double w52 = Err*w3 - w48 - w49 - w51;
    const double w53 = at.r_0*w29;
    const double w54 = 2*w23;
    const double w55 = w44*w54;
    const double w56 = at.r_0*w1;
    const double w57 = w55*w56;
    const double w58 = (1.0/((w2)*(w2)*sqrt(w2)));
    const double w59 = 3*w58;
    const double w60 = w1*w41*w59;
    const double w61 = kz*w4;
    const double w62 = w54*w61;
    const double w63 = kz*w25;
    const double w64 = w24*w63;
    const double w65 = w35*w64;
    const double w66 = K*kz*w17/((w2)*(w2)*(w2));
    const double w67 = at.r_1*w0;
    const double w68 = w34*w67;
    const double w69 = w23*w61;
    const double w70 = w50*w69;
    const double w71 = at.cq*w63;
    const double w72 = 5*w58*w71;
    const double w73 = w35*w62 - w65 - w66*w68 + w68*w70 + w68*w72;
    const double w74 = w1*w52*w53 + w46 - w47 - w57 + w60 + w73;
    const double w75 = at.r_0*w34;
    const double w76 = at.r_1*w29;
    const double w77 = w36*w75 - w40*w76;
    const double w78 = 6*w69;
    const double w79 = w21*w40;
    const double w80 = w64*w75;
    const double w81 = ((at.r_0)*(at.r_0)*(at.r_0));
    const double w82 = w34*w81;
    const double w83 = w40*w59;
    const double w84 = w67*w83;
    const double w85 = w55*w67;
    const double w86 = -at.Er*w29 - c3*at.cq*kr*w29*w37;
    const double w87 = Err*w0*w3 - w0*w48 - w0*w49 - w0*w51 - w86;
    const double w88 = w19*w29;
    const double w89 = w19*w20;
    const double w90 = w10*w15*w4*w54;
    const double w91 = w15*w27;
    const double w92 = w1*w14 - w1*w89 + w88 - w90 + w91;
    const double w93 = Err*w1*w3 - w1*w48 - w1*w49 - w1*w51 - w86;
    const double w94 = ((at.r_1)*(at.r_1)*(at.r_1));
    const double w95 = at.r_1*w45;
    const double w96 = w34*w56;
    const double w97 = w62*w75 - w66*w96 + w70*w96 + w72*w96 - w80;
    const double w98 = Ezz - w18;
    const double w99 = ((kz)*(kz)*(kz))*w5*(at.cp + 4*w33);
    const double w100 = w0*w14 - w0*w89 + w88 + w90 - w91;
    const double w101 = -w0*w52*w76 + w79 - w84 + w85 - w95 + w97;
    const double w102 = -2*K*at.sp*at.sq*w1*w23*w6*w8 + w1*w27 + w11 - w16 + w22;
    const double w103 = w34*w94;
    const double w104 = kr*at.sq;
    const double w105 = w104*w3 + w24;
    const double w106 = w105*w34;
    const double w107 = w106*w71;
    const double w108 = at.cq*w59;
    const double w109 = 3*w104*w23;
    const double w110 = w24*w50;
    const double w111 = w15*(w108 + w109 - w110);
    out <<
      -w28*w32 - w42*w74 - w77*(-w66*w82 + w70*w82 + w72*w82 + w75*w78 + w76*w87 - w79 - 3*w80 + w84 - w85),
      K*at.cq*kr*at.sp*w29*w30*w92 - w42*(-w55*w94 + w76*w93 - 3*w79 + w83*w94 + 2*w95 + w97) - w74*w77,
      w28*w77 + w32*(at.r_0*w99 + w76*w98) - w42*w92,
      -w100*w32 - w101*w42 - w77*(-2*w46 + 3*w47 - w53*w87 + w55*w81 + w73 - w81*w83),
      -w101*w77 + w102*w32 - w42*(-w103*w66 + w103*w70 + w103*w72 + w35*w78 + w47 - w53*w93 + w57 - w60 - 3*w65),
      w100*w77 - w102*w42 + w32*(at.r_1*w99 - w53*w98),
      w31*(-w107*w53 - w111*w42 + w77*(-w0*w108 - w0*w109 + w0*w110 + w105)),
      w31*(-w107*w76 - w111*w77 + w42*(-w1*w108 - w1*w109 + w1*w110 + w105)),
      w63*(-((at.cq)*(at.cq))*((at.sp)*(at.sp))*w3*w30*w63*w8 + at.r_0*w106*w77 + at.r_1*w106*w42);
  };
  double R, H, K, a, c1, c2, c3;
  double kr, kz;
  double p_inf;
  double rho_inf;
  Vector3d x0;
};

#endif
