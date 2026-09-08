"""The diffusive strip and sheet methods, on their own benchmark flows.

Meunier & Villermaux, J. Fluid Mech. 662 (2010) 134-172, track a material line
through a two-dimensional flow and read the striation thickness off its
elongation rho = s/s0.  Martinez-Ruiz, Meunier, Favier, Duchemin & Villermaux,
J. Fluid Mech. 837 (2018) 230-257, do the same in three dimensions with a
material sheet, where rho is the areal elongation A/A0.  Both rest on the
advected mesh reproducing the kinematics exactly, which is what is checked here.

Between them the two papers use four flows, all driven here: the point vortex of
the 2010 paper (section 2.6), its random sine flow (section 3), the Batchelor
vortex the 2018 paper validates against, and that paper's Taylor-Couette cell.
The last is an experiment, but the field handed to the DSM is the analytic fit
of its equation (3.2), which is `Expr_TaylorCouette`.

The random sine flow of the first paper is `data_example/sine_flow`.  Inside one
`tau` interval the driving coordinate is frozen, so the shear map is exact for
any timestep -- but only if `dt` divides `tau` in binary, otherwise floor(t/tau)
flips the shear a step early.  Hence dt = 0.0625, not 0.05.
"""

import os
import re
import subprocess

import numpy as np
import pytest

from paths import REPO, app

PARTRAC = app("partrac")
SINE = os.path.join(REPO, "data_example", "sine_flow", "expr_params.dat")
ABC = os.path.join(REPO, "data_example", "abc_flow", "expr_params.dat")
BATCHELOR = os.path.join(REPO, "data_example", "batchelor_vortex",
                         "expr_params.dat")
TAYLOR_COUETTE = os.path.join(REPO, "data_example", "taylor_couette",
                              "expr_params.dat")

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")

SINE_BASE = ("mode=analytic x0=0.5 y0=0.5 z0=0.5 Nrw=200 Nrw_max=1000000 "
             "ds_min=1e-9 coarsen=false Dm=0 int_order=1 dt=0.0625 "
             "dump_intv=1e9 checkpoint_intv=1e9 random=false seed=1").split()

ABC_BASE = ("mode=analytic x0=6.3 y0=5.7 z0=7.1 Nrw=100 Nrw_max=200000 "
            "ds_max=1e9 ds_min=1e-9 refine=false coarsen=false Dm=0 "
            "int_order=2 dt=0.002 dump_intv=1e9 checkpoint_intv=1e9 "
            "random=false seed=1").split()

# the vortex sits at (5, 5, 5); a strip along x from 0.6 to 1.8 out from the
# axis is the scalar of Meunier & Villermaux figure 5
R_IN, R_OUT = 0.6, 1.8
# The cell is 1 <= r <= 2.5, |z| <= 1.5.  Ekman pumping off the end plates makes
# two counter-rotating tori whose centres are at r = 1.75, z = +-0.595; z = 0 is
# the separatrix between them, where u_z vanishes identically, so a blob must be
# injected off it or it never goes round.  The paper injects at r = 2.
TC_BASE = ("mode=analytic x0=2.0 y0=0.0 z0=0.6 La=0.02 Lb=0.02 Nrw=100 "
           "Nrw_max=400000 ds_max=1e9 ds_min=1e-12 refine=false coarsen=false "
           "Dm=0 int_order=2 dt=0.01 T=60.0 stat_intv=0.5 dump_intv=1e9 "
           "checkpoint_intv=1e9 random=false seed=1").split()

VORTEX_BASE = ("mode=analytic y0=5.0 z0=5.0 Nrw=4000 Nrw_max=8000000 "
               "ds_min=1e-9 coarsen=false Dm=0 int_order=2 dt=0.002 "
               "x0=%g La=%g dump_intv=1e9 checkpoint_intv=1e9 random=false "
               "seed=1" % (5.0 + (R_IN + R_OUT) / 2, R_OUT - R_IN)).split()


def expr_params(path):
    d = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def run(tmp_path, example, base, extra, expr=None):
    """One partrac run; returns every tdata row as a dict of column arrays.

    `expr` overrides keys of the copied expr_params.dat: u0, R1 and the rest
    belong to the expression, and partrac rejects them on its own command line.
    """
    tmp_path.mkdir(parents=True, exist_ok=True)
    text = open(example).read()
    for key, value in (expr or {}).items():
        text, n = re.subn(r"(?m)^%s=.*$" % key, "%s=%.12g" % (key, value), text)
        assert n == 1, "%s is not a key of %s" % (key, example)
    (tmp_path / "expr_params.dat").write_text(text)
    keys = {a.split("=")[0] for a in extra}
    argv = [a for a in base if a.split("=")[0] not in keys] + extra
    r = subprocess.run([PARTRAC, str(tmp_path / "expr_params.dat")] + argv,
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr
    f = list(tmp_path.rglob("tdata_from_t*.dat"))
    assert len(f) == 1
    rows = [l for l in f[0].read_text().splitlines() if l.strip()]
    head = [h.strip() for h in rows[0].lstrip("# ").split("\t") if h.strip()]
    data = np.array([[float(v) for v in l.split("\t") if v.strip()]
                     for l in rows[1:]])
    assert data.shape[1] == len(head), "tdata columns do not match its header"
    return {h: data[:, i] for i, h in enumerate(head)}


def elongation(st):
    """rho = s/s0 for a strip, A/A0 for a sheet."""
    num, den = ("A", "A0") if "A" in st else ("s", "s0")
    return st[num] / st[den]


def sine_flow_map(n, La, T):
    """Length of a line initially along x, under the exact sine-flow map."""
    return sine_flow_lengths(n, La, T)[1][-1]


def sine_flow_lengths(n, La, T):
    """That length at the end of every half period."""
    p = expr_params(SINE)
    u, tau, L = float(p["u_inf"]), float(p["tau"]), float(p["Lx"])
    chi = [float(c) for c in p["chi"].split(",")]
    flowdir = [int(c) for c in p["flowdir"].split(",")]
    depdir = [int(c) for c in p["depdir"].split(",")]
    x = np.full((n, 3), 0.5)
    x[:, 0] = np.linspace(0.5 - La / 2, 0.5 + La / 2, n)
    out = []
    for i in range(int(round(T / tau))):
        j, k = flowdir[i % len(flowdir)], depdir[i % len(depdir)]
        x[:, j] += u * np.sin(2 * np.pi * x[:, k] / L + chi[i % len(chi)]) * tau
        out.append(((i + 1) * tau,
                    np.linalg.norm(np.diff(x, axis=0), axis=1).sum()))
    return np.array([o[0] for o in out]), np.array([o[1] for o in out])


def vortex_length(u0, R1, R2, q, T, sheet=False):
    """Exact size of a radial material line in a steady axisymmetric vortex.

    Nothing moves radially, so a point at radius s only turns through
    omega(s)*T and slides along the axis by u_z(s)*T; the elongation follows
    from differentiating that in s.  A sheet ruled along z picks up no axial
    term, since crossing dX/ds with z_hat drops it.
    """
    s = np.linspace(R_IN, R_OUT, 400001)
    a = s ** 2 / R1 ** 2
    s_domega = 2 * u0 * R1 * ((1 + a) * np.exp(-a) - 1) / s ** 2
    duz = -2 * s * q * u0 * np.exp(-s ** 2 / R2 ** 2) / R2 ** 2
    rho = 1. + (s_domega * T) ** 2
    if not sheet:
        rho = rho + (duz * T) ** 2
    return np.trapz(np.sqrt(rho), s)


def abc_areal_elongation(x_start, normal, T, dt=1e-3):
    """|cof(F) n| along one ABC trajectory, from dF/dt = grad(u) F by RK4."""
    p = expr_params(ABC)
    A, B, C = float(p["A"]), float(p["B"]), float(p["C"])
    k = 2 * np.pi / float(p["L"])
    x0 = np.array([float(p["x0"]), float(p["y0"]), float(p["z0"])])

    def rhs(y):
        r = y[:3] - x0
        sx, sy, sz = B * np.sin(k * r[0]), C * np.sin(k * r[1]), A * np.sin(k * r[2])
        cx, cy, cz = B * np.cos(k * r[0]), C * np.cos(k * r[1]), A * np.cos(k * r[2])
        gradu = np.array([[0., -k * sy, k * cz],
                          [k * cx, 0., -k * sz],
                          [-k * sx, k * cy, 0.]])
        return np.concatenate([[sz + cy, sx + cz, sy + cx],
                               (gradu @ y[3:].reshape(3, 3)).ravel()])

    y = np.concatenate([x_start, np.eye(3).ravel()])
    for _ in range(int(round(T / dt))):
        a = rhs(y)
        b = rhs(y + dt / 2 * a)
        c = rhs(y + dt / 2 * b)
        d = rhs(y + dt * c)
        y += dt / 6 * (a + 2 * b + 2 * c + d)
    F = y[3:].reshape(3, 3)
    assert abs(np.linalg.det(F) - 1) < 1e-9, "the ABC flow is incompressible"
    return np.linalg.norm(np.linalg.det(F) * np.linalg.inv(F).T @ normal)


# --- Meunier & Villermaux (2010): the strip in the random sine flow ---------

@needs_partrac
def test_a_strip_follows_the_exact_sine_flow_map(tmp_path):
    st = run(tmp_path, SINE, SINE_BASE,
             ["init_mode=strip_x", "La=0.5", "Nrw=2000", "T=1.0",
              "stat_intv=1.0", "ds_max=1e9", "refine=false"])
    assert elongation(st)[-1] == pytest.approx(sine_flow_map(200001, 0.5, 1.0) / 0.5,
                                               rel=1e-4)


@needs_partrac
def test_a_refined_strip_survives_a_large_elongation(tmp_path):
    # rho ~ 65 after five time units; without refinement the polyline would cut
    # the corners off the folds
    st = run(tmp_path, SINE, SINE_BASE,
             ["init_mode=strip_x", "La=0.2", "T=5.0", "stat_intv=5.0",
              "ds_max=0.002", "refine=true", "refine_intv=0.0625"])
    assert st["Nrw"][-1] > 5000                       # it really refined
    assert st["s0"][-1] == pytest.approx(0.2, rel=1e-12)   # and s0 is conserved
    assert elongation(st)[-1] == pytest.approx(sine_flow_map(400001, 0.2, 5.0) / 0.2,
                                               rel=1e-3)


@needs_partrac
def test_the_elongation_is_log_normal(tmp_path):
    # a multiplicative process: <log rho> and Var(log rho) both grow linearly
    st = run(tmp_path, SINE, SINE_BASE,
             ["init_mode=strip_x", "La=0.2", "T=6.0", "stat_intv=0.5",
              "ds_max=0.002", "refine=true", "refine_intv=0.0625"])
    late = st["t"] > 1.0
    for key, r2 in (("logelong_wmean", 0.99), ("logelong_wvar", 0.90)):
        slope = np.polyfit(st["t"][late], st[key][late], 1)[0]
        assert slope > 0
        assert np.corrcoef(st["t"][late], st[key][late])[0, 1] ** 2 > r2
    lyapunov = np.polyfit(st["t"][late], st["logelong_wmean"][late], 1)[0]
    assert 0.5 < lyapunov < 2.0


@needs_partrac
def test_a_two_dimensional_flow_preserves_the_area_of_a_patch(tmp_path):
    # the striation thickness is s0/rho only because the map is area-preserving
    st = run(tmp_path, SINE, SINE_BASE,
             ["init_mode=sheet_xy", "La=0.2", "Lb=0.2", "ds_init=0.01",
              "T=1.0", "stat_intv=1.0", "ds_max=1e9", "refine=false"])
    assert elongation(st)[-1] == pytest.approx(1.0, abs=1e-3)


@needs_partrac
def test_the_sine_flow_stretches_at_the_published_rate(tmp_path):
    # section 4.1 reads gamma = 0.91 +- 2% off figure 9.  The exact map for the
    # phase table of their table 1 gives 0.87, and partrac agrees with the map,
    # so the paper's figure is the odd one out; the band below holds both.
    st = run(tmp_path, SINE, SINE_BASE,
             ["init_mode=strip_x", "La=1.0", "Nrw=1000", "T=7.0",
              "stat_intv=0.25", "ds_max=0.02", "refine=true",
              "refine_intv=0.0625"])
    t, length = sine_flow_lengths(200001, 1.0, 7.0)
    fit = lambda x, y: np.polyfit(x[x >= 1.0], y[x >= 1.0], 1)[0]
    gamma = fit(st["t"], np.log(st["s"] / st["s0"]))
    assert gamma == pytest.approx(fit(t, np.log(length / 1.0)), rel=1e-2)
    assert elongation(st)[-1] == pytest.approx(length[-1] / 1.0, rel=5e-3)
    assert 0.85 < gamma < 0.95


@needs_partrac
def test_a_radial_strip_winds_into_the_point_vortex_spiral(tmp_path):
    # section 2.6 and figure 5: circulation 14.2, strip from 0.6 to 1.8 out from
    # the axis.  A Batchelor vortex whose core sits well inside R_IN is a point
    # vortex there, with circulation 2 pi u0 R1.  Their t = 10 s costs 20 s to
    # run for the same 5e-5; two units already wind the spiral eight-fold.
    circulation, R1, T = 14.2, 0.1, 2.0
    u0 = circulation / (2 * np.pi * R1)
    st = run(tmp_path, BATCHELOR, VORTEX_BASE,
             ["init_mode=strip_x", "T=%g" % T, "stat_intv=%g" % T,
              "ds_max=0.0005", "refine=true", "refine_intv=0.01"],
             expr=dict(u0=u0, R1=R1, R2=1.0, q=0.0))
    assert st["s"][-1] == pytest.approx(vortex_length(u0, R1, 1.0, 0.0, T),
                                        rel=1e-3)
    assert elongation(st)[-1] > 8


# --- Martinez-Ruiz et al. (2018): the sheet ---------------------------------

@needs_partrac
def test_the_sheet_reduces_to_the_strip_in_a_two_dimensional_flow(tmp_path):
    # the sine flow leaves z alone, so an x-z sheet is the x-strip extruded and
    # A/A0 must converge on s/s0
    exact = sine_flow_map(200001, 0.5, 1.0) / 0.5
    args = ["init_mode=sheet_xz", "La=0.5", "Lb=0.5", "T=1.0",
            "stat_intv=1.0", "ds_max=1e9", "refine=false"]
    err = [abs(elongation(run(tmp_path / d, SINE, SINE_BASE,
                              args + ["ds_init=%g" % ds]))[-1] / exact - 1)
           for d, ds in (("coarse", 0.02), ("fine", 0.01))]
    assert err[1] < 5e-4
    assert err[0] / err[1] > 3            # second order in the node spacing


@needs_partrac
def test_sheet_area_follows_the_deformation_gradient(tmp_path):
    # the areal elongation of a shrinking patch tends to |cof(F) n|, which is
    # what the diffusive sheet method integrates along each trajectory
    exact = abc_areal_elongation(np.array([6.3, 5.7, 7.1]), np.array([0., 0., 1.]),
                                 T=10.0)
    args = ["init_mode=sheet_xy", "T=10.0", "stat_intv=10.0"]
    err = [abs(elongation(run(tmp_path / d, ABC, ABC_BASE, args
                              + ["La=%g" % La, "Lb=%g" % La,
                                 "ds_init=%g" % (La / 12)]))[-1] / exact - 1)
           for d, La in (("coarse", 0.1), ("fine", 0.05))]
    assert err[1] < 1.5e-3
    assert err[0] / err[1] > 3            # second order in the patch size


@needs_partrac
def test_a_radial_strip_in_the_batchelor_vortex_matches_the_closed_form(tmp_path):
    # the flow the paper validates against.  q = 1 turns on the axial jet, so
    # the strip leaves its plane and the third direction is exercised
    u0, R1, R2, q, T = 2.0, 0.5, 1.0, 1.0, 1.0
    st = run(tmp_path, BATCHELOR, VORTEX_BASE,
             ["init_mode=strip_x", "Nrw=2000", "T=%g" % T, "stat_intv=%g" % T,
              "ds_max=0.0002", "refine=true", "refine_intv=0.005", "dt=0.0005"],
             expr=dict(u0=u0, R1=R1, R2=R2, q=q))
    assert st["s"][-1] == pytest.approx(vortex_length(u0, R1, R2, q, T), rel=1e-4)


@needs_partrac
def test_a_radial_sheet_in_the_batchelor_vortex_matches_the_closed_form(tmp_path):
    u0, R1, R2, q, T, Lb = 2.0, 0.5, 1.0, 1.0, 1.0, 0.4
    st = run(tmp_path, BATCHELOR, VORTEX_BASE,
             ["init_mode=sheet_xz", "Lb=%g" % Lb, "ds_init=0.01", "Nrw=100",
              "Nrw_max=2000000", "T=%g" % T, "stat_intv=%g" % T,
              "ds_max=1e9", "refine=false", "dt=0.001"],
             expr=dict(u0=u0, R1=R1, R2=R2, q=q))
    exact = vortex_length(u0, R1, R2, q, T, sheet=True) * Lb
    assert st["A"][-1] == pytest.approx(exact, rel=1e-4)
    # the axial jet stretches the line but not the sheet, whose ruling stays z
    assert exact < vortex_length(u0, R1, R2, q, T) * Lb


# --- Martinez-Ruiz et al. (2018) section 3: the Taylor-Couette cell -----------
#
# Stationary and axisymmetric, so it is the opposite of the sine flow: there is
# no succession of stretchings and foldings, and section 3.5 finds the sheet
# growing linearly rather than exponentially, with a bounded p(rho).

@needs_partrac
def test_only_the_radial_direction_stretches_in_the_cell(tmp_path):
    # figure 13(a).  Two points at the same radius sit on one stream torus and
    # share a mean angular velocity, so the segments joining them stay bounded;
    # a radial segment spans two tori, which is a shear and grows linearly.
    rho = {}
    for d in "xyz":
        st = run(tmp_path / d, TAYLOR_COUETTE, TC_BASE, ["init_mode=strip_%s" % d])
        rho[d] = elongation(st)
        t = st["t"]

    late = t >= 5.0
    half = rho["x"][np.argmin(abs(t - 30.0))]
    assert rho["x"][-1] > 20                              # radial
    assert np.corrcoef(t[late], rho["x"][late])[0, 1] ** 2 > 0.9
    assert rho["x"][-1] / half < 4                        # e^t would give 14
    for d in "yz":                                        # azimuthal, axial
        assert rho[d].min() > 0.5 and rho[d].max() < 2.0


@needs_partrac
def test_a_sheet_tangent_to_the_stream_torus_does_not_grow(tmp_path):
    # section 3.5: the triangle spanned by the azimuthal and axial segments is
    # tangent to the torus and is only modulated periodically, while any sheet
    # carrying the radial direction grows
    tangent = elongation(run(tmp_path / "yz", TAYLOR_COUETTE, TC_BASE,
                             ["init_mode=sheet_yz", "ds_init=0.002"]))
    radial = elongation(run(tmp_path / "xy", TAYLOR_COUETTE, TC_BASE,
                            ["init_mode=sheet_xy", "ds_init=0.002"]))
    assert tangent.min() > 0.5 and tangent.max() < 1.5
    assert radial[-1] > 20


@needs_partrac
def test_the_sheet_grows_linearly_and_stays_narrowly_distributed(tmp_path):
    # figure 14: the area is linear in time, not exponential, and p(rho) keeps a
    # bounded support -- the contrast the paper draws with the sine flow, whose
    # log-normal Var(log rho) is two orders of magnitude wider by this elongation
    st = run(tmp_path, TAYLOR_COUETTE, TC_BASE,
             ["init_mode=sheet_xy", "ds_init=0.002"])
    rho, t = elongation(st), st["t"]
    half = rho[np.argmin(abs(t - 30.0))]
    assert rho[-1] / half == pytest.approx(2.0, rel=0.1)   # linear; e^t gives 13
    assert np.corrcoef(t[t >= 5], rho[t >= 5])[0, 1] ** 2 > 0.95
    assert st["logelong_wmean"][-1] > 3.0
    assert st["logelong_wvar"][-1] < 0.05
