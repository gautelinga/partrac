"""The compressed time tau and the scalar reconstruction of diffusive strips and sheets.

Meunier & Villermaux (2010) reduce advection-diffusion across a strip of
striation thickness s = s0/rho to one variable, time in units of the current
diffusion time (their eq. 2.6),

    dtau/dt = D/s(t)^2,   so   tau = (D/s0^2) int_0^t rho^2 dt',

and the transverse profile follows from it (eq. 2.8),

    c(n, t) = c0 (1 + 4 tau)^(-1/2) exp[-(n/s)^2 / (1 + 4 tau)],

so the maximum concentration is c0/sqrt(1 + 4 tau). partrac integrates the bare
int rho^2 dt' per edge (per face for a sheet) and dumps it as `tau` beside
`dl`/`dl0`; D and s0 are applied afterwards, so a finished run can be rescaled
in s0.

Wherever the elongation has a closed form, so does tau. For a line laid along x
in plane Poiseuille the velocity is constant along each path, so an edge at x
has rho^2 = 1 + (a t)^2 with a = 3 u_inf x / R^2, and tau = t + a^2 t^3 / 3.
The Batchelor vortex and Taylor-Couette flow give further closed forms and
scaling laws.
"""

import os
import re

import numpy as np
import pytest

from dumps import all_dumps
from paths import REPO, app
from runs import run_app

PARTRAC = app("partrac")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
BATCHELOR = os.path.join(REPO, "data_example", "batchelor_vortex", "expr_params.dat")
TAYLOR_COUETTE = os.path.join(REPO, "data_example", "taylor_couette", "expr_params.dat")

U_INF, R = 1.0, 1.0

BASE = ("mode=analytic ds_min=1e-12 refine=false coarsen=false Dm=0 "
        "int_order=1 checkpoint_intv=1e9 random=false seed=1 "
        "integrate_tau=true tau_max=0").split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")


def run(tmp_path, example, extra, expr=None):
    """Run partrac with tau integration, optionally overriding expression keys; return time -> datasets as written."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    with open(example) as f:
        text = f.read()
    for key, value in (expr or {}).items():
        text, n = re.subn(r"(?m)^%s=.*$" % key, "%s=%.12g" % (key, value), text)
        assert n == 1, "%s is not a key of %s" % (key, example)
    (tmp_path / "expr_params.dat").write_text(text)
    run_app(PARTRAC, tmp_path / "expr_params.dat", BASE, extra)
    assert len(list(tmp_path.rglob("data_from_t*.h5"))) == 1
    return all_dumps(tmp_path, raw=True)


def c_max(tau, D, s0):
    """Meunier & Villermaux (2010) eq. (2.8) at n = 0, in units of c0."""
    return 1.0 / np.sqrt(1.0 + 4.0 * D * np.asarray(tau) / s0 ** 2)


def poiseuille_tau(x, t):
    """Exact int_0^t rho^2 dt' in plane Poiseuille for an edge at x of a line along x."""
    a = 3 * U_INF * np.asarray(x) / R ** 2
    return t + a ** 2 * t ** 3 / 3


def edge_midpoints_x(dump):
    """Initial x of each edge midpoint, which labels the edge for the whole run."""
    first = sorted(dump)[0]
    p = np.array(dump[first]["points"])
    e = np.array(dump[first]["edges"])
    return 0.5 * (p[e[:, 0], 0] + p[e[:, 1], 0])


@needs_partrac
def test_tau_matches_the_exact_integral(tmp_path):
    """In plane Poiseuille the dumped tau matches t + a^2 t^3/3 per edge at
    every dump, and is zero at t = 0. Every scalar concentration the user
    reconstructs from a run is a function of this value, so the maximum
    concentration c0/sqrt(1 + 4 D tau/s0^2) follows its closed form and
    decays monotonically, fastest where the strip stretches most."""
    T = 0.5
    d = run(tmp_path, POISEUILLE,
            ["init_mode=strip_x", "La=1.0", "x0=0", "y0=0", "z0=0", "Nrw=200",
             "Nrw_max=5000", "ds_max=1e9", "dt=0.0025", "T=%g" % T,
             "stat_intv=%g" % T, "dump_intv=0.05", "tau_intv=0.0025"])
    x = edge_midpoints_x(d)
    assert max(d) == T and len(d) == 11
    for t in sorted(d):
        tau = np.array(d[t]["tau"]).ravel()
        if t == 0:
            # nothing has diffused before the run begins
            assert tau == pytest.approx(0.0, abs=1e-15)
        else:
            assert tau == pytest.approx(poiseuille_tau(x, t), rel=1e-5), t


@needs_partrac
def test_tau_converges_at_second_order(tmp_path):
    """The tau error falls by 4 each time tau_intv halves, and its size is the
    trapezoid-rule remainder. This confirms tau is integrated by the trapezoid
    rule on tau_intv, so users can choose tau_intv from a known error bound."""
    # the integrand 1 + (a t)^2 is quadratic in t, so the trapezoid error is
    # exactly -a^2 h^2 T/6
    T, err = 0.5, []
    for intv in ("0.02", "0.01", "0.005"):
        d = run(tmp_path / intv, POISEUILLE,
                ["init_mode=strip_x", "La=1.0", "x0=0", "y0=0", "z0=0",
                 "Nrw=200", "Nrw_max=5000", "ds_max=1e9", "dt=0.0025",
                 "T=%g" % T, "stat_intv=%g" % T, "dump_intv=%g" % T,
                 "tau_intv=" + intv])
        tau = np.array(d[max(d)]["tau"]).ravel()
        exact = poiseuille_tau(edge_midpoints_x(d), T)
        err.append(np.max(np.abs(tau / exact - 1)))
    assert err[0] / err[1] == pytest.approx(4.0, rel=0.1)
    assert err[1] / err[2] == pytest.approx(4.0, rel=0.1)

    # the largest relative error is at the strip's end, x = 0.5
    a = 3 * U_INF * 0.5 / R ** 2
    predicted = (a ** 2 * 0.02 ** 2 * T / 6) / poiseuille_tau(0.5, T)
    assert err[0] == pytest.approx(predicted, rel=0.05)


@needs_partrac
def test_tau_max_drops_strips_that_have_finished_mixing(tmp_path):
    """With tau_max > 0, edges whose tau exceeds it are removed, and a lower
    tau_max removes more. This lets a run stop spending work on strips that are
    already mixed."""
    # at T = 0.5, tau ranges from 0.5 at x = 0 to 0.59375 at |x| = 0.5, so both
    # thresholds cut part of the strip
    T = 0.5
    args = ["init_mode=strip_x", "La=1.0", "x0=0", "y0=0", "z0=0", "Nrw=200",
            "Nrw_max=5000", "ds_max=1e9", "dt=0.0025", "T=%g" % T,
            "stat_intv=%g" % T, "dump_intv=%g" % T, "tau_intv=0.0025"]
    kept = {}
    for tau_max in ("0", "0.55", "0.52"):
        d = run(tmp_path / tau_max, POISEUILLE, args + ["tau_max=" + tau_max])
        tau = np.array(d[max(d)]["tau"]).ravel()
        kept[tau_max] = len(tau)
        if float(tau_max) > 0:
            assert tau.max() <= float(tau_max)
    assert kept["0"] > kept["0.55"] > kept["0.52"] > 0


@needs_partrac
def test_refinement_carries_tau_across_a_split(tmp_path):
    """A split edge inherits its parent's tau and rho_prev, so a refined strip
    still matches the closed form. If splitting reset tau, refined regions,
    which are the most stretched, would report the least mixing."""
    # x is invariant in this flow, so the closed form applies to child edges
    T = 0.5
    common = ["init_mode=strip_x", "La=1.0", "x0=0", "y0=0", "z0=0",
              "Nrw_max=20000", "dt=0.0025", "T=%g" % T, "stat_intv=%g" % T,
              "dump_intv=%g" % T, "tau_intv=0.0025"]
    err = {}
    for name, extra in (("plain", ["Nrw=200", "ds_max=1e9", "refine=false"]),
                        ("refined", ["Nrw=40", "ds_max=0.02", "refine=true",
                                     "refine_intv=0.01"])):
        d = run(tmp_path / name, POISEUILLE, common + extra)
        last = max(d)
        p = np.array(d[last]["points"])
        e = np.array(d[last]["edges"])
        tau = np.array(d[last]["tau"]).ravel()
        x = 0.5 * (p[e[:, 0], 0] + p[e[:, 1], 0])
        err[name] = np.max(np.abs(tau / poiseuille_tau(x, T) - 1))
        assert len(tau) > 0
    assert err["refined"] < 1e-4        # a reset tau would miss by far more
    assert err["refined"] < 2 * err["plain"]


@needs_partrac
def test_a_sheet_carries_the_same_tau_as_the_strip(tmp_path):
    """Faces of an x-y sheet integrate the same tau as the edges of an x strip.
    The flow depends on x alone, so the sheet is the strip extruded, and the
    sheet method must reproduce the strip's closed form."""
    T, ds_init = 0.5, 0.02
    d = run(tmp_path, POISEUILLE,
            ["init_mode=sheet_xy", "La=1.0", "Lb=0.4", "ds_init=%g" % ds_init,
             "x0=0", "y0=0", "z0=0", "Nrw=100", "Nrw_max=20000", "ds_max=1e9",
             "dt=0.0025", "T=%g" % T, "stat_intv=%g" % T, "dump_intv=%g" % T,
             "tau_intv=0.0025"])
    tau = np.array(d[max(d)]["tau"]).ravel()
    assert tau.min() == pytest.approx(poiseuille_tau(0.0, T), rel=1e-3)
    # the outermost face centre lies within ds_init of the sheet edge x = 0.5
    assert poiseuille_tau(0.5 - ds_init, T) < tau.max() < poiseuille_tau(0.5, T)


@needs_partrac
def test_tau_in_a_vortex_matches_the_exact_integral(tmp_path):
    """In a Batchelor vortex with axial jet, tau per edge matches
    T + ((s dOmega/ds)^2 + (du_z/ds)^2) T^3/3. This checks tau in a flow with
    two shear components rather than one."""
    # nothing moves radially, so each edge keeps its radius s from the axis at
    # x = 5 and the shear rates are constant along its path; q = 1 adds the jet
    u0, R1, R2, q, T = 2.0, 0.5, 1.0, 1.0, 1.0
    d = run(tmp_path, BATCHELOR,
            ["init_mode=strip_x", "x0=6.2", "y0=5.0", "z0=5.0", "La=1.2",
             "Nrw=400", "Nrw_max=20000", "ds_max=1e9", "dt=0.0005",
             "T=%g" % T, "stat_intv=%g" % T, "dump_intv=%g" % T,
             "tau_intv=0.0005", "int_order=2"],
            expr=dict(u0=u0, R1=R1, R2=R2, q=q))
    first, last = min(d), max(d)
    p0 = np.array(d[first]["points"])
    e = np.array(d[last]["edges"])
    s = np.abs(0.5 * (p0[e[:, 0], 0] + p0[e[:, 1], 0]) - 5.0)

    a = s ** 2 / R1 ** 2
    s_domega = 2 * u0 * R1 * ((1 + a) * np.exp(-a) - 1) / s ** 2
    duz = -2 * s * q * u0 * np.exp(-s ** 2 / R2 ** 2) / R2 ** 2
    exact = T + (s_domega ** 2 + duz ** 2) * T ** 3 / 3

    tau = np.array(d[last]["tau"]).ravel()
    assert tau == pytest.approx(exact, rel=1e-4)


@needs_partrac
def test_the_taylor_couette_concentration_decays_as_the_paper_predicts(tmp_path):
    """In Taylor-Couette flow (Martinez-Ruiz et al., section 3.5) a sheet tangent
    to the stream torus keeps its thickness, so tau ~ t and c_max ~ t^-1/2,
    while one carrying the radial direction thins as t^-1, so tau ~ t^3 and
    c_max ~ t^-3/2. The sheet method must reproduce these published decay laws."""
    # tangency holds for an infinitesimal triangle; a finite patch drifts off it
    # at later times (around t ~ 100) because of its size, so the fit stops at 60
    T = 60.0
    args = ["x0=2.0", "y0=0.0", "z0=0.6", "La=0.02", "Lb=0.02", "ds_init=0.002",
            "Nrw=100", "Nrw_max=400000", "ds_max=1e9", "int_order=2", "dt=0.02",
            "T=%g" % T, "stat_intv=%g" % T, "dump_intv=2.0", "tau_intv=0.02"]

    def tau_series(name, mode):
        d = run(tmp_path / name, TAYLOR_COUETTE, ["init_mode=" + mode] + args)
        ks = sorted(d)
        t = np.array([float(k) for k in ks])
        tau = np.array([np.mean(np.array(d[k]["tau"]).ravel()) for k in ks])
        return t[t > 0], tau[t > 0]

    t, tangent = tau_series("tangent", "sheet_yz")
    _, radial = tau_series("radial", "sheet_xy")

    assert np.corrcoef(t, tangent)[0, 1] ** 2 > 0.95          # tau ~ t
    late = t >= 20
    assert np.polyfit(np.log(t[late]), np.log(tangent[late]), 1)[0] < 1.5

    slope = np.polyfit(np.log(t[late]), np.log(radial[late]), 1)[0]
    assert slope == pytest.approx(3.0, rel=0.05)              # tau ~ t^3

    # so the concentrations decay as t^-1/2 and t^-3/2
    for series, expected in ((tangent, -0.5), (radial, -1.5)):
        c = c_max(series, D=1e4, s0=1.0)     # well past the mixing time
        got = np.polyfit(np.log(t[late]), np.log(c[late]), 1)[0]
        assert got == pytest.approx(expected, abs=0.15)


@needs_partrac
def test_tau_max_drops_faces_as_well_as_edges(tmp_path):
    """tau_max removes sheet faces whose tau exceeds it, as it does strip edges,
    and a lower tau_max removes more. The sheet method must honour the same cutoff
    as the strip method."""
    T = 0.5
    args = ["init_mode=sheet_xy", "La=1.0", "Lb=0.4", "ds_init=0.05", "x0=0",
            "y0=0", "z0=0", "Nrw=100", "Nrw_max=20000", "ds_max=1e9",
            "dt=0.0025", "T=%g" % T, "stat_intv=%g" % T, "dump_intv=%g" % T,
            "tau_intv=0.0025"]
    kept = {}
    for tau_max in ("0", "0.55", "0.52"):
        d = run(tmp_path / tau_max, POISEUILLE, args + ["tau_max=" + tau_max])
        tau = np.array(d[max(d)]["tau"]).ravel()
        kept[tau_max] = len(tau)
        if float(tau_max) > 0:
            assert tau.max() <= float(tau_max)
    assert kept["0"] > kept["0.55"] > kept["0.52"] > 0
