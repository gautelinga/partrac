"""The scalar reconstruction of the diffusive strip and sheet methods.

Meunier & Villermaux (2010) reduce the advection-diffusion problem across a
strip to one variable. With the striation thickness s = s0/rho, their equation
(2.6) counts time in units of the current diffusion time,

    dtau/dt = D/s(t)^2,   so   tau = (D/s0^2) int_0^t rho^2 dt',

and equation (2.8) gives the whole transverse profile from it,

    c(n, t) = c0 (1 + 4 tau)^(-1/2) exp[-(n/s)^2 / (1 + 4 tau)],

so the maximum concentration is c0/sqrt(1 + 4 tau). partrac integrates the bare
int rho^2 dt' per edge, or per face for a sheet, and writes it as `tau` beside
`dl`/`dl0` in the dump; D and s0 are applied afterwards, which is what lets a
finished run be rescaled in s0 (the paper makes the same point).

Wherever the elongation has a closed form, so does tau. For a line laid along x
in plane Poiseuille the velocity is constant along each path, so an edge at x
has rho^2 = 1 + (a t)^2 with a = 3 u_inf x / R^2, and tau = t + a^2 t^3 / 3.
"""

import os
import re
import shutil
import subprocess

import h5py
import numpy as np
import pytest

from paths import REPO, app

PARTRAC = app("partrac")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
SINE = os.path.join(REPO, "data_example", "sine_flow", "expr_params.dat")
BATCHELOR = os.path.join(REPO, "data_example", "batchelor_vortex", "expr_params.dat")
TAYLOR_COUETTE = os.path.join(REPO, "data_example", "taylor_couette", "expr_params.dat")

U_INF, R = 1.0, 1.0

BASE = ("mode=analytic ds_min=1e-12 refine=false coarsen=false Dm=0 "
        "int_order=1 checkpoint_intv=1e9 random=false seed=1 "
        "integrate_tau=true tau_max=0").split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")


def run(tmp_path, example, extra, expr=None):
    """One partrac run with tau integration on; returns the opened dump."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    text = open(example).read()
    for key, value in (expr or {}).items():
        text, n = re.subn(r"(?m)^%s=.*$" % key, "%s=%.12g" % (key, value), text)
        assert n == 1, "%s is not a key of %s" % (key, example)
    (tmp_path / "expr_params.dat").write_text(text)
    keys = {a.split("=")[0] for a in extra}
    argv = [a for a in BASE if a.split("=")[0] not in keys] + extra
    r = subprocess.run([PARTRAC, str(tmp_path / "expr_params.dat")] + argv,
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr
    dumps = list(tmp_path.rglob("data_from_t*.h5"))
    assert len(dumps) == 1
    return h5py.File(dumps[0], "r")


def times(dump):
    return sorted(dump.keys(), key=float)


def c_max(tau, D, s0):
    """Meunier & Villermaux equation (2.8) at n = 0, in units of c0."""
    return 1.0 / np.sqrt(1.0 + 4.0 * D * np.asarray(tau) / s0 ** 2)


def poiseuille_tau(x, t):
    """int_0^t rho^2 dt' for an edge of a line laid along x."""
    a = 3 * U_INF * np.asarray(x) / R ** 2
    return t + a ** 2 * t ** 3 / 3


def edge_midpoints_x(dump):
    """Initial x of each edge, which labels it for the whole run."""
    first = times(dump)[0]
    p = np.array(dump[first + "/points"])
    e = np.array(dump[first + "/edges"])
    return 0.5 * (p[e[:, 0], 0] + p[e[:, 1], 0])


@needs_partrac
def test_tau_matches_the_exact_integral(tmp_path):
    T = 0.5
    d = run(tmp_path, POISEUILLE,
            ["init_mode=strip_x", "La=1.0", "x0=0", "y0=0", "z0=0", "Nrw=200",
             "Nrw_max=5000", "ds_max=1e9", "dt=0.0025", "T=%g" % T,
             "stat_intv=%g" % T, "dump_intv=%g" % T, "tau_intv=0.0025"])
    tau = np.array(d[times(d)[-1] + "/tau"]).ravel()
    assert tau == pytest.approx(poiseuille_tau(edge_midpoints_x(d), T), rel=1e-5)
    # tau starts at zero: nothing has diffused before the run begins
    assert np.array(d[times(d)[0] + "/tau"]).ravel() == pytest.approx(0.0, abs=1e-15)


@needs_partrac
def test_tau_converges_at_second_order(tmp_path):
    # the integrand is quadratic in t, so the trapezoid error is -a^2 h^2 T/6
    T, err = 0.5, []
    for intv in ("0.02", "0.01", "0.005"):
        d = run(tmp_path / intv, POISEUILLE,
                ["init_mode=strip_x", "La=1.0", "x0=0", "y0=0", "z0=0",
                 "Nrw=200", "Nrw_max=5000", "ds_max=1e9", "dt=0.0025",
                 "T=%g" % T, "stat_intv=%g" % T, "dump_intv=%g" % T,
                 "tau_intv=" + intv])
        tau = np.array(d[times(d)[-1] + "/tau"]).ravel()
        exact = poiseuille_tau(edge_midpoints_x(d), T)
        err.append(np.max(np.abs(tau / exact - 1)))
    assert err[0] / err[1] == pytest.approx(4.0, rel=0.1)
    assert err[1] / err[2] == pytest.approx(4.0, rel=0.1)

    # and the size of it is the trapezoid remainder, not something else
    a = 3 * U_INF * 0.5 / R ** 2
    predicted = (a ** 2 * 0.02 ** 2 * T / 6) / poiseuille_tau(0.5, T)
    assert err[0] == pytest.approx(predicted, rel=0.05)


@needs_partrac
def test_the_maximum_concentration_follows_the_ranz_form(tmp_path):
    T = 0.5
    d = run(tmp_path, POISEUILLE,
            ["init_mode=strip_x", "La=1.0", "x0=0", "y0=0", "z0=0", "Nrw=200",
             "Nrw_max=5000", "ds_max=1e9", "dt=0.0025", "T=%g" % T,
             "stat_intv=%g" % T, "dump_intv=0.05", "tau_intv=0.0025"])
    x = edge_midpoints_x(d)
    ts = times(d)
    c_end = c_max(np.array(d[ts[-1] + "/tau"]).ravel(), D=1e-3, s0=0.01)
    assert c_end == pytest.approx(c_max(poiseuille_tau(x, T), 1e-3, 0.01), rel=1e-5)

    # it starts at c0 and only ever decays, fastest where the strip stretches most
    series = np.array([c_max(np.array(d[t + "/tau"]).ravel(), 1e-3, 0.01) for t in ts])
    assert series[0] == pytest.approx(1.0)
    assert np.all(np.diff(series, axis=0) <= 0)
    assert series[-1][np.argmax(np.abs(x))] < series[-1][np.argmin(np.abs(x))]


@needs_partrac
def test_the_reconstruction_rescales_in_the_initial_thickness(tmp_path):
    # tau is proportional to s0^-2, so a finished run can be re-read at another
    # initial thickness without advecting anything again
    T = 0.5
    d = run(tmp_path, POISEUILLE,
            ["init_mode=strip_x", "La=1.0", "x0=0", "y0=0", "z0=0", "Nrw=200",
             "Nrw_max=5000", "ds_max=1e9", "dt=0.0025", "T=%g" % T,
             "stat_intv=%g" % T, "dump_intv=%g" % T, "tau_intv=0.0025"])
    tau = np.array(d[times(d)[-1] + "/tau"]).ravel()
    thin, thick = c_max(tau, 1e-3, 0.005), c_max(tau, 1e-3, 0.01)
    assert np.all(thin < thick)                      # thinner mixes sooner
    assert thin == pytest.approx(c_max(4 * tau, 1e-3, 0.01))


@needs_partrac
def test_tau_max_drops_strips_that_have_finished_mixing(tmp_path):
    # the branch that does this is marked "untested!" in the source
    T = 0.5
    args = ["init_mode=strip_x", "La=1.0", "x0=0", "y0=0", "z0=0", "Nrw=200",
            "Nrw_max=5000", "ds_max=1e9", "dt=0.0025", "T=%g" % T,
            "stat_intv=%g" % T, "dump_intv=%g" % T, "tau_intv=0.0025"]
    kept = {}
    for tau_max in ("0", "0.55", "0.52"):
        d = run(tmp_path / tau_max, POISEUILLE, args + ["tau_max=" + tau_max])
        tau = np.array(d[times(d)[-1] + "/tau"]).ravel()
        kept[tau_max] = len(tau)
        if float(tau_max) > 0:
            assert tau.max() <= float(tau_max)
    assert kept["0"] > kept["0.55"] > kept["0.52"] > 0


@needs_partrac
def test_refinement_carries_tau_across_a_split(tmp_path):
    # a split edge takes its parent's tau and rho_prev, so its history is not
    # restarted; x is invariant in this flow, so the closed form still applies
    T = 0.5
    common = ["init_mode=strip_x", "La=1.0", "x0=0", "y0=0", "z0=0",
              "Nrw_max=20000", "dt=0.0025", "T=%g" % T, "stat_intv=%g" % T,
              "dump_intv=%g" % T, "tau_intv=0.0025"]
    err = {}
    for name, extra in (("plain", ["Nrw=200", "ds_max=1e9", "refine=false"]),
                        ("refined", ["Nrw=40", "ds_max=0.02", "refine=true",
                                     "refine_intv=0.01"])):
        d = run(tmp_path / name, POISEUILLE, common + extra)
        last = times(d)[-1]
        p = np.array(d[last + "/points"])
        e = np.array(d[last + "/edges"])
        tau = np.array(d[last + "/tau"]).ravel()
        x = 0.5 * (p[e[:, 0], 0] + p[e[:, 1], 0])
        err[name] = np.max(np.abs(tau / poiseuille_tau(x, T) - 1))
        assert len(tau) > 0
    assert err["refined"] < 1e-4        # a reset tau would be far below this
    assert err["refined"] < 2 * err["plain"]


@needs_partrac
def test_a_sheet_carries_the_same_tau_as_the_strip(tmp_path):
    # the flow depends on x alone, so a sheet spanning x and y is the strip
    # extruded, and its faces integrate the elongation the edges do
    T, ds_init = 0.5, 0.02
    d = run(tmp_path, POISEUILLE,
            ["init_mode=sheet_xy", "La=1.0", "Lb=0.4", "ds_init=%g" % ds_init,
             "x0=0", "y0=0", "z0=0", "Nrw=100", "Nrw_max=20000", "ds_max=1e9",
             "dt=0.0025", "T=%g" % T, "stat_intv=%g" % T, "dump_intv=%g" % T,
             "tau_intv=0.0025"])
    tau = np.array(d[times(d)[-1] + "/tau"]).ravel()
    assert tau.min() == pytest.approx(poiseuille_tau(0.0, T), rel=1e-3)
    # the outermost face centre sits inside the strip's last node
    assert poiseuille_tau(0.5 - ds_init, T) < tau.max() < poiseuille_tau(0.5, T)


@needs_partrac
def test_tau_in_a_vortex_matches_the_exact_integral(tmp_path):
    # nothing moves radially, so an edge keeps its radius and rho^2 has the
    # closed form of the deformation tests; q = 1 adds the axial jet
    u0, R1, R2, q, T = 2.0, 0.5, 1.0, 1.0, 1.0
    d = run(tmp_path, BATCHELOR,
            ["init_mode=strip_x", "x0=6.2", "y0=5.0", "z0=5.0", "La=1.2",
             "Nrw=400", "Nrw_max=20000", "ds_max=1e9", "dt=0.0005",
             "T=%g" % T, "stat_intv=%g" % T, "dump_intv=%g" % T,
             "tau_intv=0.0005", "int_order=2"],
            expr=dict(u0=u0, R1=R1, R2=R2, q=q))
    first, last = times(d)[0], times(d)[-1]
    p0 = np.array(d[first + "/points"])
    e = np.array(d[last + "/edges"])
    s = np.abs(0.5 * (p0[e[:, 0], 0] + p0[e[:, 1], 0]) - 5.0)

    a = s ** 2 / R1 ** 2
    s_domega = 2 * u0 * R1 * ((1 + a) * np.exp(-a) - 1) / s ** 2
    duz = -2 * s * q * u0 * np.exp(-s ** 2 / R2 ** 2) / R2 ** 2
    exact = T + (s_domega ** 2 + duz ** 2) * T ** 3 / 3

    tau = np.array(d[last + "/tau"]).ravel()
    assert tau == pytest.approx(exact, rel=1e-4)


@needs_partrac
def test_the_taylor_couette_concentration_decays_as_the_paper_predicts(tmp_path):
    # Martinez-Ruiz et al. section 3.5. A sheet tangent to the stream torus keeps
    # its thickness, so tau ~ t and c_max ~ t^-1/2; one carrying the radial
    # direction thins as t^-1, so tau ~ t^3 and c_max ~ t^-3/2. The tangency is
    # a statement about an infinitesimal triangle, and a finite patch drifts off
    # it after t ~ 100 -- unchanged when dt is quartered, so it is the patch size
    # and not the integration -- which is why the linear fit stops at 60.
    T = 60.0
    args = ["x0=2.0", "y0=0.0", "z0=0.6", "La=0.02", "Lb=0.02", "ds_init=0.002",
            "Nrw=100", "Nrw_max=400000", "ds_max=1e9", "int_order=2", "dt=0.02",
            "T=%g" % T, "stat_intv=%g" % T, "dump_intv=2.0", "tau_intv=0.02"]

    def tau_series(name, mode):
        d = run(tmp_path / name, TAYLOR_COUETTE, ["init_mode=" + mode] + args)
        ks = times(d)
        t = np.array([float(k) for k in ks])
        tau = np.array([np.mean(np.array(d[k + "/tau"]).ravel()) for k in ks])
        return t[t > 0], tau[t > 0]

    t, tangent = tau_series("tangent", "sheet_yz")
    _, radial = tau_series("radial", "sheet_xy")

    assert np.corrcoef(t, tangent)[0, 1] ** 2 > 0.95          # tau ~ t
    late = t >= 20
    assert np.polyfit(np.log(t[late]), np.log(tangent[late]), 1)[0] < 1.5

    slope = np.polyfit(np.log(t[late]), np.log(radial[late]), 1)[0]
    assert slope == pytest.approx(3.0, rel=0.05)              # tau ~ t^3

    # and so the concentrations decay as t^-1/2 and t^-3/2
    for series, expected in ((tangent, -0.5), (radial, -1.5)):
        c = c_max(series, D=1e4, s0=1.0)     # well past the mixing time
        got = np.polyfit(np.log(t[late]), np.log(c[late]), 1)[0]
        assert got == pytest.approx(expected, abs=0.15)


@needs_partrac
def test_tau_max_drops_faces_as_well_as_edges(tmp_path):
    T = 0.5
    args = ["init_mode=sheet_xy", "La=1.0", "Lb=0.4", "ds_init=0.05", "x0=0",
            "y0=0", "z0=0", "Nrw=100", "Nrw_max=20000", "ds_max=1e9",
            "dt=0.0025", "T=%g" % T, "stat_intv=%g" % T, "dump_intv=%g" % T,
            "tau_intv=0.0025"]
    kept = {}
    for tau_max in ("0", "0.55", "0.52"):
        d = run(tmp_path / tau_max, POISEUILLE, args + ["tau_max=" + tau_max])
        tau = np.array(d[times(d)[-1] + "/tau"]).ravel()
        kept[tau_max] = len(tau)
        if float(tau_max) > 0:
            assert tau.max() <= float(tau_max)
    assert kept["0"] > kept["0.55"] > kept["0.52"] > 0


@needs_partrac
def test_a_restart_carries_the_compressed_time(tmp_path):
    # the edge and face checkpoints carry tau and rho_prev, so a restarted run
    # resumes the history rather than starting from an unmixed element
    T1 = 0.25
    common = ["init_mode=strip_x", "La=1.0", "x0=0", "y0=0", "z0=0", "Nrw=100",
              "Nrw_max=5000", "ds_max=1e9", "dt=0.0025", "stat_intv=1e9",
              "tau_intv=0.0025"]
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(POISEUILLE, tmp_path / "expr_params.dat")

    def call(extra):
        argv = [a for a in BASE if a.split("=")[0]
                not in {b.split("=")[0] for b in common + extra}] + common + extra
        r = subprocess.run([PARTRAC, str(tmp_path / "expr_params.dat")] + argv,
                           capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stdout + r.stderr

    call(["T=%g" % T1, "dump_intv=1e9", "checkpoint_intv=%g" % T1])
    checkpoint = list(tmp_path.rglob("edges.edge"))
    assert len(checkpoint) == 1
    # a widened checkpoint: node, node, ds0, tau, rho_prev
    assert len(checkpoint[0].read_text().split("\n")[0].split()) == 5
    call(["T=0.5", "dump_intv=0.2", "checkpoint_intv=1e9",
          "restart_folder=" + str(checkpoint[0].parent.parent)])

    resumed = [f for f in sorted(tmp_path.rglob("data_from_t*.h5"))
               if "t0.000000" not in f.name]
    assert len(resumed) == 1
    d = h5py.File(resumed[0], "r")
    first = times(d)[0]
    tau = np.array(d[first + "/tau"]).ravel()
    x = edge_midpoints_x(d)
    # the checkpoint holds tau belonging to the positions it stores, so the
    # resumed tau matches the time the dump is labelled with
    assert tau == pytest.approx(poiseuille_tau(x, float(first)), rel=1e-4)
    assert tau.min() > 0.2                      # emphatically not reset to zero
