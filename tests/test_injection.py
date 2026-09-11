"""Injection of a sheet into Hagen-Poiseuille flow, against the swept area.

`inject_edges` advects an inlet curve and stitches each generation to the last,
so the mesh is the streak surface that curve traces out. For a steady parallel
flow u = u(x,y) z-hat that surface is z = u(x,y) tau, and both of the measures
partrac reports have a closed form.

The reference area is fixed where material is created and conserved afterwards,
so it accumulates at the rate the curve sweeps,

    d/dt sum dA0 = int_C u |dl| ,

which along a diameter of the pipe, where u_z = 2 u_inf (1 - r^2/R^2), is

    int_-R^R 2 u_inf (1 - x^2/R^2) dx = (8/3) u_inf R .

The inlet is stitched with straight edges, so what is actually integrated is
the trapezoidal rule over the inlet nodes. The integrand is quadratic, so that
rule's error is exactly (b-a) h^2 |f''| / 12 = (2/3) u_inf h^2 with no higher
term, and the discrete rate is closed form as well:

    d/dt sum dA0 = (8/3) u_inf R - (2/3) u_inf h^2 .

The current area is the other measure, and in a parallel flow it equals the
first. Each triangle of the streak surface has an edge along z joining two
generations of the same inlet node, and the two vertices sharing (x, y) share
their velocity, so the shear slides the triangle without changing its area:
dA/dA0 = 1 at every face, for all time. Remeshing is what breaks that, since a
split puts the new node on the chord rather than on the surface -- the area it
carries still does not move.
"""

import os
import shutil
import subprocess

import h5py
import numpy as np
import pytest

from paths import REPO, app

PARTRAC = app("partrac")
HAGEN = os.path.join(REPO, "data_example", "hagen_poiseuille", "expr_params.dat")
ABC = os.path.join(REPO, "data_example", "abc_flow_unsteady", "expr_params.dat")
SINE3D = os.path.join(REPO, "data_example", "sine_flow_3d", "expr_params.dat")
PI = 3.14159265358979

U_INF, R, DT, INTV = 1.0, 1.0, 0.005, 0.05
BASE = ("mode=analytic init_mode=uniform_x x0=0 y0=0 z0=0 Nrw_max=200000 "
        "inject=true inject_edges=true inject_intv=%g T_inject=1e10 "
        "ds_max=1e9 ds_min=1e-9 refine=false coarsen=false Dm=0 int_order=1 "
        "dt=%g stat_intv=1e9 checkpoint_intv=1e9 random=false seed=1"
        % (INTV, DT)).split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")


def sweep_rate(n):
    """The rate the inlet sweeps area, as its straight edges measure it."""
    h = 2 * R / (n - 1)
    return 8. / 3. * U_INF * R - 2. / 3. * U_INF * h ** 2


def series(tmp_path, extra, example=HAGEN):
    """Every dump that has faces, as (t, dA, dA0)."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(example, tmp_path / "expr_params.dat")
    keys = {a.split("=")[0] for a in extra}
    argv = [a for a in BASE if a.split("=")[0] not in keys] + extra
    r = subprocess.run([PARTRAC, str(tmp_path / "expr_params.dat")] + argv,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    dump = list(tmp_path.rglob("data_from_t*.h5"))
    assert len(dump) == 1
    h = h5py.File(dump[0], "r")
    out = []
    for k in sorted(h.keys(), key=float):
        if "dA0" not in h[k]:
            continue                   # t = 0, before the first injection
        out.append((float(k),
                    np.array(h[k + "/dA"]).ravel(),
                    np.array(h[k + "/dA0"]).ravel()))
    return out


def lines(tmp_path, extra):
    """The last dump of a run that stays one-dimensional, as (t, group)."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(HAGEN, tmp_path / "expr_params.dat")
    keys = {a.split("=")[0] for a in extra}
    argv = [a for a in BASE if a.split("=")[0] not in keys] + extra
    r = subprocess.run([PARTRAC, str(tmp_path / "expr_params.dat")] + argv,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    h = h5py.File(list(tmp_path.rglob("data_from_t*.h5"))[0], "r")
    k = sorted(h.keys(), key=float)[-1]
    return float(k), h[k]


def rake_rate(n):
    """Sum of the speed over the inlet nodes: no quadrature, just a sum."""
    x = np.linspace(-R, R, n)
    return (2 * U_INF * (1 - x ** 2 / R ** 2)).sum()


@needs_partrac
def test_the_swept_area_grows_at_the_rate_the_inlet_sweeps(tmp_path):
    n = 41
    s = series(tmp_path, ["Nrw=%d" % n, "T=0.3", "dump_intv=%g" % INTV])
    assert len(s) == 6                                  # one dump per injection
    for t, _, dA0 in s:
        assert dA0.sum() == pytest.approx(sweep_rate(n) * t, rel=1e-12)
        assert dA0.min() > 0                            # no face without area
    # two faces per inlet edge, less the one each end loses to the no-slip wall
    assert [len(dA0) for _, _, dA0 in s] == [78 * k for k in range(1, 7)]


@needs_partrac
def test_the_streak_surface_of_a_parallel_flow_does_not_stretch(tmp_path):
    s = series(tmp_path, ["Nrw=41", "T=0.3", "dump_intv=%g" % INTV])
    for t, dA, dA0 in s:
        assert np.max(np.abs(dA / dA0 - 1)) < 1e-12     # every face keeps its area
        assert dA.sum() == pytest.approx(dA0.sum(), rel=1e-12)


@needs_partrac
@pytest.mark.parametrize("n", [11, 21, 41, 81])
def test_the_swept_area_converges_on_the_exact_integral(tmp_path, n):
    t, _, dA0 = series(tmp_path, ["Nrw=%d" % n, "T=0.2", "dump_intv=0.2"])[-1]
    h = 2 * R / (n - 1)
    shortfall = 8. / 3. * U_INF * R - dA0.sum() / t
    assert shortfall == pytest.approx(2. / 3. * U_INF * h ** 2, rel=1e-9)


def chord_rate(La, x0, n):
    """The same rate over a chord of length La centred on x0, not the diameter."""
    a, b = x0 - La / 2, x0 + La / 2
    exact = 2 * U_INF * ((b - a) - (b ** 3 - a ** 3) / (3 * R ** 2))
    h = La / (n - 1)
    return exact - La * h ** 2 * U_INF / (3 * R ** 2)


@needs_partrac
@pytest.mark.parametrize("La,x0", [(2.0, 0.0), (1.0, 0.0), (0.5, 0.25)])
def test_a_strip_inlet_sweeps_the_chord_it_covers(tmp_path, La, x0):
    # strip_x is the inlet that takes an extent: uniform_x always spans the
    # domain, and La is what lets a chord be shorter than the diameter
    n = 41
    s = series(tmp_path, ["init_mode=strip_x", "La=%g" % La, "x0=%g" % x0,
                          "Nrw=%d" % n, "T=0.2", "dump_intv=0.1"])
    assert len(s) == 2
    for t, _, dA0 in s:
        assert dA0.sum() == pytest.approx(chord_rate(La, x0, n) * t, rel=1e-12)
        assert dA0.min() > 0


@needs_partrac
def test_a_strip_across_the_whole_pipe_is_the_uniform_inlet(tmp_path):
    # the two initializers lay the same 41 nodes down when La spans the domain,
    # so this is the same sheet twice -- to round-off, since strip runs from
    # +La/2 back to -La/2 and the faces come out in the other order
    args = ["Nrw=41", "T=0.2", "dump_intv=0.2"]
    _, dA_u, dA0_u = series(tmp_path / "uniform",
                            ["init_mode=uniform_x"] + args)[-1]
    _, dA_s, dA0_s = series(tmp_path / "strip",
                            ["init_mode=strip_x", "La=2.0"] + args)[-1]
    assert len(dA0_s) == len(dA0_u)
    assert dA0_s.sum() == pytest.approx(dA0_u.sum(), rel=1e-14)
    assert np.allclose(np.sort(dA0_s), np.sort(dA0_u), rtol=1e-13, atol=0)
    assert np.allclose(np.sort(dA_s), np.sort(dA_u), rtol=1e-13, atol=0)


@needs_partrac
def test_a_point_inlet_traces_a_chain_not_a_sheet(tmp_path):
    # a point has no edges to sweep, so injection stitches each generation to
    # the last with one rung apiece and the mesh stays one-dimensional
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(HAGEN, tmp_path / "expr_params.dat")
    extra = ["init_mode=point", "Nrw=20", "z0=-0.5", "T=0.2", "dump_intv=0.2"]
    keys = {a.split("=")[0] for a in extra}
    argv = [a for a in BASE if a.split("=")[0] not in keys] + extra
    r = subprocess.run([PARTRAC, str(tmp_path / "expr_params.dat")] + argv,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    h = h5py.File(list(tmp_path.rglob("data_from_t*.h5"))[0], "r")
    key = sorted(h.keys(), key=float)[-1]
    assert "faces" not in h[key]                      # never became a sheet
    assert len(np.array(h[key + "/points"])) == 20 * 5    # four injections
    assert len(np.array(h[key + "/edges"])) == 20 * 4     # a rung each time


@needs_partrac
@pytest.mark.parametrize("init_mode,extra", [
    ("sheet_xy", ["La=0.5", "Lb=0.5", "ds_init=0.1"]),
    ("ellipsoid_xy", ["La=0.5", "Lb=0.5"]),
])
def test_an_inlet_with_faces_is_refused(tmp_path, init_mode, extra):
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(HAGEN, tmp_path / "expr_params.dat")
    ex = ["init_mode=" + init_mode, "Nrw=41", "T=0.1", "dump_intv=0.1"] + extra
    keys = {a.split("=")[0] for a in ex}
    argv = [a for a in BASE if a.split("=")[0] not in keys] + ex
    r = subprocess.run([PARTRAC, str(tmp_path / "expr_params.dat")] + argv,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode != 0
    assert "would sweep a volume" in r.stdout + r.stderr


@needs_partrac
def test_a_cleared_inlet_draws_streaklines(tmp_path):
    # clear_initial_edges leaves a cloud, and injection raises what an inlet
    # traces by one dimension, so a cleared curve draws lines rather than a
    # surface. Each is straight and grows at its own node's speed, so the total
    # is a plain sum over the inlet -- no quadrature, and no error term.
    n = 41
    t, g = lines(tmp_path, ["init_mode=uniform_x", "Nrw=%d" % n, "T=0.2",
                            "dump_intv=0.2", "clear_initial_edges=true"])
    assert "faces" not in g                            # a line, not a surface
    dl0 = np.array(g["dl0"]).ravel()
    assert dl0.sum() == pytest.approx(rake_rate(n) * t, rel=1e-12)
    assert dl0.min() > 0            # the wall nodes are reused, not doubled
    # two of the 41 never leave the wall, so they lay no rung
    assert len(dl0) == 4 * (n - 2)
    assert len(np.array(g["points"])) == n + 4 * (n - 2)


@needs_partrac
def test_an_unstitched_inlet_repeats_itself(tmp_path):
    # without inject_edges the generations are not joined to each other, so the
    # inlet reappears as it is and the dimension does not rise
    n = 41
    t, g = lines(tmp_path, ["init_mode=uniform_x", "Nrw=%d" % n, "T=0.2",
                            "dump_intv=0.2", "inject_edges=false"])
    assert "faces" not in g
    dl0 = np.array(g["dl0"]).ravel()
    assert len(dl0) == 5 * (n - 1)                     # five copies of the inlet
    assert len(np.array(g["points"])) == 5 * n
    assert dl0.sum() == pytest.approx(5 * 2 * R, rel=1e-12)   # each spans the pipe


@needs_partrac
def test_a_cleared_sheet_inlet_draws_lines_too(tmp_path):
    # the rule is about what survives the clear, not about the init_mode: a
    # sheet is refused as an inlet, but a cleared one is a cloud like any other
    t, g = lines(tmp_path, ["init_mode=sheet_xy", "La=0.5", "Lb=0.5",
                            "ds_init=0.1", "Nrw=41", "T=0.2", "dump_intv=0.2",
                            "clear_initial_edges=true"])
    assert "faces" not in g
    assert np.array(g["dl0"]).min() > 0


@needs_partrac
def test_refining_faster_than_injecting(tmp_path):
    # everything before the first injection is a strip, and every one of its
    # edges is the inlet. Refining one leaves edges_inlet naming half a template
    # edge, and edge2faces one row short of edges for every split.
    #
    # The injection interval here is long enough that a generation stretches
    # past ds_max before the next one arrives -- the inlet spans the pipe, so
    # the shear pulls its ends 0.1 downstream in 0.5. Splitting such an edge
    # splits the template with it, so later generations are laid down on more
    # nodes than the first and the closed form for a fixed inlet no longer
    # holds. What is still true is that the rate only ever improves: the
    # trapezoidal rule under-reads a concave profile, and refining the template
    # walks it toward the exact integral without ever passing it.
    n = 41
    s = series(tmp_path, ["Nrw=%d" % n, "T=1.5", "dump_intv=0.5",
                          "inject_intv=0.5", "ds_max=0.1", "ds_min=0.02",
                          "refine=true", "refine_intv=0.05",
                          "coarsen=true", "coarsen_intv=0.05"])
    assert len(s) == 3
    assert len(s[-1][2]) > 3000                    # it really refined
    exact = 8. / 3. * U_INF * R
    for t, _, dA0 in s:
        # the first generation is still on the inlet as laid down, so the
        # lower bound is met exactly there and only exceeded later
        assert dA0.sum() >= sweep_rate(n) * t * (1 - 1e-12)
        assert dA0.sum() <= exact * t              # never past the integral
        assert dA0.min() > 0
    # the last generation is swept on a finer inlet than the first
    assert s[-1][2].sum() / s[-1][0] > s[0][2].sum() / s[0][0]


@needs_partrac
def test_an_exit_plane_does_not_cut_the_inlet_away(tmp_path):
    # the inlet spans the pipe, so half of it is beyond the plane from the
    # start. It has to outlive the cull by a generation or the next injection
    # has nothing to stitch to.
    s = series(tmp_path, ["Nrw=41", "T=0.5", "dump_intv=0.1",
                          "inject_intv=0.1", "exit_plane=x", "Ln=0.5",
                          "filter_intv=0.1"])
    assert len(s) == 5
    total = [dA0.sum() for _, _, dA0 in s]
    assert all(b > a for a, b in zip(total, total[1:]))   # still being fed
    assert all(dA0.min() > 0 for _, _, dA0 in s)


@needs_partrac
def test_remeshing_conserves_the_swept_area(tmp_path):
    s = series(tmp_path, ["Nrw=41", "T=0.3", "dump_intv=0.1", "ds_max=0.06",
                          "ds_min=0.02", "refine=true", "refine_intv=%g" % INTV,
                          "coarsen=true", "coarsen_intv=%g" % INTV])
    assert len(s) == 3
    assert len(s[-1][2]) > 1000                         # it really refined
    for t, dA, dA0 in s:
        assert dA0.sum() == pytest.approx(sweep_rate(41) * t, rel=1e-12)
        assert dA0.min() > 0
        assert np.max(np.abs(dA / dA0 - 1)) < 0.02      # the chord, not the mass


# --- a flow that is three-dimensional and does not stand still ------------------
#
# Blazevski & Haller, Physica D 273-274 (2014) 46-62, eq. (22): the ABC flow with
# A forced in time, A = sqrt(3), B = sqrt(2), C = 1 on [0, 2pi]^3, forcing period
# 2pi. Divergence free, and its steady limit is an exact solution of Euler's
# equation. Nothing here is closed form -- the sheet folds -- so what is checked
# is the invariant that does not care about the flow: reference area is created
# at the inlet and moved by nothing afterwards.

ABC_ARGS = ["init_mode=uniform_x", "Nrw=41", "Nrw_max=400000",
            "x0=%.14f" % PI, "y0=%.14f" % PI, "z0=%.14f" % PI,
            "inject_intv=0.25", "int_order=2", "dt=0.02", "T=1.5",
            "dump_intv=0.5", "integrate_tau=true", "tau_intv=0.02", "tau_max=0"]
ABC_REMESH = ["ds_max=0.4", "ds_min=0.1", "refine=true", "refine_intv=0.25",
              "coarsen=true", "coarsen_intv=0.25"]


@needs_partrac
def test_remeshing_a_sheet_in_an_unsteady_flow_moves_no_area(tmp_path):
    # the swept-area law needs a steady flow; this one does not. What survives
    # a folding sheet in a chaotic flow is that dA0 is laid down at the inlet
    # and untouched after, so refining and coarsening cannot shift the total.
    plain = series(tmp_path / "plain", ABC_ARGS, example=ABC)
    remeshed = series(tmp_path / "remeshed", ABC_ARGS + ABC_REMESH, example=ABC)
    assert len(plain) == len(remeshed) == 3
    assert len(remeshed[-1][2]) > 2 * len(plain[-1][2])      # it really remeshed
    for (t, _, a), (t2, _, b) in zip(plain, remeshed):
        assert t == t2
        assert b.sum() == pytest.approx(a.sum(), rel=1e-12)
        assert a.min() > 0 and b.min() > 0


@needs_partrac
def test_an_unsteady_three_dimensional_flow_stretches_the_sheet(tmp_path):
    # the counterpart of the parallel-flow case, where dA/dA0 is 1 to 1e-12:
    # here the sheet is stretched and compressed by an order of magnitude, and
    # tau stays finite through all of it
    s = series(tmp_path, ABC_ARGS + ABC_REMESH, example=ABC)
    t, dA, dA0 = s[-1]
    assert (dA / dA0).max() > 2                     # stretched
    assert (dA / dA0).min() < 0.6                   # and compressed
    assert (dA / dA0).max() / (dA / dA0).min() > 5  # measured 7.5 at T = 1.5
    tmp = list(tmp_path.rglob("data_from_t*.h5"))
    h = h5py.File(tmp[0], "r")
    key = sorted(h.keys(), key=float)[-1]
    for name in ("dA", "dA0", "tau"):
        assert np.isfinite(np.array(h[key + "/" + name]).astype(float)).all(), name


@needs_partrac
@pytest.mark.parametrize("remesh", [
    [],
    ["ds_max=0.2", "ds_min=0.02", "refine=true", "refine_intv=0.05",
     "coarsen=true", "coarsen_intv=0.05"],
])
def test_an_inlet_laid_along_the_flow_sweeps_nothing(tmp_path, remesh):
    # uniform_z runs down the pipe, so the inlet moves along itself and there is
    # no area between one generation and the next. The corners are laid down all
    # the same, since the generation after has to have something to stitch to,
    # and go as soon as they are not at the inlet any more -- which is driven by
    # the injection, so it holds whether or not anything is being remeshed.
    s = series(tmp_path, ["init_mode=uniform_z", "Nrw=41", "T=0.3",
                          "dump_intv=0.1", "stat_intv=0.05",
                          "integrate_tau=true", "tau_intv=0.005",
                          "tau_max=0"] + remesh)
    assert len(s) == 3
    for t, dA, dA0 in s:
        assert len(dA0) == 40              # one generation's worth, never more
        assert (dA0 == 0).all()            # because nothing was swept
    h = h5py.File(list(tmp_path.rglob("data_from_t*.h5"))[0], "r")
    key = sorted(h.keys(), key=float)[-1]
    assert len(np.array(h[key + "/points"])) == 81      # nor do the nodes pile up
    assert np.isfinite(np.array(h[key + "/tau"]).astype(float)).all()
    stats = list(tmp_path.rglob("tdata_from_t*.dat"))
    assert len(stats) == 1
    text = stats[0].read_text().lower()
    assert "nan" not in text and "inf" not in text


# --- and one where the inlet is parallel to the flow every third step ----------
#
# Meunier & Villermaux, J. Fluid Mech. 951 (2022) A33, section 4.1, eqs (4.1)-(4.3):
# the three-dimensional random sine flow, wavelength 1, its period cut into three
# steps of a third -- flow along y depending on x, then along z depending on y,
# then along x depending on z, with the phases of their table 1 and U = 0.3. A
# uniform_x inlet lies along the third step's flow, so one step in three sweeps
# no area at all, which no steady flow does unless it is pointed at deliberately.

SINE_ARGS = ["init_mode=uniform_x", "Nrw=41", "Nrw_max=200000",
             "x0=0.5", "y0=0.5", "z0=0.5", "inject_intv=%.10f" % (1. / 15),
             "dt=%.10f" % (1. / 300), "T=2.0", "dump_intv=%.10f" % (1. / 3),
             "stat_intv=%.10f" % (1. / 15), "integrate_tau=true",
             "tau_intv=%.10f" % (1. / 300), "tau_max=0"]
SINE_REMESH = ["ds_max=0.05", "ds_min=0.01", "refine=true",
               "refine_intv=%.10f" % (1. / 15),
               "coarsen=true", "coarsen_intv=%.10f" % (1. / 15)]


@needs_partrac
def test_a_step_that_sweeps_nothing_neither_adds_nor_accumulates(tmp_path):
    s = series(tmp_path, SINE_ARGS, example=SINE3D)
    assert len(s) == 6
    flat = [int((dA0 <= 0).sum()) for _, _, dA0 in s]
    assert set(flat) <= {0, 40}                  # a generation's worth, or none
    assert flat.count(40) == 2                   # the two parallel steps
    total = [dA0.sum() for _, _, dA0 in s]
    assert all(b >= a for a, b in zip(total, total[1:]))         # never lost
    gain = [b - a for a, b in zip([0.] + total, total)]
    assert min(g for g, f in zip(gain, flat) if f == 0) > 20 * max(
        g for g, f in zip(gain, flat) if f == 40)   # the parallel steps add ~nothing
    h = h5py.File(list(tmp_path.rglob("data_from_t*.h5"))[0], "r")
    key = sorted(h.keys(), key=float)[-1]
    for name in ("dA", "dA0", "tau"):
        assert np.isfinite(np.array(h[key + "/" + name]).astype(float)).all(), name
    p = np.array(h[key + "/points"])
    assert min(p[:, i].ptp() for i in range(3)) > 0.1     # it really is 3-D
    stats = list(tmp_path.rglob("tdata_from_t*.dat"))
    text = stats[0].read_text().lower()
    assert "nan" not in text and "inf" not in text


@needs_partrac
def test_remeshing_moves_no_area_across_a_parallel_step(tmp_path):
    # the same invariant as in the ABC flow, but over a sheet that stops growing
    # and starts again, with flat faces created and culled in between
    plain = series(tmp_path / "plain", SINE_ARGS, example=SINE3D)
    remeshed = series(tmp_path / "remeshed", SINE_ARGS + SINE_REMESH,
                      example=SINE3D)
    for (t, _, a), (t2, _, b) in zip(plain, remeshed):
        assert t == t2
        assert b.sum() == pytest.approx(a.sum(), rel=1e-9)
    # it really remeshed -- here coarsening wins, so the count falls
    assert abs(len(remeshed[-1][2]) - len(plain[-1][2])) > 0.2 * len(plain[-1][2])


@needs_partrac
def test_a_run_names_its_columns_for_the_dimension_it_settles_into(tmp_path):
    # The header is written once, before the loop, when an injecting run is
    # still the inlet curve it started from. Naming the columns for that put
    # the swept area under `s` and `s0` for the whole run: the column labelled
    # `s` read 0.13325 at t = 0.05, which is an area. Both guards missed it --
    # test_apps and test_degenerate check that the counts agree, and they did.
    def stats(case, extra):
        case.mkdir(parents=True, exist_ok=True)
        shutil.copy(HAGEN, case / "expr_params.dat")
        keys = {a.split("=")[0] for a in extra}
        argv = [a for a in BASE if a.split("=")[0] not in keys] + extra
        r = subprocess.run([PARTRAC, str(case / "expr_params.dat")] + argv,
                           capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stdout + r.stderr
        f = list(case.rglob("tdata_from_t*.dat"))
        assert len(f) == 1
        rows = [l for l in f[0].read_text().splitlines() if l.strip()]
        names = [h for h in rows[0].lstrip("# ").rstrip().split("\t") if h.strip()]
        return names, [dict(zip(names, [v for v in r_.rstrip().split("\t")
                                        if v.strip()])) for r_ in rows[1:]]

    names, rows = stats(tmp_path / "injecting",
                        ["Nrw=41", "T=0.15", "stat_intv=0.025", "dump_intv=1e9"])
    assert "A" in names and "A0" in names
    assert "s" not in names and "s0" not in names
    # nothing has been swept before the first injection, and saying so is not
    # the same as reporting the length of the curve about to sweep it
    assert float(rows[0]["A0"]) == 0.0
    assert float(rows[2]["A0"]) == pytest.approx(sweep_rate(41) * INTV, rel=1e-12)

    # the same run without injection stays a strip, and still says so
    names, _ = stats(tmp_path / "strip",
                     ["Nrw=41", "T=0.15", "stat_intv=0.025", "dump_intv=1e9",
                      "inject=false"])
    assert "s" in names and "s0" in names
    assert "A" not in names and "A0" not in names
