"""Injection of a sheet into Hagen-Poiseuille flow, against the swept area.

`inject_edges` advects an inlet curve and stitches each generation to the last,
so the mesh is the streak surface that curve traces out. For a steady parallel
flow u = u(x,y) z-hat that surface is z = u(x,y) tau, and both area measures
partrac reports have a closed form.

The reference area dA0 is fixed where material is created and conserved
afterwards, so it accumulates at the rate the inlet sweeps, d/dt sum dA0 =
int_C u |dl|. Along a diameter of the pipe, where u_z = 2 u_inf (1 - r^2/R^2),
that is (8/3) u_inf R. The inlet is stitched with straight edges, so what is
integrated is the trapezoidal rule over the inlet nodes; the integrand is
quadratic, so the rule's error is exactly (2/3) u_inf h^2 and

    d/dt sum dA0 = (8/3) u_inf R - (2/3) u_inf h^2 .

The current area dA equals dA0 in a parallel flow: each triangle has an edge
along z joining two generations of one inlet node, whose vertices share (x, y)
and so their velocity, and the shear slides the triangle without changing its
area. Remeshing puts split nodes on the chord rather than the surface, which
moves dA slightly but never dA0.

The later sections use the unsteady ABC flow and the 3-D random sine flow, where
nothing is closed form and the checked invariant is that dA0 is created at the
inlet and moved by nothing afterwards.
"""

import os

import numpy as np
import pytest

from dumps import all_dumps, read_stats
from paths import REPO, app
from runs import copy_example, run_app

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
    """Swept-area rate of an n-node inlet across the diameter, as its straight edges measure it."""
    h = 2 * R / (n - 1)
    return 8. / 3. * U_INF * R - 2. / 3. * U_INF * h ** 2


def run(tmp_path, extra, example=HAGEN):
    """Run partrac on a copy of example with extra overriding BASE; return time -> datasets as written."""
    run_app(PARTRAC, copy_example(example, tmp_path), BASE, extra, timeout=600)
    assert len(list(tmp_path.rglob("data_from_t*.h5"))) == 1
    return all_dumps(tmp_path, raw=True)


def faces_of(dumps):
    """Every dump that has faces, as (t, dA, dA0); t = 0 comes before the first injection."""
    return [(t, g["dA"].ravel(), g["dA0"].ravel())
            for t, g in sorted(dumps.items()) if "dA0" in g]


def series(tmp_path, extra, example=HAGEN):
    """Every dump of a run that has faces, as (t, dA, dA0)."""
    return faces_of(run(tmp_path, extra, example))


def lines(tmp_path, extra):
    """The last dump of a run that stays one-dimensional, as (t, datasets)."""
    dumps = run(tmp_path, extra)
    return max(dumps), dumps[max(dumps)]


def module_run(tmp_path_factory, name, extra, example=HAGEN):
    """(folder, run()) of a run in a folder of its own, for a run that several tests read."""
    if not os.path.exists(PARTRAC):
        pytest.skip("partrac is not built")
    d = tmp_path_factory.mktemp(name)
    return d, run(d, extra, example)


def rake_rate(n):
    """Sum of the speed over the n inlet nodes: the growth rate of their streaklines."""
    x = np.linspace(-R, R, n)
    return (2 * U_INF * (1 - x ** 2 / R ** 2)).sum()


@pytest.fixture(scope="module")
def hagen(tmp_path_factory):
    """A 41-node inlet across the pipe, injected and dumped every INTV to T = 0.3."""
    return module_run(tmp_path_factory, "hagen", ["Nrw=41", "T=0.3", "dump_intv=%g" % INTV])[1]


def test_the_swept_area_grows_at_the_rate_the_inlet_sweeps(hagen):
    """At every dump sum dA0 equals the discrete sweep rate times t to
    round-off, every face has positive area and three distinct nodes, and each
    injection adds the same number of faces. dA0 is the reference every A/A0
    statistic divides by, so a wrong rate would bias the elongation of the
    whole sheet."""
    n = 41
    s = faces_of(hagen)
    assert len(s) == 6                                  # one dump per injection
    for t, _, dA0 in s:
        assert dA0.sum() == pytest.approx(sweep_rate(n) * t, rel=1e-12)
        assert dA0.min() > 0
    # two faces per inlet edge, less the one each end loses to the no-slip wall
    assert [len(dA0) for _, _, dA0 in s] == [78 * k for k in range(1, 7)]
    tri = hagen[max(hagen)]["faces"]
    assert all(len(set(row.tolist())) == 3 for row in tri)   # no repeated node


def test_the_streak_surface_of_a_parallel_flow_does_not_stretch(hagen):
    """In a parallel flow the shear slides each face without changing its area,
    so dA/dA0 = 1 at every face to round-off. A departure means the current
    area is measured from the wrong vertices or the faces are stitched wrongly."""
    for t, dA, dA0 in faces_of(hagen):
        assert np.max(np.abs(dA / dA0 - 1)) < 1e-12
        assert dA.sum() == pytest.approx(dA0.sum(), rel=1e-12)


@needs_partrac
@pytest.mark.parametrize("n", [11, 21, 41, 81])
def test_the_swept_area_converges_on_the_exact_integral(tmp_path, n):
    """The shortfall from the exact rate (8/3) u_inf R is exactly the
    trapezoidal error (2/3) u_inf h^2 at every inlet resolution, so the swept
    area converges on the true integral at second order and has no other error."""
    t, _, dA0 = series(tmp_path, ["Nrw=%d" % n, "T=0.2", "dump_intv=0.2"])[-1]
    h = 2 * R / (n - 1)
    shortfall = 8. / 3. * U_INF * R - dA0.sum() / t
    assert shortfall == pytest.approx(2. / 3. * U_INF * h ** 2, rel=1e-9)


def chord_rate(La, x0, n):
    """Swept-area rate of an n-node inlet along a chord of length La centred on x0."""
    a, b = x0 - La / 2, x0 + La / 2
    exact = 2 * U_INF * ((b - a) - (b ** 3 - a ** 3) / (3 * R ** 2))
    h = La / (n - 1)
    return exact - La * h ** 2 * U_INF / (3 * R ** 2)


@needs_partrac
@pytest.mark.parametrize("La,x0", [(2.0, 0.0), (1.0, 0.0), (0.5, 0.25)])
def test_a_strip_inlet_sweeps_the_chord_it_covers(tmp_path, La, x0):
    """A strip_x inlet of length La centred on x0 sweeps area at the closed-form
    rate for that chord, including off-centre chords. strip_x is the inlet that
    takes an extent (uniform_x always spans the domain), so this is how a user
    injects over part of the pipe."""
    n = 41
    s = series(tmp_path, ["init_mode=strip_x", "La=%g" % La, "x0=%g" % x0,
                          "Nrw=%d" % n, "T=0.2", "dump_intv=0.1"])
    assert len(s) == 2
    for t, _, dA0 in s:
        assert dA0.sum() == pytest.approx(chord_rate(La, x0, n) * t, rel=1e-12)
        assert dA0.min() > 0


@needs_partrac
def test_a_strip_across_the_whole_pipe_is_the_uniform_inlet(tmp_path):
    """With La equal to the diameter, strip_x and uniform_x lay down the same
    nodes, so they must produce the same sheet: same face count, same dA and dA0."""
    # equal only to round-off: strip_x runs from +La/2 back to -La/2, so the
    # faces come out in the other order
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
    """A point inlet has no edges to sweep, so each generation is joined to the
    last by one rung per particle and the mesh stays one-dimensional, with no faces."""
    _, g = lines(tmp_path, ["init_mode=point", "Nrw=20", "z0=-0.5", "T=0.2", "dump_intv=0.2"])
    assert "faces" not in g
    assert len(g["points"]) == 20 * 5    # the inlet and four injections
    assert len(g["edges"]) == 20 * 4     # a rung per particle per injection


@needs_partrac
@pytest.mark.parametrize("init_mode,extra", [
    ("sheet_xy", ["La=0.5", "Lb=0.5", "ds_init=0.1"]),
    ("ellipsoid_xy", ["La=0.5", "Lb=0.5"]),
])
def test_an_inlet_with_faces_is_refused(tmp_path, init_mode, extra):
    """An inlet that already has faces would sweep a volume, which the surface
    mesh cannot represent, so partrac must stop with a clear error instead of
    producing a meaningless mesh."""
    r = run_app(PARTRAC, copy_example(HAGEN, tmp_path), BASE,
                ["init_mode=" + init_mode, "Nrw=41", "T=0.1", "dump_intv=0.1"], extra,
                check=False, timeout=600)
    assert r.returncode != 0
    assert "would sweep a volume" in r.stdout + r.stderr


@needs_partrac
def test_a_cleared_inlet_draws_streaklines(tmp_path):
    """clear_initial_edges leaves a cloud, and injection raises what an inlet
    traces by one dimension, so a cleared curve draws one streakline per node
    rather than a surface. Each line is straight and grows at its own node's
    speed, so sum dl0 is a plain sum over the inlet with no quadrature error."""
    n = 41
    t, g = lines(tmp_path, ["init_mode=uniform_x", "Nrw=%d" % n, "T=0.2",
                            "dump_intv=0.2", "clear_initial_edges=true"])
    assert "faces" not in g
    dl0 = g["dl0"].ravel()
    assert dl0.sum() == pytest.approx(rake_rate(n) * t, rel=1e-12)
    assert dl0.min() > 0            # the wall nodes are reused, not doubled
    # the two wall nodes never move, so they lay no rung
    assert len(dl0) == 4 * (n - 2)
    assert len(g["points"]) == n + 4 * (n - 2)


@needs_partrac
def test_an_unstitched_inlet_repeats_itself(tmp_path):
    """Without inject_edges the generations are not joined to each other, so
    the inlet reappears as it is at each injection and the mesh stays a set of
    separate curves, each spanning the pipe."""
    n = 41
    t, g = lines(tmp_path, ["init_mode=uniform_x", "Nrw=%d" % n, "T=0.2",
                            "dump_intv=0.2", "inject_edges=false"])
    assert "faces" not in g
    dl0 = g["dl0"].ravel()
    assert len(dl0) == 5 * (n - 1)                     # five copies of the inlet
    assert len(g["points"]) == 5 * n
    assert dl0.sum() == pytest.approx(5 * 2 * R, rel=1e-12)


@needs_partrac
def test_a_cleared_sheet_inlet_draws_lines_too(tmp_path):
    """What an inlet traces depends on what survives the clear, not on the
    init_mode: a sheet is refused as an inlet, but a cleared sheet is a cloud
    like any other and draws streaklines."""
    t, g = lines(tmp_path, ["init_mode=sheet_xy", "La=0.5", "Lb=0.5",
                            "ds_init=0.1", "Nrw=41", "T=0.2", "dump_intv=0.2",
                            "clear_initial_edges=true"])
    assert "faces" not in g
    assert g["dl0"].min() > 0


@needs_partrac
def test_refining_faster_than_injecting(tmp_path):
    """Before the first injection every edge of the strip is an inlet edge, and
    refining one must split the inlet template with it, so later generations are
    laid down on more nodes. The swept area then stays between the trapezoidal
    rate of the original inlet and the exact integral, and the rate rises: the
    trapezoidal rule under-reads a concave profile, and refining the template
    walks it toward the exact integral without passing it."""
    # inject_intv = 0.5 is long enough for a generation to stretch past
    # ds_max = 0.1 before the next arrives: the shear pulls the inlet's outer
    # nodes about 0.1 downstream in 0.5
    n = 41
    s = series(tmp_path, ["Nrw=%d" % n, "T=1.5", "dump_intv=0.5",
                          "inject_intv=0.5", "ds_max=0.1", "ds_min=0.02",
                          "refine=true", "refine_intv=0.05",
                          "coarsen=true", "coarsen_intv=0.05"])
    assert len(s) == 3
    assert len(s[-1][2]) > 3000                    # it really refined
    exact = 8. / 3. * U_INF * R
    for t, _, dA0 in s:
        # the first generation is swept on the inlet as laid down, so the lower
        # bound is met exactly there and only exceeded later
        assert dA0.sum() >= sweep_rate(n) * t * (1 - 1e-12)
        assert dA0.sum() <= exact * t
        assert dA0.min() > 0
    # the last generation is swept on a finer inlet than the first
    assert s[-1][2].sum() / s[-1][0] > s[0][2].sum() / s[0][0]


@needs_partrac
def test_an_exit_plane_does_not_cut_the_inlet_away(tmp_path):
    """An exit plane culls material beyond it, but the inlet must outlive the
    cull by a generation or the next injection has nothing to stitch to. With
    half the inlet beyond the plane from the start, the sheet must keep growing."""
    s = series(tmp_path, ["Nrw=41", "T=0.5", "dump_intv=0.1",
                          "inject_intv=0.1", "exit_plane=x", "Ln=0.5",
                          "filter_intv=0.1"])
    assert len(s) == 5
    total = [dA0.sum() for _, _, dA0 in s]
    assert all(b > a for a, b in zip(total, total[1:]))   # still being fed
    assert all(dA0.min() > 0 for _, _, dA0 in s)


@needs_partrac
def test_remeshing_conserves_the_swept_area(tmp_path):
    """Refining and coarsening an injected sheet redistribute dA0 among faces
    but never create or destroy it, so sum dA0 stays on the closed-form sweep.
    The current area moves off dA0 only by the chord error of split nodes."""
    s = series(tmp_path, ["Nrw=41", "T=0.3", "dump_intv=0.1", "ds_max=0.06",
                          "ds_min=0.02", "refine=true", "refine_intv=%g" % INTV,
                          "coarsen=true", "coarsen_intv=%g" % INTV])
    assert len(s) == 3
    assert len(s[-1][2]) > 1000                         # it really refined
    for t, dA, dA0 in s:
        assert dA0.sum() == pytest.approx(sweep_rate(41) * t, rel=1e-12)
        assert dA0.min() > 0
        # split nodes lie on the chord, not the surface: dA moves, dA0 does not
        assert np.max(np.abs(dA / dA0 - 1)) < 0.02


# --- a flow that is three-dimensional and does not stand still ------------------
#
# Blazevski & Haller, Physica D 273-274 (2014) 46-62, eq. (22): the ABC flow with
# A forced in time, A = sqrt(3), B = sqrt(2), C = 1 on [0, 2pi]^3, forcing period
# 2pi. Divergence free, and its steady limit is an exact solution of Euler's
# equation. Nothing here is closed form -- the sheet folds -- so what is checked
# is the invariant that does not depend on the flow: reference area is created
# at the inlet and moved by nothing afterwards.

ABC_ARGS = ["init_mode=uniform_x", "Nrw=41", "Nrw_max=400000",
            "x0=%.14f" % PI, "y0=%.14f" % PI, "z0=%.14f" % PI,
            "inject_intv=0.25", "int_order=2", "dt=0.02", "T=1.5",
            "dump_intv=0.5", "integrate_tau=true", "tau_intv=0.02", "tau_max=0"]
ABC_REMESH = ["ds_max=0.4", "ds_min=0.1", "refine=true", "refine_intv=0.25",
              "coarsen=true", "coarsen_intv=0.25"]


@pytest.fixture(scope="module")
def abc_remeshed(tmp_path_factory):
    """The injected sheet in the unsteady ABC flow, remeshed."""
    return module_run(tmp_path_factory, "abc_remeshed", ABC_ARGS + ABC_REMESH, ABC)[1]


def test_remeshing_a_sheet_in_an_unsteady_flow_moves_no_area(tmp_path, abc_remeshed):
    """In a folding sheet in an unsteady 3-D flow the swept-area law does not
    apply, but dA0 is laid down at the inlet and untouched afterwards, so a
    remeshed run must carry the same total reference area as an unremeshed one."""
    plain = series(tmp_path / "plain", ABC_ARGS, example=ABC)
    remeshed = faces_of(abc_remeshed)
    assert len(plain) == len(remeshed) == 3
    assert len(remeshed[-1][2]) > 2 * len(plain[-1][2])      # it really remeshed
    for (t, _, a), (t2, _, b) in zip(plain, remeshed):
        assert t == t2
        assert b.sum() == pytest.approx(a.sum(), rel=1e-12)
        assert a.min() > 0 and b.min() > 0


def test_an_unsteady_three_dimensional_flow_stretches_the_sheet(abc_remeshed):
    """The counterpart of the parallel-flow case, where dA/dA0 = 1: here faces
    are both stretched and compressed, with dA/dA0 spread over more than a
    factor five, and dA, dA0 and tau stay finite throughout."""
    t, dA, dA0 = faces_of(abc_remeshed)[-1]
    assert (dA / dA0).max() > 2                     # stretched
    assert (dA / dA0).min() < 0.6                   # and compressed
    assert (dA / dA0).max() / (dA / dA0).min() > 5
    last = abc_remeshed[max(abc_remeshed)]
    for name in ("dA", "dA0", "tau"):
        assert np.isfinite(last[name].astype(float)).all(), name


@needs_partrac
@pytest.mark.parametrize("remesh", [
    [],
    ["ds_max=0.2", "ds_min=0.02", "refine=true", "refine_intv=0.05",
     "coarsen=true", "coarsen_intv=0.05"],
])
def test_an_inlet_laid_along_the_flow_sweeps_nothing(tmp_path, remesh):
    """uniform_z runs down the pipe, so the inlet moves along itself and sweeps
    no area. The zero-area faces of one generation are laid down so the next has
    something to stitch to, and removed once they leave the inlet, with or without
    remeshing: faces and nodes never pile up, and tau and the statistics stay
    free of NaN and inf despite the degenerate faces."""
    dumps = run(tmp_path, ["init_mode=uniform_z", "Nrw=41", "T=0.3",
                           "dump_intv=0.1", "stat_intv=0.05",
                           "integrate_tau=true", "tau_intv=0.005",
                           "tau_max=0"] + remesh)
    s = faces_of(dumps)
    assert len(s) == 3
    for t, dA, dA0 in s:
        assert len(dA0) == 40              # one generation's worth, never more
        assert (dA0 == 0).all()            # because nothing was swept
    last = dumps[max(dumps)]
    assert len(last["points"]) == 81
    assert np.isfinite(last["tau"].astype(float)).all()
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
# no area at all. All intervals below are whole fractions of the step length 1/3,
# so steps, injections, dumps and statistics line up.

SINE_ARGS = ["init_mode=uniform_x", "Nrw=41", "Nrw_max=200000",
             "x0=0.5", "y0=0.5", "z0=0.5", "inject_intv=%.10f" % (1. / 15),
             "dt=%.10f" % (1. / 300), "T=2.0", "dump_intv=%.10f" % (1. / 3),
             "stat_intv=%.10f" % (1. / 15), "integrate_tau=true",
             "tau_intv=%.10f" % (1. / 300), "tau_max=0"]
SINE_REMESH = ["ds_max=0.05", "ds_min=0.01", "refine=true",
               "refine_intv=%.10f" % (1. / 15),
               "coarsen=true", "coarsen_intv=%.10f" % (1. / 15)]


@pytest.fixture(scope="module")
def sine_plain(tmp_path_factory):
    """(folder, dumps) of the injected sheet in the 3-D random sine flow, not remeshed."""
    return module_run(tmp_path_factory, "sine_plain", SINE_ARGS, SINE3D)


def test_a_step_that_sweeps_nothing_neither_adds_nor_accumulates(sine_plain):
    """During the steps whose flow is parallel to the inlet, injection lays
    down one generation of flat faces that adds essentially no area and does not
    accumulate; the total dA0 never decreases. The sheet still becomes fully
    three-dimensional, and dA, dA0, tau and the statistics stay finite."""
    d, dumps = sine_plain
    s = faces_of(dumps)
    assert len(s) == 6
    flat = [int((dA0 <= 0).sum()) for _, _, dA0 in s]
    assert set(flat) <= {0, 40}                  # a generation's worth, or none
    assert flat.count(40) == 2                   # the two parallel steps
    total = [dA0.sum() for _, _, dA0 in s]
    assert all(b >= a for a, b in zip(total, total[1:]))
    gain = [b - a for a, b in zip([0.] + total, total)]
    assert min(g for g, f in zip(gain, flat) if f == 0) > 20 * max(
        g for g, f in zip(gain, flat) if f == 40)   # the parallel steps add ~nothing
    last = dumps[max(dumps)]
    for name in ("dA", "dA0", "tau"):
        assert np.isfinite(last[name].astype(float)).all(), name
    p = last["points"]
    assert min(np.ptp(p[:, i]) for i in range(3)) > 0.1     # it really is 3-D
    stats = list(d.rglob("tdata_from_t*.dat"))
    text = stats[0].read_text().lower()
    assert "nan" not in text and "inf" not in text


def test_remeshing_moves_no_area_across_a_parallel_step(tmp_path, sine_plain):
    """The same invariant as in the ABC flow, over a sheet that stops growing
    and starts again, with flat faces created and culled in between: remeshing
    leaves the total dA0 unchanged at every dump."""
    plain = faces_of(sine_plain[1])
    remeshed = series(tmp_path / "remeshed", SINE_ARGS + SINE_REMESH,
                      example=SINE3D)
    for (t, _, a), (t2, _, b) in zip(plain, remeshed):
        assert t == t2
        assert b.sum() == pytest.approx(a.sum(), rel=1e-9)
    # it really remeshed; here coarsening wins, so the face count falls
    assert abs(len(remeshed[-1][2]) - len(plain[-1][2])) > 0.2 * len(plain[-1][2])


@needs_partrac
def test_a_run_names_its_columns_for_the_dimension_it_settles_into(tmp_path):
    """The tdata header is written before the loop, while an injecting run is
    still its inlet curve, but must name the columns for the surface the run
    becomes: A and A0, not s and s0. Otherwise every row reports an area under a
    length label. A0 is 0 before the first injection and follows the sweep rate
    after it; without injection the same inlet stays a strip and reports s, s0."""
    def stats(case, extra):
        run_app(PARTRAC, copy_example(HAGEN, case), BASE, extra, timeout=600)
        return read_stats(case)

    st = stats(tmp_path / "injecting",
               ["Nrw=41", "T=0.15", "stat_intv=0.025", "dump_intv=1e9"])
    assert "A" in st and "A0" in st
    assert "s" not in st and "s0" not in st
    # nothing has been swept before the first injection; the length of the
    # curve about to sweep is not an area
    assert st["A0"][0] == 0.0
    assert st["A0"][2] == pytest.approx(sweep_rate(41) * INTV, rel=1e-12)

    # the same run without injection stays a strip, and says so
    st = stats(tmp_path / "strip",
               ["Nrw=41", "T=0.15", "stat_intv=0.025", "dump_intv=1e9",
                "inject=false"])
    assert "s" in st and "s0" in st
    assert "A" not in st and "A0" not in st
