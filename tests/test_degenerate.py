"""What partrac does at the edges of its domain of validity: meshes coarsened
to almost nothing, point clouds without edges, exit planes that cut a mesh,
and runs with no particles.

The flow is plane Poiseuille, u_z = 1.5 (1 - x^2), which only shears along z.
Several cases set u_inf = 0 so nothing moves: then every face keeps
dA/dA0 = 1 and every face accumulates the same tau, and coarsening must
preserve both exactly. Coarsening hands the reference length or area of a
collapsed element to its neighbours, so the total s0 or dA0 is conserved.
A run whose mesh has lost its dimension, or that has no particles left, is
stopped with a message rather than left to write rows nobody can interpret.
"""

import os

import numpy as np
import pytest

from dumps import all_dumps, read_stats
from paths import REPO, app
from runs import copy_example, run_app

PARTRAC = app("partrac")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
TAYLOR_COUETTE = os.path.join(REPO, "data_example", "taylor_couette", "expr_params.dat")

BASE = ("mode=analytic ds_max=1e9 ds_min=1e-9 refine=false coarsen=false Dm=0 "
        "int_order=1 dt=0.01 stat_intv=0.01 dump_intv=1e9 checkpoint_intv=1e9 "
        "random=false seed=1").split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")


def run(tmp_path, example, extra):
    """Run partrac on a copy of `example` with `extra` overriding BASE; return the process result."""
    return run_app(PARTRAC, copy_example(example, tmp_path), BASE, extra, check=False, timeout=600)


def last_row(tmp_path):
    """The last statistics row under tmp_path, as column name -> value."""
    return {k: v[-1] for k, v in read_stats(tmp_path).items()}


def dump(tmp_path, first=True):
    """The first or the last group of the one dump file under tmp_path, as written."""
    assert len(list(tmp_path.rglob("data_from_t*.h5"))) == 1
    groups = all_dumps(tmp_path, raw=True)
    return groups[min(groups) if first else max(groups)]


@needs_partrac
def test_a_strip_coarsened_to_nothing_keeps_its_end_edges(tmp_path):
    """Coarsening with ds_min longer than the strip collapses every interior edge
    but never an end edge: a tip node carries a single edge, and collapsing it
    would shorten the strip with nowhere to put its ds0. The two end edges
    survive and keep the whole reference length."""
    r = run(tmp_path, POISEUILLE,
            ["init_mode=strip_x", "La=0.5", "x0=0", "y0=0", "z0=0", "Nrw=10",
             "Nrw_max=2000", "ds_min=100", "coarsen=true", "coarsen_intv=0.01",
             "T=0.1"])
    assert r.returncode == 0, r.stdout + r.stderr
    st = last_row(tmp_path)
    assert int(st["Nrw"]) == 3                        # the two end edges survive
    assert float(st["s0"]) == pytest.approx(0.5)      # and keep the whole s0


@needs_partrac
@pytest.mark.parametrize("ds_min", ["0.02", "0.06", "0.2"])
def test_coarsening_conserves_the_reference_length(tmp_path, ds_min):
    """Every edge collapse hands its ds0 to the neighbours it had, so the total
    reference length of the strip stays La however many edges go. Without that
    the stretch s/s0 of a coarsened line would drift with the remeshing."""
    d = tmp_path / ds_min
    r = run(d, POISEUILLE,
            ["init_mode=strip_x", "La=0.5", "x0=0", "y0=0", "z0=0", "Nrw=50",
             "Nrw_max=2000", "ds_min=" + ds_min, "coarsen=true",
             "coarsen_intv=0.01", "T=0.1"])
    assert r.returncode == 0, r.stdout + r.stderr
    st = last_row(d)
    assert float(st["s0"]) == pytest.approx(0.5, rel=1e-9)
    assert 3 <= int(st["Nrw"]) < 50            # it really coarsened


@needs_partrac
@pytest.mark.parametrize("ds_init", ["0.1", "0.05", "0.025"])
def test_coarsening_a_flat_sheet_keeps_every_elongation(tmp_path, ds_init):
    """In a sheet that does not move, every face has dA/dA0 = 1, and coarsening
    must keep that on every surviving face while conserving the total dA0. A face
    whose dA0 went wrong would report a stretch the flow never produced."""
    # u_inf = 0: nothing moves
    d = tmp_path / ds_init
    d.mkdir(parents=True)
    with open(POISEUILLE) as f:
        (d / "expr_params.dat").write_text(f.read().replace("u_inf=1.0", "u_inf=0.0"))
    ds = float(ds_init)

    def areas(ds_min):
        """dA and dA0 of every face in the first dump, coarsened with this ds_min."""
        case = d / ds_min
        run_app(PARTRAC, copy_example(d / "expr_params.dat", case), BASE,
                ["init_mode=sheet_xy", "La=0.5", "Lb=0.5", "ds_init=" + ds_init,
                 "x0=0", "y0=0", "z0=0", "Nrw=100", "Nrw_max=200000",
                 "ds_min=" + ds_min, "coarsen=true", "coarsen_intv=1e9",
                 "T=0.01", "dump_intv=0.01"], timeout=600)
        g = dump(case)
        return g["dA"].ravel(), g["dA0"].ravel()

    dA_ref, _ = areas("1e-12")                      # nothing collapses
    # just above the initial edge length, so many edges collapse
    dA, dA0 = areas("%g" % (1.2 * ds))

    assert len(dA) < 0.6 * len(dA_ref)              # it really coarsened
    assert dA0.min() > 0                            # and never through zero
    assert np.max(np.abs(dA / dA0 - 1)) < 1e-12     # every elongation kept
    assert dA0.sum() == pytest.approx(dA_ref.sum(), rel=1e-9)   # mass conserved


@needs_partrac
def test_coarsening_a_sheet_of_uniform_tau_keeps_it_uniform(tmp_path):
    """In a sheet that does not move, tau is the same on every face, and
    coarsening must redistribute it so it stays the same, including at t = 0
    where every tau is zero. A non-uniform result would be an artefact of the
    remeshing, not of the flow."""
    # u_inf = 0: nothing moves
    d = tmp_path
    with open(POISEUILLE) as f:
        (d / "expr_params.dat").write_text(f.read().replace("u_inf=1.0", "u_inf=0.0"))

    def taus(ds_min):
        """tau of every face in the last dump, coarsened with this ds_min."""
        case = d / ds_min
        run_app(PARTRAC, copy_example(d / "expr_params.dat", case), BASE,
                ["init_mode=sheet_xy", "La=0.5", "Lb=0.5", "ds_init=0.05",
                 "x0=0", "y0=0", "z0=0", "Nrw=100", "Nrw_max=200000",
                 "integrate_tau=true", "tau_intv=0.001", "tau_max=0",
                 "dt=0.001", "T=0.01", "dump_intv=0.01",
                 "ds_min=" + ds_min, "coarsen=true", "coarsen_intv=0.005"], timeout=600)
        return dump(case, first=False)["tau"].ravel()

    tau_ref = taus("1e-12")                         # nothing collapses
    tau = taus("0.06")

    assert len(tau) < 0.6 * len(tau_ref)            # it really coarsened
    assert tau_ref.min() > 0                        # and tau really ran
    assert tau_ref == pytest.approx(tau_ref[0], rel=1e-12)
    assert tau == pytest.approx(tau_ref[0], rel=1e-12)


@needs_partrac
def test_a_run_that_loses_every_particle_stops(tmp_path):
    """When the exit plane removes every particle, the run stops with a non-zero
    exit and says "no particles left", instead of carrying on and writing empty
    or NaN statistics."""
    # at Ln = 0.05 along z the whole strip crosses the plane long before T
    r = run(tmp_path, POISEUILLE,
            ["init_mode=strip_x", "La=0.5", "x0=0", "y0=0", "z0=0", "Nrw=50",
             "Nrw_max=2000", "T=2.0", "stat_intv=0.1", "exit_plane=z",
             "Ln=0.05", "filter_intv=0.05"])
    assert r.returncode != 0
    assert "no particles left" in r.stdout + r.stderr


@needs_partrac
def test_an_exit_plane_on_a_point_cloud_removes_only_what_crosses(tmp_path):
    """A point cloud has no edges, so the rule that a node no surviving edge
    reaches is dead must not apply to it. With the plane out of reach nothing
    crosses, and all 50 particles must survive the filter."""
    r = run(tmp_path, POISEUILLE,
            ["init_mode=uniform_x", "clear_initial_edges=true", "La=0.5",
             "x0=0", "y0=0", "z0=0", "Nrw=50", "Nrw_max=2000", "T=0.05",
             "exit_plane=z", "Ln=100", "filter_intv=0.01"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert int(last_row(tmp_path)["Nrw"]) == 50


@needs_partrac
def test_coarsening_a_point_cloud_keeps_it(tmp_path):
    """A point cloud has no edges to collapse, so a coarsening pass must leave
    every particle in place. Deriving which nodes are dead from the (absent)
    edges would delete the whole cloud."""
    r = run(tmp_path, POISEUILLE,
            ["init_mode=point", "x0=0", "y0=0.1", "z0=0", "Nrw=20",
             "Nrw_max=2000", "coarsen=true", "coarsen_intv=0.01", "T=0.05"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert int(last_row(tmp_path)["Nrw"]) == 20


@needs_partrac
def test_an_exit_plane_on_a_strip_removes_exactly_the_nodes_beyond(tmp_path):
    """On a strip the exit plane removes the nodes beyond it plus whatever they
    orphan, which for a strip is nothing: the node nearest the cut keeps its
    inward edge. The result is still one open strip with consistent indices."""
    def points(case, Ln):
        """Points and edges of the first dump with the exit plane at x = Ln."""
        r = run(tmp_path / case, POISEUILLE,
                ["init_mode=strip_x", "La=0.5", "x0=0", "y0=0", "z0=0",
                 "Nrw=51", "Nrw_max=2000", "T=0", "exit_plane=x",
                 "Ln=%g" % Ln, "filter_intv=0.01"])
        assert r.returncode == 0, r.stdout + r.stderr
        g = dump(tmp_path / case)
        return g["points"], g["edges"]

    whole, _ = points("whole", 1e9)               # strip_x is centred on x0
    # between two nodes, which are 0.01 apart
    kept, edges = points("cut", 0.105)
    assert len(kept) < len(whole)                     # it really cut
    assert np.array_equal(kept, whole[whole[:, 0] <= 0.105])
    assert len(edges) == len(kept) - 1                # still one open strip
    assert edges.max() == len(kept) - 1               # and no dangling index


@needs_partrac
def test_an_exit_plane_on_a_sheet_cuts_it(tmp_path):
    """On a sheet the exit plane removes the nodes beyond it and every face and
    edge that reached one. What remains is a non-empty sheet on the near side,
    with valid face indices and no face that lost its area."""
    def state(case, Ln):
        """Points, faces and dA0 of the first dump with the exit plane at x = Ln."""
        r = run(tmp_path / case, POISEUILLE,
                ["init_mode=sheet_xy", "La=0.5", "Lb=0.5", "ds_init=0.05",
                 "x0=0", "y0=0", "z0=0", "Nrw=100", "Nrw_max=200000", "T=0",
                 "exit_plane=x", "Ln=%g" % Ln, "filter_intv=0.01"])
        assert r.returncode == 0, r.stdout + r.stderr
        g = dump(tmp_path / case)
        return tuple(g[n] for n in ("points", "faces", "dA0"))

    whole, faces_whole, dA0_whole = state("whole", 1e9)
    kept, faces, dA0 = state("cut", 0.1)
    assert (kept[:, 0] <= 0.1).all()                  # nothing beyond survived
    assert 0 < len(kept) < len(whole)                 # and it is not empty
    assert 0 < len(faces) < len(faces_whole)
    assert faces.max() < len(kept) * 3                # no index left dangling
    assert dA0.min() > 0                              # no face lost its area
    assert dA0.sum() < dA0_whole.sum()


@needs_partrac
def test_an_initializer_that_places_nothing_stops(tmp_path):
    """An initializer that places no particle stops the run with "no particles
    left" and a non-zero exit, rather than running an empty simulation."""
    # the axis is inside the solid inner cylinder, so every particle is rejected
    r = run(tmp_path, TAYLOR_COUETTE,
            ["init_mode=point", "x0=0", "y0=0", "z0=0", "Nrw=1", "Nrw_max=1",
             "T=0.02"])
    assert r.returncode != 0
    assert "no particles left" in r.stdout + r.stderr


@needs_partrac
def test_a_single_particle_reports_zero_spread_not_nan(tmp_path):
    """A run with one particle writes finite statistics and a spread of exactly
    zero. The sample variance divides by Nrw - 1, which would otherwise make a
    single-particle run report NaN."""
    r = run(tmp_path, POISEUILLE,
            ["init_mode=point", "x0=0", "y0=0.2", "z0=0", "Nrw=1", "Nrw_max=1",
             "T=0.02"])
    assert r.returncode == 0, r.stdout + r.stderr
    st = last_row(tmp_path)
    for key, value in st.items():
        assert np.isfinite(float(value)), key
    assert float(st["dx2_mean"]) == 0.0
