"""Refinement of injected sheets, their inlet template, and zero-length edge cleanup.

partrac injects a uniform_x line of Nrw nodes across a Hagen-Poiseuille pipe
(radius R = 1, so inlet spacing a = 2R/(Nrw - 1)) and refines the sheet it
sweeps out with longest-edge splitting. Splitting the longest edge of a face
always shortens something, since the median to the longest side is at most
sqrt(3)/2 of it. If one edge `a` of a face could not be split, the sweep would
split the second longest edge `L` instead, and with third edge `s = L/2` the
median `m^2 = (2a^2 + 2s^2 - L^2)/4` equals `L` when

    4 L^2 = 2a^2 - L^2/2   ->   L* = (2/3) a ,

so the face would regenerate itself forever for any ds_max below (2/3) a.
Refinement therefore splits the inlet template as well, and the tests check
that it terminates for any ds_max, that no edge is left longer than ds_max,
and that the refined template is kept and checkpointed.

A split also raises a median from the new node to the opposite vertex, and
ds_max does not bound its length; on a nearly collinear face it can be zero,
leaving two vertices at one point. Zero-length edges are therefore collapsed
on the refinement interval whether or not coarsening is enabled. Healthy runs
do not produce them, so the tests weld two nodes of a checkpoint by hand and
resume from it.
"""

import os
import shutil
import subprocess

import h5py
import numpy as np
import pytest

from paths import REPO, app

PARTRAC = app("partrac")
SSS = app("static_space_stepper")
HAGEN = os.path.join(REPO, "data_example", "hagen_poiseuille", "expr_params.dat")
ABC = os.path.join(REPO, "data_example", "abc_flow_unsteady", "expr_params.dat")
PLANE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
PI = 3.14159265358979

R, DT, NRW_MAX = 1.0, 0.005, 400000

BASE = ("mode=analytic init_mode=uniform_x x0=0 y0=0 z0=0 Nrw_max=%d "
        "inject=true inject_edges=true T_inject=1e10 inject_intv=0.05 "
        "refine=true refine_intv=0.05 coarsen=false ds_min=1e-9 Dm=0 "
        "int_order=1 dt=%g random=false seed=1 stat_intv=0.05 dump_intv=0.1 "
        "checkpoint_intv=1e9 T=0.3" % (NRW_MAX, DT)).split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")
needs_sss = pytest.mark.skipif(not os.path.exists(SSS),
                               reason="static_space_stepper is not built")


def spacing(n):
    """Edge length of a uniform_x inlet of n nodes spanning the pipe diameter."""
    return 2 * R / (n - 1)


# a sheet with no inlet, so no node is protected from collapsing
SHEET = ("mode=analytic init_mode=sheet_xy La=1.0 Lb=1.0 ds_init=0.15 "
         "x0=%.14f y0=%.14f z0=%.14f Nrw=100 Nrw_max=200000 inject=false "
         "refine=true ds_max=0.25 coarsen=false ds_min=1e-9 Dm=0 int_order=2 "
         "dt=0.05 stat_intv=1e9 random=false seed=1" % (PI, PI, PI)).split()


# the stepper advances in position, not time: Ln bounds the loop, T is the
# per-node local-time cutoff, and stat_intv=1e9 turns statistics off
SHEET_SSS = ("mode=analytic init_mode=sheet_xy La=0.5 Lb=0.5 ds_init=0.1 "
             "Nrw=100 Nrw_max=20000 refine=true ds_max=0.4 coarsen=false "
             "ds_min=1e-9 Dm=0 int_order=1 dt=0.005 T=1e9 dx_max=0.1 "
             "dxn=0.05 stat_intv=1e9 random=false seed=1").split()


def run_sss(tmp_path, extra):
    """Run static_space_stepper on plane Poiseuille with SHEET_SSS plus overrides."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(PLANE, tmp_path / "expr_params.dat")
    argv = {}
    for a in SHEET_SSS + extra:
        argv[a.split("=")[0]] = a
    return subprocess.run([SSS, str(tmp_path / "expr_params.dat")]
                          + list(argv.values()),
                          capture_output=True, text=True, timeout=900)


def run(tmp_path, extra, example=HAGEN, base=None):
    """Run partrac on `example` with `base` (default BASE) plus overrides."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(example, tmp_path / "expr_params.dat")
    argv = {}
    for a in (BASE if base is None else base) + extra:
        argv[a.split("=")[0]] = a
    return subprocess.run([PARTRAC, str(tmp_path / "expr_params.dat")]
                          + list(argv.values()),
                          capture_output=True, text=True, timeout=900)


def dumps(tmp_path):
    """Every dump that has faces, as (t, points, faces, dA, dA0)."""
    f = list(tmp_path.rglob("data_from_t*.h5"))
    assert len(f) == 1
    h = h5py.File(f[0], "r")
    out = []
    for k in sorted(h.keys(), key=float):
        if "faces" not in h[k]:
            continue
        out.append((float(k),
                    np.array(h[k + "/points"]),
                    np.array(h[k + "/faces"])[:, :3],
                    np.array(h[k + "/dA"]).ravel(),
                    np.array(h[k + "/dA0"]).ravel()))
    return out


def edge_lengths(points, faces):
    """Length of every distinct edge of the triangle mesh."""
    e = np.unique(np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]],
                                     faces[:, [2, 0]]]), axis=1), axis=0)
    return np.linalg.norm(points[e[:, 0]] - points[e[:, 1]], axis=1)


# --- termination ----------------------------------------------------------------

@needs_partrac
@pytest.mark.parametrize("n", [21, 41, 81])
@pytest.mark.parametrize("mult", [0.95, 0.5])
def test_refinement_terminates_below_the_self_similar_floor(tmp_path, n, mult):
    """With ds_max below the self-similar length (2/3) a, refinement settles well
    short of Nrw_max with finite, positive face areas. A non-terminating sweep
    would fill the node budget with slivers, halve dA0 into denormals and write
    NaN statistics while still exiting with code zero."""
    ds_max = mult * 2. / 3. * spacing(n)
    r = run(tmp_path, ["Nrw=%d" % n, "ds_max=%.10f" % ds_max])
    assert r.returncode == 0, r.stdout + r.stderr
    s = dumps(tmp_path)
    assert s, "nothing was ever injected"
    for t, points, faces, dA, dA0 in s:
        assert len(points) < NRW_MAX / 4      # settled, far from the budget
        assert dA0.min() > 1e-300             # no dA0 halved into a denormal
        assert (dA > 0).all()                 # no collapsed face
    assert np.isfinite(s[-1][3]).all()


@needs_partrac
def test_the_mesh_does_not_depend_on_the_node_budget(tmp_path):
    """The final mesh is the same size for Nrw_max = 5000 and 400000. A sweep
    that stops because nothing is left to split is independent of the budget; one
    that stops because it ran out would return the budget itself."""
    ds_max = 0.95 * 2. / 3. * spacing(41)          # just below (2/3) a
    sizes = []
    for cap in (5000, 400000):
        case = tmp_path / str(cap)
        r = run(case, ["Nrw=41", "ds_max=%.10f" % ds_max, "Nrw_max=%d" % cap])
        assert r.returncode == 0, r.stdout + r.stderr
        sizes.append(len(dumps(case)[-1][1]))
    assert sizes[0] == sizes[1]
    assert sizes[0] < 5000                          # neither budget was reached


@needs_partrac
def test_the_node_count_grows_smoothly_across_the_old_threshold(tmp_path):
    """Lowering ds_max across (2/3) a raises the node count monotonically and
    gradually. A jump at that threshold would mean the self-similar faces are
    still regenerating."""
    # (2/3) a = 0.0333 for a 41-node inlet; the values bracket it
    counts = []
    for ds_max in (0.0345, 0.0335, 0.0325, 0.0315):
        r = run(tmp_path / ("ds%g" % ds_max), ["Nrw=41", "ds_max=%g" % ds_max])
        assert r.returncode == 0, r.stdout + r.stderr
        counts.append(len(dumps(tmp_path / ("ds%g" % ds_max))[-1][1]))
    assert all(b >= a for a, b in zip(counts, counts[1:]))
    assert counts[-1] < 2 * counts[0]


@needs_partrac
def test_refinement_leaves_no_edge_longer_than_ds_max(tmp_path):
    """Straight after a sweep no edge, the inlet edges included, is longer than
    ds_max. This is the resolution guarantee the user asks for with ds_max."""
    # dumping on the refinement interval sees the mesh right after each sweep
    ds_max = 0.04                                    # below the 0.05 inlet spacing
    r = run(tmp_path, ["Nrw=41", "ds_max=%g" % ds_max, "dump_intv=0.05"])
    assert r.returncode == 0, r.stdout + r.stderr
    s = dumps(tmp_path)
    assert len(s) > 1
    for t, points, faces, _, _ in s:
        if t == 0:
            continue                                 # no refinement at step 0
        assert edge_lengths(points, faces).max() <= ds_max * (1 + 1e-9), t


# --- inlet template -------------------------------------------------------------

def template_size(tmp_path):
    """Number of nodes in the checkpointed injection template."""
    f = list(tmp_path.rglob("positions_inj.pos"))
    assert len(f) == 1
    return len([l for l in f[0].read_text().splitlines() if l.strip()])


@needs_partrac
def test_the_template_refines_only_when_the_inlet_is_too_coarse(tmp_path):
    """The template keeps its 41 nodes when ds_max exceeds the inlet spacing and
    grows when it does not. Leaving a fine enough template untouched keeps the
    injected sheet, and the closed-form sweep rates other tests rely on, exact."""
    coarse = tmp_path / "coarse"
    r = run(coarse, ["Nrw=41", "ds_max=0.1", "checkpoint_intv=0.3"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert template_size(coarse) == 41

    fine = tmp_path / "fine"
    r = run(fine, ["Nrw=41", "ds_max=0.02", "checkpoint_intv=0.3"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert template_size(fine) > 41          # the inlet was split and kept


@needs_partrac
def test_a_template_too_coarse_for_ds_max_is_reported(tmp_path):
    """A warning is printed when ds_max is below the injection template spacing,
    and not otherwise. The run still succeeds, but the injected curve is not at
    the resolution the user gave, so the user should be told."""
    r = run(tmp_path / "warned", ["Nrw=41", "ds_max=0.02"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert "below the injection template" in r.stdout + r.stderr

    r = run(tmp_path / "quiet", ["Nrw=41", "ds_max=0.1"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert "below the injection template" not in r.stdout + r.stderr


@needs_partrac
def test_a_refined_template_survives_a_restart(tmp_path):
    """pos_inj and edges_inj are checkpointed, so a resumed run injects the curve
    as refined, not as it started. Otherwise every injection after a restart
    would stitch quads over edges the mesh has already split."""
    stop = 0.15
    case = tmp_path / "case"
    r = run(case, ["Nrw=41", "ds_max=0.02", "T=%g" % stop,
                   "checkpoint_intv=%g" % stop])
    assert r.returncode == 0, r.stdout + r.stderr
    grown = template_size(case)
    assert grown > 41
    checkpoint = list(case.rglob("edges_inj.edge"))
    assert len(checkpoint) == 1
    r = run(case, ["Nrw=41", "ds_max=0.02", "T=0.3", "checkpoint_intv=0.3",
                   "restart_folder=" + str(checkpoint[0].parent.parent)])
    assert r.returncode == 0, r.stdout + r.stderr
    assert template_size(case) >= grown      # kept what it had, and may add


# --- zero-length edges ----------------------------------------------------------

def weld_two_nodes(checkpoint):
    """Move one end of the checkpoint's first edge onto the other; return (i, j)."""
    pos, edges = checkpoint / "positions.pos", checkpoint / "edges.edge"
    x = [[float(v) for v in line.split()]
         for line in pos.read_text().splitlines() if line.strip()]
    first = edges.read_text().splitlines()[0].split()
    i, j = int(first[0]), int(first[1])
    x[j] = x[i]
    pos.write_text("".join("%.17g %.17g %.17g\n" % tuple(row) for row in x))
    return i, j


def degeneracies(tmp_path):
    """(duplicate points, zero-length edges, zero-area faces) in the last dump."""
    t, points, faces, dA, dA0 = dumps(tmp_path)[-1]
    duplicates = len(points) - len(np.unique(points, axis=0))
    zero_edges = int((edge_lengths(points, faces) == 0).sum())
    dead = int(((dA == 0) & (dA0 > 0)).sum())
    return duplicates, zero_edges, dead


@needs_partrac
def test_an_edge_of_no_length_is_collapsed_with_coarsening_off(tmp_path):
    """In partrac, a zero-length edge is collapsed on the refinement interval even
    with coarsen=false. Left in place, its faces have two vertices at one point
    and zero area, which corrupts the area and stretching statistics."""
    # the control resumes the same welded mesh with refine_intv=0, which also
    # turns the cleanup off; it proves the weld is really there to be removed
    out = {}
    for label, refine_intv in (("cleaned", "0.1"), ("control", "0")):
        case = tmp_path / label
        # dump_intv=0 disables dumps (any positive value dumps at step 0), so
        # the resumed run writes the only dump in the folder
        r = run(case, ["refine_intv=0.1", "T=0.2", "checkpoint_intv=0.2",
                       "dump_intv=0"], example=ABC, base=SHEET)
        assert r.returncode == 0, r.stdout + r.stderr
        checkpoint = list(case.rglob("edges.edge"))
        assert len(checkpoint) == 1
        weld_two_nodes(checkpoint[0].parent)
        r = run(case, ["refine_intv=" + refine_intv, "T=0.4",
                       "checkpoint_intv=1e9", "dump_intv=0.4",
                       "restart_folder=" + str(checkpoint[0].parent.parent)],
                example=ABC, base=SHEET)
        assert r.returncode == 0, r.stdout + r.stderr
        out[label] = degeneracies(case)

    assert out["control"] == (1, 1, 1)      # the weld survives when nothing runs
    assert out["cleaned"] == (0, 0, 0)      # and is removed when the cleanup runs


@needs_sss
def test_an_edge_of_no_length_is_collapsed_in_the_space_stepper(tmp_path):
    """The same zero-length edge cleanup holds in static_space_stepper, the other
    app that refines, with coarsening off. A refining march must not carry
    coincident vertices to its end."""
    out = {}
    for label, refine_intv in (("cleaned", "0.1"), ("control", "0")):
        case = tmp_path / label
        r = run_sss(case, ["refine_intv=0.1", "Ln=0.1", "checkpoint_intv=0.1",
                           "dump_intv=0"])
        assert r.returncode == 0, r.stdout + r.stderr
        checkpoint = list(case.rglob("edges.edge"))
        assert len(checkpoint) == 1
        weld_two_nodes(checkpoint[0].parent)
        r = run_sss(case, ["refine_intv=" + refine_intv, "Ln=0.3",
                           "checkpoint_intv=1e9", "dump_intv=0.3",
                           "restart_folder=" + str(checkpoint[0].parent.parent)])
        assert r.returncode == 0, r.stdout + r.stderr
        out[label] = degeneracies(case)

    assert out["control"] == (1, 1, 1)
    assert out["cleaned"] == (0, 0, 0)
