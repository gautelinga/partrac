"""Refinement against an inlet it is not allowed to leave alone.

Splitting the longest edge of a triangle always shortens something: the median
to a side is at most sqrt(3)/2 of it when that side is the longest, so the
sweep runs out of work. Exempt one edge from splitting and that stops being
true. The sweep then takes the second longest, and the median it raises can
come back as long as the edge it replaced.

For a face with an exempt edge `a`, a split edge `L` and a third edge `s`, the
split reproduces the same face when `s = L/2` and the new median equals `L`.
With the median relation `m^2 = (2a^2 + 2s^2 - L^2)/4` that is

    4 L^2 = 2a^2 - L^2/2   ->   L* = (2/3) a ,

a triangle that regenerates itself forever, halving dA0 each pass until it
denormalises. So exempting the inlet made any `ds_max` below two thirds of an
inlet edge refine until Nrw_max stopped it, leaving a mesh of slivers and a
statistics file of NaN -- with an exit code of zero.

The inlet is exempt for a reason: injection rebuilds `edges_inlet` from a
template each time, so a split the template does not follow is forgotten at the
next injection, and the quad stitched over it spans half an edge. The way out
is to split the template too, which is what these check -- refinement now
terminates for any ds_max, and the inlet it injects gets finer to match.

The other half of the same geometry: a split also raises a median from the new
node to the opposite vertex of each face it divides, and nothing bounds that
edge. `ds_max` governs the edge being divided; the median's length is the
face's business. On a nearly collinear face -- which is what a tight fold looks
like here -- it comes out at nothing, and once it is exactly nothing the faces
on it have two vertices at one point. So coarsening runs whether or not it was
asked for: with it off the threshold drops to what is numerically zero, on the
refinement interval, because refinement is what raises those edges.

Healthy runs never produce one, which is the difficulty in testing it. The
mesh here is welded by hand instead -- a checkpoint is written, one node moved
onto its neighbour, and the run resumed from it.

One consequence has no test here because it is a property of the design rather
than of a run: the template only ever grows. Coarsening refuses to collapse any
edge touching an inlet node, so nothing gives the nodes back.
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
    """The uniform_x inlet spans the pipe, so its edges are this long."""
    return 2 * R / (n - 1)


# a sheet with no inlet at all, so nothing is protected from collapsing
SHEET = ("mode=analytic init_mode=sheet_xy La=1.0 Lb=1.0 ds_init=0.15 "
         "x0=%.14f y0=%.14f z0=%.14f Nrw=100 Nrw_max=200000 inject=false "
         "refine=true ds_max=0.25 coarsen=false ds_min=1e-9 Dm=0 int_order=2 "
         "dt=0.05 stat_intv=1e9 random=false seed=1" % (PI, PI, PI)).split()


# the stepper advances a position rather than a time: Ln bounds the loop, T is
# the local-time cutoff per node, and stat_intv=1e9 keeps it off
SHEET_SSS = ("mode=analytic init_mode=sheet_xy La=0.5 Lb=0.5 ds_init=0.1 "
             "Nrw=100 Nrw_max=20000 refine=true ds_max=0.4 coarsen=false "
             "ds_min=1e-9 Dm=0 int_order=1 dt=0.005 T=1e9 dx_max=0.1 "
             "dxn=0.05 stat_intv=1e9 random=false seed=1").split()


def run_sss(tmp_path, extra):
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(PLANE, tmp_path / "expr_params.dat")
    argv = {}
    for a in SHEET_SSS + extra:
        argv[a.split("=")[0]] = a
    return subprocess.run([SSS, str(tmp_path / "expr_params.dat")]
                          + list(argv.values()),
                          capture_output=True, text=True, timeout=900)


def run(tmp_path, extra, example=HAGEN, base=None):
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
    e = np.unique(np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]],
                                     faces[:, [2, 0]]]), axis=1), axis=0)
    return np.linalg.norm(points[e[:, 0]] - points[e[:, 1]], axis=1)


# --- the runaway itself ---------------------------------------------------------

@needs_partrac
@pytest.mark.parametrize("n", [21, 41, 81])
@pytest.mark.parametrize("mult", [0.95, 0.5])
def test_refinement_terminates_below_the_self_similar_floor(tmp_path, n, mult):
    # every one of these used to run to Nrw_max, whatever Nrw_max was
    ds_max = mult * 2. / 3. * spacing(n)
    r = run(tmp_path, ["Nrw=%d" % n, "ds_max=%.10f" % ds_max])
    assert r.returncode == 0, r.stdout + r.stderr
    s = dumps(tmp_path)
    assert s, "nothing was ever injected"
    for t, points, faces, dA, dA0 in s:
        assert len(points) < NRW_MAX / 4      # it settled, it did not run out
        assert dA0.min() > 1e-300             # no dA0 halved into a denormal
        assert (dA > 0).all()                 # and no collapsed face
    assert np.isfinite(s[-1][3]).all()


@needs_partrac
def test_the_mesh_does_not_depend_on_the_node_budget(tmp_path):
    # the sharpest way to say the sweep terminates. A refinement that stops
    # because it has nothing left to split gives the same mesh whatever the
    # budget; one that stops because it ran out gives back the budget, which
    # is what both of these used to do.
    ds_max = 0.95 * 2. / 3. * spacing(41)          # just inside the old cliff
    sizes = []
    for cap in (5000, 400000):
        case = tmp_path / str(cap)
        r = run(case, ["Nrw=41", "ds_max=%.10f" % ds_max, "Nrw_max=%d" % cap])
        assert r.returncode == 0, r.stdout + r.stderr
        sizes.append(len(dumps(case)[-1][1]))
    assert sizes[0] == sizes[1]
    assert sizes[0] < 5000                          # neither budget was met


@needs_partrac
def test_the_node_count_grows_smoothly_across_the_old_threshold(tmp_path):
    # the floor for a 41-node inlet is at ds_max = 0.0333, and the old cliff
    # was vertical: 445 nodes one side, Nrw_max the other
    counts = []
    for ds_max in (0.0345, 0.0335, 0.0325, 0.0315):
        r = run(tmp_path / ("ds%g" % ds_max), ["Nrw=41", "ds_max=%g" % ds_max])
        assert r.returncode == 0, r.stdout + r.stderr
        counts.append(len(dumps(tmp_path / ("ds%g" % ds_max))[-1][1]))
    assert all(b >= a for a, b in zip(counts, counts[1:]))   # monotone
    assert counts[-1] < 2 * counts[0]                        # and no cliff


@needs_partrac
def test_refinement_leaves_no_edge_longer_than_ds_max(tmp_path):
    # dumping on the refinement interval catches the mesh straight after a
    # sweep, and nothing is exempt from it any more. The inlet edges used to
    # be, so this is the assertion the old code could not meet.
    ds_max = 0.04                                    # below the 0.05 spacing
    r = run(tmp_path, ["Nrw=41", "ds_max=%g" % ds_max, "dump_intv=0.05"])
    assert r.returncode == 0, r.stdout + r.stderr
    s = dumps(tmp_path)
    assert len(s) > 1
    for t, points, faces, _, _ in s:
        if t == 0:
            continue                                 # refinement skips step 0
        assert edge_lengths(points, faces).max() <= ds_max * (1 + 1e-9), t


# --- the template that makes it possible ----------------------------------------

def template_size(tmp_path):
    f = list(tmp_path.rglob("positions_inj.pos"))
    assert len(f) == 1
    return len([l for l in f[0].read_text().splitlines() if l.strip()])


@needs_partrac
def test_the_template_refines_only_when_the_inlet_is_too_coarse(tmp_path):
    # ds_max above the inlet spacing leaves the template alone, which is what
    # keeps the closed-form sweep rate exact for the runs that assert it
    coarse = tmp_path / "coarse"
    r = run(coarse, ["Nrw=41", "ds_max=0.1", "checkpoint_intv=0.3"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert template_size(coarse) == 41

    fine = tmp_path / "fine"
    r = run(fine, ["Nrw=41", "ds_max=0.02", "checkpoint_intv=0.3"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert template_size(fine) > 41          # the inlet was split, and kept


@needs_partrac
def test_a_template_too_coarse_for_ds_max_is_reported(tmp_path):
    # the run is fine, but it is not the run that was asked for: the injected
    # curve will not stay at the resolution it was given
    r = run(tmp_path / "warned", ["Nrw=41", "ds_max=0.02"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert "below the injection template" in r.stdout + r.stderr

    r = run(tmp_path / "quiet", ["Nrw=41", "ds_max=0.1"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert "below the injection template" not in r.stdout + r.stderr


@needs_partrac
def test_a_refined_template_survives_a_restart(tmp_path):
    # pos_inj and edges_inj ride in the checkpoint, so a resumed run injects
    # the curve as it had been refined, not as it started
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


# --- the edges refinement raises and does not control --------------------------

def weld_two_nodes(checkpoint):
    """Move one end of the first edge onto the other, so it has zero length."""
    pos, edges = checkpoint / "positions.pos", checkpoint / "edges.edge"
    x = [[float(v) for v in line.split()]
         for line in pos.read_text().splitlines() if line.strip()]
    first = edges.read_text().splitlines()[0].split()
    i, j = int(first[0]), int(first[1])
    x[j] = x[i]
    pos.write_text("".join("%.17g %.17g %.17g\n" % tuple(row) for row in x))
    return i, j


def degeneracies(tmp_path):
    """(duplicate points, zero-length edges, faces with no area) at the end."""
    t, points, faces, dA, dA0 = dumps(tmp_path)[-1]
    duplicates = len(points) - len(np.unique(points, axis=0))
    zero_edges = int((edge_lengths(points, faces) == 0).sum())
    dead = int(((dA == 0) & (dA0 > 0)).sum())
    return duplicates, zero_edges, dead


@needs_partrac
def test_an_edge_of_no_length_is_collapsed_with_coarsening_off(tmp_path):
    # No run reaches this state any more, so the mesh is welded by hand: take a
    # checkpoint, move a node onto its neighbour, resume. The control resumes
    # the same welded mesh with refine_intv = 0, which turns the cleanup off
    # with it -- without that, this test would pass on a mesh that was never
    # degenerate in the first place.
    out = {}
    for label, refine_intv in (("cleaned", "0.1"), ("control", "0")):
        case = tmp_path / label
        # dump_intv = 0 is off entirely; anything positive still dumps at
        # step 0, and the resumed run needs to be the only dump in the folder
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
    assert out["cleaned"] == (0, 0, 0)      # and does not when the cleanup does


@needs_sss
def test_an_edge_of_no_length_is_collapsed_in_the_space_stepper(tmp_path):
    # the same claim as above for the other app that refines. It kept the
    # cleanup behind the coarsen flag, so a refining run with coarsening off
    # carried the weld to the end
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
