"""The filaments app: its two edge resizes, a resumed run, and reinjection of
stuck edges.

resize=doublings: an edge longer than the target is cut by 2^n, the least n
that brings it to the target or below; its reference length dl0 stays and n is
counted, so logelong = log(dl/dl0) + doublings log 2 is the total stretch.
resize=rescale brings the edge to the target and scales its reference with it,
so dl/dl0 carries the stretch instead. Pairs across plane Poiseuille,
u_z = 1.5 (1 - x^2), are sheared steadily, so every resize step has edges to
shorten.

outside=reinject is for edges stuck because the field is underresolved, not a
boundary condition: a node whose step the field carries into a solid is
declined, and on a steady field declined again at every step. The case is a
voxel channel along the x-y diagonal, five lattice nodes across, whose wall
nodes still carry flow; at dt = 2 a step covers more of it than the lattice
resolves. A diffusive move out of the domain is not stuck, and the schema
refuses to reinject it.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

from paths import REPO, app

FILAMENTS = app("filaments")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

# the resize target, and the initial edge length, so edges exceed it as soon as they stretch
TARGET = 0.05
BASE = ("mode=analytic init_mode=pairs_xz_x x0=0 y0=0 z0=0 Nrw=200 Nrw_max=400 int_order=1 "
        "ds_init=%g Dm=0 scheme=RK4 dt=0.01 dump_intv=0.1 stat_intv=0.1 checkpoint_intv=1e9 "
        "resize=doublings resize_target=ds_init resize_intv=0.01 random=false seed=1" % TARGET).split()

needs_filaments = pytest.mark.skipif(not os.path.exists(FILAMENTS), reason="filaments is not built")


def run(d, extra):
    """Run filaments on plane Poiseuille in d with `extra` overriding BASE; return d."""
    d.mkdir(parents=True, exist_ok=True)
    shutil.copy(POISEUILLE, d / "expr_params.dat")
    argv = {}
    for a in BASE + extra.split():
        argv[a.split("=")[0]] = a
    r = subprocess.run([FILAMENTS, str(d / "expr_params.dat")] + list(argv.values()),
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    return d


def groups(d):
    """time -> {dataset: array}, over every dump file under d."""
    h5py = pytest.importorskip("h5py")
    out = {}
    for f in sorted(d.rglob("data_from_t*.h5")):
        with h5py.File(f, "r") as h:
            for g in h:
                out[round(float(g), 9)] = {k: np.array(h[g][k]) for k in h[g]}
    return out


def statistics(d):
    """The statistics file's column names, and its rows as column -> value."""
    lines = [l for l in next(d.rglob("tdata_from_t*.dat")).read_text().splitlines() if l.strip()]
    head = [h for h in lines[0].lstrip("# ").rstrip().split("\t") if h.strip()]
    return head, [dict(zip(head, map(float, l.split()))) for l in lines[1:]]


@needs_filaments
def test_resizing_by_halvings_keeps_the_stretch(tmp_path):
    """With resize=doublings every dumped edge is at or below the target, dl0 is
    unchanged, the halving count never decreases, and logelong equals
    log(dl/dl0) + doublings log 2. The statistics report the mean of that
    logelong, so the stretch is not lost when edges are shortened."""
    d = run(tmp_path / "run", "T=0.5")
    g = groups(d)
    times = sorted(g)
    dl0 = g[times[0]]["dl0"][:, 0]
    prev = np.zeros(len(dl0))
    for t in times:
        e = g[t]
        dl, n, logelong = e["dl"][:, 0], e["doublings"][:, 0], e["logelong"][:, 0]
        assert np.array_equal(e["dl0"][:, 0], dl0)                 # the reference stays
        assert np.all(dl <= TARGET * (1 + 1e-12))                  # resized before the dump
        assert np.all(n >= prev)                                   # halvings only accumulate
        prev = n
        assert np.allclose(logelong, np.log(dl / dl0) + n * np.log(2), rtol=0, atol=1e-12)
    assert prev.max() >= 1                                         # something was halved
    head, rows = statistics(d)
    assert head[-5:] == ["n_edges", "elong_mean", "elong2_mean", "logelong_mean", "logelong_var"]
    assert rows[-1]["n_edges"] == len(dl0)
    # the statistics are printed to six significant digits
    assert np.isclose(rows[-1]["logelong_mean"], logelong.mean(), rtol=5e-6, atol=1e-9)


@needs_filaments
def test_resizing_by_rescaling_scales_the_reference(tmp_path):
    """With resize=rescale edges stay at or below the target, no halving count or
    logelong is written, and the reference length shrinks with the edge so that
    dl/dl0 carries the stretch."""
    d = run(tmp_path / "run", "T=0.5 resize=rescale")
    g = groups(d)
    times = sorted(g)
    for t in times:
        e = g[t]
        assert "doublings" not in e and "logelong" not in e
        assert np.all(e["dl"][:, 0] <= TARGET * (1 + 1e-12))
    assert np.any(g[times[-1]]["dl0"][:, 0] < g[times[0]]["dl0"][:, 0])


@needs_filaments
def test_a_resumed_run_keeps_its_halvings(tmp_path):
    """A run stopped and resumed from its checkpoint dumps exactly what an
    uninterrupted run dumps. The checkpoint carries each edge's halvings; without
    them a resumed run's logelong would lose the stretch halved away before the
    stop."""
    cont = run(tmp_path / "cont", "T=0.4")
    split = run(tmp_path / "split", "T=0.19")
    # the final checkpoint is written one step past T: t = 0.2, step 20
    folder = os.path.dirname(os.path.dirname(next(split.rglob("Checkpoints/positions.pos"))))
    run(tmp_path / "split", "T=0.4 restart_folder=" + folder)
    a, b = groups(cont), groups(split)
    for t in (0.3, 0.4):
        assert t in b, "the resumed run has no dump at t = %g (it dumped at %s)" % (t, sorted(b))
        assert set(a[t]) == set(b[t])
        for k in a[t]:
            assert np.array_equal(a[t][k], b[t][k]), "%s differs at t = %g after a restart" % (k, t)


N = 32          # lattice nodes along each axis, spacing 1, periodic
HALF = 2.0      # the channel's half-width in lattice units: five nodes across


def voxel_channel(d):
    """Write a felbm case in d, a channel along the x-y diagonal, uniform in z,
    with velocity parallel to it; return the fluid mask, indexed (z, y, x).

    A point is inside where its nearest lattice node is fluid, so the wall is a
    staircase; the parabola vanishes beyond that staircase, so the fluid nodes
    next to it still carry flow."""
    h5py = pytest.importorskip("h5py")
    d.mkdir(parents=True, exist_ok=True)
    iz, iy, ix = np.meshgrid(np.arange(N), np.arange(N), np.arange(N), indexing="ij")
    dist = np.abs((ix - iy + N // 2) % N - N // 2) / np.sqrt(2.)
    fluid = dist < HALF
    # the parabola reaches zero at 1.5 HALF, past the last fluid node
    u = np.where(fluid, 1. - (dist / (1.5 * HALF)) ** 2, 0.) / np.sqrt(2.)
    # written (z, y, x) on a cube, as the loader reads x fastest
    with h5py.File(d / "output_is_solid.h5", "w") as f:
        f.create_dataset("is_solid", data=(~fluid).astype(np.int32))
    # two identical snapshots: a steady field
    for k in range(2):
        with h5py.File(d / ("output_%d.h5" % k), "w") as f:
            for name, a in (("u_x", u), ("u_y", u), ("u_z", np.zeros_like(u)),
                            ("density", np.ones_like(u)), ("pressure", np.zeros_like(u))):
                f.create_dataset(name, data=a)
    (d / "timestamps.dat").write_text("0\toutput_0.h5\n1000\toutput_1.h5\n")
    (d / "felbm_params.dat").write_text("timestamps=timestamps.dat\nis_solid_file=output_is_solid.h5\n")
    return fluid


# at dt = 2 a step covers more of the channel than the lattice resolves, so nodes stick
CHANNEL = ("mode=felbm scheme=RK4 int_order=1 Dm=0 dt=2 T=100 Nrw=2000 Nrw_max=4000 "
           "init_mode=pairs_xy_xy ds_init=0.5 x0=16 y0=16 z0=16 dump_intv=10 stat_intv=10 "
           "checkpoint_intv=1e9 resize=doublings resize_target=ds_init random=false seed=1")


def channel_args(extra):
    """CHANNEL's arguments with extra's keys taking precedence."""
    argv = {}
    for a in (CHANNEL + " " + extra).split():
        argv[a.split("=")[0]] = a
    return list(argv.values())


def run_channel(d, extra):
    """Run filaments on the voxel channel; return the fluid mask, stdout and the dumps."""
    fluid = voxel_channel(d)
    r = subprocess.run([FILAMENTS, str(d / "felbm_params.dat")] + channel_args(extra),
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    return fluid, r.stdout, groups(d)


def by_id(group):
    """A dump group's points, ordered by particle id."""
    order = np.argsort(group["id"][:, 0])
    return group["points"][order]


@needs_filaments
def test_an_edge_stuck_in_an_underresolved_channel_is_reinjected(tmp_path):
    """With outside=reinject no edge is dropped, no node ever sits in the solid,
    nodes stay in their initial z plane, and no node stays stuck between the
    last two dumps. Left alone, the same case has nodes declined at every step,
    which would freeze those edges for the rest of the run."""
    # control: without reinjection the case does stick
    _, out, g = run_channel(tmp_path / "ignore", "outside=ignore")
    t = sorted(g)
    stuck = np.all(by_id(g[t[-2]]) == by_id(g[t[-1]]), axis=1)
    assert "Some nodes are outside." in out
    assert stuck.sum() > 0, "no node stuck: the case no longer tests reinjection"

    fluid, out, g = run_channel(tmp_path / "reinject", "outside=reinject")
    t = sorted(g)
    assert "Some nodes are outside." in out
    z0 = by_id(g[t[0]])[:, 2]
    for s in t:
        p = by_id(g[s])
        assert len(p) == len(z0)                                   # edges are moved, never dropped
        idx = np.mod(np.round(p), N).astype(int)
        assert fluid[idx[:, 2], idx[:, 1], idx[:, 0]].all(), "a node sits in the solid at t = %g" % s
        # the flow has no z, and reinjection draws in init_mode's third token, xy
        assert np.array_equal(p[:, 2], z0), "a node moved in z by t = %g" % s
    # nothing stays stuck: every edge that stuck was moved on
    assert not np.any(np.all(by_id(g[t[-2]]) == by_id(g[t[-1]]), axis=1))


@needs_filaments
def test_reinjection_is_refused_with_diffusion(tmp_path):
    """outside=reinject together with diffusion is refused before the run starts.
    A diffusive move out of the domain is declined, not stuck, so reinjecting it
    would relocate particles that were never trapped."""
    voxel_channel(tmp_path)
    r = subprocess.run([FILAMENTS, str(tmp_path / "felbm_params.dat")] +
                       channel_args("scheme=explicit Dm=0.01 outside=reinject"),
                       capture_output=True, text=True, timeout=600)
    assert r.returncode != 0
    assert "not for diffusion" in r.stdout + r.stderr
