"""scheme=RK4cells on the divergence-free split (divfree=true).

The split field is P2 on each cell's barycentric split, so its gradient jumps
on the planes inside a cell as well as on the cell's facets: RK4 loses its
order there as on any mesh, and RK4cells cuts its steps at every plane. A
cleaned channel, periodic in x and y with walls at z = 0 and 1, two stamps:
the order of F and x, restarts with the cells in the checkpoint, the thread
count, and the walls, where the field is at rest and tracers crowd.

The kernel's own properties on the split (a step in one sub-cell is RK4's bit
for bit; planes crossed; stale cells relocated) are in test_cells.cpp.
"""

import os
import shutil

import numpy as np
import pytest

from divfree_cases import D, channel3d, smooth_noslip, write_case
from dumps import all_dumps, deformation_gradient, dump_at, read_stats
from paths import app
from runs import (cells_counts, checkpoint_file, checkpoint_folder, halving_ratios, read_checkpoint, run_app,
                  same_checkpoint, same_dumps)

TRACERS = app("tracers")
TENSORS = app("tracertensors")

needs_apps = pytest.mark.skipif(not all(os.path.exists(a) for a in (TRACERS, TENSORS)),
                                reason="the tracer apps are not built")

# 200 particles over the channel
ARGS = ("mode=tet init_mode=points_xyz x0=0.5 y0=0.5 z0=0.5 Nrw=200 Nrw_max=200 Dm=0 int_order=1 "
        "stat_intv=1e9 checkpoint_intv=1e9 random=false seed=1")


@pytest.fixture(scope="module")
def channel(tmp_path_factory):
    """The cleaned channel: 4^3 boxes of six tets, the flow at rest on the
    walls and moving towards them, at t = 0 and half of it at t = 1."""
    root = tmp_path_factory.mktemp("split_channel")
    X, cells = channel3d(4)
    per = [True, True, False]
    U = smooth_noslip(D.Topo(X, cells, per).node_x)
    cfg = write_case(root / "in", X, cells, [U, 0.5 * U], per)
    D.clean_case(cfg, root / "out", verbose=False)
    assert "divfree=true" in (root / "out" / "dolfin_params.dat").read_text()
    return root / "out"


def run(channel, root, name, *args, binary=TENSORS):
    """binary on a copy of the channel; (folder, stdout)."""
    d = root / name
    shutil.copytree(channel, d)
    r = run_app(binary, d / "dolfin_params.dat", ARGS, *args)
    return d, r.stdout


@needs_apps
def test_rk4cells_converges_at_fourth_order_on_the_split(channel, tmp_path):
    """Halving dt divides the change of F and x by about 16 under RK4cells
    (p90 over particles, T = 0.3); under RK4 F's change about halves."""
    T = 0.3
    out = {}
    for scheme in ("RK4cells", "RK4"):
        runs = []
        for k in range(4):
            d, _ = run(channel, tmp_path, "%s%d" % (scheme, k), "scheme=" + scheme,
                       "dt=%r T=%r dump_intv=%r" % (0.025 / 2 ** k, T, T * (1 + 1e-9)))
            g = dump_at(d, T)
            runs.append((g["points"], deformation_gradient(g)))
        out[scheme] = halving_ratios(runs)
    for rF, rx in out["RK4cells"]:
        assert rF > 12 and rx > 12, out
    assert all(rF < 3 for rF, _ in out["RK4"]), out


@needs_apps
def test_a_resumed_run_continues_bit_for_bit(channel, tmp_path):
    """A run resumed from its checkpoint, cells in it, steps on bit for bit;
    without the cells it is located afresh and resumes to round-off."""
    common = ["scheme=RK4cells", "dt=0.0625", "dump_intv=0.25"]
    end = ["T=1.0", "checkpoint_intv=1e9"]
    cont, _ = run(channel, tmp_path, "cont", common, end)
    split, _ = run(channel, tmp_path, "split", common, ["T=0.4375", "checkpoint_intv=0.4375"])
    ck = read_checkpoint(split)
    assert (ck["cell_id"] >= 0).all()
    run_app(TENSORS, split / "dolfin_params.dat", ARGS, common, end,
            ["restart_folder=%s" % checkpoint_folder(split)])
    assert len(same_dumps(cont, split, after=0.5)) == 2
    old, _ = run(channel, tmp_path, "old", common, ["T=0.4375", "checkpoint_intv=0.4375"])
    h5py = pytest.importorskip("h5py")
    with h5py.File(checkpoint_file(old), "r+") as h:
        del h["cell_id"]
    run_app(TENSORS, old / "dolfin_params.dat", ARGS, common, end,
            ["restart_folder=%s" % checkpoint_folder(old)])
    ga, go = dump_at(cont, 1.0), dump_at(old, 1.0)
    assert np.abs(ga["points"] - go["points"]).max() < 1e-9
    assert np.abs(deformation_gradient(ga) - deformation_gradient(go)).max() < 1e-7


@needs_apps
def test_one_and_eight_threads_agree_bit_for_bit(channel, tmp_path):
    """Every dump, the cells in the checkpoint and the counts printed are the
    same at 1 and 8 threads."""
    out = {}
    for n in (1, 8):
        d, stdout = run(channel, tmp_path, "t%d" % n, "scheme=RK4cells", "dt=0.05", "T=0.5",
                        "dump_intv=0.25", "checkpoint_intv=0.25", "num_threads=%d" % n)
        out[n] = (d, [l for l in stdout.splitlines() if "RK4cells:" in l])
    (a, sa), (b, sb) = out[1], out[8]
    assert sa and sa == sb
    same_dumps(a, b)
    same_checkpoint(a, b)


@needs_apps
@pytest.mark.parametrize("dt", [0.05, 0.2, 0.5])
def test_at_the_walls_rk4cells_declines_no_more_than_rk4_and_stays_in_the_fluid(channel, tmp_path, dt):
    """Tracers crowd at the walls, where the flow comes to rest: RK4cells
    declines no more of them than RK4, falls back to RK4 on under one step in
    a thousand, and stores no position outside the channel."""
    T = 2.0
    declined = {}
    for scheme in ("RK4cells", "RK4"):
        d, stdout = run(channel, tmp_path, scheme, "scheme=" + scheme, "Nrw=1000", "Nrw_max=1000",
                        "dt=%r T=%r dump_intv=%r stat_intv=%r" % (dt, T, T / 4, dt), binary=TRACERS)
        declined[scheme] = read_stats(d)["n_declined"].sum()
        if scheme == "RK4cells":
            c = cells_counts(stdout)
            assert c and c["fallbacks"] < 1e-3, stdout
            for t, g in all_dumps(d).items():
                z = g["points"][:, 2]
                assert (z >= 0).all() and (z <= 1).all(), t
    assert declined["RK4cells"] <= declined["RK4"], declined


@needs_apps
@pytest.mark.parametrize("case", ["in", "out"])
@pytest.mark.parametrize("scheme", ["explicit", "RK4", "RK4cells"])
def test_past_the_last_stamp_the_field_is_held(channel, tmp_path, case, scheme):
    """T beyond the stamps is cut to the last one, t = 1; the loop's last step,
    0.9 to 1.2, is cut at that stamp and goes on past it on the last stamp's
    field held, for every scheme, on the plain field ("in") and the split one
    ("out"): its checkpoint at 1.2 is that of a run to T = 0.9 on the same
    case with the last stamp repeated at t = 2. A Debug build also checks
    that every evaluation lies in its bracket."""
    h5py = pytest.importorskip("h5py")
    runs = {}
    for name, extra, T in (("held", None, 2.0), ("repeated", "2 ", 0.9)):
        d = tmp_path / name
        shutil.copytree(channel.parent / case, d)
        stamps = d / "timestamps.dat"
        if extra:
            last = stamps.read_text().split("\n")[-2].split()[1]
            stamps.write_text(stamps.read_text() + extra + last + "\n")
        run_app(TRACERS, d / "dolfin_params.dat", ARGS, "scheme=" + scheme, "dt=0.3", "T=%r" % T,
                "dump_intv=0.3")
        runs[name] = d
    assert max(all_dumps(runs["held"])) == pytest.approx(0.9)
    with h5py.File(checkpoint_file(runs["held"]), "r") as h:
        assert h.attrs["t"] == pytest.approx(1.2)
    held, repeated = (read_checkpoint(runs[k])["points"] for k in ("held", "repeated"))
    assert np.abs(held - repeated).max() < 1e-12


@needs_apps
def test_frozen_fields_on_the_split_do_not_depend_on_the_time(channel, tmp_path):
    """Fields frozen between the two stamps, interior values too, are their
    blend at t_frozen at every time: a run from t0 = 0 and the same run from
    t0 = 2, past the last stamp, end at the same points bit for bit."""
    ends = []
    for t0 in (0.0, 2.0):
        d, _ = run(channel, tmp_path, "t0_%g" % t0, "scheme=RK4 dt=0.05 frozen_fields=true t_frozen=0.5",
                   "t0=%r T=%r dump_intv=0.4" % (t0, t0 + 0.4), binary=TRACERS)
        start, end = dump_at(d, t0)["points"], dump_at(d, t0 + 0.4)["points"]
        assert np.abs(end - start).max() > 1e-3
        ends.append(end)
    assert np.array_equal(ends[0], ends[1])

