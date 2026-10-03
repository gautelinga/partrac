"""scheme=RK4cells on the felbm lattice (mode=felbm, interpolation=linear).

The lattice's velocity is trilinear in each floor cell whose nodes are all
fluid and, next to a solid node, in each of the cell's eight sub-cubes; its
gradient jumps on every face between them, so RK4 loses its order as on a
mesh, and RK4cells cuts its steps at the faces. A sub-cube whose nearest node
is solid is a wall. On wall_felbm (walls at z = 0 and 15, a pillar): the order
of F and x, restarts, the thread count, the walls, and what it refuses:
interpolation=constant, which has no regions, and mode=fenics.

The kernel's own properties on the lattice (a step in one region is RK4's bit
for bit; faces crossed; resting at walls) are in test_cells.cpp.
"""

import os

import numpy as np
import pytest

from dumps import all_dumps, deformation_gradient, dump_at, read_stats
from paths import app
from runs import (cells_counts, checkpoint_folder, copy_case, halving_ratios, read_checkpoint, run_app,
                  same_checkpoint, same_dumps)
from test_structured import wall_felbm

TRACERS = app("tracers")
TENSORS = app("tracertensors")

needs_apps = pytest.mark.skipif(not all(os.path.exists(a) for a in (TRACERS, TENSORS)),
                                reason="the tracer apps are not built")

# 400 particles over the lattice's fluid
ARGS = ("mode=felbm init_mode=points_xyz x0=8 y0=8 z0=8 Nrw=400 Nrw_max=400 Dm=0 int_order=1 "
        "stat_intv=1e9 checkpoint_intv=1e9 random=false seed=1")


def run(root, name, *args, binary=TENSORS):
    """binary on a fresh wall_felbm case under root/name; (folder, stdout)."""
    d = root / name
    d.mkdir()
    wall_felbm(d)
    r = run_app(binary, d / "felbm_params.dat", ARGS, *args)
    return d, r.stdout


@needs_apps
def test_rk4cells_converges_at_fourth_order_on_the_lattice(tmp_path):
    """Halving dt divides the change of F and x by about 16 under RK4cells
    (p90 over particles, T = 2), from dt |u| near a cell down; under RK4 F's
    change about halves."""
    T = 2.0
    out = {}
    for scheme in ("RK4cells", "RK4"):
        runs = []
        for k in range(4):
            d, _ = run(tmp_path, "%s%d" % (scheme, k), "scheme=" + scheme,
                       "dt=%r T=%r dump_intv=%r" % (0.25 / 2 ** k, T, T * (1 + 1e-9)))
            g = dump_at(d, T)
            runs.append((g["points"], deformation_gradient(g)))
        out[scheme] = halving_ratios(runs)
    for rF, rx in out["RK4cells"]:
        assert rF > 12 and rx > 12, out
    assert all(rF < 3 for rF, _ in out["RK4"]), out


@needs_apps
def test_a_resumed_run_continues_bit_for_bit(tmp_path):
    """A run resumed from its checkpoint steps on as the run never stopped
    would, the region found again from the stored node and the position."""
    common = ["scheme=RK4cells", "dt=0.25", "dump_intv=1"]
    end = ["T=4.0", "checkpoint_intv=1e9"]
    cont, _ = run(tmp_path, "cont", common, end)
    split, _ = run(tmp_path, "split", common, ["T=1.75", "checkpoint_intv=1.75"])
    ck = read_checkpoint(split)
    assert (ck["cell_id"] >= 0).all()
    run_app(TENSORS, split / "felbm_params.dat", ARGS, common, end,
            ["restart_folder=%s" % checkpoint_folder(split)])
    assert len(same_dumps(cont, split, after=2)) == 2


@needs_apps
def test_one_and_eight_threads_agree_bit_for_bit(tmp_path):
    """Every dump, the nodes in the checkpoint and the counts printed are the
    same at 1 and 8 threads."""
    out = {}
    for n in (1, 8):
        d, stdout = run(tmp_path, "t%d" % n, "scheme=RK4cells", "dt=0.25", "T=2",
                        "dump_intv=1", "checkpoint_intv=1", "num_threads=%d" % n)
        out[n] = (d, [l for l in stdout.splitlines() if "RK4cells:" in l])
    (a, sa), (b, sb) = out[1], out[8]
    assert sa and sa == sb
    same_dumps(a, b)
    same_checkpoint(a, b)


@needs_apps
@pytest.mark.parametrize("dt", [0.25, 1.0, 4.0])
def test_at_the_walls_rk4cells_declines_no_more_than_rk4_and_stays_in_the_fluid(tmp_path, dt):
    """Tracers crowd at the walls and the pillar, where the flow comes to
    rest: RK4cells declines no more of them than RK4, falls back to RK4 on
    under one step in a thousand, and stores no position in a solid (every
    dumped point's nearest node is fluid)."""
    T = 16.0
    declined = {}
    for scheme in ("RK4cells", "RK4"):
        d, stdout = run(tmp_path, scheme, "scheme=" + scheme, "Nrw=1000", "Nrw_max=1000",
                        "dt=%r T=%r dump_intv=%r stat_intv=%r" % (dt, T, T / 4, dt), binary=TRACERS)
        declined[scheme] = read_stats(d)["n_declined"].sum()
        if scheme == "RK4cells":
            c = cells_counts(stdout)
            assert c and c["fallbacks"] < 1e-3, stdout
            for t, g in all_dumps(d).items():
                node = np.mod(np.floor(g["points"] + 0.5).astype(int), 16)
                # wall_felbm's solids: z = 0 and 15, the pillar 6 <= x, y <= 9
                solid = (node[:, 2] == 0) | (node[:, 2] == 15) | (
                    (node[:, 0] >= 6) & (node[:, 0] <= 9) & (node[:, 1] >= 6) & (node[:, 1] <= 9))
                assert not solid.any(), t
    assert declined["RK4cells"] <= declined["RK4"], declined


@needs_apps
def test_rk4cells_refuses_the_constant_lattice_and_mode_fenics(tmp_path, mesh_dir):
    """interpolation=constant has no regions; mode=fenics (dolfin's own
    evaluation) none either: both refused, naming the modes RK4cells takes."""
    d = tmp_path / "const"
    d.mkdir()
    wall_felbm(d)
    with open(d / "felbm_params.dat", "a") as f:
        f.write("interpolation=constant\n")
    r = run_app(TRACERS, d / "felbm_params.dat", ARGS, "scheme=RK4cells dt=0.1 T=0.1", check=False)
    assert r.returncode == 2 and "scheme=RK4cells steps in the cells of a mesh" in r.stderr, r.stderr
    assert "felbm (interpolation=linear)" in r.stderr and "use RK4" in r.stderr, r.stderr
    tet = copy_case(mesh_dir("tet"), tmp_path / "fenics")
    r = run_app(TRACERS, tet / "dolfin_params.dat",
                "mode=fenics init_mode=points_xyz x0=0.5 y0=0.5 z0=0.5 Nrw=10 Nrw_max=10 Dm=0 int_order=1 "
                "dt=0.1 T=0.1 random=false seed=1 scheme=RK4cells", check=False)
    assert r.returncode == 2 and "scheme=RK4cells steps in the cells of a mesh (not mode=fenics)" in r.stderr, r.stderr
    assert "use RK4" in r.stderr, r.stderr
