"""init_mode=file:<path>: a line read from an HDF5 file's `nodes` rows.

The points outside the fluid are dropped and consecutive kept points are
joined by edges, so a line through an obstacle starts as the pieces outside
it; a file with no point in the fluid is refused. The fluid is the Stokes
flow around the sphere of radius 1 centred at (0, 0, 1), whose interior is
outside.
"""

import os

import numpy as np
import pytest

from dumps import dump_at
from paths import REPO, app
from runs import copy_example, run_app

PARTRAC = app("partrac")
SPHERE = os.path.join(REPO, "data_example", "stokes_sphere", "expr_params.dat")
BASE = ("mode=analytic Nrw=1 Nrw_max=1000 ds_max=1e9 ds_min=1e-9 refine=false coarsen=false "
        "Dm=0 int_order=1 dt=0.01 T=0.01 dump_intv=0.01 stat_intv=1e9 checkpoint_intv=1e9 "
        "random=false seed=1").split()

needs = pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")


def nodes_file(path, X):
    h5py = pytest.importorskip("h5py")
    with h5py.File(path, "w") as f:
        f["nodes"] = np.asarray(X, float)
    return path


@needs
def test_a_line_through_the_sphere_starts_as_its_pieces_outside(tmp_path):
    """41 points on the z axis from -2 to 2: the 19 with 0 < z < 2 lie in the
    sphere; the other 22 start where the file puts them, in its order."""
    z = np.linspace(-2, 2, 41)
    X = np.c_[np.zeros_like(z), np.zeros_like(z), z]
    f = nodes_file(tmp_path / "line_of_points.h5", X)
    params = copy_example(SPHERE, tmp_path / "c")
    run_app(PARTRAC, params, BASE, "init_mode=file:%s" % f)
    p = dump_at(tmp_path / "c", 0.0)["points"]
    kept = X[np.abs(z - 1) >= 1]
    assert len(p) == len(kept) == 22
    assert np.array_equal(p, kept)


@needs
def test_a_file_with_no_point_in_the_fluid_is_refused(tmp_path):
    """Every point inside the sphere: nothing to start from, said so."""
    f = nodes_file(tmp_path / "inside.h5", [[0., 0., 1.], [0.1, 0., 1.], [0., 0.2, 1.3]])
    params = copy_example(SPHERE, tmp_path / "c")
    r = run_app(PARTRAC, params, BASE, "init_mode=file:%s" % f, check=False)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr[-1000:]
    assert "no points inside the domain" in r.stdout + r.stderr


@needs
def test_a_node_reinjected_from_a_file_start_moves_along_every_axis(tmp_path):
    """outside=reinject moves a declined node by a random offset along
    init_mode's directions, and a file start has none: all three, whatever
    the path holds (its second '_' piece here, 'file', names no axis). Five
    nodes on a line below the sphere; one is put into it through the
    checkpoint, so its step is declined; reinjected, it ends in the fluid and
    off the plane y = 0, which the flow never leaves (u_y = 0 there). The
    timeout catches a node that can never be moved."""
    from runs import checkpoint_folder, read_checkpoint, write_checkpoint
    f = nodes_file(tmp_path / "line.h5", [[x, 0., -2.] for x in (-1., -0.5, 0., 0.5, 1.)])
    params = copy_example(SPHERE, tmp_path / "c")
    args = BASE + ["init_mode=file:%s" % f, "outside=reinject"]
    run_app(PARTRAC, params, args)
    ck = checkpoint_folder(tmp_path / "c")
    pos = read_checkpoint(ck)["points"]
    pos[2] = [0., 0., 1.]                     # the sphere's centre
    write_checkpoint(ck, points=pos)
    run_app(PARTRAC, params, args, "T=0.03", "restart_folder=%s" % ck, timeout=60)
    p = dump_at(tmp_path / "c", 0.03)["points"]
    assert len(p) == 5
    assert np.all(np.linalg.norm(p - [0., 0., 1.], axis=1) >= 1.)
    moved = np.abs(p[:, 1]) > 0
    assert moved.sum() == 1, p


@needs
@pytest.mark.parametrize("nodes,message", [
    (np.zeros((3, 4)), "columns, not one to three coordinates"),
    (np.zeros(6), "not a two-dimensional array"),
])
def test_a_nodes_array_of_the_wrong_shape_is_refused(tmp_path, nodes, message):
    """nodes holds one point a row, one to three coordinates: four columns
    would write past a point, and a flat array has no rows; both are said so."""
    f = nodes_file(tmp_path / "bad.h5", nodes)
    params = copy_example(SPHERE, tmp_path / "c")
    r = run_app(PARTRAC, params, BASE, "init_mode=file:%s" % f, check=False)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr[-1000:]
    assert message in r.stderr, r.stderr
