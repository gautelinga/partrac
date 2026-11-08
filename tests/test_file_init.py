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
from paths import REPO, app, needs
from runs import copy_example, run_app

PARTRAC = app("partrac")
SPHERE = os.path.join(REPO, "data_example", "stokes_sphere", "expr_params.dat")
BASE = ("mode=analytic Nrw=1 Nrw_max=1000 ds_max=1e9 ds_min=1e-9 refine=false coarsen=false "
        "Dm=0 int_order=1 dt=0.01 T=0.01 dump_intv=0.01 stat_intv=1e9 checkpoint_intv=1e9 "
        "random=false seed=1").split()

needs_partrac = needs(PARTRAC)


def nodes_file(path, X):
    h5py = pytest.importorskip("h5py")
    with h5py.File(path, "w") as f:
        f["nodes"] = np.asarray(X, float)
    return path


@needs_partrac
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


@needs_partrac
def test_a_file_with_no_point_in_the_fluid_is_refused(tmp_path):
    """Every point inside the sphere: nothing to start from, said so."""
    f = nodes_file(tmp_path / "inside.h5", [[0., 0., 1.], [0.1, 0., 1.], [0., 0.2, 1.3]])
    params = copy_example(SPHERE, tmp_path / "c")
    r = run_app(PARTRAC, params, BASE, "init_mode=file:%s" % f, check=False)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr[-1000:]
    assert "no points inside the domain" in r.stdout + r.stderr


@needs_partrac
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


@needs_partrac
def test_a_file_of_more_points_than_nrw_max_is_refused(tmp_path):
    """30 points in the fluid with Nrw_max=20: the particle arrays hold 20,
    so the run stops naming both counts instead of writing past them."""
    z = np.linspace(-2, -1.5, 30)
    f = nodes_file(tmp_path / "many.h5", np.c_[np.zeros_like(z), np.zeros_like(z), z])
    params = copy_example(SPHERE, tmp_path / "c")
    r = run_app(PARTRAC, params, BASE, "Nrw_max=20", "init_mode=file:%s" % f, check=False)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr[-1000:]
    assert "the initial state has 30 nodes, more than Nrw_max = 20" in r.stderr, r.stderr
