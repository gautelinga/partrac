"""The mesh interpolators, on meshes generated at test time.

The data_example folders ship generate_up.py but not the mesh.h5 it writes, so
nothing exercised mode=triangle or mode=tet. Generating needs dolfin, which the
dolfin-off build does not have either, so the whole module skips without it.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

from paths import REPO, app, built_with_dolfin

PARTRAC = app("partrac")
DATA = os.path.join(REPO, "data_example")

pytestmark = pytest.mark.skipif(not built_with_dolfin(),
                                reason="partrac was built without dolfin")

# mode -> the mesh kind that conftest generates
CASES = [("triangle", "triangle"), ("tet", "tet")]

BASE = ("init_mode=uniform_x Nrw=200 Nrw_max=5000 ds_max=0.4 ds_min=0.1 "
        "Dm=0 int_order=1 dt=0.005 T=0.05 dump_intv=0.05 stat_intv=0.05 "
        "checkpoint_intv=0.05 random=false seed=1").split()

FILES = ("dolfin_params.dat", "mesh.h5", "up_0.h5", "timestamps.dat")


@pytest.fixture(params=CASES, ids=[c[0] for c in CASES])
def mesh_case(request, mesh_dir):
    mode, kind = request.param
    return mode, mesh_dir(kind)


def run(case_dir, args):
    r = subprocess.run([PARTRAC, str(case_dir / "dolfin_params.dat")] + args,
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr
    return r


def make_run(mesh_case, dest, extra=()):
    mode, src = mesh_case
    dest.mkdir()
    for f in FILES:
        shutil.copy(src / f, dest / f)
    run(dest, BASE + ["mode=" + mode] + list(extra))
    return dest


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_particles_move_and_stay_in_the_periodic_domain(mesh_case, tmp_path):
    d = make_run(mesh_case, tmp_path / "run")
    h5py = pytest.importorskip("h5py")
    dumps = sorted(d.glob("**/data_from_t*.h5"))
    assert dumps
    with h5py.File(dumps[0]) as h:
        keys = sorted(h.keys(), key=float)
        first = np.array(h[keys[0]]["points"])
        last = np.array(h[keys[-1]]["points"])
    assert np.abs(last - first).max() > 1e-3          # it advected
    assert np.isfinite(last).all()
    # |u| <= sqrt(3) on both meshes, so over T=0.05 nothing can travel 0.2.
    # Positions are not wrapped: a periodic run has to keep the dispersion.
    assert np.abs(last - first).max() < 0.2


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_mesh_run_is_reproducible(mesh_case, tmp_path):
    a = make_run(mesh_case, tmp_path / "a")
    b = make_run(mesh_case, tmp_path / "b")
    xa = np.loadtxt(sorted(a.glob("**/Checkpoints/positions.pos"))[0])
    xb = np.loadtxt(sorted(b.glob("**/Checkpoints/positions.pos"))[0])
    assert np.array_equal(xa, xb)


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_plane_poiseuille_on_a_mesh_conserves_x(mesh_case, tmp_path):
    if mesh_case[0] != "triangle":
        pytest.skip("only the plane Poiseuille case has an invariant axis")
    d = make_run(mesh_case, tmp_path / "axis")
    h5py = pytest.importorskip("h5py")
    with h5py.File(sorted(d.glob("**/data_from_t*.h5"))[0]) as h:
        keys = sorted(h.keys(), key=float)
        first = np.array(h[keys[0]]["points"])
        last = np.array(h[keys[-1]]["points"])
    # u is along y and depends only on x
    assert np.abs(last[:, 0] - first[:, 0]).max() < 1e-12
    assert np.abs(last[:, 1] - first[:, 1]).max() > 1e-3
