"""Renumbering cells by their dofs changes nothing a user can see.

`mode=fenics` renumbers its cells at load: `renumber_cells=auto` does so when
consecutive cells rarely share a dof, which is when the file's cell order
carries no locality. Cell ids are internal, so a run must not depend on the
choice: with `always` and with `never`, every dump must agree particle by
particle in id order. This is also where a point lying exactly on a shared face
would show up, since the neighbour walk then tests the cells around it in a
different order and could settle in the other one.

It is the only loader that reads the key. The rest load their cells in Morton
order of the centroids, which is the order their tree wants, so their schemas
do not declare `renumber_cells` at all and a file that sets it is refused at
startup rather than silently ignored; the second test here checks that refusal
and its message.

What those loaders do renumber, by themselves, is the nodes: a mesh file whose
vertex numbering carries no locality makes every gather a scatter over the
whole field, so when a cell's vertex ids
span more than a quarter of the vertices on average the nodes are put along
the same Morton curve as the cells. Node ids are internal too, so the last test
takes a P1 case, relabels its vertices at random, and requires the same dumps
and the announcement only for the relabelled file.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

from cases import KINDS
from dumps import dump_at
from paths import REPO, app, built_with_dolfin
from runs import run_app

TRACERS = app("tracers")
DT, T = 0.001, 0.01
BASE = ("init_mode=points_xyz Nrw=2000 Nrw_max=2000 int_order=1 Dm=0 scheme=RK4 "
        "dt=%g T=%g dump_intv=%g stat_intv=1e9 checkpoint_intv=1e9 random=false "
        "seed=1" % (DT, T, DT))

# the mesh mode=fenics reads, and the kind of input it is generated from
KIND, MODE = "tet", "fenics"

needs_dolfin = pytest.mark.skipif(not built_with_dolfin(),
                                  reason="mode=fenics needs a build with dolfin")


def run(case, tmp_path, renumber):
    """Run tracers on a copy of the case with renumber_cells=renumber; return its folder."""
    d = tmp_path / (MODE + "_" + renumber)
    shutil.copytree(case, d)
    with open(d / "dolfin_params.dat", "a") as f:
        f.write("\nrenumber_cells=%s\n" % renumber)
    r = run_app(TRACERS, d / "dolfin_params.dat", BASE, "mode=" + MODE, KINDS[KIND][1])
    # only `always` renumbers a mesh whose cells are already local
    assert ("renumbering cells by dofs" in r.stdout) == (renumber == "always"), r.stdout[-800:]
    return d


@pytest.mark.fenics
@needs_dolfin
@pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")
def test_renumbered_cells_give_identical_dumps(mesh_dir, tmp_path):
    """Renumbered and unrenumbered runs agree dataset by dataset, in particle-id
    order, at every step."""
    a = run(mesh_dir(KIND), tmp_path, "never")
    b = run(mesh_dir(KIND), tmp_path, "always")
    for k in range(int(round(T / DT)) + 1):
        t = k * DT
        da, db = dump_at(a, t), dump_at(b, t)
        assert set(da) == set(db)
        for name in da:
            assert np.array_equal(da[name], db[name]), (t, name)


# the keys each schema requires, so that the only thing wrong with the file
# below is the key under test
REQUIRED = {
    "tet": "velocity_space=P1\npressure_space=P1\ntimestamps=timestamps.dat\nmesh=mesh.h5\n",
    "triangle": "velocity_space=P1\npressure_space=P1\ntimestamps=timestamps.dat\nmesh=mesh.h5\n",
    "trianglefreq": ("velocity_space=P1\npressure_space=P1\nfreqstamps=freqstamps.dat\n"
                     "mesh=mesh.h5\ntau=1.0\nt_min=0\nt_max=1\n"),
    "xdmftriangle": "u=u.xdmf\np=p.xdmf\n",
}


@pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")
@pytest.mark.parametrize("mode", sorted(REQUIRED))
def test_a_loader_that_does_not_renumber_refuses_the_key(tmp_path, mode):
    """A mode whose schema has no renumber_cells says so instead of ignoring it.

    The schema is read before any mesh, so the case need not be a working one:
    what is checked is that the key is refused by name, and that the same file
    without it gets past the schema and fails further in, on its missing files.
    """
    args = (BASE + " mode=" + mode + " x0=0.5 y0=0.5 z0=0.5").split()

    def start(name, text):
        f = tmp_path / (name + ".dat")
        f.write_text(text)
        return run_app(TRACERS, f, args, check=False)

    r = start("with_key", REQUIRED[mode] + "renumber_cells=always\n")
    assert r.returncode != 0
    assert "unknown parameter 'renumber_cells'" in r.stdout + r.stderr, r.stdout + r.stderr
    # without it the file passes the schema and stops on its missing mesh
    r = start("without_key", REQUIRED[mode])
    assert r.returncode != 0
    assert "unknown parameter" not in r.stdout + r.stderr, r.stdout + r.stderr


P1_TET = os.path.join(REPO, "data_example", "test_tet_p1")


@pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")
def test_a_mesh_numbered_without_locality_gives_identical_dumps(tmp_path):
    """P1 only: a P2 file's edge dofs follow the vertex ids within each cell, so
    relabelling its vertices would make the checkpoint inconsistent."""
    pytest.importorskip("dolfin", reason="generating the case needs dolfin")
    h5py = pytest.importorskip("h5py")
    local = tmp_path / "local"
    shutil.copytree(P1_TET, local)
    r = subprocess.run(["python3", "generate_up.py"], cwd=local, capture_output=True,
                       text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr
    shuffled = tmp_path / "shuffled"
    shutil.copytree(local, shuffled)
    with h5py.File(shuffled / "mesh.h5", "r+") as f:
        x, topo = f["mesh/coordinates"][()], f["mesh/topology"][()]
        label = np.random.default_rng(2).permutation(len(x))   # new id of old vertex
        moved = np.empty_like(x)
        moved[label] = x
        f["mesh/coordinates"][...] = moved
        f["mesh/topology"][...] = label[topo].astype(topo.dtype)   # local order kept

    out = {}
    for d in (local, shuffled):
        out[d] = run_app(TRACERS, d / "dolfin_params.dat", BASE, "mode=tet").stdout
    assert "renumbering nodes" not in out[local], out[local][-600:]
    assert "renumbering nodes" in out[shuffled], out[shuffled][-600:]
    for k in range(int(round(T / DT)) + 1):
        da, db = dump_at(local, k * DT), dump_at(shuffled, k * DT)
        assert set(da) == set(db)
        for name in da:
            assert np.array_equal(da[name], db[name]), (k * DT, name)
