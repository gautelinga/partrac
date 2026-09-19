"""Renumbering cells by their dofs changes nothing a user can see.

Every mesh interpolator may renumber its cells at load: `renumber_cells=auto`
does so when consecutive cells rarely share a dof, which is when the file's
cell order carries no locality. Cell ids are internal, so a run must not depend
on the choice: with `always` and with `never`, every dump must agree particle
by particle in id order. This is also where a point lying exactly on a shared
face would show up, since the neighbour walk then tests the cells around it in
a different order and could settle in the other one.

The cases cover one interpolator each: TetInterpol, TriangleInterpol,
TriangleFreqInterpol, XDMFTriangleInterpol and DolfInterpol (mode=fenics).
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

from dumps import dump_at
from paths import app, built_with_dolfin
from test_apps import KINDS

TRACERS = app("tracers")
DT, T = 0.001, 0.01
BASE = ("init_mode=points_xyz Nrw=2000 Nrw_max=2000 int_order=1 Dm=0 scheme=RK4 "
        "dt=%g T=%g dump_intv=%g stat_intv=1e9 checkpoint_intv=1e9 random=false "
        "seed=1" % (DT, T, DT))

# (input kind, mode): one per mesh interpolator
CASES = [("tet", "tet"), ("triangle", "triangle"), ("trianglefreq", "trianglefreq"),
         ("xdmf", "xdmftriangle"), ("tet", "fenics")]


def run(case, tmp_path, kind, mode, renumber):
    """Run tracers on a copy of the case with renumber_cells=renumber; return its folder."""
    d = tmp_path / (mode + "_" + renumber)
    shutil.copytree(case, d)
    with open(d / "dolfin_params.dat", "a") as f:
        f.write("\nrenumber_cells=%s\n" % renumber)
    args = (BASE + " mode=" + mode + " " + KINDS[kind][1]).split()
    r = subprocess.run([TRACERS, str(d / "dolfin_params.dat")] + args,
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr
    # only `always` renumbers a mesh whose cells are already local
    assert ("renumbering cells by dofs" in r.stdout) == (renumber == "always"), r.stdout[-800:]
    return d


@pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")
@pytest.mark.skipif(not built_with_dolfin(), reason="partrac was built without dolfin")
@pytest.mark.parametrize("kind,mode", CASES, ids=["%s-%s" % c for c in CASES])
def test_renumbered_cells_give_identical_dumps(mesh_dir, xdmf_dir, tmp_path, kind, mode):
    """Renumbered and unrenumbered runs agree dataset by dataset, in particle-id
    order, at every step."""
    case = xdmf_dir if kind == "xdmf" else mesh_dir(kind)
    a = run(case, tmp_path, kind, mode, "never")
    b = run(case, tmp_path, kind, mode, "always")
    for k in range(int(round(T / DT)) + 1):
        t = k * DT
        da, db = dump_at(a, t), dump_at(b, t)
        assert set(da) == set(db)
        for name in da:
            assert np.array_equal(da[name], db[name]), (t, name)
