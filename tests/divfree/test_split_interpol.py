"""SplitInterpol on cases the cleaner writes, read through the apps.

A dataset solved on several ranks labels its rows by global cell id, so the
loader composes them through `mesh/cell_indices`, and the cleaner's output has
to carry labels that survive that. Each cell's interior values are built on
their own, in parallel, so the thread count must not change the field.
"""

import glob
import os
import shutil
import subprocess

import numpy as np
import pytest

from divfree_cases import (D, channel2d, channel3d, interpol_probe,
                           relabel_by_global_id, smooth_noslip, write_case)
from paths import app

PER = [True, True, False]


def cleaned_channel(folder, dim):
    """A cleaned channel of two stamps, periodic but across the walls: its
    parameter file."""
    X, cells = channel2d(6) if dim == 2 else channel3d(3)
    per = PER[3 - dim:]
    U = smooth_noslip(D.Topo(X, cells, per).node_x)
    cfg = write_case(folder / "in", X, cells, [U, 0.5 * U], per)
    D.clean_case(cfg, folder / "out", verbose=False)
    params = folder / "out" / "dolfin_params.dat"
    assert "divfree=true" in params.read_text()
    return params


@pytest.mark.skipif(not os.path.exists(app("tracers")), reason="tracers is not built")
@pytest.mark.parametrize("dim,mode", [(2, "triangle"), (3, "tet")], ids=["triangle", "tet"])
def test_thread_count_does_not_change_the_field(tmp_path, dim, mode):
    """Each cell's interior values are built on their own, so a run at 1 and at
    8 threads writes the same values to the byte.

    `random=false` so that the initializer places the same particles at both
    counts; the dump then differs only if the field does.
    """
    h5py = pytest.importorskip("h5py")
    params = cleaned_channel(tmp_path, dim)
    out = []
    for threads in (1, 8):
        folder = tmp_path / ("t%d" % threads)
        shutil.copytree(params.parent, folder)
        r = subprocess.run([app("tracers"), str(folder / params.name), "mode=" + mode,
                            "init_mode=points_" + "xyz"[:dim], "x0=0.5", "y0=0.5", "z0=0.5",
                            "random=false", "seed=1",
                            "Nrw=2000", "Nrw_max=2000", "int_order=1", "Dm=0", "dt=0.002",
                            "scheme=RK4", "T=0.2", "dump_intv=0.2", "stat_intv=1000",
                            "checkpoint_intv=1e9", "num_threads=%d" % threads],
                           capture_output=True, text=True, timeout=600,
                           env=dict(os.environ, OMP_NUM_THREADS=str(threads)))
        assert r.returncode == 0, r.stdout + r.stderr
        [dump] = glob.glob(str(folder / "Tracers" / "**" / "data_*.h5"), recursive=True)
        # Every dataset as its own bytes: the files themselves differ whenever
        # the two runs straddle a second, since HDF5 stamps its object headers
        sets = {}
        with h5py.File(dump, "r") as f:
            f.visititems(lambda n, o: sets.__setitem__(n, np.array(o).tobytes())
                         if isinstance(o, h5py.Dataset) else None)
        out.append(sets)
    a, b = out
    assert sorted(a) == sorted(b)
    assert "0.200000/u" in a, "the run did not step"
    for k in a:
        assert a[k] == b[k], k


@pytest.mark.slow
@pytest.mark.skipif(not os.path.exists(app("interpol")), reason="interpol is not built")
def test_a_case_labelled_by_global_cell_id_reads_as_the_same_field(tmp_path):
    """A cleaned case whose mesh's `cell_indices` is a permutation -- what
    `mpirun -n 8` writes -- read through `divfree=true` gives the field the same
    case labelled the identity way gives, at the same points.

    The two meshes differ in nothing but those labels, so the loader either
    composes the checkpoint's rows onto the right cells or refuses the file.
    """
    pytest.importorskip("h5py")
    X, cells = channel3d(3)
    U = smooth_noslip(D.Topo(X, cells, PER).node_x)
    plain = write_case(tmp_path / "plain", X, cells, [U], PER)
    labelled = write_case(tmp_path / "labelled", X, cells, [U], PER)
    relabel_by_global_id(tmp_path / "labelled")
    for cfg, out in ((plain, "plain_out"), (labelled, "labelled_out")):
        D.clean_case(cfg, tmp_path / out, verbose=False)
        assert "divfree=true" in (tmp_path / out / "dolfin_params.dat").read_text()

    a = interpol_probe(tmp_path / "plain_out" / "dolfin_params.dat", "tet", 4000, int_order=2)
    b = interpol_probe(tmp_path / "labelled_out" / "dolfin_params.dat", "tet", 4000, int_order=2)
    assert len(a["x"]) == len(b["x"]) > 0
    for c in "xyz":
        assert np.abs(a[c] - b[c]).max() == 0.0, c
    scale = max(np.abs(a["u" + c]).max() for c in "xyz")
    for c in "xyz":
        assert np.abs(a["u" + c] - b["u" + c]).max() < 1e-12 * scale, c


@pytest.mark.skipif(not os.path.exists(app("interpol")), reason="interpol is not built")
@pytest.mark.parametrize("dim,mode", [(2, "triangle"), (3, "tet")], ids=["triangle", "tet"])
def test_the_loader_gives_the_field_the_cleaner_writes_on_the_split_mesh(tmp_path, dim, mode):
    """`divfree=true` on the cleaned macro field and a plain P2 read of the
    cleaner's `--split` output are one field at the same points, value and
    gradient: loader and tool build the same interior values."""
    X, cells = channel2d(6) if dim == 2 else channel3d(3)
    per = PER[3 - dim:]
    cfg = write_case(tmp_path / "in", X, cells, [smooth_noslip(D.Topo(X, cells, per).node_x)], per)
    D.clean_case(cfg, tmp_path / "macro", verbose=False)
    D.clean_case(cfg, tmp_path / "split", split=True, verbose=False)
    a, b = (interpol_probe(tmp_path / f / "dolfin_params.dat", mode, 2000, int_order=2)
            for f in ("macro", "split"))
    assert sorted(a) == sorted(b) and len(a["x"]) > 0
    for k in a:
        assert np.abs(a[k] - b[k]).max() <= 1e-12 * max(1.0, np.abs(a[k]).max()), k
