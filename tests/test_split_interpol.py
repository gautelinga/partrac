"""SplitInterpol against the offline reconstruction of the same field.

`harness/divfree_research/d1/split/<case>/` holds, for each research case, the
cleaned P2 velocity on the original mesh (`clean_p2/`) and the same
divergence-free field written out on the barycentrically split mesh
(`split_p2/`), which `mode=triangle` and `mode=tet` read as an ordinary P2
field. Reading `clean_p2` with `divfree=true` has to give that field back: the
loader builds in C++ what the Python tool wrote to disk.

The `interpol` app probes both at the same random points -- the point set comes
from the seed and the bounding box, which the two meshes share -- and the test
requires the velocities and their gradients to agree to round-off. The cases
are 16 GB of research output that is not in git, so the test skips when they
are absent. One thread: the app draws its points per thread.

A case built here covers the other half of the loader's input: a dataset solved
on several ranks labels its rows by global cell id, so the loader composes them
through `mesh/cell_indices`, and the cleaner's output has to carry labels that
survive that. Another cleans one case on one rank and on three and reads both
through the loader, which is where a solve shared out over ranks has to end up.
"""

import glob
import os
import subprocess

import numpy as np
import pytest

from paths import REPO, app
# the cleaner and the fixtures that write a case for it
from test_divfree_clean import (D, channel3d, mpi_clean, relabel_by_global_id,
                                smooth_noslip, write_case)

CASES = os.path.join(REPO, "harness", "divfree_research", "d1", "split")

# case folder -> the loader's mode
MODES = {"cyl20_h0.08": "triangle", "cyl20_h0.04": "triangle",
         "obstacles_h0.04": "triangle", "obstacles_h0.02": "triangle",
         "sphere_h0.1": "tet", "sphere_h0.06": "tet"}

# One 2D case and the sphere; the rest are the same construction on more cells
CHECKED = ["obstacles_h0.04", "sphere_h0.06"]

NPTS = 4000


def link_case(src, dst, divfree):
    """A run folder whose inputs are symlinks into the research cases, so that
    the probe's output does not land in 16 GB of cached results."""
    os.makedirs(dst, exist_ok=True)
    for f in glob.glob(os.path.join(src, "*.h5")) + [os.path.join(src, "timestamps.dat")]:
        link = os.path.join(dst, os.path.basename(f))
        if not os.path.exists(link):
            os.symlink(f, link)
    with open(os.path.join(src, "dolfin_params.dat")) as f:
        prm = f.read()
    if divfree:
        prm += "divfree=true\n"
    with open(os.path.join(dst, "dolfin_params.dat"), "w") as f:
        f.write(prm)
    return os.path.join(dst, "dolfin_params.dat")


def probe(params, mode):
    """Velocity and gradient at NPTS points, sorted by position."""
    env = dict(os.environ, OMP_NUM_THREADS="1")
    r = subprocess.run([app("interpol"), params, "mode=" + mode, "Nrw=%d" % NPTS,
                        "int_order=2", "random=false", "seed=1", "t0=0"],
                       capture_output=True, text=True, timeout=1800, env=env)
    assert r.returncode == 0, r.stdout + r.stderr
    import h5py
    [out] = glob.glob(os.path.join(os.path.dirname(params), "Interpolation", "**",
                                   "interpolation.h5part"), recursive=True)
    with h5py.File(out, "r") as h:
        g = h[list(h.keys())[0]]
        d = {k: np.array(g[k]) for k in g}
    order = np.lexsort((d["z"], d["y"], d["x"]))
    return {k: v[order] for k, v in d.items()}


@pytest.mark.slow
@pytest.mark.parametrize("case", CHECKED)
def test_matches_the_offline_split_field(case, tmp_path):
    """The velocity and its gradient through SplitInterpol are the offline field."""
    src = os.path.join(CASES, case)
    if not os.path.isdir(os.path.join(src, "clean_p2")):
        pytest.skip("the research cases are not in this tree")
    pytest.importorskip("h5py")
    mode = MODES[case]
    clean = probe(link_case(os.path.join(src, "clean_p2"), str(tmp_path / "clean"), True), mode)
    split = probe(link_case(os.path.join(src, "split_p2"), str(tmp_path / "split"), False), mode)

    assert len(clean["x"]) == len(split["x"])
    for c in "xyz":
        assert np.allclose(clean[c], split[c], rtol=0, atol=1e-13)
    comps = ["ux", "uy", "uz"]
    grads = ["uxx", "uxy", "uxz", "uyx", "uyy", "uyz", "uzx", "uzy", "uzz"]
    scale = max(np.abs(split[c]).max() for c in comps)
    gscale = max(np.abs(split[c]).max() for c in grads)
    for c in comps:
        assert np.abs(clean[c] - split[c]).max() < 1e-12 * scale, c
    for c in grads:
        assert np.abs(clean[c] - split[c]).max() < 1e-12 * gscale, c
    # What the construction is for: the divergence is zero in both
    assert np.abs(clean["divu"]).max() < 1e-10 * gscale


@pytest.mark.slow
def test_thread_count_does_not_change_the_field(tmp_path):
    """Each cell's interior values are built on their own, so a run at 1 and at
    8 threads writes the same values to the byte.

    `random=false` so that the initializer places the same particles at both
    counts; the dump then differs only if the field does.
    """
    src = os.path.join(CASES, "obstacles_h0.04", "clean_p2")
    if not os.path.isdir(src):
        pytest.skip("the research cases are not in this tree")
    h5py = pytest.importorskip("h5py")
    out = []
    for threads in (1, 8):
        folder = str(tmp_path / ("t%d" % threads))
        params = link_case(src, folder, True)
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
        r = subprocess.run([app("tracers"), params, "mode=triangle", "init_mode=points_xy",
                            "x0=0.5", "y0=0.5", "random=false", "seed=1",
                            "Nrw=2000", "Nrw_max=2000", "int_order=1", "Dm=0", "dt=0.002",
                            "scheme=RK4", "T=0.2", "dump_intv=0.2", "stat_intv=1000",
                            "checkpoint_intv=1e9", "num_threads=%d" % threads],
                           capture_output=True, text=True, timeout=1800, env=env)
        assert r.returncode == 0, r.stdout + r.stderr
        [dump] = glob.glob(os.path.join(folder, "Tracers", "**", "data_*.h5"), recursive=True)
        # Every dataset as its own bytes: the files themselves differ whenever
        # the two runs straddle a second, since HDF5 stamps its object headers
        sets = {}
        with h5py.File(dump, "r") as f:
            f.visititems(lambda n, o: sets.__setitem__(n, np.array(o).tobytes())
                         if isinstance(o, h5py.Dataset) else None)
        out.append(sets)
    assert sorted(out[0]) == sorted(out[1])
    for k in out[0]:
        assert out[0][k] == out[1][k], k


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
    per = [True, True, False]
    U = smooth_noslip(D.Topo(X, cells, per).node_x)
    plain = write_case(tmp_path / "plain", X, cells, [U], per)
    labelled = write_case(tmp_path / "labelled", X, cells, [U], per)
    relabel_by_global_id(tmp_path / "labelled")
    for cfg, out in ((plain, "plain_out"), (labelled, "labelled_out")):
        D.clean_case(cfg, tmp_path / out, verbose=False)
        assert "divfree=true" in (tmp_path / out / "dolfin_params.dat").read_text()

    a = probe(str(tmp_path / "plain_out" / "dolfin_params.dat"), "tet")
    b = probe(str(tmp_path / "labelled_out" / "dolfin_params.dat"), "tet")
    assert len(a["x"]) == len(b["x"]) > 0
    for c in "xyz":
        assert np.abs(a[c] - b[c]).max() == 0.0, c
    scale = max(np.abs(a["u" + c]).max() for c in "xyz")
    for c in "xyz":
        assert np.abs(a["u" + c] - b["u" + c]).max() < 1e-12 * scale, c


@pytest.mark.slow
@pytest.mark.skipif(not os.path.exists(app("interpol")), reason="interpol is not built")
def test_a_case_cleaned_on_three_ranks_reads_as_the_one_rank_field(tmp_path):
    """The whole path under MPI: a case cleaned by `mpirun -n 3` is read through
    `divfree=true` and gives, at the same points, what the one-rank clean of the
    same case gives.

    The ranks assemble their own cell blocks, so the matrices differ in their
    last bits and the solves are not the same iteration; the field they agree on
    is what a reader sees, and the loader's own flux criterion is what either
    output has to pass to be read at all."""
    pytest.importorskip("h5py")
    X, cells = channel3d(3)
    per = [True, True, False]
    U = smooth_noslip(D.Topo(X, cells, per).node_x)
    cfg = write_case(tmp_path / "in", X, cells, [U], per)
    got = {}
    for n in (1, 3):
        out = tmp_path / ("out%d" % n)
        mpi_clean(cfg, out, n)
        assert "divfree=true" in (out / "dolfin_params.dat").read_text()
        got[n] = probe(str(out / "dolfin_params.dat"), "tet")
    a, b = got[1], got[3]
    assert len(a["x"]) == len(b["x"]) > 0
    scale = max(np.abs(a["u" + c]).max() for c in "xyz")
    for c in "xyz":
        assert np.abs(a[c] - b[c]).max() == 0.0, c
        assert np.abs(a["u" + c] - b["u" + c]).max() < 1e-8 * scale, c
