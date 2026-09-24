"""python/divfree/divfree_clean.py as a program: a cleaned case read back by the apps,
and the tool run on several ranks.

A case cleaned on three ranks meets the criterion the one-rank run meets, gives
the same field, read through the loader too, and reports the same numbers; a
write only the first rank can fail takes the job down rather than leaving the
others waiting for it. The MPI tests skip where no launcher starts a job of
several ranks for mpi4py, and fail there with PARTRAC_REQUIRE_MPI set.
"""

import os
import re
import subprocess
import sys

import numpy as np
import pytest

from divfree_cases import (CASES, D, IDS, INTERPOL, channel2d, channel3d,
                           interpol_probe, mpi_clean, mpi_run, pocket_case, same_log,
                           same_output, smooth_noslip, write_case)
from paths import REPO


# ------------------------------------------------------- read back by the apps


@pytest.mark.parametrize("dim,mode", [(2, "triangle"), (3, "tet")], ids=["triangle", "tet"])
@pytest.mark.skipif(not os.path.exists(INTERPOL), reason="interpol is not built")
def test_a_cleaned_case_reads_back_as_the_field_it_holds(tmp_path, dim, mode):
    """The whole path: a case written and cleaned by the tool, read by the app as
    a plain P2 field with --no-key, gives at its own probe points the velocity
    the Python evaluation of the same file gives. A dof table written by hand is
    exactly where a silent permutation would hide."""
    pytest.importorskip("h5py")
    X, cells = channel2d(6) if dim == 2 else channel3d(3)
    per = [True] * (dim - 1) + [False]
    t = D.Topo(X, cells, per)
    cfg = write_case(tmp_path / "in", X, cells, [smooth_noslip(t.node_x)], per)
    out = tmp_path / "out"
    r = subprocess.run([sys.executable, os.path.join(REPO, "python", "divfree", "divfree_clean.py"),
                        str(cfg), "--out", str(out), "--no-key"],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    g = interpol_probe(out / "dolfin_params.dat", mode, 60)
    pts = np.stack([g[a] for a in "xyz"[:dim]], axis=1)
    got = np.stack([g["u" + a] for a in "xyz"[:dim]], axis=1)
    case = D.read_case(out / "dolfin_params.dat")
    ct = D.Topo(case["X"], case["cells"], case["periodic"])
    U = D.DofTable(out / "u_0000.h5", "u", ct, case["cell_indices"]).values(out / "u_0000.h5")
    cell = np.full(len(pts), -1)
    for k, c in enumerate(ct.cells):
        lam = np.linalg.solve((ct.X[c[1:]] - ct.X[c[0]]).T, (pts - ct.X[c[0]]).T)
        inside = (lam >= -1e-12).all(axis=0) & (lam.sum(axis=0) <= 1 + 1e-12)
        cell[inside] = k
    assert (cell >= 0).all(), "a probe point is in no cell"
    want = D.p2_eval(ct, U, pts, cell)
    assert np.abs(got - want).max() < 1e-12 * np.abs(U).max()


# ------------------------------------------------------------ on several ranks


@pytest.mark.slow
def test_a_failed_write_of_the_first_rank_takes_the_job_down(tmp_path):
    """The mesh is linked into the output by the first rank alone, so a full
    disk or a permission there is a failure the other ranks cannot have: they
    would wait for it in the next collective forever and the job would look
    like one that is still running. It stops instead, naming what could not be
    written."""
    pytest.importorskip("h5py")
    X, cells = channel2d(6)
    t = D.Topo(X, cells, [True, False])
    cfg = write_case(tmp_path / "in", X, cells, [smooth_noslip(t.node_x)], [True, False])
    out = tmp_path / "out"
    out.mkdir()
    out.chmod(0o555)
    if os.access(out, os.W_OK):
        out.chmod(0o755)
        pytest.skip("this user writes where the mode says it cannot")
    try:
        r = mpi_clean(cfg, out, 2, ok=False, timeout=120)
    finally:
        out.chmod(0o755)
    assert r.returncode != 0
    assert "rank 0 could not" in r.stdout + r.stderr


@pytest.mark.slow
@pytest.mark.parametrize("dim,mesh", CASES, ids=IDS)
@pytest.mark.skipif(not os.path.exists(INTERPOL), reason="interpol is not built")
def test_three_ranks_clean_the_case_one_rank_cleans(tmp_path, dim, mesh):
    """Every rank holds, solves and writes its share. A case cleaned on three
    ranks balances to what a loader accepts, leaves every held node at the data's value, and
    reports each stamp once -- every printed line is the first rank's -- and the
    field is the one-rank run's, as the tool reads it and as the divergence-free
    loader does. Every number it prints is the one-rank run's too.

    Bit-for-bit agreement is not asked for and is not there: the ranks add their
    cell blocks in another order, so GAMG aggregates differently and the
    iteration counts differ. What is compared is the criterion the loader
    applies, and the fields to well inside it."""
    pytest.importorskip("h5py")
    X, cells = mesh
    per = [True] * (dim - 1) + [False]
    t = D.Topo(X, cells, per)
    fields = [smooth_noslip(t.node_x), 0.5 * smooth_noslip(t.node_x)]
    cfg = write_case(tmp_path / "in", X, cells, fields, per, stamps=["0", "1"])
    runs = {n: mpi_clean(cfg, tmp_path / ("out%d" % n), n) for n in (1, 3)}

    same_log(runs[1].stdout, runs[3].stdout)
    got = {}
    for n, r in runs.items():
        lines = [l for l in r.stdout.split("\n") if l.strip()]
        assert sum(l.startswith("  up_") for l in lines) == 2, "one line a stamp"
        assert sum("wrote 2 stamps" in l for l in lines) == 1
        out = tmp_path / ("out%d" % n)
        assert D.check_case(out / "dolfin_params.dat", quiet=True) <= D.FLUX_TOL
        case = D.read_case(out / "dolfin_params.dat")
        ct = D.Topo(case["X"], case["cells"], case["periodic"])
        got[n] = [D.DofTable(out / "u_0000.h5", "u", ct, case["cell_indices"]).values(
            out / name) for name in ("u_0000.h5", "u_0001.h5")]
    held = D.Topo(X, cells, per)
    held.set_held(D.at_rest_nodes(np.maximum(*[np.linalg.norm(f, axis=1) for f in fields]),
                                  max(float(np.abs(f).max()) for f in fields)))
    assert held.held_node.sum() > 0
    scale = max(float(np.abs(f).max()) for f in fields)
    for k, (a, b) in enumerate(zip(got[1], got[3])):
        # a held midpoint has no unknown, so no rank can write one
        for u in (a, b):
            assert np.abs(u[held.held_node] - fields[k][held.held_node]).max() == 0.0
        assert np.abs(a - b).max() < 1e-12 * scale
    # the loader reads what the first rank wrote under MPI
    mode = "triangle" if dim == 2 else "tet"
    assert "divfree=true" in (tmp_path / "out3" / "dolfin_params.dat").read_text()
    a, b = (interpol_probe(tmp_path / ("out%d" % n) / "dolfin_params.dat", mode, 200)
            for n in (1, 3))
    for c in "xyz"[:dim]:
        assert np.array_equal(a[c], b[c]), c
        assert np.abs(a["u" + c] - b["u" + c]).max() < 1e-8 * scale, c


@pytest.mark.slow
@pytest.mark.parametrize("partition", ["blocks", "ptscotch"])
@pytest.mark.parametrize("kind", ["closed", "moving wall"])
@pytest.mark.parametrize("dim", [2, 3])
def test_three_ranks_report_every_number_one_rank_reports(tmp_path, dim, kind, partition):
    """What the tool prints -- the counts of cells, edges, held nodes, free and
    boundary-free midpoints and changed edges, the fluxes, the volume means,
    the drift, the divergence left, and under --check the Taylor-Hood moments
    -- is the whole mesh's, each a reduction over the ranks, and a reduction
    over some of them moves a number by far more than round-off. So every line
    of the three-rank run is the one-rank run's, the numbers to round-off.

    The fixtures reach what the smooth channel does not: a pocket whose cells
    are measured against the floor of the stamp's largest facet flux, which on
    blocks one rank holds alone, so that rank's own largest is not the
    stamp's; a closed domain whose data leaves a net flux the step has to
    share out over every cell; a moving wall with free midpoints on it."""
    pytest.importorskip("h5py")
    cfg = pocket_case(tmp_path / "in", dim, kind)
    tool = os.path.join(REPO, "python", "divfree", "divfree_clean.py")
    runs, checks = {}, {}
    for n in (1, 3):
        runs[n] = mpi_clean(cfg, tmp_path / ("out%d" % n), n, extra=["--partition", partition],
                            timeout=120).stdout
        r = mpi_run(n, [tool, cfg, "--check", "--partition", partition], timeout=120)
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
        checks[n] = r.stdout
    free_on_wall = int(re.search(r"(\d+) of them on the boundary", runs[1]).group(1))
    assert (free_on_wall > 0) == (kind == "moving wall")
    assert ("singular" in runs[1]) == (kind == "closed")
    same_log(runs[1], runs[3])
    same_log(checks[1], checks[3])
    same_output(tmp_path / "out1", tmp_path / "out3", 1.25)
