"""python/divfree/divfree_clean.py on several ranks, each holding its share:
every refusal fires under MPI as it does on one rank, naming the same thing,
the two write routes give the same files, and a rank without a cell or
without an unknown of its own takes its part.

A refusal on several ranks has to be collective -- a rank that raised alone
would leave the others inside the next exchange -- so each case is cleaned as
an MPI job that catches the refusal on every rank and prints the first rank's,
and the message is compared with the one-rank job's. What is refused is decided
from the whole mesh: a ragged or unknown row, a table longer than the mesh or
with fewer labels than rows, a cell no row names, a dof past the values, a node
no cell gives a value, a node two cells give two values -- within a rank or
across two, the owner deciding --, a non-finite value by its dof, a mesh whose
every midpoint is held, and a step that leaves a cell unbalanced, named by its
serial id. `--split` stays serial and is refused on several ranks.

Run as a script, this file is the MPI job. The jobs skip where no launcher
starts a job of several ranks for mpi4py, and fail there with
PARTRAC_REQUIRE_MPI set.
"""

import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import numpy as np                 # noqa: E402
import pytest                      # noqa: E402

from divfree_cases import (D, channel2d, channel3d, mpi_clean, mpi_run,  # noqa: E402
                           read_back, same_log, same_output, smooth_noslip,
                           two_cells, write_case)
from paths import REPO             # noqa: E402
from test_divfree_files import poke_checkpoint   # noqa: E402

from mpi4py import MPI             # noqa: E402


# ------------------------------------------------------------------ the jobs


def job(cfg, out, partition, unbalanced):
    """Clean the case; print the refusal every rank raised, or that none did.
    With unbalanced the step is taken to be zero, so the data's own worst cell
    is the one the balance refusal has to name."""
    comm = MPI.COMM_WORLD
    if unbalanced == "1":
        D.Equil._step = lambda self, U, rhs, who, facet: (np.zeros(len(self.gcol)), 0, "")
    try:
        D.clean_case(cfg, out, verbose=False, partition=partition)
        said = "CLEANED"
    except (ValueError, RuntimeError) as exc:
        said = "REFUSED %s" % exc
    every = comm.gather(said, root=0)
    if comm.rank == 0:
        assert len(set(every)) == 1, every
        print(said)


def fault_job(cfg, out):
    """The command line with a read that fails on the last rank alone, while
    the others go on into the exchange that reads it with them."""
    comm = MPI.COMM_WORLD
    real = D.DofTable.values

    def values(self, path):
        if comm.rank == comm.size - 1:
            raise OSError("a read only the last rank fails")
        return real(self, path)
    D.DofTable.values = values
    return D.main([cfg, "--out", out])


def csr_job():
    """A matrix block past an 8-bit PETSc's index range on the first rank
    alone, handed to PETSc on every rank; print what every rank said."""
    import types
    import scipy.sparse as sp
    comm = MPI.COMM_WORLD
    n = 20 if comm.rank == 0 else 5
    A = sp.csr_matrix(np.ones((n, n)))
    small = types.SimpleNamespace(IntType=np.int8)
    try:
        D._petsc_csr(small, A, comm)
        said = "ACCEPTED"
    except ValueError as exc:
        said = "REFUSED %s" % exc
    every = comm.gather(said, root=0)
    if comm.rank == 0:
        print("\n".join(every))


def cache_job(cfg, out):
    """Clean the case with a stamp cache whose budget a rank lies between the
    ranks' shares of a stamp: the smallest share fits it and the largest does
    not."""
    comm = MPI.COMM_WORLD
    case = D.read_case(cfg, mesh=False)
    topo = D._read_topology(case, comm, "ptscotch")
    share = np.array(comm.allgather(topo.nnodes * topo.dim * 8))
    assert share.min() < share.max(), share
    D.READ_CACHE_BYTES = int((share.min() + share.max()) // 2) * comm.size
    D.clean_case(cfg, out, verbose=False, partition="ptscotch")


JOBS = {"refuse": job, "csr": csr_job, "cache": cache_job}


if __name__ == "__main__":
    if sys.argv[1] == "fault":
        # nothing here catches it: the tool itself has to end the job
        sys.exit(fault_job(*sys.argv[2:]))
    try:
        JOBS[sys.argv[1]](*sys.argv[2:])
    except BaseException:
        sys.stderr.write("rank %d:\n%s" % (MPI.COMM_WORLD.rank, traceback.format_exc()))
        sys.stderr.flush()
        MPI.COMM_WORLD.Abort(1)


# --------------------------------------------------------------- the cases


def run_job(cfg, out, ranks, partition="blocks", unbalanced=False):
    """This file's refusal job; what it printed."""
    r = mpi_run(ranks, [os.path.abspath(__file__), "refuse", cfg, out, partition,
                        "1" if unbalanced else "0"], timeout=120)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return r.stdout.strip().split("\n")[-1]


def two_rank_disagreement(folder):
    """Cells that disagree across two ranks and nowhere else: every cell of the
    second half of the rows names another dof for one vertex it shares with the
    first half, so each rank's cells agree among themselves and only the
    vertex's owner sees the two values."""
    import h5py
    X, cells, _ = D.read_mesh_h5(folder / "mesh.h5")
    half = len(cells) // 2
    v = np.intersect1d(cells[:half].ravel(), cells[half:].ravel())[0]
    other = [w for w in range(len(X)) if w not in cells[half:].ravel()][0]
    with h5py.File(folder / "up_0.h5", "r+") as f:
        cd = np.array(f["u/cell_dofs"]).reshape(len(cells), -1)
        vec = np.array(f["u/vector_0"])
        for k in range(half, len(cells)):
            for j in np.nonzero(cells[k] == v)[0]:
                # the first component of vertex j, the columns component by component
                assert vec[cd[k, j]] != vec[other * 2]
                cd[k, j] = other * 2
        f["u/cell_dofs"][...] = cd.ravel()


def all_held(folder):
    """Three triangles that share nothing, at rest: every midpoint is on an
    exterior facet and held."""
    X = np.array([[0, 0], [1, 0], [0, 1], [2, 0], [3, 0], [2, 1], [4, 0], [5, 0], [4, 1]],
                 float)
    cells = np.arange(9).reshape(3, 3)
    t = D.Topo(X, cells, [False, False])
    return write_case(folder, X, cells, [np.zeros((t.nnodes, 2))], [False, False])


def broken(tmp_path, kind):
    """The case of one refusal, written into tmp_path/in."""
    folder = tmp_path / "in"
    if kind == "every midpoint held":
        return all_held(folder)
    X, cells = channel2d(6)
    t = D.Topo(X, cells, [True, False])
    U = smooth_noslip(t.node_x)
    if kind.startswith("cells disagree"):
        # no two nodes alike, as a wall at rest would make them
        U = np.random.default_rng(12).normal(size=U.shape)
    fields = [U, 0.5 * U]
    if kind == "non-finite":
        fields[1] = fields[1].copy()
        fields[1][t.nverts + 17, 1] = np.inf
    # poke_checkpoint adds its vertex in no cell to a closed square
    cfg = write_case(folder, X, cells, fields, [True, False] if kind != "node without a value"
                     else [False, False])
    if kind in ("ragged", "unknown cell", "cells disagree", "node without a value",
                "a row past the mesh", "short cells", "a cell named twice",
                "a dof past the values"):
        poke_checkpoint(folder, kind, len(cells))
    if kind == "cells disagree across ranks":
        two_rank_disagreement(folder)
    return cfg


REFUSED = {
    "ragged": "'u/x_cell_dofs' is ragged",
    "unknown cell": "'u/cells' names a cell the mesh does not",
    "node without a value": "up_0.h5: 'u' leaves a node without a value",
    "cells disagree": "up_0.h5: cells disagree on a node of 'u'",
    "cells disagree across ranks": "up_0.h5: cells disagree on a node of 'u'",
    "non-finite": "up_1.h5: 'u' holds a non-finite value at dof",
    "every midpoint held": "every midpoint of the 3 cells is held",
    "a row past the mesh": "'u/cells' names a cell the mesh does not",
    "short cells": "'u/x_cell_dofs' is ragged",
    "a cell named twice": "'u/cells' leaves a cell of the mesh without a row",
    "a dof past the values": "up_0.h5: 'u' holds 337 values, fewer than its dof table names",
}


@pytest.mark.slow
@pytest.mark.parametrize("ranks", [2, 3])
@pytest.mark.parametrize("kind", sorted(REFUSED))
def test_a_refusal_fires_on_every_rank_and_names_what_one_rank_names(tmp_path, kind, ranks):
    """Each rank reads only its rows and values, so what is wrong may be in
    another rank's share: the verdict is reduced over the ranks and every rank
    refuses together, with the one-rank run's message -- the same file, the same
    dataset, the same dof."""
    pytest.importorskip("h5py")
    cfg = broken(tmp_path, kind)
    one = run_job(cfg, tmp_path / "out1", 1)
    many = run_job(cfg, tmp_path / "outn", ranks)
    assert one.startswith("REFUSED") and REFUSED[kind] in one, one
    assert many == one


@pytest.mark.slow
@pytest.mark.parametrize("partition", ["blocks", "ptscotch"])
@pytest.mark.parametrize("field", ["smooth", "tied"])
def test_the_balance_refusal_names_the_serial_cell(tmp_path, field, partition):
    """A step that leaves a cell above what a loader accepts is refused naming
    the cell; on several ranks the cells are numbered by rank, so the name has
    to be the cell's serial id, the row of mesh/topology, and the worst cell
    the whole mesh's. With no step at all the worst cell is the data's, the
    same at every rank count. Where several cells are exactly as bad -- a field
    of z alone on a grid whose coordinates are exact binary fractions, so every
    translated copy of a cell has the same bits -- the one named is the one
    with the smallest serial id, whichever ranks hold the others."""
    pytest.importorskip("h5py")
    X, cells = channel3d(3 if field == "smooth" else 4)
    per = [True, True, False]
    t = D.Topo(X, cells, per)
    U = smooth_noslip(t.node_x)
    if field == "tied":
        U = np.zeros_like(U)
        U[:, 2] = t.node_x[:, 2] * (1 - t.node_x[:, 2])
    ratio, _ = D.flux_ratios(t, U)
    worst = np.nonzero(ratio == ratio.max())[0]
    assert len(worst) > 1 or field == "smooth"
    cfg = write_case(tmp_path / "in", X, cells, [U], per)
    one = run_job(cfg, tmp_path / "out1", 1, unbalanced=True)
    assert one.startswith("REFUSED up_0.h5: the flux solve left cell %d " % worst.min()), one
    for ranks in (2, 3):
        assert run_job(cfg, tmp_path / ("out%d" % ranks), ranks, partition,
                       unbalanced=True) == one


@pytest.mark.slow
def test_split_is_refused_on_several_ranks(tmp_path):
    """--split writes the split mesh and its field from one rank; on several it
    is refused, saying so, before anything is read."""
    pytest.importorskip("h5py")
    X, cells = channel2d(4)
    t = D.Topo(X, cells, [True, False])
    cfg = write_case(tmp_path / "in", X, cells, [smooth_noslip(t.node_x)], [True, False])
    r = mpi_clean(cfg, tmp_path / "out", 2, extra=["--split"], ok=False, timeout=120)
    assert r.returncode != 0
    assert "--split writes the split mesh and its field from one rank, and this job has 2" \
        in r.stdout + r.stderr


@pytest.mark.slow
def test_both_write_routes_write_the_same_case_and_say_which(tmp_path, monkeypatch):
    """The rows go to their blocks and are written by every rank through
    parallel HDF5, or gathered to the first rank where PARTRAC_HDF5_GATHER is
    set; the log says which, once. Both routes give the same datasets, the
    integers identical; the solve on several ranks is not bit-reproducible run
    to run, so the floats agree to round-off. Both routes put the rows in
    place the same way, so each is also read back by node through the tool's
    own reader and compared with the one-rank run's field: values written at
    another node's row would agree between the routes and read back wrong."""
    h5py = pytest.importorskip("h5py")
    X, cells = channel3d(3)
    t = D.Topo(X, cells, [True, True, False])
    phi = [np.ones((t.nnodes, 1)), 2 * np.ones((t.nnodes, 1))]
    fields = [smooth_noslip(t.node_x), 0.5 * smooth_noslip(t.node_x)]
    cfg = write_case(tmp_path / "in", X, cells, fields, [True, True, False], phi=phi)
    one = mpi_clean(cfg, tmp_path / "one", 1)
    assert one.stdout.count("wrote 2 stamps") == 1
    routes = {}
    for gather in ("", "1"):
        monkeypatch.setenv("PARTRAC_HDF5_GATHER", gather)
        r = mpi_clean(cfg, tmp_path / ("out" + gather), 3)
        said = [l for l in r.stdout.split("\n") if l.startswith("  wrote 2 stamps")]
        assert len(said) == 1, r.stdout
        routes[gather] = said[0]
    assert routes["1"].endswith("gathered to rank 0")
    if routes[""].endswith("gathered to rank 0"):
        pytest.skip("the h5py of the jobs has no MPI, so both runs were gathered")
    assert routes[""].endswith("through parallel HDF5")
    scale = max(float(np.abs(f).max()) for f in fields)
    for name in ("u_0000.h5", "u_0001.h5"):
        with h5py.File(tmp_path / "out" / name, "r") as a, \
                h5py.File(tmp_path / "out1" / name, "r") as b:
            got = []
            a.visititems(lambda n, o: got.append(n) if isinstance(o, h5py.Dataset) else None)
            assert "phi/vector_0" in got and "p/vector_0" in got
            for n in got:
                x, y = a[n][()], b[n][()]
                assert x.dtype == y.dtype and x.shape == y.shape, n
                if x.dtype.kind == "f":
                    assert np.abs(x - y).max() < 1e-12 * scale, n
                else:
                    assert np.array_equal(x, y), n
                assert dict(a[n].attrs).keys() == dict(b[n].attrs).keys(), n
    names = ["u_0000.h5", "u_0001.h5"]
    want = read_back(tmp_path / "one", names)
    for out in ("out", "out1"):
        for a, b in zip(want, read_back(tmp_path / out, names)):
            assert np.abs(a - b).max() < 1e-12 * scale, out


@pytest.mark.slow
def test_the_gather_is_forced_by_one_or_true_and_by_nothing_else(tmp_path, monkeypatch):
    """PARTRAC_HDF5_GATHER=0 or =false reads as asking for no gather, so only 1
    and true, in any case, force it; the log says which route was taken."""
    pytest.importorskip("h5py")
    X, cells = channel2d(4)
    t = D.Topo(X, cells, [True, False])
    cfg = write_case(tmp_path / "in", X, cells, [smooth_noslip(t.node_x)], [True, False])
    said = {}
    for value in ("", "0", "false", "1", "True"):
        monkeypatch.setenv("PARTRAC_HDF5_GATHER", value)
        r = mpi_clean(cfg, tmp_path / ("out" + value), 2)
        said[value] = r.stdout.rstrip().split("\n")[-1]
    if said[""].endswith("gathered to rank 0"):
        pytest.skip("the h5py of the jobs has no MPI, so every run was gathered")
    for value, route in (("0", "through parallel HDF5"), ("false", "through parallel HDF5"),
                         ("1", "gathered to rank 0"), ("True", "gathered to rank 0")):
        assert said[value].endswith(route), (value, said[value])


@pytest.mark.slow
def test_the_check_on_several_ranks_reports_what_one_rank_reports(tmp_path):
    """--check reads and measures its share on every rank and reduces what it
    reports: the worst cell, which is the whole mesh's, is the one-rank run's
    to the bit, each stamp is reported once, and every line is the one-rank
    run's, the numbers to round-off."""
    pytest.importorskip("h5py")
    X, cells = channel3d(3)
    t = D.Topo(X, cells, [True, True, False])
    U = smooth_noslip(t.node_x)
    cfg = write_case(tmp_path / "in", X, cells, [U, 0.5 * U], [True, True, False],
                     stamps=["0", "1"])
    said = {}
    for n in (1, 3):
        r = mpi_run(n, [os.path.join(REPO, "python", "divfree", "divfree_clean.py"), cfg,
                        "--check"], timeout=120)
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
        lines = r.stdout.rstrip().split("\n")
        assert sum(l.startswith("  0 up_0.h5:") for l in lines) == 1, r.stdout
        assert sum(l.startswith("  1 up_1.h5:") for l in lines) == 1, r.stdout
        said[n] = r.stdout
    assert said[1].rstrip().split("\n")[-1].startswith("worst cell balance")
    assert said[3].rstrip().split("\n")[-1] == said[1].rstrip().split("\n")[-1]
    same_log(said[1], said[3])


@pytest.mark.slow
def test_a_failure_on_one_rank_alone_takes_the_job_down(tmp_path):
    """Whatever raises on one rank and not on the others -- a read, a bug, a
    full disk outside the writer -- leaves the others waiting in the next
    exchange for a rank that has gone, and the job would sit there looking
    busy. The command line ends it instead, the failing rank saying what it
    was."""
    pytest.importorskip("h5py")
    X, cells = channel2d(4)
    t = D.Topo(X, cells, [True, False])
    cfg = write_case(tmp_path / "in", X, cells, [smooth_noslip(t.node_x)], [True, False])
    r = mpi_run(2, [os.path.abspath(__file__), "fault", cfg, tmp_path / "out"], timeout=60)
    assert r.returncode != 0
    assert "rank 1: OSError: a read only the last rank fails" in r.stderr, r.stderr[-3000:]


@pytest.mark.slow
@pytest.mark.parametrize("partition", ["blocks", "ptscotch"])
@pytest.mark.parametrize("ranks", [2, 3])
def test_ranks_with_no_cell_or_no_unknown_clean_what_one_rank_cleans(tmp_path, ranks,
                                                                    partition):
    """A rank may own no unknown -- its cells reach only unknowns whose first
    cell is another rank's -- or hold no cell at all, where there are more
    ranks than cells. It still takes its part in every exchange and in the
    assembly, with no rows of its own, and the case comes out as one rank
    cleans it. PETSc takes an empty preallocation array for something else
    than no rows, and the assembly then waits forever."""
    pytest.importorskip("h5py")
    cfg = two_cells(tmp_path / "in")
    mpi_clean(cfg, tmp_path / "out1", 1)
    mpi_clean(cfg, tmp_path / "outn", ranks, extra=["--partition", partition], timeout=120)
    same_output(tmp_path / "out1", tmp_path / "outn", 1.0)


@pytest.mark.slow
def test_a_matrix_past_the_index_range_on_one_rank_is_refused_on_every_rank():
    """Each rank hands PETSc its own block of K and C, and the blocks differ in
    size, so a block past PETSc's integer range can be one rank's alone. A rank
    that refused alone would leave the others in the assembly forever: every
    rank refuses, with the same message, naming the largest block."""
    r = mpi_run(3, [os.path.abspath(__file__), "csr"], timeout=60)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    said = r.stdout.strip().split("\n")
    assert len(said) == 3 and len(set(said)) == 1, said
    assert said[0].startswith("REFUSED") and "400 nonzeros" in said[0], said[0]
    assert "past the 127 that this PETSc's 8-bit indices hold" in said[0], said[0]


@pytest.mark.slow
def test_a_stamp_is_kept_for_the_second_pass_by_every_rank_or_by_none(tmp_path):
    """The values read in the wall pass are kept for the cleaning pass where
    they fit the rank's part of the budget, and the shares differ by rank. A
    rank that kept a stamp another rank did not would skip the collective read
    the other makes, and the two would pair up the wrong exchanges. With the
    budget between the smallest and the largest share the case still cleans,
    and to what one rank writes."""
    pytest.importorskip("h5py")
    X, cells = channel3d(3)
    t = D.Topo(X, cells, [True, True, False])
    fields = [smooth_noslip(t.node_x), 0.5 * smooth_noslip(t.node_x)]
    cfg = write_case(tmp_path / "in", X, cells, fields, [True, True, False])
    mpi_clean(cfg, tmp_path / "one", 1)
    r = mpi_run(3, [os.path.abspath(__file__), "cache", cfg, tmp_path / "three"], timeout=120)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    same_output(tmp_path / "one", tmp_path / "three", float(np.abs(fields[0]).max()))
