"""python/divfree/divfree_clean.py on several ranks, each holding its share:
every refusal fires under MPI as it does on one rank, naming the same thing,
the two write routes give the same files, and a rank without a cell or
without an unknown of its own takes its part.

A refusal on several ranks has to be collective -- a rank that raised alone
would leave the others inside the next exchange -- so the cases of one rank
count are cleaned in turn by one MPI job that catches each refusal on every
rank and prints the first rank's, and each message is compared with the one
this process gives, where COMM_WORLD is one rank. What is refused is decided
from the whole mesh: a ragged or unknown row, a table longer than the mesh or
with fewer labels than rows, a cell no row names, a dof past the values, a node
no cell gives a value, a node two cells give two values -- within a rank or
across two, the owner deciding --, a non-finite value by its dof, a mesh whose
every midpoint is held, and a step that leaves a cell unbalanced, named by its
serial id. `--split` stays serial and is refused on several ranks.

The jobs run this file's cases through divfree_cases.mpi_batch; run as a
script, it is the job of a failure on one rank alone. The jobs skip where no
launcher starts a job of several ranks for mpi4py, and fail there with
PARTRAC_REQUIRE_MPI set.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import numpy as np                 # noqa: E402
import pytest                      # noqa: E402

from divfree_cases import (D, channel2d, channel3d, cleaner_main, each_case,  # noqa: E402
                           launcher_or_skip, mpi_batch, mpi_clean, mpi_or_skip, mpi_run, one_rank_main,
                           read_back, same_output, smooth_noslip, two_cells,
                           two_stamp_channel, write_case)
from test_divfree_files import poke_checkpoint   # noqa: E402

from mpi4py import MPI             # noqa: E402


# ------------------------------------------------------------------ the jobs


def refused(cfg, out, partition="blocks"):
    """Clean the case; the refusal it raised, or that there was none."""
    try:
        D.clean_case(cfg, out, verbose=False, partition=partition)
        return "CLEANED"
    except (ValueError, RuntimeError) as exc:
        return "REFUSED %s" % exc


def no_step(self, U, rhs, who, facet):
    """Equil._step taken to be zero, so the data's own worst cell is the one
    the balance refusal has to name."""
    return np.zeros(len(self.gcol)), 0, ""


STEP = D.Equil._step


def refusal_case(cfg, out, partition, unbalanced):
    """A case of mpi_batch: clean the case, with no step where unbalanced, and
    print the refusal it raised, or that there was none."""
    D.Equil._step = no_step if unbalanced else STEP
    print(refused(cfg, out, partition))


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


def csr_case():
    """A case of mpi_batch: a matrix block past an 8-bit PETSc's index range
    on the first rank alone, handed to PETSc on every rank; print what the
    rank said."""
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
    print(said)


def cache_case(cfg, out):
    """A case of mpi_batch: clean the case with a stamp cache whose budget a
    rank lies between the ranks' shares of a stamp: the smallest share fits it
    and the largest does not."""
    comm = MPI.COMM_WORLD
    case = D.read_case(cfg, mesh=False)
    topo = D._read_topology(case, comm, "ptscotch")
    share = np.array(comm.allgather(topo.nnodes * topo.dim * 8))
    assert share.min() < share.max(), share
    D.READ_CACHE_BYTES = int((share.min() + share.max()) // 2) * comm.size
    D.clean_case(cfg, out, verbose=False, partition="ptscotch")


if __name__ == "__main__" and sys.argv[1:2] == ["fault"]:
    # nothing here catches it: the tool itself has to end the job
    sys.exit(fault_job(*sys.argv[2:]))


# --------------------------------------------------------------- the cases


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


# each refusal at two or three ranks, in turn; the case built for two at two
KIND_RANKS = [(k, 2 if k == "cells disagree across ranks" else 2 + i % 2)
              for i, k in enumerate(sorted(REFUSED))]


@pytest.mark.slow
@pytest.mark.parametrize("ranks", [2, 3])
def test_a_refusal_fires_on_every_rank_and_names_what_one_rank_names(tmp_path, ranks):
    """Each rank reads only its rows and values, so what is wrong may be in
    another rank's share: the verdict is reduced over the ranks and every rank
    refuses together, with the one-rank run's message -- the same file, the same
    dataset, the same dof. Every refusal of this rank count, one job."""
    pytest.importorskip("h5py")
    launcher_or_skip()
    kinds = [k for k, n in KIND_RANKS if n == ranks]
    folder = {k: tmp_path / k.replace(" ", "_") for k in kinds}
    cfg = {k: broken(folder[k], k) for k in kinds}
    one = {k: refused(cfg[k], folder[k] / "out1") for k in kinds}
    many = mpi_batch(ranks, refusal_case,
                     {k: [cfg[k], folder[k] / "outn", "blocks", False] for k in kinds},
                     tmp_path / "jobs", timeout=300)

    def check(kind):
        assert one[kind].startswith("REFUSED") and REFUSED[kind] in one[kind], one[kind]
        # every rank refused, with the one-rank run's message
        assert many[kind] == [one[kind] + "\n"] * ranks, many[kind]
    each_case(kinds, check)


@pytest.mark.slow
def test_the_balance_refusal_names_the_serial_cell(tmp_path, monkeypatch):
    """A step that leaves a cell above what a loader accepts is refused naming
    the cell; on several ranks the cells are numbered by rank, so the name has
    to be the cell's serial id, the row of mesh/topology, and the worst cell
    the whole mesh's. With no step at all the worst cell is the data's, the
    same at every rank count. Where several cells are exactly as bad -- a field
    of z alone on a grid whose coordinates are exact binary fractions, so every
    translated copy of a cell has the same bits -- the one named is the one
    with the smallest serial id, whichever ranks hold the others. Both fields
    under both partitions, one job of three ranks."""
    pytest.importorskip("h5py")
    launcher_or_skip()
    monkeypatch.setattr(D.Equil, "_step", no_step)
    cfg, one = {}, {}
    for field in ("smooth", "tied"):
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
        cfg[field] = write_case(tmp_path / field / "in", X, cells, [U], per)
        one[field] = refused(cfg[field], tmp_path / field / "out1")
        assert one[field].startswith(
            "REFUSED up_0.h5: the flux solve left cell %d " % worst.min()), one[field]
    cases = {"%s-%s" % (f, p): (f, p) for f in ("smooth", "tied") for p in ("blocks", "ptscotch")}
    many = mpi_batch(3, refusal_case,
                     {n: [cfg[f], tmp_path / n, p, True] for n, (f, p) in cases.items()},
                     tmp_path / "jobs", timeout=300)

    def check(name):
        assert many[name] == [one[cases[name][0]] + "\n"] * 3, many[name]
    each_case(cases, check)


@pytest.mark.slow
def test_split_is_refused_on_several_ranks(tmp_path):
    """--split writes the split mesh and its field from one rank; on several it
    is refused, saying so, before anything is read."""
    pytest.importorskip("h5py")
    launcher_or_skip()
    X, cells = channel2d(4)
    t = D.Topo(X, cells, [True, False])
    cfg = write_case(tmp_path / "in", X, cells, [smooth_noslip(t.node_x)], [True, False])
    r = mpi_clean(cfg, tmp_path / "out", 2, extra=["--split"], ok=False, timeout=120)
    assert r.returncode != 0
    assert "--split writes the split mesh and its field from one rank, and this job has 2" \
        in r.stdout + r.stderr


@pytest.mark.slow
def test_both_write_routes_write_the_same_case_and_say_which(tmp_path, monkeypatch,
                                                              channel_runs):
    """The rows go to their blocks and are written by every rank through
    parallel HDF5, or gathered to the first rank where PARTRAC_HDF5_GATHER is
    set; the log says which, once. Both routes give the same datasets, the
    integers identical; the solve on several ranks is not bit-reproducible run
    to run, so the floats agree to round-off. Both routes put the rows in
    place the same way, so each is also read back by node through the tool's
    own reader and compared with the one-rank run's field: values written at
    another node's row would agree between the routes and read back wrong.
    The one-rank run and the writer's own route's are channel_runs'."""
    h5py = pytest.importorskip("h5py")
    d = channel_runs(3)
    fields = two_stamp_channel(3)[3]
    assert (d / "one.log").read_text().count("wrote 2 stamps") == 1
    monkeypatch.setenv("PARTRAC_HDF5_GATHER", "1")
    out = {"": d / "three", "1": tmp_path / "out1"}
    logs = {"": (d / "three.log").read_text(),
            "1": mpi_clean(d / "in" / "dolfin_params.dat", out["1"], 3).stdout}
    routes = {}
    for gather, log in logs.items():
        said = [l for l in log.split("\n") if l.startswith("  wrote 2 stamps")]
        assert len(said) == 1, log
        routes[gather] = said[0]
    assert routes["1"].endswith("gathered to rank 0")
    if routes[""].endswith("gathered to rank 0"):
        mpi_or_skip("the h5py of the jobs has no MPI, so both runs were gathered")
    assert routes[""].endswith("through parallel HDF5")
    scale = max(float(np.abs(f).max()) for f in fields)
    for name in ("u_0000.h5", "u_0001.h5"):
        with h5py.File(out[""] / name, "r") as a, h5py.File(out["1"] / name, "r") as b:
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
    want = read_back(d / "one", names)
    for gather in out:
        for a, b in zip(want, read_back(out[gather], names)):
            assert np.abs(a - b).max() < 1e-12 * scale, gather


def test_the_gather_is_forced_by_one_or_true_and_by_nothing_else(monkeypatch):
    """PARTRAC_HDF5_GATHER=0 or =false reads as asking for no gather, so only 1
    and true, in any case, force it. Which route the log names is the write
    routes' test's."""
    from divfree_write import gather_forced
    for value, forced in (("", False), ("0", False), ("false", False), ("yes", False),
                          ("1", True), ("True", True), (" TRUE ", True)):
        monkeypatch.setenv("PARTRAC_HDF5_GATHER", value)
        assert gather_forced() == forced, value


@pytest.mark.slow
def test_a_failure_on_one_rank_alone_takes_the_job_down(tmp_path):
    """Whatever raises on one rank and not on the others -- a read, a bug, a
    full disk outside the writer -- leaves the others waiting in the next
    exchange for a rank that has gone, and the job would sit there looking
    busy. The command line ends it instead, the failing rank saying what it
    was."""
    pytest.importorskip("h5py")
    launcher_or_skip()
    X, cells = channel2d(4)
    t = D.Topo(X, cells, [True, False])
    cfg = write_case(tmp_path / "in", X, cells, [smooth_noslip(t.node_x)], [True, False])
    r = mpi_run(2, [os.path.abspath(__file__), "fault", cfg, tmp_path / "out"], timeout=300)
    assert r.returncode != 0
    assert "rank 1: OSError: a read only the last rank fails" in r.stderr, r.stderr[-3000:]


@pytest.mark.slow
@pytest.mark.parametrize("ranks", [2, 3])
def test_ranks_with_no_cell_or_no_unknown_clean_what_one_rank_cleans(tmp_path, ranks):
    """A rank may own no unknown -- its cells reach only unknowns whose first
    cell is another rank's -- or hold no cell at all, where there are more
    ranks than cells. It still takes its part in every exchange and in the
    assembly, with no rows of its own, and the case comes out as one rank
    cleans it. PETSc takes an empty preallocation array for something else
    than no rows, and the assembly then waits forever. Both partitions, one
    job."""
    pytest.importorskip("h5py")
    launcher_or_skip()
    cfg = two_cells(tmp_path / "in")
    one_rank_main([cfg, "--out", tmp_path / "out1"])
    partitions = ["blocks", "ptscotch"]
    mpi_batch(ranks, cleaner_main,
              {p: [cfg, "--out", tmp_path / p, "--partition", p] for p in partitions},
              tmp_path / "jobs")
    each_case(partitions, lambda p: same_output(tmp_path / "out1", tmp_path / p, 1.0))


@pytest.mark.slow
def test_a_matrix_past_the_index_range_on_one_rank_is_refused_on_every_rank(tmp_path):
    """Each rank hands PETSc its own block of K and C, and the blocks differ in
    size, so a block past PETSc's integer range can be one rank's alone. A rank
    that refused alone would leave the others in the assembly forever: every
    rank refuses, with the same message, naming the largest block."""
    said = [s.strip() for s in mpi_batch(3, csr_case, {"csr": []}, tmp_path, timeout=300)["csr"]]
    assert len(said) == 3 and len(set(said)) == 1, said
    assert said[0].startswith("REFUSED") and "400 nonzeros" in said[0], said[0]
    assert "past the 127 that this PETSc's 8-bit indices hold" in said[0], said[0]


@pytest.mark.slow
def test_a_stamp_is_kept_for_the_second_pass_by_every_rank_or_by_none(tmp_path, channel_runs):
    """The values read in the wall pass are kept for the cleaning pass where
    they fit the rank's part of the budget, and the shares differ by rank. A
    rank that kept a stamp another rank did not would skip the collective read
    the other makes, and the two would pair up the wrong exchanges. With the
    budget between the smallest and the largest share the case still cleans,
    and to what one rank writes, channel_runs' run."""
    pytest.importorskip("h5py")
    d = channel_runs(3)
    fields = two_stamp_channel(3)[3]
    mpi_batch(3, cache_case, {"cache": [d / "in" / "dolfin_params.dat", tmp_path / "three"]},
              tmp_path / "jobs", timeout=120)
    same_output(d / "one", tmp_path / "three", float(np.abs(fields[0]).max()))
