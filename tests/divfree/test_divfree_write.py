"""python/divfree/divfree_write.py: the cleaner's output written by every rank.

Each rank hands the writer its contiguous block of the cells' dof rows and of
the node values, and the files must be the ones the serial writers give,
dataset by dataset and attribute by attribute, whichever route the writer
takes: parallel HDF5, or everything gathered to the first rank. The blocks may
be uneven and a rank may hold none; a field appended to a file and a group
copied through are written the same way; a write that fails on one rank takes
the job down instead of leaving the others inside a collective; what the
writer wrote reads back through the tool's reader and the loader.

Files are compared by their datasets and attributes, never byte for byte: HDF5
stores modification times in its object headers. The MPI tests run as jobs of
this module (`python3 test_divfree_write.py job ...`), skip where no launcher
starts a job of several ranks for mpi4py, and fail there with
PARTRAC_REQUIRE_MPI set.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from divfree_cases import (D, INTERPOL, channel2d, channel3d, interpol_probe,
                           mpi_launcher, relabel_by_global_id, smooth_noslip,
                           write_case)

import divfree_write as W          # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


# ------------------------------------------------------------------ the cases


def make_case(folder, dim):
    """A two-stamp case with a P2 velocity, a P1 pressure and a P2 phase field,
    its cells labelled by a permuted global id as a run on several ranks writes
    them."""
    X, cells = channel2d(6) if dim == 2 else channel3d(3)
    per = [True] * (dim - 1) + [False]
    t = D.Topo(X, cells, per)
    U = smooth_noslip(t.node_x)
    phi = [np.sin(3 * t.node_x[:, :1]), np.cos(2 * t.node_x[:, -1:])]
    cfg = write_case(folder, X, cells, [U, 0.5 * U], per, stamps=["0", "1"], phi=phi)
    relabel_by_global_id(folder)
    return cfg


def _read(cfg):
    """(case, topology, the fields of every stamp file by name)."""
    case = D.read_case(cfg)
    t = D.Topo(case["X"], case["cells"], case["periodic"])
    first = case["folder"] / case["stamps"][0][1]
    names = [n for _, n in case["stamps"]]
    tables = {f: D.DofTable(first, f, t, case["cell_indices"]) for f in ("u", "phi")}
    vals = {n: {f: d.values(case["folder"] / n) for f, d in tables.items()} for n in names}
    return case, t, vals


def serial_output(cfg, out):
    """The serial writers' files of the case: the velocity and the phase field
    written, the phase field appended to the velocity's file, the pressure
    copied through."""
    case, t, vals = _read(cfg)
    out.mkdir(parents=True, exist_ok=True)
    dim = t.dim
    for i, (name, v) in enumerate(vals.items()):
        if i == 0:
            D.write_checkpoint(out / name, "u", v["u"], t, 2, cell_indices=case["cell_indices"])
            D.write_checkpoint(out / name, "phi", v["phi"], t, 2, mode="a",
                               cell_indices=case["cell_indices"])
        else:
            D.write_vector(out / name, "u", v["u"], 2, dim, dim)
            D.write_vector(out / name, "phi", v["phi"], 2, 1, dim, mode="a")
        assert D.copy_group(case["folder"] / name, out / name, "p")
    (out / "timestamps.dat").write_text((case["folder"] / "timestamps.dat").read_text())
    return list(vals)


def block_output(cfg, out, writer, cells, nodes):
    """The same files from this rank's blocks: cells [c0, c1) and nodes
    [n0, n1) in serial order. Returns what copy_group said of a group the
    input does not hold, stamp by stamp."""
    case, t, vals = _read(cfg)
    out.mkdir(parents=True, exist_ok=True)
    dim = t.dim
    (c0, c1), (n0, n1) = cells(t.ncells), nodes(t.nnodes)
    gid = case["cell_indices"].astype(np.uint64)[c0:c1]
    missing = []
    for i, (name, v) in enumerate(vals.items()):
        if i == 0:
            writer.write_checkpoint(out / name, "u", D._cell_dofs(t, 2, dim)[c0:c1], gid,
                                    v["u"][n0:n1], 2, dim, dim)
            writer.write_checkpoint(out / name, "phi", D._cell_dofs(t, 2, 1)[c0:c1], gid,
                                    v["phi"][n0:n1], 2, 1, dim, mode="a")
        else:
            writer.write_vector(out / name, "u", v["u"][n0:n1], 2, dim, dim)
            writer.write_vector(out / name, "phi", v["phi"][n0:n1], 2, 1, dim, mode="a")
        assert writer.copy_group(case["folder"] / name, out / name, "p")
        missing.append(writer.copy_group(case["folder"] / name, out / name, "nothing"))
    # every rank returns once the first has written it
    writer.write_text(out / "timestamps.dat", (case["folder"] / "timestamps.dat").read_text())
    assert (out / "timestamps.dat").read_text() == (case["folder"] / "timestamps.dat").read_text()
    return list(vals), missing


def differences(a, b):
    """What differs between two HDF5 files: the objects each holds, and for
    every one both hold its kind, a dataset's shape, type and values, and every
    attribute's type, shape and value. Empty when they are equal."""
    import h5py

    def objects(f):
        out = {"/": f}
        f.visititems(lambda n, o: out.__setitem__(n, o))
        return out

    out = []
    with h5py.File(str(a), "r") as f, h5py.File(str(b), "r") as g:
        fo, go = objects(f), objects(g)
        if set(fo) != set(go):
            out.append("objects %s against %s" % (sorted(fo), sorted(go)))
        for n in sorted(set(fo) & set(go)):
            x, y = fo[n], go[n]
            if isinstance(x, h5py.Dataset) != isinstance(y, h5py.Dataset):
                out.append("%s: a dataset in one file only" % n)
                continue
            if isinstance(x, h5py.Dataset):
                if x.shape != y.shape or x.dtype != y.dtype:
                    out.append("%s: %s %s against %s %s" % (n, x.shape, x.dtype, y.shape, y.dtype))
                elif not np.array_equal(x[()], y[()]):
                    out.append("%s: the values differ" % n)
            if set(x.attrs) != set(y.attrs):
                out.append("%s: attributes %s against %s" % (n, sorted(x.attrs), sorted(y.attrs)))
                continue
            for k in x.attrs:
                p, q = x.attrs.get_id(k), y.attrs.get_id(k)
                if p.dtype != q.dtype or p.shape != q.shape \
                        or not np.array_equal(x.attrs[k], y.attrs[k]):
                    out.append("%s@%s: %r %s against %r %s"
                               % (n, k, x.attrs[k], p.dtype, y.attrs[k], q.dtype))
    return out


# ------------------------------------------------------------------ one rank


@pytest.mark.parametrize("dim", [2, 3], ids=["2d", "3d"])
def test_one_rank_writes_what_the_serial_writers_write(tmp_path, dim):
    """On one rank the block is the whole and the route is the gather to
    itself: the files must be the serial writers' own, the appended phase
    field and the copied pressure included, and a missing group is reported
    rather than refused."""
    pytest.importorskip("h5py")
    from mpi4py import MPI
    cfg = make_case(tmp_path / "in", dim)
    names = serial_output(cfg, tmp_path / "ref")
    writer = W.BlockWriter(MPI.COMM_SELF)
    assert writer.route == "gathered to rank 0"
    got, missing = block_output(cfg, tmp_path / "out", writer,
                                lambda n: (0, n), lambda n: (0, n))
    assert got == names and missing == [False] * len(names)
    for name in names:
        assert differences(tmp_path / "out" / name, tmp_path / "ref" / name) == [], name


def test_rows_of_the_wrong_width_are_refused_on_one_rank(tmp_path):
    """A block whose rows do not match the element is a caller's mistake that
    would write a dof table the loaders misread; on one rank it is an
    exception naming the widths."""
    pytest.importorskip("h5py")
    from mpi4py import MPI
    writer = W.BlockWriter(MPI.COMM_SELF)
    with pytest.raises(ValueError, match="3 components"):
        writer.write_vector(tmp_path / "u.h5", "u", np.zeros((4, 2)), 2, 3, 3)


# ------------------------------------------------------------ several ranks


def run_job(args, ranks, timeout=300):
    """This module as a job of that many ranks. One OpenMP thread a rank, set
    in the job's environment alone."""
    launcher = mpi_launcher()
    if launcher is None:
        if os.environ.get("PARTRAC_REQUIRE_MPI"):
            pytest.fail("no launcher here starts an MPI job for mpi4py")
        pytest.skip("no launcher here starts an MPI job for mpi4py")
    env = dict(os.environ, OMP_NUM_THREADS="1",
               PYTHONPATH=os.pathsep.join([os.path.dirname(HERE), HERE,
                                           os.environ.get("PYTHONPATH", "")]))
    env.pop("PARTRAC_HDF5_GATHER", None)
    if args[0] == "gather":
        env["PARTRAC_HDF5_GATHER"] = "1"
    cmd = [launcher, "-n", str(ranks), sys.executable, os.path.abspath(__file__), "job"]
    return subprocess.run(cmd + [str(a) for a in args[1:]], capture_output=True, text=True,
                          timeout=timeout, env=env)


def expect_route(r, route):
    """The route the job took is the one asked for; the parallel one is
    skipped only where the job's h5py has no MPI driver, never assumed."""
    text = r.stdout + r.stderr
    if route == "parallel":
        if "h5py mpi: False" in text:
            pytest.skip("this h5py has no MPI driver, so there is no parallel route")
        assert "h5py mpi: True" in text, text[-3000:]
        assert "route: parallel HDF5" in text, text[-3000:]
    else:
        assert "route: gathered to rank 0" in text, text[-3000:]


@pytest.mark.slow
@pytest.mark.parametrize("route", ["parallel", "gather"])
@pytest.mark.parametrize("ranks", [2, 3, 4])
@pytest.mark.parametrize("dim", [2, 3], ids=["2d", "3d"])
def test_every_rank_writes_what_the_serial_writers_write(tmp_path, dim, ranks, route):
    """Each rank writes its own uneven block of the cells and of the nodes --
    the two split differently, rank 1 holding no cells and at three ranks and
    more one rank no nodes -- in pieces far smaller than a block, so every
    loop runs. The files must equal the serial writers' dataset by dataset and
    attribute by attribute on either route, and a group the input does not
    hold is reported missing by every rank."""
    pytest.importorskip("h5py")
    cfg = make_case(tmp_path / "in", dim)
    serial_output(cfg, tmp_path / "ref")
    r = run_job([route, cfg, tmp_path / "out", tmp_path / "ref"], ranks)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    expect_route(r, route)
    assert "equal: up_0.h5 up_1.h5" in r.stdout, r.stdout[-3000:]
    assert "missing on every rank" in r.stdout


@pytest.mark.slow
@pytest.mark.parametrize("route,rank,kind", [
    ("parallel", 1, "path"), ("parallel", 0, "path"), ("gather", 0, "path"),
    ("parallel", 1, "rows"), ("gather", 1, "rows")],
    ids=["parallel-rank1-path", "parallel-rank0-path", "gather-rank0-path",
         "parallel-rank1-rows", "gather-rank1-rows"])
def test_a_write_that_fails_on_one_rank_takes_the_job_down(tmp_path, route, rank, kind):
    """A write that fails on one rank leaves the others inside a collective --
    the open, a gather, the close -- where they would wait forever and the job
    would look like one still running. It stops instead, within the timeout,
    naming the rank and what it could not write: a file where the rank's
    output directory should be (on the gathered route only the first rank
    opens a file), or rows that do not fit the element."""
    pytest.importorskip("h5py")
    cfg = make_case(tmp_path / "in", 2)
    (tmp_path / "blocker").write_text("a file where a directory is expected\n")
    r = run_job([route, cfg, tmp_path / "out", "-", "fail", rank, kind, tmp_path / "blocker"],
                2, timeout=120)
    text = r.stdout + r.stderr
    assert r.returncode != 0, text[-3000:]
    expect_route(r, route)
    assert "rank %d could not write" % rank in text, text[-3000:]
    if kind == "path":
        assert str(tmp_path / "blocker") in text
    else:
        assert "ValueError" in text


@pytest.mark.slow
@pytest.mark.parametrize("dim,mode", [(2, "triangle"), (3, "tet")], ids=["triangle", "tet"])
def test_three_ranks_write_a_case_that_reads_back(tmp_path, dim, mode):
    """A case written by three ranks on the route this machine takes reads
    back through the tool's reader as the fields that went in, and through
    the loader as the P2 field the Python evaluation of the same file gives:
    a dof table assembled from blocks is exactly where a silent permutation
    would hide."""
    pytest.importorskip("h5py")
    cfg = make_case(tmp_path / "in", dim)
    out = tmp_path / "out"
    r = run_job(["default", cfg, out, "-"], 3)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    for name in ("mesh.h5", "dolfin_params.dat"):
        shutil.copy(tmp_path / "in" / name, out / name)
    case, t, want = _read(cfg)
    got = _read(out / "dolfin_params.dat")[2]
    for name in want:
        for f in ("u", "phi"):
            assert np.array_equal(got[name][f], want[name][f]), (name, f)
    if not os.path.exists(INTERPOL):
        pytest.skip("interpol is not built")
    g = interpol_probe(out / "dolfin_params.dat", mode, 60)
    pts = np.stack([g[a] for a in "xyz"[:dim]], axis=1)
    val = np.stack([g["u" + a] for a in "xyz"[:dim]], axis=1)
    cell = np.full(len(pts), -1)
    for k, c in enumerate(t.cells):
        lam = np.linalg.solve((t.X[c[1:]] - t.X[c[0]]).T, (pts - t.X[c[0]]).T)
        inside = (lam >= -1e-12).all(axis=0) & (lam.sum(axis=0) <= 1 + 1e-12)
        cell[inside] = k
    assert (cell >= 0).all(), "a probe point is in no cell"
    U = got["up_0.h5"]["u"]
    assert np.abs(val - D.p2_eval(t, U, pts, cell)).max() < 1e-12 * np.abs(U).max()


# -------------------------------------------------------------- the MPI job


def _uneven(n, weights, rank):
    """Rank's block of n rows when the ranks take them in these proportions."""
    b = np.floor(n * np.concatenate([[0], np.cumsum(weights)]) / sum(weights)).astype(int)
    b[-1] = n
    return int(b[rank]), int(b[rank + 1])


def _job(argv):
    """route cfg out ref [fail rank kind blocker]: the case written from this
    rank's blocks and, on the first rank, compared with the serial files in
    ref; or, with fail, one rank's write made to fail."""
    import h5py
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    cfg, out, ref = argv[0], Path(argv[1]), argv[2]
    W.PIECE_BYTES = 40          # every block written and gathered in many pieces
    writer = W.BlockWriter(comm)
    if comm.rank == 0:
        print("h5py mpi: %s" % h5py.get_config().mpi)
        print("route: %s" % writer.route, flush=True)
    size, rank = comm.size, comm.rank
    wc = [(r + 1) ** 2 for r in range(size)]
    wc[1] = 0
    wn = [(size - r) ** 2 for r in range(size)]
    if size >= 3:
        wn[-1] = 0
    cells = lambda n: _uneven(n, wc, rank)
    nodes = lambda n: _uneven(n, wn, rank)
    if len(argv) > 3:
        bad, kind, blocker = int(argv[4]), argv[5], Path(argv[6])
        case, t, vals = _read(cfg)
        (c0, c1), (n0, n1) = cells(t.ncells), nodes(t.nnodes)
        U = vals["up_0.h5"]["u"][n0:n1]
        if rank == bad and kind == "path":
            out = blocker
        else:
            out.mkdir(parents=True, exist_ok=True)
        if rank == bad and kind == "rows":
            U = U[:, :1]
        comm.Barrier()
        writer.write_checkpoint(out / "up_0.h5", "u", D._cell_dofs(t, 2, t.dim)[c0:c1],
                                case["cell_indices"].astype(np.uint64)[c0:c1], U, 2,
                                t.dim, t.dim)
        comm.Barrier()
        print("rank %d wrote without failing" % rank)
        return 0
    names, missing = block_output(cfg, out, writer, cells, nodes)
    missing = comm.gather(missing)
    if rank == 0 and ref != "-":
        bad = {n: differences(out / n, Path(ref) / n) for n in names}
        bad = {n: d for n, d in bad.items() if d}
        if bad:
            print("different: %s" % bad)
            return 1
        print("equal: %s" % " ".join(names))
    if rank == 0:
        if any(any(m) for m in missing):
            print("copy_group found a group that is not there: %s" % missing)
            return 1
        print("missing on every rank")
    return 0


if __name__ == "__main__" and sys.argv[1:2] == ["job"]:
    sys.exit(_job(sys.argv[2:]))
