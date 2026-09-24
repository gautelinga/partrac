"""python/divfree/divfree_dist.py: the cleaner's mesh tables shared out over
the ranks.

Every table a rank builds for its own cells is the serial `Topo`'s restricted
to them, and gathered they are `Topo`'s array for array: the numbering is the
serial one whatever the rank count and the partition, which is what makes the
output the same at every rank count. The building blocks the later steps are
built from -- the keyed exchange, the sample sort with its way back, the
exclusive scan, the reductions and the exact median -- are checked on their
own.

The one-rank path runs in process. The tables at one to four ranks run as MPI
jobs, one a mesh, each running every check and reporting by its exit code;
run as a script, this file is that job. The jobs skip where no launcher starts
a job of several ranks for mpi4py, and fail there with PARTRAC_REQUIRE_MPI set.
"""

import os
import sys
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import pytest                      # noqa: E402

from divfree_cases import (D, channel3d, mpi_run, relabel_by_global_id,  # noqa: E402
                           topo_cases)

import divfree_dist as DD          # noqa: E402
from mpi4py import MPI             # noqa: E402

TABLES = ["cells", "cell_edge", "edges", "master", "boundary_node", "boundary_edge",
          "held_node", "held_edge", "facet_exterior", "facet_periodic", "facet_boundary",
          "cell_boundary", "facet_verts", "vol", "facet_n", "J", "Jinv", "detJ", "node_x", "X"]


# ------------------------------------------------------------------ the checks


def rest_rule(t):
    """An at_rest set of the serial nodes that holds periodic images whose
    masters are not at rest: every node on the high face of the first periodic
    direction, and a scatter of others."""
    rest = (np.arange(t.nnodes) * 2654435761) % 7 < 3
    per = np.nonzero(t.periodic)[0]
    if len(per):
        a = per[0]
        rest |= t.node_x[:, a] > t.x_max[a] - 1e-9
        rest[t.master[rest]] = False
        rest |= t.node_x[:, a] > t.x_max[a] - 1e-9
    return rest


def serial_first(t):
    """(the smallest cell holding each node, the smallest holding any image of
    each edge's master) of a serial topology, as `free_unknowns` finds them."""
    nodes = np.concatenate([t.cells, t.nverts + t.cell_edge], axis=1)
    first = np.full(t.nnodes, t.ncells, np.int64)
    first[nodes[::-1].ravel()] = np.repeat(np.arange(t.ncells - 1, -1, -1), nodes.shape[1])
    m = t.master[t.nverts:] - t.nverts
    fm = np.full(t.nedges, t.ncells, np.int64)
    fm[m[t.cell_edge[::-1]].ravel()] = np.repeat(np.arange(t.ncells - 1, -1, -1),
                                                t.cell_edge.shape[1])
    return first, fm[m]


def check_tables(mesh, per, partition, comm, ref, ci, label=None):
    """Every check of one partition of one mesh; returns their count."""
    n = 0
    t = DD.DistTopo.read(mesh, per, comm, partition=partition)
    rank, size = comm.rank, comm.size
    block = np.arange(*D._split(ref.ncells, size, rank))
    moved = comm.allreduce(not np.array_equal(t.cell_gid, block), op=MPI.LOR)
    assert not (moved and partition == "blocks")
    if rank == 0 and moved:
        print("%s moved cells" % (label or partition))
    if label:
        assert comm.allreduce(t.ncells == 0, op=MPI.LOR), "no rank is without a cell"
    # the globals, the same on every rank and the serial ones
    for name in ("dim", "nv"):
        assert getattr(t, name) == getattr(ref, name), name
    for name in ("ncells", "nverts", "nedges", "nnodes"):
        assert getattr(t, name + "_global") == getattr(ref, name), name
        assert len(set(comm.allgather(getattr(t, name + "_global")))) == 1, name
    assert np.array_equal(t.periodic, ref.periodic)
    assert np.array_equal(t.x_min, ref.x_min) and np.array_equal(t.x_max, ref.x_max)
    assert t.mesh_h == D.mesh_h(ref), (t.mesh_h, D.mesh_h(ref))
    assert abs(t.vol_total - ref.vol.sum()) <= 1e-14 * ref.vol.sum()
    assert tuple(t.node_block) == D._split(ref.nnodes, size, rank)
    n += 10
    # the rank's own tables are the serial ones restricted to its cells
    c = t.cell_gid
    assert np.all(np.diff(c) > 0) and np.all(np.diff(t.vert_gid) > 0)
    assert np.all(np.diff(t.edge_gid) > 0)
    assert np.array_equal(t.vert_gid[t.cells], ref.cells[c])
    assert np.array_equal(t.edge_gid[t.cell_edge], ref.cell_edge[c])
    assert np.array_equal(t.vert_gid[t.edges], ref.edges[t.edge_gid])
    assert np.array_equal(t.vert_gid[t.facet_verts], ref.facet_verts[c])
    for name in ("vol", "detJ", "J", "Jinv", "facet_n", "facet_exterior", "facet_periodic",
                 "facet_boundary", "cell_boundary"):
        assert np.array_equal(getattr(t, name), getattr(ref, name)[c]), name
    g = t.node_gid
    for name in ("node_x", "master", "boundary_node"):
        assert np.array_equal(getattr(t, name), getattr(ref, name)[g]), name
    assert np.array_equal(t.boundary_edge, ref.boundary_edge[t.edge_gid])
    assert np.array_equal(t.X, ref.X[t.vert_gid])
    if ci is not None:
        assert np.array_equal(t.cell_index, ci[c])
    n += 22
    # the cells are shared out once and whole
    allc = np.concatenate(comm.allgather(c))
    assert len(allc) == ref.ncells and np.array_equal(np.sort(allc), np.arange(ref.ncells))
    rank_of = np.zeros(ref.ncells, np.int64)
    for r, cr in enumerate(comm.allgather(c)):
        rank_of[cr] = r
    # ownership by the smallest cell holding the node, as free_unknowns numbers them
    first, first_master = serial_first(ref)
    assert np.array_equal(t.first_cell, first[g])
    assert np.array_equal(t.node_owner, rank_of[first[g]])
    # the local view: what the rank holds under Topo's names
    assert (t.ncells, t.nverts, t.nedges) == (len(c), len(t.vert_gid), len(t.edge_gid))
    assert t.nnodes == t.nverts + t.nedges == len(g)
    lv = t.nverts
    got = DD.reduce_keyed(t.master[lv:], t.first_cell[lv:], "min", comm, t.nnodes_global)
    assert np.array_equal(got, first_master[t.edge_gid])
    n += 5
    # reduce_nodes: every op against the copies every rank holds
    every = comm.allgather(g)
    copies = np.zeros(ref.nnodes, np.int64)
    top = np.full(ref.nnodes, -1)
    low = np.full(ref.nnodes, size)
    for r, gr in enumerate(every):
        copies[gr] += 1
        top[gr] = np.maximum(top[gr], r)
        low[gr] = np.minimum(low[gr], r)
    assert np.array_equal(t.reduce_nodes(np.ones(len(g), np.int64), "sum"), copies[g])
    two = np.stack([np.ones(len(g)), g.astype(float)], axis=1)
    assert np.array_equal(t.reduce_nodes(two, "sum"),
                          np.stack([copies[g], copies[g] * g], axis=1).astype(float))
    assert np.array_equal(t.reduce_nodes(np.full(len(g), rank), "max"), top[g])
    assert np.array_equal(t.reduce_nodes(np.full(len(g), rank), "min"), low[g])
    assert np.array_equal(t.reduce_nodes(np.full(len(g), rank == size - 1), "or"),
                          top[g] == size - 1)
    n += 5
    # the held set, some periodic images held and their masters not at rest
    rest = rest_rule(ref)
    ref.set_held(rest)
    if ref.dim == 3 and ref.boundary_node.any() and ref.periodic.any():
        assert (ref.held_node & ~rest).any(), "no master is held through an image alone"
    t.set_held(rest[g])
    assert np.array_equal(t.held_node, ref.held_node[g])
    assert np.array_equal(t.held_edge, ref.held_edge[t.edge_gid])
    n += 2
    # the unknowns: the serial ones, numbered in the serial order under blocks
    # and shared out evenly there; the first cell's rank's under a partitioner
    u, su = D.free_unknowns(t), D.free_unknowns(ref)
    assert u.total == su.total
    reached = t.master[lv:][~t.held_edge] - t.nverts_global
    assert np.array_equal(np.unique(reached), np.unique(u.master[u.col[~t.held_edge]]))
    own = (u.gid >= u.start) & (u.gid < u.start + u.count)
    assert own.sum() == u.count and np.all(own[u.master < 0])
    if partition == "blocks":
        assert (u.start, u.start + u.count) == D._split(u.total, size, rank)
        hit = u.master >= 0
        assert np.array_equal(su.master[u.gid[hit]], u.master[hit])
    else:
        m = ref.master[ref.nverts:] - ref.nverts
        mine = np.unique(m[~ref.held_edge][rank_of[first_master[~ref.held_edge]] == rank])
        assert np.array_equal(np.sort(u.master[own]), mine)
    n += 4
    # gathered, array for array
    out = t.gather(0, mesh)
    if rank == 0:
        for name in TABLES:
            a, b = getattr(out, name), getattr(ref, name)
            assert a.shape == b.shape and np.array_equal(a, b), name
        if ci is not None:
            assert np.array_equal(out.cell_index, ci)
        n += len(TABLES)
    else:
        assert out is None
    ref.set_held(None)
    return n


def last_rank_empty(gcells, gid, ncells, nverts, comm, kind):
    """A partition that leaves the last rank without a cell: every rank's cells
    go to itself, the last rank's to the one before it."""
    return np.full(len(gcells), min(comm.rank, comm.size - 2), np.int64)


def job(mesh, per):
    """The MPI job of one mesh: blocks and every partitioner this PETSc has,
    the refusal of each it lacks, a partition that leaves a rank without a
    cell, and every check of each."""
    comm = MPI.COMM_WORLD
    X, cells, ci = D.read_mesh_h5(mesh)
    X = X[:, :cells.shape[1] - 1]
    ref = D.Topo(X, cells, per)
    n = 0
    for partition in ("blocks",) + DD.PARTITIONERS:
        if partition == "blocks" or D.PETSc.Sys.hasExternalPackage(partition):
            n += check_tables(mesh, per, partition, comm, ref, ci)
            continue
        try:
            DD.DistTopo.read(mesh, per, comm, partition=partition)
            raise AssertionError("%s was not refused" % partition)
        except ValueError as e:
            want = "partition=%s needs a PETSc built with %s" % (partition, partition)
            assert want in str(e), str(e)
        n += 1
    if comm.size > 1:
        real = DD._partition
        DD._partition = last_rank_empty
        try:
            n += check_tables(mesh, per, "ptscotch", comm, ref, ci, "a rank without a cell")
        finally:
            DD._partition = real
    return n


def blocks_job():
    """The building blocks across the ranks against their serial answers."""
    comm = MPI.COMM_WORLD
    rank, size = comm.rank, comm.size
    rng = np.random.default_rng(10 + rank)
    n = 0
    # the sample sort: every rank's keys, many repeated, the second rank with none
    k = rng.integers(0, 50 if rank % 2 else 5000, size=(0 if rank == 1 else 3000))
    if rank != 1:
        k = np.concatenate([k, np.full(40, 7)])
    assert len(k) == 0 or rank != 1
    allk = np.concatenate(comm.allgather(k))
    ids, total, cnt, s = DD.number_keys(k, comm)
    u, inv, uc = np.unique(allk, return_inverse=True, return_counts=True)
    off = sum(len(a) for a in comm.allgather(k)[:rank])
    assert total == len(u)
    assert np.array_equal(ids, inv[off:off + len(k)])
    assert np.array_equal(cnt, uc[inv[off:off + len(k)]])
    shares = comm.allgather(s.keys[:, 0])
    assert np.array_equal(np.concatenate(shares), np.sort(allk))
    assert np.array_equal(s.back(s.keys[:, 0]), k)
    n += 5
    # two-column keys order lexicographically, routed by the first
    m2 = 0 if rank == 1 else 500
    k2 = np.stack([rng.integers(0, 30, m2), rng.integers(0, 4, m2)], axis=1)
    all2 = np.concatenate(comm.allgather(k2))
    ids2, total2, _, _ = DD.number_keys(k2, comm)
    u2, inv2 = np.unique(all2, axis=0, return_inverse=True)
    off2 = sum(len(a) for a in comm.allgather(k2)[:rank])
    assert total2 == len(u2) and np.array_equal(ids2, inv2.ravel()[off2:off2 + m2])
    n += 1
    # the exclusive scan
    assert DD.exscan(rank + 1, comm) == (rank * (rank + 1) // 2, size * (size + 1) // 2)
    n += 1
    # the exact median: odd and even counts, many equal values, negative ones
    for m in (1001 + rank, 1000 + 2 * rank):
        for kind in ("wide", "equal", "signed"):
            v = rng.random(m)
            if kind == "equal":
                v = np.round(v * 3) / 3
            if kind == "signed":
                v = rng.normal(size=m) * 1e-3
            allv = np.concatenate(comm.allgather(v))
            assert DD.exact_median(v, comm) == float(np.median(allv)), (m, kind)
            n += 1
    # a reduction keyed by any id, each value returned to every copy
    keys = rng.integers(0, 97, 300)
    vals = rng.integers(0, 1000, 300)
    allk, allv = np.concatenate(comm.allgather(keys)), np.concatenate(comm.allgather(vals))
    for op, f in (("max", np.maximum), ("min", np.minimum), ("sum", np.add)):
        want = {int(q): int(f.reduce(allv[allk == q])) for q in np.unique(allk)}
        got = DD.reduce_keyed(keys, vals, op, comm, 97)
        assert [want[int(q)] for q in keys] == list(got), op
        n += 1
    # an exchange there and back keeps the rows and their order
    dest = rng.integers(0, size, 200)
    rows = np.stack([np.arange(200) + 1000 * rank, dest], axis=1)
    ex = DD.Exchange(comm, dest)
    got = ex.forward(rows)
    assert np.all(got[:, 1] == rank)
    assert np.array_equal(ex.sources(), got[:, 0] // 1000)
    assert np.array_equal(ex.back(got), rows)
    n += 3
    return n


def main(argv):
    comm = MPI.COMM_WORLD
    try:
        if argv[0] == "blocks":
            n = blocks_job()
        else:
            n = job(argv[0], [c == "T" for c in argv[1]])
    except BaseException:
        sys.stderr.write("rank %d:\n%s" % (comm.rank, traceback.format_exc()))
        sys.stderr.flush()
        comm.Abort(1)
    n = comm.reduce(n, root=0)
    if comm.rank == 0:
        print("PASSED %d checks on %d ranks" % (n, comm.size))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


# ------------------------------------------------------------------ the cases


def mesh_cases():
    """(X, cells, periodic, labelled) by name: every mesh of `topo_cases()`,
    and the 2D and 3D channels labelled as a run on several ranks writes them."""
    out = {k: v + (False,) for k, v in topo_cases().items()}
    for base in ("2d", "3d"):
        out[base + "_labelled"] = topo_cases()[base] + (True,)
    return out


CASES = mesh_cases()


def write_mesh(folder, name):
    """The case's mesh as the tool reads it; (its path, the periodic flags)."""
    X, cells, per, labelled = CASES[name]
    folder.mkdir(parents=True, exist_ok=True)
    D.write_mesh_h5(folder / "mesh.h5", X, cells)
    if labelled:
        relabel_by_global_id(folder)
    return folder / "mesh.h5", per


def run_job(ranks, *args):
    """This file as an MPI job; its output when it passed."""
    r = mpi_run(ranks, [os.path.abspath(__file__)] + list(args), timeout=300)
    assert r.returncode == 0 and "PASSED" in r.stdout, r.stdout[-3000:] + r.stderr[-5000:]
    return r.stdout


# -------------------------------------------------------------- in process


@pytest.mark.parametrize("name", sorted(CASES))
def test_one_rank_builds_the_serial_tables(tmp_path, name):
    """On one rank nothing is shared out, and the tables are `Topo`'s as they
    are: the same numbering, the same bits, the same held set."""
    pytest.importorskip("h5py")
    mesh, per = write_mesh(tmp_path, name)
    X, cells, ci = D.read_mesh_h5(mesh)
    ref = D.Topo(X[:, :cells.shape[1] - 1], cells, per)
    assert check_tables(mesh, per, "ptscotch", MPI.COMM_SELF, ref, ci) > 40


def test_the_sample_sort_numbers_keys_by_their_rank_among_the_distinct_ones():
    """An edge's id is its rank among the sorted vertex pairs of the whole mesh,
    so the numbering is the serial one: repeated keys get one id, the ids count
    up from zero in key order, the multiplicities are the keys' counts, and the
    way back returns every answer to the key it belongs to."""
    rng = np.random.default_rng(4)
    k = rng.integers(0, 200, 5000)
    ids, total, cnt, s = DD.number_keys(k, MPI.COMM_SELF, payload=(np.arange(5000),))
    u, inv, uc = np.unique(k, return_inverse=True, return_counts=True)
    assert total == len(u)
    assert np.array_equal(ids, inv) and np.array_equal(cnt, uc[inv])
    assert np.array_equal(s.keys[:, 0], np.sort(k))
    assert np.array_equal(k[s.payload[0]], s.keys[:, 0])
    assert np.array_equal(s.back(s.payload[0]), np.arange(5000))


def test_an_exchange_takes_only_rows_of_its_own_addressing():
    """An exchange is built for one set of rows and moves any array of them;
    handed another number of rows it would send the first ones by the counts,
    or pick a subset by its order, and the rows would land beside cells they
    do not belong to. It refuses them, on one rank too, where it would
    otherwise hand them back untouched."""
    for dest in ([0, 0, 0], [0, 0, 0][::-1]):
        ex = DD.Exchange(MPI.COMM_SELF, dest)
        assert np.array_equal(ex.forward(np.arange(3)), np.arange(3))
        for n in (2, 4):
            with pytest.raises(AssertionError, match="rows for another exchange"):
                ex.forward(np.arange(n))


def test_the_exclusive_scan_of_one_rank_is_zero():
    """The first rank's ids start at zero, and the total is its own count."""
    assert DD.exscan(17, MPI.COMM_SELF) == (0, 17)


@pytest.mark.parametrize("m", [1, 2, 999, 1000])
@pytest.mark.parametrize("kind", ["wide", "equal", "signed"])
def test_the_exact_median_is_numpys(m, kind):
    """gamma = g h^2 enters every output value, and h is the median edge
    length, so the distributed median must be np.median bit for bit: the
    middle value of an odd count, the mean of the two middle ones of an even
    count, and among many equal values or across zero too."""
    rng = np.random.default_rng(m)
    v = rng.random(m)
    if kind == "equal":
        v = np.round(v * 3) / 3
    if kind == "signed":
        v = rng.normal(size=m)
    assert DD.exact_median(v, MPI.COMM_SELF) == float(np.median(v))


def test_a_partitioner_this_petsc_lacks_is_refused_by_name():
    """A partition option names a PETSc external package; one this PETSc was
    not built with is refused, naming the option and the ones it has, before
    anything is read."""
    for p in DD.PARTITIONERS:
        if D.PETSc.Sys.hasExternalPackage(p):
            continue
        with pytest.raises(ValueError, match="partition=%s needs a PETSc built with %s" % (p, p)):
            DD.DistTopo.read("no such file", [False] * 3, MPI.COMM_SELF, partition=p)
    with pytest.raises(ValueError, match="partition 'metis' is none of blocks"):
        DD.check_partition("metis")


# ----------------------------------------------------------- on several ranks


@pytest.mark.slow
@pytest.mark.parametrize("ranks", [2, 3, 4])
def test_the_building_blocks_across_ranks_give_the_serial_answers(ranks):
    """The sample sort with its way back, the scan, the keyed reductions and
    the exact median, each on keys and values spread unevenly over the ranks --
    the second rank with no key at all, many repeated -- against the same
    computation on all of them together."""
    assert "PASSED" in run_job(ranks, "blocks")


@pytest.mark.slow
@pytest.mark.parametrize("ranks", [1, 2, 3, 4])
@pytest.mark.parametrize("name", sorted(CASES))
def test_the_distributed_tables_are_the_serial_ones(tmp_path, name, ranks):
    """Read in blocks, and partitioned by PT-Scotch and moved: each rank's
    tables are `Topo`'s restricted to its cells, the globals are the serial
    ones on every rank (the median edge length bit for bit), every node is
    owned by the rank holding the smallest cell that holds it, every reduction
    reaches every copy, the held set holds a master whose image is held, and
    gathered, every table is `Topo`'s array for array. A partitioner this
    PETSc lacks is refused. A rank left without a cell builds its empty share
    and the others' tables are still the serial ones: a rank that raised alone
    there would leave the rest waiting in the next exchange."""
    pytest.importorskip("h5py")
    mesh, per = write_mesh(tmp_path, name)
    out = run_job(ranks, mesh, "".join("T" if p else "F" for p in per))
    if ranks > 1 and name.endswith("shuffled"):
        assert "ptscotch moved cells" in out


def test_the_triply_periodic_fixture_has_corners_that_need_the_chain():
    """The tables are only checked on the chain of links if a fixture has one:
    the corners of the triply periodic channel are linked to a face, that face
    to an edge and that to the origin, and all eight share the one master."""
    X, cells = channel3d()
    t = D.Topo(X, cells, [True, True, True])
    corners = np.nonzero(np.all((X == 0) | (X == 1), axis=1))[0]
    assert len(corners) == 8
    origin = np.nonzero(np.all(X == 0, axis=1))[0][0]
    assert np.all(t.master[corners] == origin)


def test_a_facet_key_does_not_wrap_past_two_million_vertices():
    """Three vertex ids packed into one int64 wrap once the vertex count passes
    about 2.1M, and two different facets then share a key: each is counted as
    held by two cells, so a boundary facet reads as interior. With 2^22
    vertices, (5, b, c) and (2^20 + 5, b, c) are such a pair; the key is two
    columns, and a mesh of two tets on those vertices has every facet exterior
    but the one they share."""
    n = 1 << 22
    F = np.array([[5, (1 << 20) + 6, (1 << 20) + 7],
                  [(1 << 20) + 5, (1 << 20) + 6, (1 << 20) + 7],
                  [5, (1 << 20) + 6, (1 << 20) + 7]])
    packed = (F[:, 0] * n + F[:, 1]) * n + F[:, 2]
    assert packed[0] == packed[1], "the ids do not collide in one int64"
    assert list(D._facet_multiplicity(F, n)) == [2, 1, 2]
    X = np.zeros((n, 3))
    X[[5, 9]] = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]
    X[[(1 << 20) + 5, (1 << 20) + 6, (1 << 20) + 7]] = np.eye(3)
    cells = np.array([[5, (1 << 20) + 6, (1 << 20) + 7, 9],
                      [(1 << 20) + 5, (1 << 20) + 6, (1 << 20) + 7, 9]])
    t = D.Topo(X, cells, [False] * 3)
    assert t.facet_exterior.sum() == 6
