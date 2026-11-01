"""python/divfree/divfree_dist.py: the cleaner's mesh tables shared out over
the ranks.

Every table a rank builds for its own cells is the serial `Topo`'s restricted
to them, and gathered they are `Topo`'s array for array: the numbering is the
serial one whatever the rank count and the partition, which is what makes the
output the same at every rank count. The building blocks the later steps are
built from -- the keyed exchange, the sample sort with its way back, the
exclusive scan, the reductions and the exact median -- are checked on their
own.

The one-rank path runs in process. The tables at two to four ranks run as MPI
jobs through divfree_cases.mpi_batch, one a rank count, the rank counts shared
out over the meshes, each running every check of each of its meshes as a
case. The one-rank job runs in process too, where COMM_WORLD is one rank. The
jobs skip where no launcher starts a job of several ranks for mpi4py, and fail
there with PARTRAC_REQUIRE_MPI set.
"""

import types

import numpy as np
import pytest

from divfree_cases import (D, channel2d, channel3d, launcher_or_skip, mpi_batch,
                           read_mesh_h5, relabel_by_global_id, shuffled, topo_cases)

import divfree_dist as DD          # noqa: E402
from mpi4py import MPI             # noqa: E402

TABLES = ["cells", "cell_edge", "edges", "master", "boundary_node", "boundary_edge",
          "held_node", "held_edge", "facet_boundary", "cell_boundary", "vol", "facet_n", "J",
          "Jinv", "node_x", "X"]


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


def local_gids(t):
    """(the serial ids of the rank's vertices, of its edges)."""
    return t.node_gid[:t.nverts], t.node_gid[t.nverts:] - t.nverts_global


def gather(t, mesh_path, root=0):
    """The rank's tables assembled on root in the serial order, with `Topo`'s
    names; None on the other ranks. A vertex in no cell is taken from
    mesh_path's coordinates, as unpaired and off the boundary."""
    comm = t.comm
    lv = t.nverts
    vg, eg = local_gids(t)
    cellpart = dict(
        gid=t.cell_gid, cells=vg[t.cells], cell_edge=eg[t.cell_edge],
        facet_boundary=t.facet_boundary, cell_boundary=t.cell_boundary,
        vol=t.vol, facet_n=t.facet_n, J=t.J, Jinv=t.Jinv)
    if t.cell_index is not None:
        cellpart["cell_index"] = t.cell_index
    own = t.node_owner == comm.rank
    nodepart = dict(gid=t.node_gid[own], master=t.master[own],
                    boundary_node=t.boundary_node[own], held_node=t.held_node[own],
                    at_rest=t.at_rest[own], node_x=t.node_x[own])
    eown = own[lv:]
    edgepart = dict(gid=eg[eown], edges=vg[t.edges[eown]])
    parts = comm.gather((cellpart, nodepart, edgepart), root=root)
    if comm.rank != root:
        return None
    out = types.SimpleNamespace()

    def place(i, n, names, dtypes):
        for name in names:
            g = np.concatenate([p[i]["gid"] for p in parts])
            v = np.concatenate([p[i][name] for p in parts])
            full = np.zeros((n,) + v.shape[1:], dtypes.get(name, v.dtype))
            full[g] = v
            setattr(out, name, full)

    i32 = np.int32
    place(0, t.ncells_global, [k for k in cellpart if k != "gid"],
          dict(cells=i32, cell_edge=i32))
    place(2, t.nedges_global, ["edges"], dict(edges=i32))
    n = t.nnodes_global
    g = np.concatenate([p[1]["gid"] for p in parts])
    out.master = np.arange(n, dtype=i32)
    out.master[g] = np.concatenate([p[1]["master"] for p in parts])
    for name in ("boundary_node", "held_node", "at_rest"):
        a = np.zeros(n, bool)
        a[g] = np.concatenate([p[1][name] for p in parts])
        setattr(out, name, a)
    out.node_x = np.zeros((n, t.dim))
    out.node_x[g] = np.vstack([p[1]["node_x"] for p in parts])
    missing = np.setdiff1d(np.arange(t.nverts_global), g[g < t.nverts_global])
    if len(missing):
        X, _, _ = read_mesh_h5(mesh_path)
        out.node_x[missing] = X[missing, :t.dim]
    out.X = out.node_x[:t.nverts_global]
    out.boundary_edge = out.boundary_node[t.nverts_global:]
    out.held_edge = out.held_node[t.nverts_global:]
    return out


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
    assert D.mesh_h(t) == D.mesh_h(ref), (D.mesh_h(t), D.mesh_h(ref))
    assert abs(t.vol_total - ref.vol.sum()) <= 1e-14 * ref.vol.sum()
    n += 9
    # the rank's own tables are the serial ones restricted to its cells
    c = t.cell_gid
    vg, eg = local_gids(t)
    assert np.all(np.diff(c) > 0) and np.all(np.diff(t.node_gid) > 0)
    assert np.all(vg < t.nverts_global) and np.all(eg >= 0)
    assert np.array_equal(vg[t.cells], ref.cells[c])
    assert np.array_equal(eg[t.cell_edge], ref.cell_edge[c])
    assert np.array_equal(vg[t.edges], ref.edges[eg])
    for name in ("vol", "J", "Jinv", "facet_n", "facet_boundary", "cell_boundary"):
        assert np.array_equal(getattr(t, name), getattr(ref, name)[c]), name
    g = t.node_gid
    for name in ("node_x", "master", "boundary_node"):
        assert np.array_equal(getattr(t, name), getattr(ref, name)[g]), name
    assert np.array_equal(t.boundary_edge, ref.boundary_edge[eg])
    assert np.array_equal(t.X, ref.X[vg])
    if ci is not None:
        assert np.array_equal(t.cell_index, ci[c])
    n += 18
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
    assert (t.ncells, t.nverts, t.nedges) == (len(c), len(vg), len(eg))
    assert t.nnodes == t.nverts + t.nedges == len(g)
    lv = t.nverts
    got = DD.reduce_keyed(t.master[lv:], t.first_cell[lv:], "min", comm, t.nnodes_global)
    assert np.array_equal(got, first_master[eg])
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
    assert np.array_equal(t.held_edge, ref.held_edge[eg])
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
    out = gather(t, mesh)
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


def last_rank_empty(nb, ncells, comm, kind):
    """A partition that leaves the last rank without a cell: every rank's cells
    go to itself, the last rank's to the one before it."""
    return np.full(len(nb), min(comm.rank, comm.size - 2), np.int64)


def box_facts(ref, per):
    """What the box's geometry and its literal flags say, not `Topo`: a node's
    master is it moved by whole periods along periodic axes and not at all along
    the others, every node on a periodic axis's high face has its master on the
    low face, and the periodic facets and the walls are the facets on the faces
    of the periodic axes and of the others. Returns the count."""
    per = np.asarray(per, bool)[:ref.dim]
    period = ref.x_max - ref.x_min
    x, m = ref.node_x, ref.master
    d = x - x[m]
    assert np.allclose(d[:, ~per], 0.0)
    k = d[:, per] / period[per]
    assert np.allclose(k, np.round(k)) and set(np.round(k).ravel()) <= {0.0, 1.0}
    assert np.array_equal(m[m], m)
    for a in np.nonzero(per)[0]:
        hi = np.isclose(x[:, a], ref.x_max[a])
        assert hi.any() and np.allclose(x[m[hi], a], ref.x_min[a]), a
    fx = ref.X[ref.facet_verts]
    on = [np.zeros(ref.facet_exterior.shape, bool) for _ in range(2)]
    for a in range(ref.dim):
        for v in (ref.x_min[a], ref.x_max[a]):
            on[int(per[a])] |= np.all(np.isclose(fx[..., a], v), axis=-1)
    assert np.array_equal(ref.facet_periodic, on[1])
    assert np.array_equal(ref.facet_exterior & ~ref.facet_periodic, on[0] & ~on[1])
    return 6


def job(mesh, per):
    """The MPI job of one mesh: blocks and every partitioner this PETSc has,
    the refusal of each it lacks, a partition that leaves a rank without a
    cell, and every check of each."""
    comm = MPI.COMM_WORLD
    X, cells, ci = read_mesh_h5(mesh)
    X = X[:, :cells.shape[1] - 1]
    ref = D.Topo(X, cells, per)
    n = box_facts(ref, per)
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
    ids, total = DD.number_keys(k, comm)
    u, inv = np.unique(allk, return_inverse=True)
    off = sum(len(a) for a in comm.allgather(k)[:rank])
    assert total == len(u)
    assert np.array_equal(ids, inv[off:off + len(k)])
    s = DD.SampleSort(k, comm, payload=(off + np.arange(len(k)),))
    shares = comm.allgather(s.keys[:, 0])
    assert np.array_equal(np.concatenate(shares), np.sort(allk))
    # equal keys in the order of their ranks, and of their positions within
    pos = np.concatenate(comm.allgather(s.payload[0]))
    assert np.array_equal(pos, np.argsort(allk, kind="stable"))
    assert np.array_equal(s.back(s.keys[:, 0]), k)
    n += 5
    # two-column keys order lexicographically, routed by the first
    m2 = 0 if rank == 1 else 500
    k2 = np.stack([rng.integers(0, 30, m2), rng.integers(0, 4, m2)], axis=1)
    all2 = np.concatenate(comm.allgather(k2))
    ids2, total2 = DD.number_keys(k2, comm)
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


def jittered(X, per, tol, seed=8):
    """X with every vertex on a periodic face moved by up to 0.2 tol along each
    axis it is inside the box on: the images then sit on either side of the
    edges of the cells of the face that the pairing shares out."""
    rng = np.random.default_rng(seed)
    side = np.isclose(X, X.min(axis=0)) | np.isclose(X, X.max(axis=0))
    move = side[:, per].any(axis=1)[:, None] & ~side
    return X + move * rng.uniform(-0.2, 0.2, X.shape) * tol


def pairing_meshes(size):
    """(X, cells, periodic, tol) by name: periodic along one to three axes, the
    corners and edges of the box reaching their masters through a chain, the
    vertices on the edges of the cells the faces are shared out in at this
    rank count and the images moved across them; a tolerance wider than half
    a midpoint's spacing, so several nodes are within it; and an image moved
    out of it."""
    tol = 1e-3
    n2 = int(DD._face_grid(0, np.zeros(2), np.ones(2), 2 * tol, size)[1][0])
    n3 = int(DD._face_grid(0, np.zeros(3), np.ones(3), 2 * tol, size)[1][0])
    out = {}
    for name, (X, cells), per in (("2d_x", channel2d(n2), [True, False]),
                                  ("2d_xy", channel2d(n2), [True, True]),
                                  ("3d_x", channel3d(n3), [True, False, False]),
                                  ("3d_xz", channel3d(n3), [True, False, True]),
                                  ("3d_xyz", channel3d(n3), [True, True, True]),
                                  ("3d_xyz_shuffled", shuffled(*channel3d(n3)),
                                   [True, True, True])):
        out[name] = (jittered(X, np.array(per), tol), cells, per, tol)
    X, cells = channel2d(16)
    out["2d_wide"] = (X, cells, [True, True], 0.7 / 16)
    X, cells = channel3d(n3)
    X = X.copy()
    k = np.flatnonzero((X[:, 0] == 1) & (X[:, 1] > 0) & (X[:, 1] < 1)
                       & (X[:, 2] > 0) & (X[:, 2] < 1))[-1]
    X[k, 1] += 1.5 * tol
    out["3d_x_off"] = (X, cells, [True, False, False], tol)
    return out


def pairing_job():
    """periodic_masters across the ranks against `Topo`'s masters: the face
    nodes given out scattered, the last rank given none, every rank asking
    for some of them and for nodes off the faces; the same refusal where an
    image is out of the tolerance."""
    comm = MPI.COMM_WORLD
    rank, size = comm.rank, comm.size
    n = 0
    for name, (X, cells, per, tol) in pairing_meshes(size).items():
        try:
            ref = D.Topo(X, cells, per, periodic_tol=tol)
            refused = None
        except ValueError as e:
            ref, refused = D.Topo(X, cells, [False] * X.shape[1], periodic_tol=tol), str(e)
            ref.periodic = np.asarray(per)
        x, p = ref.node_x, np.flatnonzero(ref.periodic)
        face = np.flatnonzero(np.any((x[:, p] < ref.x_min[p] + tol)
                                     | (x[:, p] > ref.x_max[p] - tol), axis=1))
        to = (face * 2654435761) % 7919 % max(size - 1, 1)
        g = face[to == rank] if rank < max(size - 1, 1) else face[:0]
        ask = np.concatenate([face[rank::size], np.flatnonzero(~np.isin(
            np.arange(ref.nnodes), face))[rank::3 * size]])
        args = (g, x[g], ask, ref.periodic, ref.x_min, ref.x_max, tol, comm)
        if refused is not None:
            with pytest.raises(ValueError) as e:
                DD.periodic_masters(*args)
            assert str(e.value) == refused, (name, str(e.value), refused)
            n += 1
            continue
        got = DD.periodic_masters(*args)
        assert np.array_equal(got, ref.master[ask]), name
        n += 1
        if name == "2d_wide":
            continue
        # an image in another cell of the face grid than its node, along each axis
        for a in p:
            grid = DD._face_grid(a, ref.x_min, ref.x_max, 2 * tol, size)
            one = (x[:, a] > ref.x_max[a] - tol) & (ref.master != np.arange(ref.nnodes))
            for b in p[p < a]:
                one &= x[:, b] <= ref.x_max[b] - tol
            i = np.flatnonzero(one)
            c = [s[:, 0] for s in DD._face_cells(x[i], grid, 0.0)]
            m = [s[:, 0] for s in DD._face_cells(x[ref.master[i]], grid, 0.0)]
            assert any(np.any(u != v) for u, v in zip(c, m)), (name, a)
            n += 1
    if size > 1:
        assert n > 0
    return n


def pairing_case():
    """A case of mpi_batch: pairing_job, printing how many checks the rank
    ran."""
    print("%d checks" % pairing_job())


def mesh_case(mesh, per):
    """A case of mpi_batch: job, printing the flags it was given and how many
    checks the rank ran."""
    print("flags %s" % [bool(p) for p in per])
    print("%d checks" % job(mesh, per))


def blocks_case():
    """A case of mpi_batch: blocks_job, printing how many checks the rank
    ran."""
    print("%d checks" % blocks_job())


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


# -------------------------------------------------------------- in process


@pytest.mark.parametrize("name", sorted(CASES))
def test_one_rank_builds_the_serial_tables(tmp_path, name):
    """On one rank nothing is shared out, and the tables are `Topo`'s as they
    are: the same numbering, the same bits, the same held set."""
    pytest.importorskip("h5py")
    mesh, per = write_mesh(tmp_path, name)
    X, cells, ci = read_mesh_h5(mesh)
    ref = D.Topo(X[:, :cells.shape[1] - 1], cells, per)
    assert check_tables(mesh, per, "ptscotch", MPI.COMM_SELF, ref, ci) > 40


def test_the_sample_sort_numbers_keys_by_their_rank_among_the_distinct_ones():
    """An edge's id is its rank among the sorted vertex pairs of the whole mesh,
    so the numbering is the serial one: repeated keys get one id, the ids count
    up from zero in key order, equal keys keep their order, and the way back
    returns every answer to the key it belongs to."""
    rng = np.random.default_rng(4)
    k = rng.integers(0, 200, 5000)
    ids, total = DD.number_keys(k, MPI.COMM_SELF)
    u, inv = np.unique(k, return_inverse=True)
    assert total == len(u) and np.array_equal(ids, inv)
    s = DD.SampleSort(k, MPI.COMM_SELF, payload=(np.arange(5000),))
    assert np.array_equal(s.keys[:, 0], np.sort(k))
    assert np.array_equal(s.payload[0], np.argsort(k, kind="stable"))
    assert np.array_equal(s.back(s.payload[0]), np.arange(5000))


def test_a_block_owner_is_the_split_block_holding_the_id():
    """Every id of range(total) goes to the `_split` block holding it, with
    fewer ids than ranks too."""
    for total in (0, 3, 7, 64, 1001):
        for size in (1, 2, 3, 8):
            want = np.concatenate([np.full(np.subtract(*D._split(total, size, r)[::-1]), r)
                                   for r in range(size)])
            assert np.array_equal(DD.block_owner(np.arange(total), total, size), want)


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


def test_the_periodic_masters_are_topos_on_one_rank():
    """The face nodes paired at the ranks of their places on the face are
    `Topo`'s masters: along one to three axes, through the chains of the box's
    corners and edges, across the edges of the cells the faces are shared out
    in, with several nodes within the tolerance, and refused as `Topo`
    refuses an image out of it."""
    meshes = pairing_meshes(1)
    # a check a mesh, and one an axis of every mesh but the wide and the refused one
    axes = sum(int(np.sum(per)) for name, (_, _, per, _) in meshes.items()
               if name not in ("2d_wide", "3d_x_off"))
    assert pairing_job() == len(meshes) + axes


def test_a_home_is_spread_over_the_ranks_whatever_the_ids():
    """The face nodes of a structured mesh have ids in strides and in runs; a
    home by the id alone would put a stride's on one rank and a run's in one
    block."""
    for ids in (np.arange(0, 64 * 4096, 64), np.arange(4096), np.arange(1 << 40, (1 << 40) + 4096)):
        for size in (3, 64, 1000):
            c = np.bincount(DD.hash_owner(ids, size), minlength=size)
            assert c.max() < 2 * len(ids) / size + 12, (size, c.max())


def test_a_dual_graph_past_petscs_integers_is_refused():
    """A partitioner counts the dual graph's cells and adjacencies in PETSc's
    integers; past 32-bit ones (a billion tets has four billion adjacencies)
    the graph is refused, naming the way out, not handed over to wrap."""
    if np.dtype(D.PETSc.IntType).itemsize > 4:
        pytest.skip("this PETSc has 64-bit indices")
    nb = np.array([[1, -1, -1, -1], [0, -1, -1, -1]])
    with pytest.raises(ValueError, match="needs a PETSc built with 64-bit indices"):
        DD._partition(nb, 1 << 31, MPI.COMM_SELF, "ptscotch")


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
@pytest.mark.parametrize("ranks", [2, 4])
def test_the_building_blocks_across_ranks_give_the_serial_answers(tmp_path, ranks):
    """The sample sort with its way back, the scan, the keyed reductions and
    the exact median, each on keys and values spread unevenly over the ranks --
    the second rank with no key at all, many repeated -- against the same
    computation on all of them together."""
    said = mpi_batch(ranks, blocks_case, {"blocks": []}, tmp_path)["blocks"]
    assert all(int(s.split()[-2]) > 0 for s in said), said


@pytest.mark.slow
@pytest.mark.parametrize("ranks", [2, 3, 4, 5])
def test_the_periodic_masters_across_ranks_are_topos(tmp_path, ranks):
    """The pairing of the periodic faces on two to five ranks, the last given
    no node: `Topo`'s masters for every node asked for, its refusal word for
    word, the face grid's cells edged on the vertices at each rank count."""
    said = mpi_batch(ranks, pairing_case, {"pairing": []}, tmp_path)["pairing"]
    assert all(int(s.split()[-2]) > 0 for s in said), said


# each mesh at one rank count, every count on 2D and 3D meshes; one mesh at one rank
MESH_RANKS = [(name, 2 + k % 3) for k, name in enumerate(sorted(CASES))] + [(sorted(CASES)[0], 1)]


@pytest.mark.slow
@pytest.mark.parametrize("ranks", [1, 2, 3, 4])
def test_the_distributed_tables_are_the_serial_ones(tmp_path, ranks):
    """Read in blocks, and partitioned by PT-Scotch and moved: each rank's
    tables are `Topo`'s restricted to its cells, the globals are the serial
    ones on every rank (the median edge length bit for bit), every node is
    owned by the rank holding the smallest cell that holds it, every reduction
    reaches every copy, the held set holds a master whose image is held, and
    gathered, every table is `Topo`'s array for array. A partitioner this
    PETSc lacks is refused. A rank left without a cell builds its empty share
    and the others' tables are still the serial ones: a rank that raised alone
    there would leave the rest waiting in the next exchange. Every mesh of
    this rank count, one job; one rank in this process."""
    pytest.importorskip("h5py")
    if ranks > 1:
        launcher_or_skip()
    names = [n for n, r in MESH_RANKS if r == ranks]
    meshes = {name: write_mesh(tmp_path / name, name) for name in names}
    if ranks == 1:
        # COMM_WORLD here is one rank
        for name, (mesh, per) in meshes.items():
            assert job(mesh, per) > 40, name
        return
    said = mpi_batch(ranks, mesh_case, {name: [str(mesh), [bool(p) for p in per]]
                                        for name, (mesh, per) in meshes.items()},
                     tmp_path / "jobs")
    for name in names:
        assert all(int(s.split()[-2]) > 0 for s in said[name]), (name, said[name])
        # the flags the job used are the case's, not lost on their way to it
        flags = "flags %s" % [bool(p) for p in CASES[name][2]]
        assert all(flags in s for s in said[name]), (name, said[name])
        if name.endswith("shuffled"):
            assert "ptscotch moved cells" in said[name][0], (name, said[name])


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
