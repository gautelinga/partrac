#!/usr/bin/env python3
"""The divergence cleaner's mesh tables, shared out over the MPI ranks.

`DistTopo` builds `divfree_clean.Topo`'s edge, facet and periodic tables for
each rank's own cells, in the serial numbering (a vertex its row of
`mesh/coordinates`, a cell its row of `mesh/topology`, an edge its rank among
the sorted vertex pairs), so the output is the same at every rank count. What
needs the whole mesh is found by exchange; the building blocks (`Exchange`,
`exscan`, `SampleSort`, `number_keys`, `Reduction`, `reduce_keyed`,
`exact_select`, `periodic_masters`) take arrays and a communicator and hold no
mesh.
"""

import os
import types

import numpy as np
from mpi4py import MPI

import divfree_clean as D
from divfree_clean import LOC_EDGES, PERIODIC_TOL, _split
from petsc4py import PETSc

PARTITIONERS = ("ptscotch", "parmetis")
SAMPLES = 64             # splitter samples a rank contributes to a sample sort

_OPS = {"max": np.maximum, "min": np.minimum, "or": np.logical_or, "sum": np.add}


# ------------------------------------------------------------ building blocks


def _mpitype(dtype):
    return MPI._typedict[np.dtype(dtype).char]


def _displs(counts):
    return np.concatenate([[0], np.cumsum(counts)[:-1]]).astype(np.int64)


class Exchange:
    """A keyed all-to-all: every row goes to the rank it is addressed to, and
    an answer to each received row goes back to where the row came from, in the
    sender's original order. Built once for an addressing and used for any
    number of arrays of rows, of any trailing shape."""

    def __init__(self, comm, dest):
        # ranks as uint16 sort by radix
        dest = np.asarray(dest).astype(np.uint16 if comm.size <= 1 << 16 else np.int64)
        self.comm = comm
        # rows already in rank order are sent as they are, without a copy
        self.order = (None if np.all(dest[1:] >= dest[:-1])
                      else np.argsort(dest, kind="stable"))
        self.sent = np.bincount(dest, minlength=comm.size).astype(np.int64)
        # a rank past the communicator, or a negative one wrapped
        assert len(self.sent) == comm.size, "an address outside the communicator"
        self.received = np.empty(comm.size, np.int64)
        comm.Alltoall(self.sent, self.received)
        self.n_in = int(self.received.sum())
        # MPI counts are C ints: refused on every rank together
        rows = comm.allreduce(max(int(self.sent.sum()), self.n_in), op=MPI.MAX)
        top = int(np.iinfo(np.int32).max)
        if rows > top:
            raise ValueError("a rank exchanges %d rows, past the %d MPI's counts hold; "
                             "more ranks reach it" % (rows, top))

    def _move(self, a, sc, rc):
        a = np.ascontiguousarray(a)
        # one rank sends to itself alone: the rows are the answer
        if self.comm.size == 1:
            return a
        if a.dtype == bool:
            return self._move(a.view(np.uint8), sc, rc).view(bool)
        rest = a.shape[1:]
        w = int(np.prod(rest, dtype=np.int64))
        # a row one datatype, the counts in rows
        out = np.empty((int(rc.sum()),) + rest, a.dtype)
        t = _mpitype(a.dtype)
        row = t.Create_contiguous(w).Commit() if w > 1 else t
        try:
            self.comm.Alltoallv([a, (sc, _displs(sc)), row], [out, (rc, _displs(rc)), row])
        finally:
            if row is not t:
                row.Free()
        return out

    def forward(self, a):
        """The rows sent to their ranks: what this rank receives, grouped by
        the sending rank in rank order and in the sender's order within."""
        a = np.asarray(a)
        assert len(a) == self.sent.sum(), "rows for another exchange"
        return self._move(a if self.order is None else a[self.order], self.sent,
                          self.received)

    def back(self, a):
        """One row for every row received, returned to the rank that sent it
        and put back in that rank's original order."""
        got = self._move(a, self.received, self.sent)
        if self.order is None:
            return got
        out = np.empty_like(got)
        out[self.order] = got
        return out

    def sources(self):
        """The rank every received row came from."""
        return np.repeat(np.arange(self.comm.size), self.received)


def exscan(n, comm):
    """(the sum of n over the lower ranks, the sum over all ranks)."""
    low = comm.exscan(int(n))
    return (0 if comm.rank == 0 else int(low)), int(comm.allreduce(int(n)))


def block_owner(ids, total, size):
    """The rank whose `_split` block of range(total) holds each id."""
    q, rem = divmod(int(total), size)
    ids = np.asarray(ids, np.int64)
    # the first rem blocks hold q + 1 ids, the others q
    big = rem * (q + 1)
    return np.where(ids < big, ids // (q + 1), rem + (ids - big) // max(q, 1))


def _lex_order(k):
    """The stable order of the rows of a 2D key array, first column first."""
    if k.shape[1] == 1:
        return np.argsort(k[:, 0], kind="stable")
    return np.lexsort(k.T[::-1])


class SampleSort:
    """Keys of all ranks sorted across the ranks: this rank's share, every key
    equal to one of them among it, and the lower ranks' keys below it.

    keys is int64, (n,) or (n, w); rows compare lexicographically and are
    routed by their first column, so equal rows meet on one rank. `keys` and
    `payload` are the received rows in sorted order; `back` returns one value a
    sorted row to the position its key came from."""

    def __init__(self, keys, comm, payload=()):
        k = np.asarray(keys, np.int64)
        k = k[:, None] if k.ndim == 1 else k
        size = comm.size
        # sorted here first, so the rows arrive as one sorted run a sender
        self.local = _lex_order(k)
        k = k[self.local]
        if size == 1:
            dest = np.zeros(len(k), np.int64)
        else:
            s = k[:, 0]
            take = s[np.linspace(0, len(s) - 1, SAMPLES).astype(np.int64)] if len(s) \
                else s
            allk = np.sort(np.concatenate(comm.allgather(take)))
            if len(allk):
                split = allk[(np.arange(1, size) * len(allk)) // size]
            else:
                split = np.zeros(size - 1, np.int64)
            dest = np.searchsorted(split, k[:, 0], side="right")
        self.ex = Exchange(comm, dest)
        got = self.ex.forward(k)
        # one run is sorted already
        self.order = _lex_order(got) if size > 1 else slice(None)
        self.keys = got[self.order]
        self.payload = [self.ex.forward(np.asarray(p)[self.local])[self.order]
                        for p in payload]

    def back(self, values):
        v = np.asarray(values)
        out = np.empty_like(v)
        out[self.order] = v
        got = self.ex.back(out)
        out = np.empty_like(got)
        out[self.local] = got
        return out


def _groups(k):
    """(the group of every row of sorted keys, (n,) or (n, w), the first row of
    every group)."""
    new = np.ones(len(k), bool)
    if len(k) > 1:
        new[1:] = k[1:] != k[:-1] if k.ndim == 1 else np.any(k[1:] != k[:-1], axis=1)
    return np.cumsum(new) - 1, np.flatnonzero(new)


def number_keys(keys, comm):
    """(the id of every key, the number of distinct keys). A key's id is its
    rank among the distinct keys of all ranks together, so it is the serial
    numbering whatever the rank count."""
    s = SampleSort(keys, comm)
    grp, first = _groups(s.keys)
    off, total = exscan(len(first), comm)
    return s.back(off + grp), total


def _reduceat(v, first, op):
    """The op over the runs of v starting at first; ops a column for 2D v."""
    if isinstance(op, str):
        return _OPS[op].reduceat(v, first) if len(first) else v[:0]
    return np.stack([_reduceat(v[:, j], first, o) for j, o in enumerate(op)], axis=1)


class Reduction:
    """The op over every row of all ranks with the same key, at the rank the
    exchange sends it to, returned to every row: built once for the keys and
    used for any number of value arrays. Each key's values are reduced in the
    order of their sending ranks, and in each sender's order within."""

    def __init__(self, ex, keys):
        self.ex = ex
        k = ex.forward(np.asarray(keys, np.int64))
        self.order = np.argsort(k, kind="stable")
        self.group, self.first = _groups(k[self.order])

    def __call__(self, values, op):
        """op 'max', 'min', 'or' or 'sum', or one a column of 2D values."""
        v = self.ex.forward(values)
        ans = np.empty_like(v)
        ans[self.order] = _reduceat(v[self.order], self.first, op)[self.group]
        return self.ex.back(ans)


def reduce_keyed(keys, values, op, comm, total):
    """The op ('max', 'min', 'or', 'sum') over every value of all ranks with the
    same key, returned to every one of them. The keys are ids in range(total),
    reduced at the rank whose `_split` block holds the id."""
    ex = Exchange(comm, block_owner(keys, total, comm.size))
    return Reduction(ex, keys)(np.asarray(values), op)


def _sortable(v):
    """float64 as int64 in the same order."""
    i = np.ascontiguousarray(v, np.float64).view(np.int64)
    return np.where(i < 0, i ^ np.int64(0x7fffffffffffffff), i)


def _unsortable(k):
    k = np.int64(k)
    return float(np.array([k ^ np.int64(0x7fffffffffffffff) if k < 0 else k],
                          np.int64).view(np.float64)[0])


def exact_select(values, ks, comm):
    """The ks-th smallest (0-based) of the values of all ranks together, each
    exactly one of the values: bisection on their bit patterns with allreduced
    counts, 64 rounds at most, whatever the values are."""
    key = np.sort(_sortable(values))
    ks = [int(k) for k in ks]
    big = np.iinfo(np.int64)
    lo0 = comm.allreduce(int(key[0]) if len(key) else big.max, op=MPI.MIN)
    hi0 = comm.allreduce(int(key[-1]) if len(key) else big.min, op=MPI.MAX)
    lo, hi = [lo0] * len(ks), [hi0] * len(ks)
    while any(a < b for a, b in zip(lo, hi)):
        mid = [a + (b - a) // 2 for a, b in zip(lo, hi)]
        c = np.searchsorted(key, np.array(mid, np.int64), side="right").astype(np.int64)
        comm.Allreduce(MPI.IN_PLACE, c, op=MPI.SUM)
        for i, k in enumerate(ks):
            if lo[i] < hi[i]:
                if c[i] >= k + 1:
                    hi[i] = mid[i]
                else:
                    lo[i] = mid[i] + 1
    return [_unsortable(v) for v in lo]


def exact_median(values, comm):
    """np.median of the values of all ranks together, bit for bit: the middle
    value, or the mean of the two middle ones."""
    n = comm.allreduce(len(values))
    if n % 2:
        return exact_select(values, [(n - 1) // 2], comm)[0]
    a, b = exact_select(values, [n // 2 - 1, n // 2], comm)
    return float(np.mean(np.array([a, b])))


# ------------------------------------------------------------ periodic pairing


# 2^64 / golden ratio: Fibonacci hashing
_MIX = np.uint64(0x9E3779B97F4A7C15)


def hash_owner(ids, size):
    """A rank for every id, spread evenly whatever pattern the ids follow."""
    h = np.asarray(ids, np.int64).astype(np.uint64) * _MIX
    return ((h >> np.uint64(32)) % np.uint64(size)).astype(np.int64)


def _face_grid(a, x_min, x_max, r, size):
    """(the axes but a, the cells along each, the scale to a cell index, the
    origin): about four cells a rank over the face of axis a, none narrower
    than 4 r."""
    axes = [b for b in range(len(x_min)) if b != a]
    want = int(np.ceil((4 * size) ** (1.0 / max(len(axes), 1))))
    ncell, scale = [], []
    for b in axes:
        ext = x_max[b] - x_min[b]
        n = 1 if ext <= 0 else want if 4 * r * want <= ext else max(1, int(ext // (4 * r)))
        ncell.append(n)
        scale.append(n / ext if ext > 0 else 0.0)
    return axes, np.array(ncell, np.int64), scale, [x_min[b] for b in axes]


def _face_cells(x, grid, w):
    """The first and last cell of the points x widened by w on either side,
    (n, 2) an axis of the grid: at most two apart by one."""
    axes, ncell, scale, x0 = grid
    return [np.stack([np.clip(np.floor((v - o) * k), 0, n - 1).astype(np.int64)
                      for v in (x[:, b] - w, x[:, b] + w)], axis=1)
            for b, n, k, o in zip(axes, ncell, scale, x0)]


def _cell_rank(c, ncell, size):
    """The rank of every cell, (n, axes) indices, round the ranks in turn."""
    lin = np.zeros(len(c), np.int64)
    for j, n in enumerate(ncell):
        lin = lin * n + c[:, j]
    return lin % size


def _pair_axis(hg, hx, link, a, x_min, x_max, tol, comm):
    """`Topo._periodic`'s links along axis a of the face nodes this rank is
    home to, (hg, hx) sorted by id: every node on the high face not linked yet
    to the nearest node on the low face to its image, found at the rank of the
    image's cell, where every low node within 2 tol of the cell is sent."""
    size, dim = comm.size, hx.shape[1]
    lo = np.flatnonzero(hx[:, a] < x_min[a] + tol)
    hi = np.flatnonzero((hx[:, a] > x_max[a] - tol) & (link == hg))
    r = 2.0 * tol
    grid = _face_grid(a, x_min, x_max, r, size)
    ncell = grid[1]
    span = _face_cells(hx[lo], grid, r)
    # a low node to every rank of a cell within r, once a rank
    copies = []
    for off in np.ndindex(*([2] * len(span))):
        ok = np.ones(len(lo), bool)
        c = np.empty((len(lo), len(span)), np.int64)
        for j, s in enumerate(span):
            c[:, j] = s[:, 0] + off[j]
            ok &= c[:, j] <= s[:, 1]
        copies.append(np.stack([lo, _cell_rank(c, ncell, size)], axis=1)[ok])
    lr = np.unique(np.concatenate(copies), axis=0) if copies else np.zeros((0, 2), np.int64)
    hspan = _face_cells(hx[hi], grid, 0.0)
    hc = (np.stack([s[:, 0] for s in hspan], axis=1) if hspan
          else np.zeros((len(hi), 0), np.int64))
    # a row: the id's bits, the low ones as they are and the high ones as -1 - id
    rows = np.empty((len(lr) + len(hi), 1 + dim))
    code = rows.view(np.int64)
    code[:len(lr), 0] = hg[lr[:, 0]]
    code[len(lr):, 0] = -1 - hg[hi]
    rows[:len(lr), 1:] = hx[lr[:, 0]]
    rows[len(lr):, 1:] = hx[hi]
    ex = Exchange(comm, np.concatenate([lr[:, 1], _cell_rank(hc, ncell, size)]))
    got = ex.forward(rows)
    del rows
    g = got.view(np.int64)[:, 0]
    isq = g < 0
    lg, lx = g[~isq], got[~isq, 1:]
    o = np.argsort(lg)
    lg, lx = lg[o], lx[o]
    qg, x = -1 - g[isq], got[isq, 1:]
    q = x.copy()
    q[:, a] -= x_max[a] - x_min[a]
    if len(lg) and len(q):
        from scipy.spatial import cKDTree
        d, j = cKDTree(lx).query(q)
        ok = d < tol
    else:
        j, ok = np.zeros(len(q), np.int64), np.zeros(len(q), bool)
    bad = int(comm.allreduce(int((~ok).sum())))
    if bad:
        k = np.flatnonzero(~ok)
        first = ((int(qg[k].min()), x[k[np.argmin(qg[k])]].tolist()) if len(k)
                 else (np.iinfo(np.int64).max, []))
        raise D._not_periodic(a, bad, tol, comm.allreduce(first, op=MPI.MIN)[1])
    ans = np.full(len(g), -1, np.int64)
    ans[isq] = lg[j]
    link = link.copy()
    link[hi] = ex.back(ans)[len(lr):]
    return link


def _chase(hg, link, steps, comm):
    """link applied steps times to every node this rank is home to, each step
    asking the next node's home for its link; a node linked to itself ends."""
    m = link.copy()
    act = np.flatnonzero(m != hg)
    for _ in range(steps - 1):
        if not comm.allreduce(len(act)):
            break
        ex = Exchange(comm, hash_owner(m[act], comm.size))
        nxt = ex.back(link[np.searchsorted(hg, ex.forward(m[act]))])
        moved = nxt != m[act]
        m[act] = nxt
        act = act[moved]
    return m


def periodic_masters(gid, x, ask, periodic, x_min, x_max, tol, comm):
    """`Topo._periodic`'s master of every id in ask, the pairing being of the
    nodes (gid, x) of all ranks together, each node on one rank only; an id
    given by no rank is its own master. No rank holds more than its share:
    each node is kept at a home rank by its id, and paired along each axis at
    the rank of its place on the face."""
    size = comm.size
    home = Exchange(comm, hash_owner(gid, size))
    hg = home.forward(np.asarray(gid, np.int64))
    hx = home.forward(np.asarray(x, float))
    del home
    o = np.argsort(hg)
    hg, hx = hg[o], hx[o]
    link = hg.copy()
    for a in np.nonzero(periodic)[0]:
        link = _pair_axis(hg, hx, link, a, x_min, x_max, tol, comm)
    master = _chase(hg, link, len(x_min) + 1, comm)
    ask = np.asarray(ask, np.int64)
    ex = Exchange(comm, hash_owner(ask, size))
    q = ex.forward(ask)
    at = np.minimum(np.searchsorted(hg, q), max(len(hg) - 1, 0))
    hit = hg[at] == q if len(hg) else np.zeros(len(q), bool)
    return ex.back(np.where(hit, master[at] if len(hg) else q, q))


# ------------------------------------------------------------------ facets


def facet_pass(gcells, cell_gid, nverts, comm, neighbours=True):
    """(whether one cell alone holds each facet; the other cell of a facet two
    cells hold, -1 where none does, or None without neighbours), both (n, nv)
    for the rank's cells given by global vertex ids: the exterior facets, and
    the dual graph the partition is found on."""
    n, nv = gcells.shape
    key = D.facet_keys(D.cell_facets(gcells), nverts)
    payload = (np.repeat(np.asarray(cell_gid, np.int64), nv),) if neighbours else ()
    s = SampleSort(key, comm, payload)
    grp, first = _groups(s.keys)
    cnt = np.diff(np.append(first, len(s.keys)))[grp]
    if not neighbours:
        return s.back(cnt == 1).reshape(n, nv), None
    cell = s.payload[0]
    pos = np.arange(len(cell))
    start = first[grp]
    partner = np.where(pos == start, pos + 1, start)
    nb = np.where(cnt == 2, cell[np.minimum(partner, len(cell) - 1)], -1) if len(cell) \
        else np.zeros(0, np.int64)
    both = s.back(np.stack([cnt, nb], axis=1))
    return (both[:, 0] == 1).reshape(n, nv), both[:, 1].reshape(n, nv)


def _partition(nb, ncells, comm, kind):
    """The rank each of the rank's cells goes to: PETSc's MatPartitioning of
    the dual graph, cells adjacent where they share a facet; nb the other cell
    of each facet, -1 where there is none."""
    has = nb >= 0
    # each row's neighbours ascending, the missing ones sorted past them
    rows = np.sort(np.where(has, nb, np.iinfo(np.int64).max), axis=1)
    cols = rows[rows < np.iinfo(np.int64).max]
    # the partitioners count the graph's ids and adjacencies in PETSc's integers
    adj = int(comm.allreduce(len(cols)))
    if max(int(ncells), adj) > np.iinfo(PETSc.IntType).max:
        raise ValueError("partition=%s: a dual graph of %d cells and %d adjacencies needs a "
                         "PETSc built with 64-bit indices; partition=blocks needs none"
                         % (kind, ncells, adj))
    indptr = np.concatenate([[0], np.cumsum(has.sum(axis=1))]).astype(PETSc.IntType)
    n = len(nb)
    A = PETSc.Mat().createAIJ(size=((n, ncells), (n, ncells)),
                              csr=(indptr, cols.astype(PETSc.IntType), np.ones(len(cols))),
                              comm=comm)
    A.assemble()
    # Scotch 7 runs a thread a core in every rank; ranks on a node wait on each other's
    os.environ.setdefault("SCOTCH_PTHREAD_NUMBER", "1")
    part = PETSc.MatPartitioning().create(comm=comm)
    part.setAdjacency(A)
    part.setType(kind)
    iset = PETSc.IS().create(comm=comm)
    part.apply(iset)
    dest = iset.getIndices().astype(np.int64)
    for o in (iset, part, A):
        o.destroy()
    return dest


def check_partition(partition):
    """The partition option refused where it names nothing or a partitioner
    this PETSc lacks, naming the ones it has."""
    if partition == "blocks":
        return
    if partition not in PARTITIONERS:
        raise ValueError("partition '%s' is none of blocks, %s"
                         % (partition, ", ".join(PARTITIONERS)))
    if not PETSc.Sys.hasExternalPackage(partition):
        have = [p for p in PARTITIONERS if PETSc.Sys.hasExternalPackage(p)]
        raise ValueError("partition=%s needs a PETSc built with %s, and this one has %s; "
                         "partition=blocks needs none"
                         % (partition, partition, " and ".join(have) or "no partitioner"))


# ------------------------------------------------------------------ topology


class DistTopo(D.Ranks):
    """`Topo`'s tables for this rank's cells, in the serial numbering, under
    `Topo`'s names: what the rank holds is a mesh of its own to every cell-wise
    function of `divfree_clean`, and what needs the others is a reduction on
    `comm`.

    Globals, the same on every rank: dim, nv (vertices a cell), ncells_global,
    nverts_global, nedges_global, nnodes_global (nverts_global + nedges_global),
    periodic, x_min, x_max, vol_total, partition (how the cells were shared
    out; 'blocks' on one rank).

    Per rank, local indices into the rank's own vertices, edges and nodes, of
    which it holds ncells, nverts, nedges and nnodes = nverts + nedges:

      cells       (nc, nv)        the rank's cells, local vertices; ascending cell_gid
      cell_gid    (nc,)           their serial ids, the rows of mesh/topology
      cell_index  (nc,)           their mesh/cell_indices labels, or None
      X           (lv, dim)       the coordinates of the vertices they touch
      edges       (le, 2)         the local vertex pairs of their edges, sorted
      cell_edge   (nc, ne)        local edges in LOC_EDGES order
      node_gid    (ln,)           serial node ids, ascending: the vertices' rows of
                                  mesh/coordinates, then nverts_global + the edges'
      node_x      (ln, dim)       vertices, then edge midpoints
      J, Jinv, vol, facet_n       the serial geometry of the rank's cells
      facet_boundary (nc, nv); cell_boundary (nc,)
      master      (ln,)           the serial node id of each node's periodic master
      boundary_node, held_node, at_rest (ln,); boundary_edge, held_edge (le,)
      node_owner  (ln,)           the rank that holds the serial-smallest cell
                                  containing the node
      first_cell  (ln,)           that cell's serial id
      node_shared (ln,)           whether another rank holds the node too
      node_turn   (ln,)           the rank a node's copies are reduced at: the
                                  first rank holding it at or after gid mod size,
                                  round the ranks

    Nothing of global size is held.
    """

    @classmethod
    def read(cls, mesh_path, periodic, comm=None, periodic_tol=PERIODIC_TOL,
             partition="ptscotch"):
        """The tables of a dolfin HDF5 mesh. Every rank reads its contiguous
        block of cells and of vertices; the cells then move to the ranks the
        partitioner gives them ('ptscotch', 'parmetis', or 'blocks' to keep
        the blocks as read), and every table is built on those."""
        import h5py
        comm = MPI.COMM_WORLD if comm is None else comm
        check_partition(partition)
        self = cls()
        self.comm = comm
        size, rank = comm.size, comm.rank
        self.partition = partition if size > 1 else "blocks"
        with h5py.File(str(mesh_path), "r") as f:
            T = f["mesh/topology"]
            self.ncells_global, self.nv = int(T.shape[0]), int(T.shape[1])
            self.dim = self.nv - 1
            c0, c1 = _split(self.ncells_global, size, rank)
            gcells = np.asarray(T[c0:c1]).astype(np.int64)
            ci = (np.asarray(f["mesh/cell_indices"][c0:c1]).astype(np.int64)
                  if "mesh/cell_indices" in f else None)
            Xd = f["mesh/coordinates"]
            self.nverts_global = int(Xd.shape[0])
            # edge and facet keys pack two vertex ids in an int64
            if self.nverts_global ** 2 >= 2 ** 63:
                raise ValueError("%s: %d vertices, past the %d that a pair of them packed in "
                                 "64 bits holds" % (mesh_path, self.nverts_global, int(2 ** 31.5)))
            v0, v1 = _split(self.nverts_global, size, rank)
            Xb = np.ascontiguousarray(np.asarray(Xd[v0:v1])[:, :self.dim], dtype=float)
        self.periodic = np.asarray(periodic, dtype=bool)[:self.dim]
        self.x_min = np.full(self.dim, np.inf) if not len(Xb) else Xb.min(axis=0)
        self.x_max = np.full(self.dim, -np.inf) if not len(Xb) else Xb.max(axis=0)
        comm.Allreduce(MPI.IN_PLACE, self.x_min, op=MPI.MIN)
        comm.Allreduce(MPI.IN_PLACE, self.x_max, op=MPI.MAX)
        gid = np.arange(c0, c1, dtype=np.int64)
        moved = self.partition != "blocks"
        ext, nb = facet_pass(gcells, gid, self.nverts_global, comm, neighbours=moved)
        if moved:
            ex = Exchange(comm, _partition(nb, self.ncells_global, comm, partition))
            del nb
            gcells, gid, ext = ex.forward(gcells), ex.forward(gid), ex.forward(ext)
            ci = ex.forward(ci) if ci is not None else None
            o = np.argsort(gid)
            gcells, gid, ext = gcells[o], gid[o], ext[o]
            ci = ci[o] if ci is not None else None
        self.cell_gid, self.cell_index = gid, ci
        self._build(gcells, ext, Xb, v0, periodic_tol)
        return self

    # ---- the tables, on the final cells
    def _build(self, gcells, ext, Xb, v0, tol):
        comm = self.comm
        self.ncells = len(gcells)
        vert_gid, inv = np.unique(gcells, return_inverse=True)
        self.cells = inv.reshape(gcells.shape).astype(np.int32)
        self.nverts = len(vert_gid)
        # the coordinates from the vertices' block owners
        ex = Exchange(comm, block_owner(vert_gid, self.nverts_global, comm.size))
        self.X = ex.back(Xb[ex.forward(vert_gid) - v0])
        self.loc = np.array(LOC_EDGES[self.nv])
        self._edges(vert_gid)
        tol = D.match_tol(self.X, self.edges, tol, comm)
        # the serial formulas, on the rank's cells
        geo = types.SimpleNamespace(X=self.X, cells=self.cells, ncells=self.ncells,
                                    nv=self.nv, dim=self.dim)
        D.Topo._geometry(geo)
        self.J, self.Jinv = geo.J, geo.Jinv
        self.vol, self.facet_n = geo.vol, geo.facet_n
        del geo
        self._facets(ext, tol)
        shared = self._owners()
        face = self._periodic(Xb, v0, tol)
        # a node no other rank holds, in no periodic class, is its master class alone
        sel = np.flatnonzero(shared | face)
        m = self.master[sel]
        self._masters = sel, Reduction(Exchange(comm, block_owner(m, self.nnodes_global,
                                                                  comm.size)), m)
        self._exterior()
        self.vol_total = float(np.sum(comm.allgather(float(np.sum(self.vol)))))
        self.set_held(None)

    def _edges(self, vert_gid):
        # local vertices are in serial order, so a sorted local pair is a sorted serial one
        pair = np.sort(self.cells[:, self.loc], axis=2).reshape(-1, 2)
        g = vert_gid[pair]
        key, at, inv = np.unique(g[:, 0] * np.int64(self.nverts_global) + g[:, 1],
                                 return_index=True, return_inverse=True)
        del g
        # the distinct keys ascend, and so do their ids
        edge_gid, self.nedges_global = number_keys(key, self.comm)
        self.nedges = len(edge_gid)
        self.cell_edge = inv.reshape(self.ncells, len(self.loc)).astype(np.int32)
        self.edges = pair[at]
        self.nnodes_global = self.nverts_global + self.nedges_global
        self.nnodes = self.nverts + self.nedges
        self.node_gid = np.concatenate([vert_gid, self.nverts_global + edge_gid])
        self.node_x = np.vstack([self.X,
                                 0.5 * (self.X[self.edges[:, 0]] + self.X[self.edges[:, 1]])])

    def _facets(self, ext, tol):
        """The boundary the rule classifies: exterior facets the mesh does not
        pair."""
        per = np.zeros((self.ncells, self.nv), bool)
        for o in range(self.nv):
            f = self.cells[:, [i for i in range(self.nv) if i != o]]
            for d in np.nonzero(self.periodic)[0]:
                xd = self.X[f, d]
                for v in (self.x_min[d], self.x_max[d]):
                    per[:, o] |= np.all(np.abs(xd - v) < tol, axis=1)
        self.facet_boundary = ext & ~per
        self.cell_boundary = self.facet_boundary.any(axis=1)

    def cell_nodes(self):
        """The local nodes of every cell: its vertices, then its edges."""
        return np.concatenate([self.cells, self.nverts + self.cell_edge], axis=1)

    def _owners(self):
        """node_owner, first_cell, node_shared (held by other ranks too),
        node_turn (the rank reducing it) and that reduction; returns node_shared."""
        # the cells ascend, so the last write of the reversed list is the smallest
        cn = self.cell_nodes()
        first = np.full(self.nnodes, self.ncells_global, np.int64)
        first[cn[::-1].ravel()] = np.repeat(self.cell_gid[::-1], cn.shape[1])
        del cn
        comm, g = self.comm, self.node_gid
        size, rank = comm.size, comm.rank
        # the first cell's rank, a rank touching the node by turns, the copies
        v = np.stack([first * size + rank, (rank - g) % size, np.ones(self.nnodes, np.int64)],
                     axis=1)
        del first
        red = reduce_keyed(g, v, ("min", "min", "sum"), comm, self.nnodes_global)
        self.first_cell = red[:, 0] // size
        self.node_owner = (red[:, 0] % size).astype(np.int32)
        shared = red[:, 2] > 1
        self.node_shared = shared
        self.node_turn = ((g + red[:, 1]) % size).astype(np.int32)
        self._nodes = np.flatnonzero(shared), self.shared_reduction()
        return shared

    @staticmethod
    def _reduce(plan, values, op):
        sel, red = plan
        out = np.array(values)
        out[sel] = red(out[sel], op)
        return out

    def reduce_nodes(self, values, op):
        """The op ('max', 'min', 'or', 'sum') over every rank's copy of each of
        the rank's nodes, returned to every rank that touches it. With 'sum'
        each copy counts once."""
        return self._reduce(self._nodes, values, op)

    def reduce_masters(self, values, op):
        """The op over every copy of every periodic image of each node's master,
        returned to each: a master class decides together. With 'sum' every
        rank's copy of every image counts once, so the callers pass each rank's
        partial sums, gathered at `class_index`, and they add up to the class's
        total as `Topo.reduce_masters` adds its one copy of each node."""
        return self._reduce(self._masters, values, op)

    def median(self, values):
        """np.median of every rank's values together, exactly."""
        return exact_median(values, self.comm)

    def _periodic(self, Xb, v0, tol):
        """The periodic faces' nodes, every rank's, paired as `Topo` pairs the
        whole mesh by `periodic_masters`: no other node is in a pairing, so
        the masters are the serial ones. Returns which of the rank's nodes are
        on those faces."""
        comm = self.comm
        self.master = self.node_gid.copy()
        if not self.periodic.any():
            return np.zeros(self.nnodes, bool)

        def on_face(x):
            m = np.zeros(len(x), bool)
            for a in np.nonzero(self.periodic)[0]:
                m |= (x[:, a] < self.x_min[a] + tol) | (x[:, a] > self.x_max[a] - tol)
            return m

        # the vertices by their block, so one in no cell is paired too; the
        # midpoints by their owner
        lv = self.nverts
        vb = np.nonzero(on_face(Xb))[0]
        ex = self.node_x[lv:]
        eo = np.nonzero(on_face(ex) & (self.node_owner[lv:] == comm.rank))[0]
        g = np.concatenate([v0 + vb, self.node_gid[lv + eo]])
        x = np.vstack([Xb[vb], ex[eo]])
        # every rank's copy of a face node, the same bits as its sender's
        hit = on_face(self.node_x)
        self.master[hit] = periodic_masters(g, x, self.node_gid[hit], self.periodic,
                                            self.x_min, self.x_max, tol, comm)
        return hit

    def _exterior(self):
        """A vertex is on the boundary if any rank's cell puts it there; an edge
        if any image of its master is."""
        ef = self.facet_boundary
        node = np.zeros(self.nnodes, bool)
        for o in range(self.nv):
            f = [i for i in range(self.nv) if i != o]
            loc = [j for j, (a, b) in enumerate(LOC_EDGES[self.nv]) if a != o and b != o]
            k = np.nonzero(ef[:, o])[0]
            if len(k):
                node[self.cells[np.ix_(k, f)].ravel()] = True
                node[self.nverts + self.cell_edge[np.ix_(k, loc)].ravel()] = True
        lv = self.nverts
        node[:lv] = self.reduce_nodes(node, "or")[:lv]
        node[lv:] = self.reduce_masters(node, "or")[lv:]
        self.boundary_edge = node[lv:].copy()
        self.boundary_node = node

    def set_held(self, at_rest):
        """The held set of `Topo.set_held`, at_rest given on the rank's nodes
        (the same on every copy of a node): held where at rest on the boundary,
        and a master with every image of it held if any image is."""
        lv = self.nverts
        self.at_rest = (np.zeros(self.nnodes, bool) if at_rest is None
                        else np.asarray(at_rest, bool))
        held = self.boundary_node & self.at_rest
        # nothing at rest holds nothing
        if at_rest is not None:
            held[lv:] = self.reduce_masters(held, "or")[lv:]
        self.held_node = held
        self.held_edge = held[lv:].copy()
