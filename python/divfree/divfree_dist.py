#!/usr/bin/env python3
"""The divergence cleaner's mesh tables, shared out over the MPI ranks.

`divfree_clean.Topo` builds the edge, facet and periodic tables of the whole
mesh on every rank, which at tens of millions of cells is more than a rank can
hold. `DistTopo` builds the same tables for each rank's own cells only, in the
serial numbering: a vertex is its row of `mesh/coordinates`, a cell its row of
`mesh/topology`, an edge its rank among the sorted vertex pairs of the whole
mesh, so everything written from these tables is the same at every rank count.

What needs the whole mesh is found by exchange rather than by holding it:

  * the edges and facets by a sample sort of their vertex keys across the
    ranks, numbered by an exclusive scan of the distinct keys each rank holds;
  * the periodic masters on the nodes of the periodic faces alone, a surface,
    which every rank pairs as `Topo` pairs the whole mesh;
  * a flag or value of a node by a reduction over the ranks that touch it (or
    over its periodic images), at a rank that owns it;
  * the median edge length, which the penalty is measured in, by bisection on
    the value with allreduced counts: exact, not gathered.

The cells are read in contiguous blocks and moved to the ranks a graph
partitioner gives them on their dual graph, so a rank's cells are a compact
region and few of its nodes are shared; or left in their blocks. A rank may be
left with no cell, and builds its empty share.

The building blocks -- `Exchange`, `exscan`, `sample_sort`, `number_keys`,
`reduce_keyed`, `exact_select` -- take arrays and a communicator and hold no
mesh, so the field, the step and the writer are built from them too.
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
        dest = np.asarray(dest, np.int64)
        self.comm = comm
        # rows already in rank order are sent as they are, without a copy
        self.order = (None if np.all(dest[1:] >= dest[:-1])
                      else np.argsort(dest, kind="stable"))
        self.sent = np.bincount(dest, minlength=comm.size).astype(np.int64)
        self.received = np.empty(comm.size, np.int64)
        comm.Alltoall(self.sent, self.received)
        self.n_in = int(self.received.sum())

    def _move(self, a, sc, rc):
        a = np.ascontiguousarray(a)
        # one rank sends to itself alone: the rows are the answer
        if self.comm.size == 1:
            return a
        if a.dtype == bool:
            return self._move(a.view(np.uint8), sc, rc).view(bool)
        rest = a.shape[1:]
        w = int(np.prod(rest, dtype=np.int64))
        out = np.empty((int(rc.sum()),) + rest, a.dtype)
        t = _mpitype(a.dtype)
        self.comm.Alltoallv([a, (sc * w, _displs(sc) * w), t],
                            [out, (rc * w, _displs(rc) * w), t])
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
    starts = np.array([_split(total, size, r)[0] for r in range(size)], np.int64)
    return np.searchsorted(starts, np.asarray(ids, np.int64), side="right") - 1


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
        if size == 1:
            dest = np.zeros(len(k), np.int64)
        else:
            s = np.sort(k[:, 0])
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
        self.order = _lex_order(got)
        self.keys = got[self.order]
        self.payload = [self.ex.forward(p)[self.order] for p in payload]

    def back(self, values):
        v = np.asarray(values)
        out = np.empty_like(v)
        out[self.order] = v
        return self.ex.back(out)


def sample_sort(keys, comm, payload=()):
    """`SampleSort` of the keys: the sorted share, its payload rows, and the
    way back to the keys' own ranks and positions."""
    return SampleSort(keys, comm, payload)


def _groups(k):
    """(the group of every row of sorted 2D keys, the first row of every group)."""
    new = np.ones(len(k), bool)
    if len(k) > 1:
        new[1:] = np.any(k[1:] != k[:-1], axis=1)
    return np.cumsum(new) - 1, np.flatnonzero(new)


def number_keys(keys, comm, payload=()):
    """(the id of every key, the number of distinct keys, the multiplicity of
    every key, the sort). A key's id is its rank among the distinct keys of all
    ranks together, so it is the serial numbering whatever the rank count."""
    s = SampleSort(keys, comm, payload)
    grp, first = _groups(s.keys)
    off, total = exscan(len(first), comm)
    cnt = np.diff(np.append(first, len(s.keys)))
    s.group, s.first, s.count, s.offset = grp, first, cnt, off
    return s.back(off + grp), total, s.back(cnt[grp]), s


def _reduce_sorted(keys, values, op):
    """(the distinct keys, the op over the values of each)."""
    o = np.argsort(keys, kind="stable")
    k = keys[o]
    grp, first = _groups(k.reshape(-1, 1))
    red = _OPS[op].reduceat(values[o], first) if len(first) else values[:0]
    return k[first], grp, o, red


def _reduce_exchange(ex, keys, values, op):
    """Every received copy of a key reduced with the others, each sender
    answered with the reduction of its key."""
    k = ex.forward(keys)
    v = ex.forward(values)
    _, grp, o, red = _reduce_sorted(k, v, op)
    ans = np.empty_like(v)
    ans[o] = red[grp]
    return ex.back(ans)


def reduce_keyed(keys, values, op, comm, total):
    """The op ('max', 'min', 'or', 'sum') over every value of all ranks with the
    same key, returned to every one of them. The keys are ids in range(total),
    reduced at the rank whose `_split` block holds the id."""
    ex = Exchange(comm, block_owner(keys, total, comm.size))
    return _reduce_exchange(ex, np.asarray(keys, np.int64), np.asarray(values), op)


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


# ------------------------------------------------------------------ facets


def _facet_keys(gcells, nverts):
    """(the sorted global vertices of every facet, (n*nv, nv-1); their int64
    keys, (n*nv, 1) in 2D and (n*nv, 2) in 3D). Facet o of a cell misses its
    vertex o."""
    nv = gcells.shape[1]
    F = np.stack([np.sort(gcells[:, [i for i in range(nv) if i != o]], axis=1)
                  for o in range(nv)], axis=1).reshape(-1, nv - 1)
    first = F[:, 0].astype(np.int64) * np.int64(nverts) + F[:, 1]
    key = first[:, None] if nv == 3 else np.stack([first, F[:, 2].astype(np.int64)], axis=1)
    return F, key


def facet_pass(gcells, cell_gid, nverts, comm):
    """(how many cells hold each facet, the other cell of a facet two cells
    hold or -1), both (n, nv) for the rank's cells given by global vertex ids.
    Run on the contiguous blocks it is the dual graph the partition is found
    on; run on the final cells it is the exterior facets."""
    n, nv = gcells.shape
    _, key = _facet_keys(gcells, nverts)
    owner = np.repeat(np.asarray(cell_gid, np.int64), nv)
    _, _, cnt, s = number_keys(key, comm, payload=(owner,))
    cell = s.payload[0]
    pos = np.arange(len(cell))
    start = s.first[s.group]
    two = s.count[s.group] == 2
    partner = np.where(pos == start, pos + 1, start)
    nb = np.where(two, cell[np.minimum(partner, len(cell) - 1)], -1) if len(cell) \
        else np.zeros(0, np.int64)
    return cnt.reshape(n, nv), s.back(nb).reshape(n, nv)


def _partition(gcells, gid, ncells, nverts, comm, kind):
    """The rank each of the rank's cells goes to: PETSc's MatPartitioning of
    the dual graph, cells adjacent where they share a facet."""
    _, nb = facet_pass(gcells, gid, nverts, comm)
    has = nb >= 0
    # each row's neighbours ascending, the missing ones sorted past them
    rows = np.sort(np.where(has, nb, np.iinfo(np.int64).max), axis=1)
    cols = rows[rows < np.iinfo(np.int64).max]
    indptr = np.concatenate([[0], np.cumsum(has.sum(axis=1))]).astype(PETSc.IntType)
    n = len(gcells)
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
    periodic, x_min, x_max, vol_total, mesh_h (the serial `mesh_h` bit for bit),
    node_block (the `_split` block of the serial node ids this rank writes),
    vert_block (its block of `mesh/coordinates`), partition (how the cells were
    shared out; 'blocks' on one rank).

    Per rank, local indices into the rank's own vertices, edges and nodes, of
    which it holds ncells, nverts, nedges and nnodes = nverts + nedges:

      cells       (nc, nv)        the rank's cells, local vertices; ascending cell_gid
      cell_gid    (nc,)           their serial ids, the rows of mesh/topology
      cell_index  (nc,)           their mesh/cell_indices labels, or None
      vert_gid    (lv,)           serial ids of the vertices they touch, ascending
      X           (lv, dim)       their coordinates
      edges       (le, 2)         the local vertex pairs of their edges, sorted
      edge_gid    (le,)           serial edge ids, ascending
      cell_edge   (nc, ne)        local edges in LOC_EDGES order
      node_gid    (ln,)           vert_gid, then nverts_global + edge_gid, ascending
      node_x      (ln, dim)       vertices, then edge midpoints
      J, Jinv, detJ, vol, facet_n the serial geometry of the rank's cells
      facet_verts (nc, nv, nv-1)  local vertices of facet o, the one missing vertex o
      facet_exterior, facet_periodic, facet_boundary (nc, nv); cell_boundary (nc,)
      master      (ln,)           the serial node id of each node's periodic master
      boundary_node, held_node, at_rest (ln,); boundary_edge, held_edge (le,)
      node_owner  (ln,)           the rank that holds the serial-smallest cell
                                  containing the node
      first_cell  (ln,)           that cell's serial id

    Nothing of global size is held but the periodic faces' nodes while they are
    paired.
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
            self.vert_block = _split(self.nverts_global, size, rank)
            v0, v1 = self.vert_block
            Xb = np.ascontiguousarray(np.asarray(Xd[v0:v1])[:, :self.dim], dtype=float)
        self.periodic = np.asarray(periodic, dtype=bool)[:self.dim]
        self.x_min = np.full(self.dim, np.inf) if not len(Xb) else Xb.min(axis=0)
        self.x_max = np.full(self.dim, -np.inf) if not len(Xb) else Xb.max(axis=0)
        comm.Allreduce(MPI.IN_PLACE, self.x_min, op=MPI.MIN)
        comm.Allreduce(MPI.IN_PLACE, self.x_max, op=MPI.MAX)
        gid = np.arange(c0, c1, dtype=np.int64)
        if size > 1 and partition != "blocks":
            dest = _partition(gcells, gid, self.ncells_global, self.nverts_global, comm, partition)
            ex = Exchange(comm, dest)
            gcells, gid = ex.forward(gcells), ex.forward(gid)
            ci = ex.forward(ci) if ci is not None else None
            o = np.argsort(gid)
            gcells, gid = gcells[o], gid[o]
            ci = ci[o] if ci is not None else None
        self.cell_gid, self.cell_index = gid, ci
        self._build(gcells, Xb, periodic_tol)
        return self

    # ---- the tables, on the final cells
    def _build(self, gcells, Xb, tol):
        comm = self.comm
        self.ncells = len(gcells)
        self.vert_gid, inv = np.unique(gcells, return_inverse=True)
        self.cells = inv.reshape(gcells.shape).astype(np.int32)
        self.nverts = len(self.vert_gid)
        # the coordinates from the vertices' block owners
        v0 = self.vert_block[0]
        ex = Exchange(comm, block_owner(self.vert_gid, self.nverts_global, comm.size))
        self.X = ex.back(Xb[ex.forward(self.vert_gid) - v0])
        self.loc = np.array(LOC_EDGES[self.nv])
        self._edges()
        # the serial formulas, on the rank's cells
        geo = types.SimpleNamespace(X=self.X, cells=self.cells, ncells=self.ncells,
                                    nv=self.nv, dim=self.dim)
        D.Topo._geometry(geo)
        self.J, self.detJ, self.Jinv = geo.J, geo.detJ, geo.Jinv
        self.vol, self.facet_n = geo.vol, geo.facet_n
        self._facets(gcells, tol)
        self._owners()
        self._periodic(Xb, tol)
        self._exterior()
        self.vol_total = float(np.sum(comm.allgather(float(np.sum(self.vol)))))
        self.mesh_h = D.mesh_h(self)
        self.node_block = _split(self.nnodes_global, comm.size, comm.rank)
        self.set_held(None)

    def _edges(self):
        # local vertices are in serial order, so a sorted local pair is a sorted serial one
        pair = np.sort(self.cells[:, self.loc], axis=2).reshape(-1, 2)
        g = self.vert_gid[pair]
        key = g[:, 0] * np.int64(self.nverts_global) + g[:, 1]
        ids, self.nedges_global, _, _ = number_keys(key, self.comm)
        self.edge_gid, inv = np.unique(ids, return_inverse=True)
        self.nedges = len(self.edge_gid)
        self.cell_edge = inv.reshape(self.ncells, len(self.loc)).astype(np.int32)
        self.edges = np.zeros((self.nedges, 2), np.int32)
        self.edges[inv] = pair
        self.nnodes_global = self.nverts_global + self.nedges_global
        self.nnodes = self.nverts + self.nedges
        self.node_gid = np.concatenate([self.vert_gid, self.nverts_global + self.edge_gid])
        self.node_x = np.vstack([self.X,
                                 0.5 * (self.X[self.edges[:, 0]] + self.X[self.edges[:, 1]])])

    def _facets(self, gcells, tol):
        cnt, _ = facet_pass(gcells, self.cell_gid, self.nverts_global, self.comm)
        self.facet_exterior = cnt == 1
        F, _ = _facet_keys(self.cells, self.nverts)
        self.facet_verts = F.reshape(self.ncells, self.nv, self.nv - 1)
        per = np.zeros((self.ncells, self.nv), bool)
        for d in range(self.dim):
            if not self.periodic[d]:
                continue
            xd = self.X[:, d][self.facet_verts]
            for v in (self.x_min[d], self.x_max[d]):
                per |= np.all(np.abs(xd - v) < tol, axis=2)
        self.facet_periodic = per

    def cell_nodes(self):
        """The local nodes of every cell: its vertices, then its edges."""
        return np.concatenate([self.cells, self.nverts + self.cell_edge], axis=1)

    def _owners(self):
        # the cells ascend, so the last write of the reversed list is the smallest
        cn = self.cell_nodes()
        first = np.full(self.nnodes, self.ncells_global, np.int64)
        first[cn[::-1].ravel()] = np.repeat(self.cell_gid[::-1], cn.shape[1])
        size = self.comm.size
        code = reduce_keyed(self.node_gid, first * size + self.comm.rank, "min",
                            self.comm, self.nnodes_global)
        self.first_cell = code // size
        self.node_owner = (code % size).astype(np.int32)
        self._node_ex = Exchange(self.comm, self.node_owner)

    def reduce_nodes(self, values, op):
        """The op ('max', 'min', 'or', 'sum') over every rank's copy of each of
        the rank's nodes, at the node's owner, returned to every rank that
        touches it. With 'sum' each copy counts once."""
        return _reduce_exchange(self._node_ex, self.node_gid, np.asarray(values), op)

    def reduce_masters(self, values, op):
        """The op over every copy of every periodic image of each node's master,
        returned to each: a master class decides together. With 'sum' every
        rank's copy of every image counts once, so the callers pass each rank's
        partial sums, gathered at `class_index`, and they add up to the class's
        total as `Topo.reduce_masters` adds its one copy of each node."""
        return reduce_keyed(self.master, np.asarray(values), op, self.comm, self.nnodes_global)

    def median(self, values):
        """np.median of every rank's values together, exactly."""
        return exact_median(values, self.comm)

    def _periodic(self, Xb, tol):
        """The periodic faces' nodes, every rank's, paired as `Topo` pairs the
        whole mesh: no other node is in a pairing, so the masters are the
        serial ones."""
        comm = self.comm
        self.master = self.node_gid.copy()
        if not self.periodic.any():
            return

        def on_face(x):
            m = np.zeros(len(x), bool)
            for a in np.nonzero(self.periodic)[0]:
                m |= (x[:, a] < self.x_min[a] + tol) | (x[:, a] > self.x_max[a] - tol)
            return m

        # the vertices by their block, so one in no cell is paired too; the
        # midpoints by their owner
        vb = np.nonzero(on_face(Xb))[0]
        ex = self.node_x[self.nverts:]
        eo = np.nonzero(on_face(ex) & (self.node_owner[self.nverts:] == comm.rank))[0]
        mine_gid = np.concatenate([self.vert_block[0] + vb,
                                   self.nverts_global + self.edge_gid[eo]]).astype(np.int64)
        mine_x = np.vstack([Xb[vb], ex[eo]])
        sg = np.concatenate(comm.allgather(mine_gid))
        sx = np.vstack(comm.allgather(mine_x))
        o = np.argsort(sg)
        sg, sx = sg[o], sx[o]
        surf = types.SimpleNamespace(nnodes=len(sg), node_x=sx, periodic=self.periodic,
                                     dim=self.dim, x_min=self.x_min, x_max=self.x_max)
        D.Topo._periodic(surf, tol)
        at = np.searchsorted(sg, self.node_gid)
        at = np.minimum(at, max(len(sg) - 1, 0))
        hit = (sg[at] == self.node_gid) if len(sg) else np.zeros(self.nnodes, bool)
        self.master[hit] = sg[surf.master[at[hit]]]

    def _exterior(self):
        """The boundary the rule classifies: exterior facets the mesh does not
        pair. A vertex is on it if any rank's cell puts it there; an edge if any
        image of its master is."""
        ef = self.facet_exterior & ~self.facet_periodic
        self.facet_boundary = ef
        self.cell_boundary = ef.any(axis=1)
        node = np.zeros(self.nnodes, bool)
        for o in range(self.nv):
            loc = [j for j, (a, b) in enumerate(LOC_EDGES[self.nv]) if a != o and b != o]
            k = np.nonzero(ef[:, o])[0]
            if len(k):
                node[self.facet_verts[k, o].ravel()] = True
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
        held[lv:] = self.reduce_masters(held, "or")[lv:]
        self.held_node = held
        self.held_edge = held[lv:].copy()

    # ---- for the tests and the comparisons
    def gather(self, root=0, mesh_path=None):
        """The tables assembled on root in the serial order, with `Topo`'s
        names; None on the other ranks. Root then holds every table of the whole
        mesh, the Jacobians and normals too, so this is for the tests and small
        meshes; what a large one needs on one rank is gathered by itself. A
        vertex in no cell is taken from mesh_path's coordinates, as unpaired and
        off the boundary."""
        comm = self.comm
        lv = self.nverts
        cellpart = dict(
            gid=self.cell_gid, cells=self.vert_gid[self.cells],
            cell_edge=self.edge_gid[self.cell_edge],
            facet_verts=self.vert_gid[self.facet_verts],
            facet_exterior=self.facet_exterior, facet_periodic=self.facet_periodic,
            facet_boundary=self.facet_boundary, cell_boundary=self.cell_boundary,
            vol=self.vol, facet_n=self.facet_n, J=self.J, Jinv=self.Jinv, detJ=self.detJ)
        if self.cell_index is not None:
            cellpart["cell_index"] = self.cell_index
        own = self.node_owner == comm.rank
        nodepart = dict(gid=self.node_gid[own], master=self.master[own],
                        boundary_node=self.boundary_node[own], held_node=self.held_node[own],
                        at_rest=self.at_rest[own], node_x=self.node_x[own])
        eown = own[lv:]
        edgepart = dict(gid=self.edge_gid[eown], edges=self.vert_gid[self.edges[eown]])
        parts = comm.gather((cellpart, nodepart, edgepart), root=root)
        if comm.rank != root:
            return None
        out = types.SimpleNamespace(dim=self.dim, nv=self.nv, ncells=self.ncells_global,
                                    nverts=self.nverts_global, nedges=self.nedges_global,
                                    nnodes=self.nnodes_global, periodic=self.periodic,
                                    x_min=self.x_min, x_max=self.x_max)

        def place(i, n, names, dtypes):
            for name in names:
                a = [p[i][name] for p in parts]
                g = np.concatenate([p[i]["gid"] for p in parts])
                v = np.concatenate(a)
                full = np.zeros((n,) + v.shape[1:], dtypes.get(name, v.dtype))
                full[g] = v
                setattr(out, name, full)

        i32 = np.int32
        place(0, self.ncells_global, [k for k in cellpart if k != "gid"],
              dict(cells=i32, cell_edge=i32, facet_verts=i32))
        place(2, self.nedges_global, ["edges"], dict(edges=i32))
        n = self.nnodes_global
        g = np.concatenate([p[1]["gid"] for p in parts])
        out.master = np.arange(n, dtype=i32)
        out.master[g] = np.concatenate([p[1]["master"] for p in parts])
        for name in ("boundary_node", "held_node", "at_rest"):
            a = np.zeros(n, bool)
            a[g] = np.concatenate([p[1][name] for p in parts])
            setattr(out, name, a)
        out.node_x = np.zeros((n, self.dim))
        out.node_x[g] = np.vstack([p[1]["node_x"] for p in parts])
        missing = np.setdiff1d(np.arange(self.nverts_global), g[g < self.nverts_global])
        if len(missing):
            X, _, _ = D.read_mesh_h5(mesh_path)
            out.node_x[missing] = X[missing, :self.dim]
        out.X = out.node_x[:self.nverts_global]
        out.boundary_edge = out.boundary_node[self.nverts_global:]
        out.held_edge = out.held_node[self.nverts_global:]
        return out
