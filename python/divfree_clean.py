#!/usr/bin/env python3
"""A pointwise divergence-free velocity from a dolfin HDF5 dataset.

At a no-slip facet the normal derivative of the normal velocity is div u, so a
P1 or P2 field whose divergence is only weakly zero leaks into the wall and
traps tracers there. A field that is divergence-free *pointwise* in every cell
does not, and it needs no wall normal, corner rule or mean-normal compromise:

  * the data as P2 on the mesh, a P1 field's edge midpoints the mean of their
    ends;
  * the edge midpoints moved by the weighted change that makes every cell's net
    flux zero and leaves as little divergence in the cells as it can; a midpoint
    at rest on an exterior facet the mesh does not pair is held, a free one on
    such a facet costs BOUNDARY_W times an interior one, and a periodic image
    and its master are one unknown;
  * in each cell the unique continuous P2 field on its barycentric (Alfeld)
    split with that boundary data and div u = 0 pointwise.

The divergence the second step leaves behind is what the third has to cancel,
and the size of that cancellation is what makes the split field's gradients
rough, so the second step penalises it: gamma = g h^2 with g = 0.3 in 2D and 10
in 3D, `--penalty 0` for the smallest change alone. Zero net flux in every cell
is the only constraint: the volume mean of the final field, and with it the
throughput of every periodic direction, is reported before and after and never
held.

This tool does the first two and writes an ordinary P2 checkpoint case on the
same mesh. The split field is one fixed reference matrix and a Piola map a
cell, cheap enough for the loader to redo at every load, so the interior values
are not stored. `--split` writes it out as P2 on the barycentric split instead,
which the plain mesh loaders read as an ordinary P2 field.

Everything is read and written with h5py alone, so a dataset is cleaned without
dolfin, and the core functions take arrays and run without files. The solve is
PETSc's and runs on as many ranks as it is given; one rank without `mpirun` is
the ordinary tool.

  divfree_clean.py CASE/dolfin_params.dat --out DIR   cleaned P2, the same mesh
  mpirun -n 8 divfree_clean.py ... --out DIR          the same, solved on 8 ranks
  divfree_clean.py --check DIR/dolfin_params.dat      net fluxes of a dataset
"""

import argparse
import contextlib
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp

import petsc4py

petsc4py.init([])                # this tool's options are not PETSc's

from mpi4py import MPI           # noqa: E402
from petsc4py import PETSc       # noqa: E402

# dolfin's local edge numbering: edge j of a triangle faces vertex j; a tet's
# six in the order the loaders' local_edges<4>() gives them
LOC_EDGES = {3: [(1, 2), (0, 2), (0, 1)],
             4: [(2, 3), (1, 3), (1, 2), (0, 3), (0, 2), (0, 1)]}

REST_TOL = 1e-8          # |u| <= this times max |u| is at rest; --rest-tol
BOUNDARY_W = 1e4         # least-change weight of a free midpoint on an exterior facet
PERIODIC_TOL = 1e-12     # the dolfin HDF5 loaders' own; the XDMF tet meshes use 1e-4
MINRES_RTOL = 1e-12      # relative tolerance of the KKT solve, before refinement
FLUX_TOL = 1e-9          # net flux a loader refuses, relative to the cell's scale
FLUX_FLOOR = 1e-6        # a cell's scale is at least this of the stamp's largest facet flux
FLUX_EXIT = 1e-2         # already this far below the refusal: no step at all
FLUX_REFINE = 1e-3       # refinement stops here, a decade below the early exit
REFINE_MAX = 4           # refinement passes on the true residual of a step
REFINE_GAIN = 3.0        # a pass that gains less than this factor is the last
RESIDUAL_FLOOR = 1e-14   # a residual this far below the right-hand side is round-off
RESIDUAL_STALL = 1e-16   # a solve asks no less of itself than this of its first right-hand side
MAXITER = 20000          # iterations of one MINRES solve
REFINE_MAXITER = MAXITER // 10   # iterations of one refinement pass
SCHUR_SHIFT = 1e-8       # diagonal shift of the Schur preconditioner, of its mean diagonal
DEFAULT_G = {2: 0.3, 3: 10.0}    # gamma = g h^2, the divergence penalty
READ_CACHE_BYTES = 1 << 31   # stamp values kept from the wall pass for the cleaning pass
READ_CHUNK_VALUES = 1 << 18  # values gathered at once when a field is scattered by cell


# ------------------------------------------------------------------ topology


class Topo:
    """Mesh arrays, the edge and facet tables, the periodic masters and the held
    set. Nodes are the vertices, then the edges."""

    def __init__(self, X, cells, periodic, at_rest=None, periodic_tol=PERIODIC_TOL):
        self.X = np.ascontiguousarray(X, dtype=float)
        self.cells = np.ascontiguousarray(cells, dtype=np.int32)
        self.dim = self.X.shape[1]
        self.nv = self.cells.shape[1]
        assert self.nv == self.dim + 1
        self.periodic = np.asarray(periodic, dtype=bool)[:self.dim]
        self.nverts = len(self.X)
        self.ncells = len(self.cells)
        self.loc = np.array(LOC_EDGES[self.nv])
        self.x_min = self.X.min(axis=0)
        self.x_max = self.X.max(axis=0)
        self._edges()
        self._geometry()
        self._facets(periodic_tol)
        self._periodic(periodic_tol)
        self._exterior()
        self.set_held(at_rest)

    def set_held(self, at_rest):
        """The held set, which the dataset's own stamps decide; everything above
        it is the mesh alone, so it is set on a topology already built."""
        self._held(np.zeros(self.nnodes, bool) if at_rest is None
                   else np.asarray(at_rest, bool))

    # ---- edges
    def _edges(self):
        pair = np.sort(self.cells[:, self.loc], axis=2).reshape(-1, 2)
        key = pair[:, 0].astype(np.int64) * np.int64(self.nverts) + pair[:, 1]
        uniq, inv = np.unique(key, return_inverse=True)
        self.nedges = len(uniq)
        self.cell_edge = inv.reshape(self.ncells, len(self.loc)).astype(np.int32)
        self.edges = np.stack([uniq // self.nverts, uniq % self.nverts],
                              axis=1).astype(np.int32)
        self.nnodes = self.nverts + self.nedges
        self.node_x = np.vstack([self.X,
                                 0.5 * (self.X[self.edges[:, 0]] + self.X[self.edges[:, 1]])])

    # ---- volumes and outward facet normals (|n| = length or area)
    def _geometry(self):
        P = self.X[self.cells]                                   # (nc, nv, d)
        self.J = np.transpose(P[:, 1:] - P[:, :1], (0, 2, 1))    # (nc, d, d)
        self.detJ = np.linalg.det(self.J)
        self.Jinv = np.linalg.inv(self.J)
        self.vol = np.abs(self.detJ) / (2.0 if self.dim == 2 else 6.0)
        n = np.zeros((self.ncells, self.nv, self.dim))
        for o in range(self.nv):
            f = [i for i in range(self.nv) if i != o]
            if self.dim == 2:
                t = P[:, f[1]] - P[:, f[0]]
                nn = np.stack([t[:, 1], -t[:, 0]], axis=1)
            else:
                nn = np.cross(P[:, f[1]] - P[:, f[0]], P[:, f[2]] - P[:, f[0]]) / 2.0
            s = np.einsum('ij,ij->i', nn, P[:, o] - P[:, f[0]])
            n[:, o] = np.where(s[:, None] > 0, -nn, nn)
        self.facet_n = n

    # ---- facets: exterior, periodic
    def _facets(self, tol):
        rows = [np.sort(self.cells[:, [i for i in range(self.nv) if i != o]], axis=1)
                for o in range(self.nv)]
        F = np.stack(rows, axis=1).reshape(-1, self.nv - 1)      # (nc*nv, d)
        key = np.zeros(len(F), np.int64)
        for c in range(self.nv - 1):
            key = key * np.int64(self.nverts) + F[:, c].astype(np.int64)
        _, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
        self.facet_exterior = (cnt[inv] == 1).reshape(self.ncells, self.nv)
        self.facet_verts = F.reshape(self.ncells, self.nv, self.nv - 1)
        per = np.zeros((self.ncells, self.nv), bool)
        for d in range(self.dim):
            if not self.periodic[d]:
                continue
            # One coordinate of the facets' vertices, not all of them
            xd = self.X[:, d][self.facet_verts]
            for v in (self.x_min[d], self.x_max[d]):
                per |= np.all(np.abs(xd - v) < tol, axis=2)
        self.facet_periodic = per

    # ---- periodic node masters, paired as the loaders pair them: the min face
    # is the master and a corner reaches its own through the chain of links
    def _periodic(self, tol):
        n = self.nnodes
        link = np.arange(n, dtype=np.int32)
        if self.periodic.any():
            from scipy.spatial import cKDTree
            x = self.node_x
            for a in range(self.dim):
                if not self.periodic[a]:
                    continue
                lo = np.nonzero(x[:, a] < self.x_min[a] + tol)[0]
                hi = np.nonzero((x[:, a] > self.x_max[a] - tol) & (link == np.arange(n)))[0]
                if len(lo) == 0 or len(hi) == 0:
                    continue
                q = x[hi].copy()
                q[:, a] -= self.x_max[a] - self.x_min[a]
                d, j = cKDTree(x[lo]).query(q)
                ok = d < tol
                link[hi[ok]] = lo[j[ok]]
        master = np.arange(n, dtype=np.int32)
        for _ in range(self.dim + 1):
            master = link[master]
        self.master = master

    # ---- the boundary the rule classifies: exterior facets the mesh does not pair
    def _exterior(self):
        ef = self.facet_exterior & ~self.facet_periodic
        self.facet_boundary = ef
        self.cell_boundary = ef.any(axis=1)
        node = np.zeros(self.nnodes, bool)
        # facet o of a cell misses vertex o, so it holds the local edges whose
        # pair avoids o
        for o in range(self.nv):
            loc = [j for j, (a, b) in enumerate(LOC_EDGES[self.nv]) if a != o and b != o]
            k = np.nonzero(ef[:, o])[0]
            if len(k):
                node[self.facet_verts[k, o].ravel()] = True
                node[self.nverts + self.cell_edge[np.ix_(k, loc)].ravel()] = True
        m = self.master[self.nverts:] - self.nverts
        on = np.zeros(self.nedges, bool)
        on[m[node[self.nverts:]]] = True
        self.boundary_edge = on[m]
        node[self.nverts:] = self.boundary_edge
        self.boundary_node = node

    # ---- held
    def _held(self, at_rest):
        self.at_rest = at_rest
        self.held_node = self.boundary_node & at_rest
        # a master is held if any of its images is
        m = self.master[self.nverts:] - self.nverts
        fixed = np.zeros(self.nedges, bool)
        fixed[m[self.held_node[self.nverts:]]] = True
        self.held_edge = fixed[m]
        self.held_node[self.nverts:] = self.held_edge


def at_rest_nodes(nmax, scale, rest_tol=None):
    """Nodes at rest in every stamp or component: nmax is the largest |u| each
    node ever holds, scale the largest the dataset holds anywhere. The threshold
    sits above an iterative solver's round-off at its no-slip nodes."""
    return nmax <= (REST_TOL if rest_tol is None else rest_tol) * scale


# ------------------------------------------------------------- equilibration


def cell_flux(topo, U):
    """Net P2 flux out of every cell, and the largest facet flux of each. 2D
    facet e: |e|(u_a + 4 u_m + u_b)/6 . n; 3D facet f: |f|/3 (sum of its three
    edge midpoints) . n."""
    Um = U[topo.nverts:]
    r = np.zeros(topo.ncells)
    big = np.zeros(topo.ncells)
    for o in range(topo.nv):
        n = topo.facet_n[:, o]
        f = [i for i in range(topo.nv) if i != o]
        if topo.dim == 2:
            e = topo.cell_edge[:, o]
            v = (U[topo.cells[:, f[0]]] + 4 * Um[e] + U[topo.cells[:, f[1]]]) / 6.0
        else:
            loc = [j for j, (a, b) in enumerate(LOC_EDGES[4]) if a != o and b != o]
            v = Um[topo.cell_edge[:, loc]].sum(axis=1) / 3.0
        q = np.einsum('ij,ij->i', v, n)
        r += q
        big = np.maximum(big, np.abs(q))
    return r, big


def flux_ratios(topo, U, parts=None):
    """(the balance of every cell, the stamp's largest facet flux). A cell is
    balanced when its ratio is at most FLUX_TOL, and the scale it is measured
    against is its own largest facet flux, floored at FLUX_FLOOR of the stamp's:
    a cell whose own fluxes are round-off -- a dead-end pore where the solver's
    velocity is at its noise floor -- would otherwise be judged by a ratio of
    round-off to round-off. The loader applies the same criterion, so what the
    tool accepts is what loads. `parts` is `cell_flux`'s answer for this field
    where the caller already has it."""
    r, big = cell_flux(topo, U) if parts is None else parts
    scale = np.maximum(big, FLUX_FLOOR * max(big.max(), 1e-300))
    return np.abs(r) / scale, float(big.max())


def free_unknowns(topo):
    """(free masters, edge -> free column). An unknown is numbered by the first
    cell it appears in, so a block partition of the cells cuts the unknowns into
    contiguous blocks too and almost every entry a rank assembles is a row it
    owns."""
    m = topo.master[topo.nverts:] - topo.nverts
    free = np.nonzero((m == np.arange(topo.nedges)) & ~topo.held_edge)[0]
    first = np.full(topo.nedges, topo.ncells, np.int64)
    # the cells backwards, so the smallest is the write that stands
    first[m[topo.cell_edge[::-1]].ravel()] = np.repeat(
        np.arange(topo.ncells - 1, -1, -1), topo.cell_edge.shape[1])
    free = free[np.argsort(first[free], kind="stable")]
    col_of = np.full(topo.nedges, -1, np.int32)
    col_of[free] = np.arange(len(free))
    return free, col_of[m]


def flux_matrix(topo, cells=None, unknowns=None):
    """(C, free masters, edge -> free column) with C the net flux of the cells in
    `cells` -- all of them by default -- in the change of the free edge
    midpoints; its rows are those cells in order and its columns the whole mesh's
    nfree*d unknowns."""
    free_master, col = free_unknowns(topo) if unknowns is None else unknowns
    c0, c1 = (0, topo.ncells) if cells is None else cells
    sl = slice(c0, c1)
    d = topo.dim
    rows, cols, vals = [], [], []
    for o in range(topo.nv):
        n = topo.facet_n[sl, o]
        if topo.dim == 2:
            locs, w = [o], 4.0 / 6.0
        else:
            locs = [j for j, (a, b) in enumerate(LOC_EDGES[4]) if a != o and b != o]
            w = 1.0 / 3.0
        for j in locs:
            c = col[topo.cell_edge[sl, j]]
            k = np.nonzero(c >= 0)[0]
            for q in range(d):
                rows.append(k)
                cols.append(d * c[k] + q)
                vals.append(w * n[k, q])
    C = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(c1 - c0, d * len(free_master)))
    C.sum_duplicates()
    return C, free_master, col


def _row_counts(A, r0, r1, c0, c1):
    """Entries of rows r0..r1 of a CSR whose column is in [c0, c1). reduceat over
    the mask, which is a byte a nonzero, rather than a row index an entry."""
    a, b = int(A.indptr[r0]), int(A.indptr[r1])
    if b <= a:                       # a block with no entry at all: reduceat has none
        return np.zeros(r1 - r0, np.int32)
    mask = ((A.indices[a:b] >= c0) & (A.indices[a:b] < c1)).view(np.int8)
    empty = np.diff(A.indptr[r0:r1 + 1]) == 0
    out = np.add.reduceat(mask, np.minimum(A.indptr[r0:r1] - a, len(mask) - 1))
    out = out.astype(np.int32)
    # reduceat reads one entry where a row is empty
    out[empty] = 0
    return out


def _split(total, size, rank):
    """PETSc's own division of a global size into contiguous blocks."""
    q, r = divmod(total, size)
    start = rank * q + min(rank, r)
    return start, start + q + (1 if rank < r else 0)


_KSP_TAG = [0]


def _petsc_csr(PETSc, A):
    """A's CSR arrays with the indices in PETSc's integer type, which the cast
    would wrap silently past its range: refused there. PETSc's AIJ kernels assume
    sorted columns, so the matrix is sorted first."""
    if not A.has_sorted_indices:
        A.sort_indices()
    top = int(np.iinfo(PETSc.IntType).max)
    if not max(A.nnz, A.shape[0], A.shape[1]) <= top:
        raise ValueError(
            "a %d x %d matrix with %d nonzeros is past the %d that this PETSc's %d-bit "
            "indices hold; a PETSc built with 64-bit indices reaches it"
            % (A.shape[0], A.shape[1], A.nnz, top, 8 * np.dtype(PETSc.IntType).itemsize))
    return (A.indptr.astype(PETSc.IntType, copy=False),
            A.indices.astype(PETSc.IntType, copy=False), A.data)


def mesh_h(topo):
    """The mesh scale the penalty is measured in: the median edge length."""
    e = topo.X[topo.edges[:, 0]] - topo.X[topo.edges[:, 1]]
    return float(np.median(np.linalg.norm(e, axis=1)))


# ------------------------------------------------- the divergence of the data


def p2_vertex_grads(d):
    """Gref[a, n, b] = d phi_n / d xhat_b at vertex a of the reference simplex,
    n over the P2 nodes in the canonical (vertices, LOC_EDGES) order."""
    V = np.vstack([np.zeros(d), np.eye(d)])
    return np.array([_p2_basis(V, V[a], grad=True)[1] for a in range(d + 1)])


def mass_p1(d):
    """Mhat with int_K f g = vol * f^T Mhat g for P1 f, g at the vertices."""
    return (np.eye(d + 1) + np.ones((d + 1, d + 1))) / ((d + 1) * (d + 2))


def cell_div(topo, U):
    """div u at the d+1 vertices of every cell, (ncells, d+1); div of a P2 field
    is P1 in a cell, so these values are the whole of it."""
    g = np.concatenate([U[topo.cells], U[topo.nverts + topo.cell_edge]], axis=1)
    return np.einsum('anb,kbc,knc->ka', p2_vertex_grads(topo.dim), topo.Jinv, g,
                     optimize=True)


def div_norms(topo, U):
    """||div u||_{L2(K)} for every cell."""
    dv = cell_div(topo, U)
    return np.sqrt(np.maximum(np.einsum('ka,ab,kb->k', dv, mass_p1(topo.dim), dv), 0.0)
                   * topo.vol)


def div_edge_blocks(topo, col, sl=slice(None)):
    """(Be, gcol) for a slice of the cells: Be[k, a, (j, c)] is the coefficient
    of the change of local edge midpoint j, component c, in div u at vertex a of
    cell k, and gcol its column in the unknowns (-1 for a fixed midpoint). They
    are the divergence operator, so they are built a chunk at a time and thrown
    away rather than kept for the mesh."""
    d, ne = topo.dim, topo.cell_edge.shape[1]
    Gref = p2_vertex_grads(d)[:, topo.nv:, :]
    Jinv = topo.Jinv[sl]
    nk = len(Jinv)
    Be = np.einsum('anb,kbc->kanc', Gref, Jinv).reshape(nk, d + 1, ne * d)
    c = col[topo.cell_edge[sl]]
    gcol = np.where(c[:, :, None] >= 0, c[:, :, None] * d + np.arange(d)[None, None, :], -1)
    return Be, gcol.reshape(nk, ne * d)


def _cell_chunks(topo, nb, chunk=1 << 22, cells=None):
    """Slices of a range of the cells small enough that a chunk's dense blocks
    fit."""
    c0, c1 = (0, topo.ncells) if cells is None else cells
    per = max(1, chunk // max(nb * nb, 1))
    return [slice(k0, min(k0 + per, c1)) for k0 in range(c0, c1, per)]


def _merge_part(parts, P):
    """Push a chunk's matrix onto a stack that keeps each entry at least twice
    the size of the one above it, so only log(chunks) partial sums are alive at
    once."""
    while parts and parts[-1].nnz <= 2 * P.nnz:
        P = parts.pop() + P
    parts.append(P)


def stiffness_matrix(topo, col, n, gamma, wd, cells=None):
    """W + gamma Q over a range of the cells, with s^T Q s = sum_K
    ||div s||^2_{L2(K)}: the cell blocks with the cell's P1 mass matrix, plus the
    diagonal weights, which are wd where the caller carries them and zero
    elsewhere. The cells are taken in chunks and each chunk merged as it is
    built, so neither the whole entry list nor the chunks themselves are ever
    alive together."""
    if not gamma:
        return sp.diags(wd).tocsr()
    M = mass_p1(topo.dim)
    nb = topo.cell_edge.shape[1] * topo.dim
    parts = [sp.diags(wd).tocsr()]
    for sl in _cell_chunks(topo, nb, cells=cells):
        Be, gc = div_edge_blocks(topo, col, sl)
        K = gamma * np.einsum('k,kan,ab,kbm->knm', topo.vol[sl], Be, M, Be, optimize=True)
        del Be
        rows = np.repeat(gc[:, :, None], nb, axis=2)
        cols = np.repeat(gc[:, None, :], nb, axis=1)
        ok = (rows >= 0) & (cols >= 0)
        part = sp.coo_matrix((K[ok], (rows[ok], cols[ok])), shape=(n, n)).tocsr()
        del K, rows, cols, ok, gc
        _merge_part(parts, part)
    while len(parts) > 1:
        parts.append(parts.pop() + parts.pop())
    K = parts[0]
    K.sum_duplicates()
    return K


def penalty_rhs(topo, U, col, n, cells=None):
    """q of the cross term 2 gamma q^T s over a range of the cells:
    q = sum_K B_K^T (vol_K Mhat) div u_K."""
    M = mass_p1(topo.dim)
    Gref = p2_vertex_grads(topo.dim)
    out = np.zeros(n)
    for sl in _cell_chunks(topo, topo.cell_edge.shape[1] * topo.dim, cells=cells):
        Be, gc = div_edge_blocks(topo, col, sl)
        g = np.concatenate([U[topo.cells[sl]], U[topo.nverts + topo.cell_edge[sl]]], axis=1)
        dv = np.einsum('anb,kbc,knc->ka', Gref, topo.Jinv[sl], g, optimize=True)
        v = np.einsum('k,kan,ab,kb->kn', topo.vol[sl], Be, M, dv, optimize=True)
        ok = gc >= 0
        out += np.bincount(gc[ok], weights=v[ok], minlength=n)
    return out


# ------------------------------------------------------------- the volume mean


def _p2_int_weights(d):
    """(vertex, midpoint) weights with int_K f = vol * (a_v sum_verts f +
    a_m sum_mids f) for P2 f."""
    return (0.0, 1.0 / 3.0) if d == 2 else (-1.0 / 20.0, 1.0 / 5.0)


def volume_mean(topo, U):
    """The volume mean of the P2 field on the macro cells."""
    a_v, a_m = _p2_int_weights(topo.dim)
    s = (a_v * U[topo.cells].sum(axis=1)
         + a_m * U[topo.nverts + topo.cell_edge].sum(axis=1))
    return (topo.vol @ s) / topo.vol.sum()


def div_moment(topo, U):
    """int x_c div u over the domain, per component. The split field is
    divergence-free in every cell and takes the macro field's trace on the
    boundary, so its volume mean is the macro field's plus this over the volume:
    the moment is what the final field's mean gains over the data's, and it
    vanishes for data that satisfies the discrete continuity equation against the
    P1 pressures."""
    M = mass_p1(topo.dim)
    return np.einsum('k,ab,kac,kb->c', topo.vol, M, topo.X[topo.cells], cell_div(topo, U),
                     optimize=True)


def split_volume_mean(topo, U):
    """The volume mean of the final split field of macro data with zero net cell
    fluxes."""
    return volume_mean(topo, U) + div_moment(topo, U) / topo.vol.sum()


class Equil:
    """The global step of one mesh: assembled and preconditioned once, solved
    for every stamp or component.

    The change s of the free edge midpoints minimises

        s^T W s + gamma sum_K ||div(u + s)||^2_{L2(K)},   gamma = g h^2,

    subject to zero net flux through every cell, and nothing else. A midpoint at
    rest on an exterior facet the mesh does not pair is held and has no unknown;
    a free midpoint on such a facet costs boundary_weight times an interior one
    in W, which is what keeps a moving wall's values and an open facet's flux
    near the data's. The penalty term buys a smoother field: what the per-cell
    reconstruction has to cancel is what makes its correction, and so the split
    field's gradients, rough.

    Whatever gamma is, the step is the one KKT system

        [[W + gamma Q, C^T], [C, 0]] [s; lambda] = [-gamma q; rhs],

    solved by PETSc's MINRES with an additive field split: a GAMG V-cycle on the
    SPD block K = W + gamma Q, and one on the Schur approximation
    C diag(K)^-1 C^T, shifted so that it is definite whatever the rank of the
    constraint rows. With gamma = 0 the block is diagonal and Jacobi is exact on
    it, and that approximation is then the Schur complement itself. Both blocks
    are fixed symmetric positive operators -- a V-cycle with the same pre- and
    post-smoothing is its own adjoint -- so MINRES stays valid.

    The mesh tables are every rank's; the matrices and the solve are shared out.
    The cells go to the ranks in contiguous blocks and the unknowns, numbered by
    the first cell they appear in, in blocks of their own, so each rank's rows of
    both blocks are contiguous and nearly every entry it assembles is a row it
    owns. Both shares are even: numbering an unknown by its first cell alone
    would give the low ranks most of them, since an early cell claims all six of
    its edges and a late one claims none.
    Everything measured by cell -- the fluxes, the refinement's criterion -- each
    rank computes from the field it holds, so nothing is gathered but the step.

    Where every boundary midpoint is held the constant over the cells is in the
    null space of the flux rows: the right-hand side is projected onto its
    complement, and what the data's own net flux leaves is spread beforehand.

    The step is refined against the true residual, since an iterative solver
    stops on its own estimate of it, until every cell's balance is FLUX_REFINE
    of what a loader accepts or a pass stops paying. The first solve is refused
    if it does not converge, naming the stamp; a refinement pass that does not
    is a pass that does not pay, and the step before it stands.
    """

    def __init__(self, topo, weights="volume", penalty=None, minres_rtol=None,
                 boundary_weight=None, comm=None):
        self.topo = topo
        self.minres_rtol = MINRES_RTOL if minres_rtol is None else minres_rtol
        self.h = mesh_h(topo)
        self.g = DEFAULT_G[topo.dim] if penalty is None else float(penalty)
        self.gamma = self.g * self.h ** 2
        self.comm = PETSc.COMM_WORLD if comm is None else comm
        self.mpi = self.comm.tompi4py()
        d = topo.dim
        free, col = free_unknowns(topo)
        self.free, self.col = free, col
        self.n = len(free) * d
        self.m = topo.ncells
        if self.n == 0:
            raise ValueError("every midpoint of the %d cells is held, so there is nothing "
                             "to solve for; --rest-tol says what counts as at rest"
                             % topo.ncells)
        self.cells = _split(self.m, self.mpi.size, self.mpi.rank)
        # the unknowns are in the order of their first cells, so an even share of
        # them is contiguous and is nearly the share the rank's own cells reach
        self.rows = tuple(d * k for k in _split(len(free), self.mpi.size, self.mpi.rank))
        self.nloc = self.rows[1] - self.rows[0]
        self.mloc = self.cells[1] - self.cells[0]
        # PETSc will not assemble a block with no rows, and the job would hang
        # there; the blocks are even, so every rank refuses on the same test
        if self.m < self.mpi.size or len(free) < self.mpi.size:
            raise ValueError("%d ranks is more than this mesh's %d cells and %d free "
                             "midpoints share out: a rank with no row of its own does not "
                             "assemble, so run it on at most %d"
                             % (self.mpi.size, self.m, len(free),
                                min(self.m, len(free))))
        self.ucounts = np.array(self.mpi.allgather(self.nloc))
        self.udispls = np.concatenate([[0], np.cumsum(self.ucounts)[:-1]])
        self.ubounds = np.concatenate([[0], np.cumsum(self.ucounts)])
        if weights == "volume":
            w = np.bincount(
                (topo.master[topo.nverts + topo.cell_edge] - topo.nverts).ravel(),
                weights=np.repeat(topo.vol, topo.cell_edge.shape[1]),
                minlength=topo.nedges)
            w = np.maximum(w[self.free], 1e-300)
        elif weights == "none":
            w = np.ones(len(self.free))
        else:
            raise ValueError("unknown weighting: " + weights)
        self.boundary_weight = (BOUNDARY_W if boundary_weight is None
                                else float(boundary_weight))
        # zero empties those rows of K, and the preconditioner is built on its diagonal
        if not self.boundary_weight > 0.0:
            raise ValueError("the boundary weight is a cost relative to an interior "
                             "midpoint's and must be positive, not %g"
                             % self.boundary_weight)
        # a periodic seam is not a boundary, so its midpoints keep the ordinary cost
        w = w * np.where(topo.boundary_edge[self.free], self.boundary_weight, 1.0)
        self.wd = np.repeat(w, d)
        self.iterations = 0
        self.refinements = 0
        t0 = time.time()
        Cl = flux_matrix(topo, cells=self.cells, unknowns=(free, col))[0]
        one = np.ones(self.mloc)
        sums = np.vstack([Cl.T @ one, np.abs(Cl).T @ one])
        self.mpi.Allreduce(MPI.IN_PLACE, sums, op=MPI.SUM)
        self.singular = (float(np.abs(sums[0]).max())
                         < 1e-10 * max(float(sums[1].max()), 1e-300))
        del sums, one
        self._setup(Cl)
        self.setup_time = time.time() - t0

    # ---- the distributed matrices
    def _counts(self, A, rowb, colb):
        """(nonzeros in the row owner's own columns, in the rest) of every row of
        a matrix whose rows are the whole system's. Which columns are the
        diagonal block's is the owner's business, not the assembling rank's, so
        the rows are counted a rank's block at a time. The one step of the layer
        that is O(the whole system) on every rank: one allreduce a matrix."""
        cnt = np.zeros((2, A.shape[0]), np.int32)
        per_row = np.diff(A.indptr)
        for r in range(len(rowb) - 1):
            r0, r1 = int(rowb[r]), int(rowb[r + 1])
            if r1 <= r0 or A.indptr[r1] <= A.indptr[r0]:
                continue
            inside = _row_counts(A, r0, r1, colb[r], colb[r + 1])
            cnt[0, r0:r1] = inside
            cnt[1, r0:r1] = per_row[r0:r1] - inside
        self.mpi.Allreduce(MPI.IN_PLACE, cnt, op=MPI.SUM)
        return cnt

    def _mpiaij(self, A, rowb, colb):
        """An MPIAIJ from every rank's own contribution to a scipy matrix whose
        rows are the whole system's and which is empty outside the rank's share.
        The preallocation is the ranks' counts added up, which is exact but for a
        row two ranks both touch: one on the partition's boundary is
        preallocated a little wide."""
        N, M = A.shape
        r0, rloc = int(rowb[self.mpi.rank]), int(rowb[self.mpi.rank + 1] - rowb[self.mpi.rank])
        cloc = int(colb[self.mpi.rank + 1] - colb[self.mpi.rank])
        cnt = self._counts(A, rowb, colb)[:, r0:r0 + rloc]
        P = PETSc.Mat().createAIJ(
            size=((rloc, N), (cloc, M)),
            nnz=(np.minimum(cnt[0], cloc).astype(PETSc.IntType),
                 np.minimum(cnt[1], M - cloc).astype(PETSc.IntType)),
            comm=self.comm)
        # the rows this rank touches, its own and its neighbours' alike; the
        # matrix is handed over whole rather than sliced, which would copy it
        indptr, indices, data = _petsc_csr(PETSc, A)
        P.setValuesIJV(indptr, indices, data, addv=PETSc.InsertMode.ADD_VALUES,
                       rowmap=np.arange(N, dtype=PETSc.IntType))
        P.assemble()
        return P

    def _mpiaij_rows(self, A, rloc, ntotal, colb):
        """An MPIAIJ from a matrix that is the rank's own rows and nothing else,
        so the count is exact and no entry leaves the rank."""
        M = A.shape[1]
        c0, c1 = int(colb[self.mpi.rank]), int(colb[self.mpi.rank + 1])
        inside = _row_counts(A, 0, rloc, c0, c1)
        P = PETSc.Mat().createAIJ(
            size=((rloc, ntotal), (c1 - c0, M)),
            nnz=(inside.astype(PETSc.IntType),
                 (np.diff(A.indptr) - inside).astype(PETSc.IntType)),
            comm=self.comm)
        P.setValuesCSR(*_petsc_csr(PETSc, A), addv=PETSc.InsertMode.ADD_VALUES)
        P.assemble()
        return P

    def _setup(self, Cl):
        n, m, d = self.n, self.m, self.topo.dim
        u0 = self.rows[0]
        diag = np.zeros(n)
        diag[u0:u0 + self.nloc] = self.wd[u0:u0 + self.nloc]
        Kl = stiffness_matrix(self.topo, self.col, n, self.gamma, diag, cells=self.cells)
        K = self._mpiaij(Kl, self.ubounds, self.ubounds)
        del Kl, diag
        K.setOption(PETSc.Mat.Option.SYMMETRIC, True)
        # a node's d components are one unknown to the aggregation
        K.setBlockSize(d)
        C = self._mpiaij_rows(Cl, self.mloc, m, self.ubounds)
        del Cl
        Ct = C.transpose(PETSc.Mat())
        dg = K.getDiagonal()
        dg.abs()
        dg.reciprocal()
        B = Ct.copy()
        B.diagonalScale(L=dg)
        S = C.matMult(B)
        B.destroy()
        # definite whatever the rank of the constraint rows
        S.shift(SCHUR_SHIFT * S.getDiagonal().sum() / max(m, 1))
        S.setOption(PETSc.Mat.Option.SYMMETRIC, True)
        self.A = PETSc.Mat().createNest([[K, Ct], [C, None]], comm=self.comm)
        self.A.assemble()
        self.P = PETSc.Mat().createNest([[K, None], [None, S]], comm=self.comm)
        self.P.assemble()
        self.blocks = (K, C, Ct, S)
        _KSP_TAG[0] += 1
        self.prefix = prefix = "divfree%d_" % _KSP_TAG[0]
        opts = PETSc.Options()
        for name, kind in (("u", "gamg" if self.gamma else "jacobi"), ("s", "gamg")):
            opts[prefix + "fieldsplit_%s_ksp_type" % name] = "preonly"
            opts[prefix + "fieldsplit_%s_pc_type" % name] = kind
            if kind == "gamg":
                # plain aggregation, not smoothed: smoothing fills the coarse operators
                opts[prefix + "fieldsplit_%s_pc_gamg_agg_nsmooths" % name] = 0
        self.ksp = PETSc.KSP().create(comm=self.comm)
        self.ksp.setOptionsPrefix(prefix)
        self.ksp.setOperators(self.A, self.P)
        self.ksp.setType("minres")
        self.ksp.setTolerances(rtol=self.minres_rtol, max_it=MAXITER)
        pc = self.ksp.getPC()
        pc.setType("fieldsplit")
        pc.setFieldSplitIS(*zip(("u", "s"), self.A.getNestISs()[0]))
        pc.setFieldSplitType(PETSc.PC.CompositeType.ADDITIVE)
        self.ksp.setFromOptions()
        self.ksp.setUp()
        self._b, self._x, self._r, self._dx, self._w = [self._vec() for _ in range(5)]
        self.solver = ("MINRES + GAMG blocks" if self.gamma
                       else "MINRES + GAMG on the flux block") \
            + (", %d ranks" % self.mpi.size if self.mpi.size > 1 else "")

    def destroy(self):
        """The PETSc objects and the options database entries of this instance.
        A caller that cleans several cases in one process would otherwise keep
        every KSP, nest and GAMG hierarchy alive."""
        for v in (self._b, self._x, self._r, self._dx, self._w):
            v.destroy()
        self.ksp.destroy()
        self.A.destroy()
        self.P.destroy()
        for b in self.blocks:
            b.destroy()
        opts = PETSc.Options()
        for name in list(opts.getAll()):
            if name.startswith(self.prefix):
                opts.delValue(name)

    def _vec(self):
        """A vector of the KKT system: the rank's unknowns, then its cells."""
        v = PETSc.Vec().create(comm=self.comm)
        v.setSizes((self.nloc + self.mloc, self.n + self.m))
        v.setUp()
        return v

    def _gather(self, x):
        """The unknown block of a KKT vector, whole, on every rank."""
        out = np.empty(self.n)
        with x as a:
            self.mpi.Allgatherv(np.ascontiguousarray(a[:self.nloc]),
                                [out, self.ucounts, self.udispls, MPI.DOUBLE])
        return out

    def _consistent(self, v):
        """A right-hand side orthogonal to the constant over the cells, which the
        closed domain's flux rows leave in their null space. What the data's own
        net fluxes leave is spread beforehand; this is round-off."""
        if not self.singular:
            return
        with v as a:
            s = float(a[self.nloc:].sum())
        s = self.mpi.allreduce(s) / max(self.m, 1)
        with v as a:
            a[self.nloc:] -= s

    def _pnorm(self, v):
        """The norm MINRES measures its residual in, sqrt(v^T M v): what an
        absolute tolerance has to be given in."""
        self.ksp.getPC().apply(v, self._w)
        return float(np.sqrt(max(v.dot(self._w), 0.0)))

    def _minres(self, b, x, who, atol=0.0, max_it=None):
        """One preconditioned MINRES solve of the KKT system, to minres_rtol or
        to atol where that is larger; who names the stamp in a failure. A
        correction whose right-hand side is a true residual at round-off would
        otherwise be chased below what the first solve could reach."""
        x.set(0.0)
        self.ksp.setTolerances(rtol=self.minres_rtol, atol=max(atol, 1e-50),
                               max_it=MAXITER if max_it is None else max_it)
        self.ksp.solve(b, x)
        self._its += self.ksp.getIterationNumber()
        reason = self.ksp.getConvergedReason()
        if reason < 0:
            raise RuntimeError("%sthe KKT solve did not converge: PETSc's MINRES stopped "
                               "with reason %d after %d iterations"
                               % (who, reason, self.ksp.getIterationNumber()))

    def _spread_net(self, net, scale):
        """The net flux a closed domain's cell rows cannot take away, shared out
        over the cells in proportion to their scale, so that every cell is left
        the same fraction of what a loader measures it against."""
        if not self.singular:
            return np.zeros(self.topo.ncells)
        return net * scale / scale.sum()

    def _step(self, U, rhs_c, who, facet):
        """(the step, the refinement passes it took, why they stopped). Passes on
        the true residual until every cell's balance is FLUX_REFINE below what a
        loader accepts, or a pass stops paying: it gains less than REFINE_GAIN,
        the residual is at round-off, or the pass runs out of its own iterations.

        The absolute tolerance every solve carries is the criterion's own: the
        smallest net flux it asks of any cell, FLUX_REFINE * FLUX_TOL of the
        smallest scale a cell is measured against, which is FLUX_FLOOR of the
        stamp's largest facet flux. Taken from the first right-hand side
        instead, it would pin every later correction -- whose own right-hand
        side is minres_rtol of the first -- a couple of decades below where it
        started. Nothing below RESIDUAL_STALL of the first right-hand side is
        asked for, since MINRES stops making progress there."""
        u0, c0 = self.rows[0], self.cells[0]
        q = (self._reduce(penalty_rhs(self.topo, U, self.col, self.n, cells=self.cells))
             if self.gamma else None)
        with self._b as a:
            a[:self.nloc] = -self.gamma * q[u0:u0 + self.nloc] if self.gamma else 0.0
            a[self.nloc:] = rhs_c[c0:c0 + self.mloc]
        del q
        self._consistent(self._b)
        bnorm = self._b.norm()
        floor = max(FLUX_REFINE * FLUX_TOL * FLUX_FLOOR * max(facet, 1e-300),
                    RESIDUAL_STALL * self._pnorm(self._b))
        self._minres(self._b, self._x, who, floor)
        s = self._gather(self._x)
        npass = 0
        worse = np.inf
        why = "%d passes, the limit" % REFINE_MAX
        while npass < REFINE_MAX:
            ratio, _ = flux_ratios(self.topo, U + self._spread(s))
            now = float(ratio.max()) / (FLUX_REFINE * FLUX_TOL)
            if now <= 1.0:
                why = "the target was reached"
                break
            if now > worse / REFINE_GAIN:
                why = "the last pass gained %.3gx, less than the %gx one must" \
                    % (worse / now, REFINE_GAIN)
                break
            worse = now
            self.A.mult(self._x, self._r)
            self._r.axpy(-1.0, self._b)
            self._consistent(self._r)
            # a residual at round-off holds a part the flux rows cannot reach
            if self._r.norm() < RESIDUAL_FLOOR * max(bnorm, 1e-300):
                why = "the residual is at round-off"
                break
            # a correction that does not converge is a pass that does not pay:
            # the step before it stands, where the first solve has none to fall to
            try:
                self._minres(self._r, self._dx, who, floor, max_it=REFINE_MAXITER)
            except RuntimeError:
                why = "a pass did not converge in its %d iterations" % REFINE_MAXITER
                break
            self._x.axpy(-1.0, self._dx)
            s = self._gather(self._x)
            npass += 1
        return s, npass, why

    def _reduce(self, v):
        """A per-cell sum every rank made a part of, added up."""
        if self.mpi.size > 1:
            self.mpi.Allreduce(MPI.IN_PLACE, v, op=MPI.SUM)
        return v

    # ---- the step itself
    def _spread(self, s):
        """The change of the free midpoints as a change of every P2 node."""
        topo = self.topo
        out = np.zeros((topo.nnodes, topo.dim))
        has = self.col >= 0
        out[topo.nverts:][has] = s.reshape(-1, topo.dim)[self.col[has]]
        return out

    def apply(self, U, name=""):
        """U with the free edge midpoints moved by the step. Returns (U, info);
        name is the stamp, for the failure. Whether to step at all is decided on
        the fluxes alone, so balanced data is returned as it is, mean and all."""
        topo = self.topo
        who = name + ": " if name else ""
        if not np.isfinite(U).all():
            raise ValueError("%sthe velocity holds a non-finite value at node %d"
                             % (who, int(np.flatnonzero(~np.isfinite(U).all(axis=1))[0])))
        r, big = cell_flux(topo, U)
        ratio, facet = flux_ratios(topo, U, (r, big))
        umax = max(float(np.abs(U).max()), 1e-300)
        # the interior facets cancel in the sum, so it is the net flux through
        # the boundary: reported, never refused
        imbalance = abs(float(r.sum())) / max(facet, 1e-300)
        target = self._spread_net(float(r.sum()),
                                  np.maximum(big, FLUX_FLOOR * max(facet, 1e-300)))
        mean_in = split_volume_mean(topo, U)
        npass = 0
        stepped = True
        self._its = 0
        if ratio.max() <= FLUX_EXIT * FLUX_TOL:
            # a solve here would only chase round-off
            s = np.zeros(self.n)
            stepped = False
            why = "the data was balanced already, so no step was taken"
        else:
            s, npass, why = self._step(U, target - r, who, facet)
        self.iterations = max(self.iterations, self._its)
        self.refinements = max(self.refinements, npass)
        ds = self._spread(s)
        out = U + ds
        r2, big2 = cell_flux(topo, out)
        ratio2, _ = flux_ratios(topo, out, (r2, big2))
        if not ratio2.max() <= FLUX_TOL:
            raise RuntimeError("%sthe flux solve left cell %d at %.2e of its scale, more "
                               "than the %.0e a loader accepts"
                               % (who, int(np.argmax(ratio2)), ratio2.max(), FLUX_TOL))
        mean_out = split_volume_mean(topo, out)
        change = ds[topo.nverts:]
        return out, dict(flux_before_max=float(np.abs(r).max()),
                         facet_flux_max=facet, imbalance=imbalance,
                         flux_after_max=float(np.abs(r2).max()),
                         balance_before=float(ratio.max()), balance_after=float(ratio2.max()),
                         mean_before=mean_in, mean_after=mean_out,
                         drift=(mean_out - mean_in) / umax,
                         div_rms=float(np.sqrt(np.mean(div_norms(topo, out) ** 2))),
                         stepped=stepped, refinements=npass, refine_stop=why,
                         change=change, n_free=len(self.free))


def _axes(mask):
    """The names of the coordinates a mask selects."""
    return "".join("xyz"[c] for c in np.nonzero(mask)[0]) or "no direction"


def throughput_drift(topo, drift):
    """The drift of the volume mean in the periodic directions. With
    impermeable walls the mean of a divergence-free field in a periodic
    direction is that direction's throughput, times its length over the volume;
    such a field also keeps a uniform tracer density uniform, so that mean is
    the mean Lagrangian velocity of uniformly seeded tracers and whatever the
    cleaning does to it goes one for one into the result."""
    return drift[topo.periodic] if topo.periodic.any() else np.zeros(0)


def change_report(topo, change, umax):
    """The size of the midpoint change relative to |u|max, in the cells with an
    exterior facet the mesh does not pair and in the rest."""
    mag = np.linalg.norm(change, axis=1)
    wall = np.zeros(topo.nedges, bool)
    k = np.nonzero(topo.cell_boundary)[0]
    if len(k):
        wall[topo.cell_edge[k].ravel()] = True
    out = {}
    for name, m in (("boundary", wall), ("bulk", ~wall)):
        v = mag[m] / max(umax, 1e-300)
        out[name] = (float(np.sqrt(np.mean(v ** 2))) if len(v) else 0.0,
                     float(v.max()) if len(v) else 0.0, int(m.sum()))
    return out


# ------------------------------------------------- the split field's matrix


def _exps(d, k):
    import itertools
    return [e for e in itertools.product(range(k + 1), repeat=d) if sum(e) <= k]


def _mono(p, E):
    return np.array([[np.prod(np.asarray(q, float) ** np.array(e)) for e in E] for q in p])


def _dmono(p, E, j):
    out = np.zeros((len(p), len(E)))
    for c, e in enumerate(E):
        if e[j] == 0:
            continue
        f = list(e)
        f[j] -= 1
        out[:, c] = e[j] * np.prod(np.asarray(p, float) ** np.array(f), axis=1)
    return out


def reference_nodes(d):
    """(vertices, boundary nodes, interior nodes) of the Alfeld split of the
    reference simplex, in the order the per-cell data uses."""
    V = np.vstack([np.zeros(d), np.eye(d)])
    B = [V[i] for i in range(d + 1)] + [(V[a] + V[b]) / 2 for a, b in LOC_EDGES[d + 1]]
    z = V.mean(axis=0)
    I = [z] + [(V[i] + z) / 2 for i in range(d + 1)]
    return V, np.array(B), np.array(I)


def reference_matrix(d):
    """R with u_int = R g on the reference simplex: the unique P2 field on the
    Alfeld split with boundary data g and div u = 0. Also the rank check."""
    V, B, I = reference_nodes(d)
    z = V.mean(axis=0)
    nb, ni = len(B), len(I)
    key = lambda p: tuple(np.round(p, 10))
    kb = {key(p): j for j, p in enumerate(B)}
    ki = {key(p): j for j, p in enumerate(I)}
    E = _exps(d, 2)
    rowsA, rowsB = [], []
    for o in range(d + 1):
        S = np.vstack([np.delete(V, o, axis=0), z])
        nodes = np.vstack([S] + [(S[a] + S[b]) / 2 for a, b in LOC_EDGES[d + 1]])
        coef = np.linalg.inv(_mono(nodes, E))
        G = [_dmono(S, E, c) @ coef for c in range(d)]   # div at d+1 points: P1 unisolvent
        A = np.zeros((d + 1, ni * d))
        Bm = np.zeros((d + 1, nb * d))
        for i, p in enumerate(nodes):
            k = key(p)
            for c in range(d):
                if k in kb:
                    Bm[:, kb[k] * d + c] += G[c][:, i]
                else:
                    A[:, ki[k] * d + c] += G[c][:, i]
        rowsA.append(A)
        rowsB.append(Bm)
    A = np.vstack(rowsA)
    Bm = np.vstack(rowsB)
    rank = np.linalg.matrix_rank(A, tol=1e-9)
    R = -np.linalg.pinv(A) @ Bm
    return R, dict(rank=int(rank), ndof=ni * d, rows=A.shape[0],
                   residual=float(np.abs(A @ R + Bm).max()))


_R_CACHE = {}


def split_interior(topo, U):
    """The split field's interior values for every cell: (ncells, d+2, d), the
    barycenter first and then the midpoints of the vertex-barycenter edges."""
    d = topo.dim
    if d not in _R_CACHE:
        _R_CACHE[d] = reference_matrix(d)[0]
    R = _R_CACHE[d]
    g = np.concatenate([U[topo.cells], U[topo.nverts + topo.cell_edge]], axis=1)
    gh = np.einsum('cij,cnj->cni', topo.Jinv, g)
    uh = (gh.reshape(topo.ncells, -1) @ R.T).reshape(topo.ncells, -1, d)
    return np.einsum('cij,cnj->cni', topo.J, uh)


def split_mesh(topo):
    """(vertices, cells) of the barycentric split: the original vertices, then
    one barycenter a cell; each cell replaced by its d+1 sub-simplices."""
    z = topo.X[topo.cells].mean(axis=1)
    Xs = np.vstack([topo.X, z])
    zid = topo.nverts + np.arange(topo.ncells)
    sub = []
    for o in range(topo.nv):
        f = [i for i in range(topo.nv) if i != o]
        sub.append(np.column_stack([topo.cells[:, f], zid]))
    return Xs, np.concatenate(sub, axis=0)


def split_nodes(topo, U, uint):
    """(positions, values) of every P2 node of the split mesh: the original
    vertices and edge midpoints, the barycenters, and the midpoints of the
    vertex-barycenter edges."""
    z = topo.X[topo.cells].mean(axis=1)
    inner_x = 0.5 * (topo.X[topo.cells] + z[:, None, :])
    pos = np.vstack([topo.node_x, z, inner_x.reshape(-1, topo.dim)])
    val = np.vstack([U, uint[:, 0], uint[:, 1:].reshape(-1, topo.dim)])
    return pos, val


def _p2_basis(S, p, grad=False):
    """The P2 basis of a simplex at a point, in barycentric form: vertices
    first, then the edge midpoints in the local order. (values) or (gradients)."""
    d = S.shape[1]
    T = np.vstack([S.T, np.ones(d + 1)])
    lam = np.linalg.solve(T, np.append(p, 1.0))
    if not grad:
        v = [l * (2 * l - 1) for l in lam]
        v += [4 * lam[a] * lam[b] for a, b in LOC_EDGES[d + 1]]
        return lam, np.array(v)
    G = np.linalg.solve(T, np.vstack([np.eye(d), np.zeros(d)]))
    g = [(4 * lam[i] - 1) * G[i] for i in range(d + 1)]
    g += [4 * (lam[b] * G[a] + lam[a] * G[b]) for a, b in LOC_EDGES[d + 1]]
    return lam, np.array(g)


def _cell_node_values(topo, U, uint, k):
    """A dict from the rounded position of every P2 node of macro cell k's split
    to its value."""
    V = topo.X[topo.cells[k]]
    z = V.mean(axis=0)
    out = {}
    vb = np.concatenate([U[topo.cells[k]], U[topo.nverts + topo.cell_edge[k]]])
    for j, p in enumerate(np.vstack([V] + [(V[a] + V[b]) / 2 for a, b in LOC_EDGES[topo.nv]])):
        out[tuple(np.round(p, 10))] = vb[j]
    for j, p in enumerate(np.vstack([z[None]] + [(V[i] + z) / 2 for i in range(topo.nv)])):
        out[tuple(np.round(p, 10))] = uint[k][j]
    return V, z, out


def eval_split(topo, U, uint, pts, cell_of, grad=False):
    """The split field at points, each given the macro cell it lies in; with
    grad, its pointwise divergence."""
    d = topo.dim
    out = np.zeros((len(pts), 1 if grad else d))
    for q, (p, k) in enumerate(zip(pts, cell_of)):
        V, z, nodemap = _cell_node_values(topo, U, uint, k)
        best, bestlam = None, -np.inf
        for o in range(d + 1):
            S = np.vstack([np.delete(V, o, axis=0), z])
            lam, B = _p2_basis(S, p, grad)
            if lam.min() > bestlam:
                bestlam, best = lam.min(), (S, B)
            if lam.min() > -1e-10:
                break
        S, B = best
        nodes = np.vstack([S] + [(S[a] + S[b]) / 2 for a, b in LOC_EDGES[d + 1]])
        vals = np.array([nodemap[tuple(np.round(n, 10))] for n in nodes])
        out[q] = np.einsum('nc,nc->', B, vals) if grad else B @ vals
    return out


def p2_eval(topo, U, pts, cell_of, grad=False):
    """The plain P2 field on the original cells at points."""
    d = topo.dim
    out = np.zeros((len(pts), 1 if grad else d))
    for q, (p, k) in enumerate(zip(pts, cell_of)):
        S = topo.X[topo.cells[k]]
        _, B = _p2_basis(S, p, grad)
        vals = np.concatenate([U[topo.cells[k]], U[topo.nverts + topo.cell_edge[k]]])
        out[q] = np.einsum('nc,nc->', B, vals) if grad else B @ vals
    return out


def p1_to_p2(topo, U1):
    """A P1 vertex field as P2 nodes: the edge midpoints the mean of their ends."""
    return np.vstack([U1, 0.5 * (U1[topo.edges[:, 0]] + U1[topo.edges[:, 1]])])


# ------------------------------------------------------------ reading a case


def _root():
    """Every printed line and every written file is the first rank's; the others
    read the same files, solve their share and keep quiet."""
    return MPI.COMM_WORLD.rank == 0


@contextlib.contextmanager
def _root_writes(what):
    """Writing only the first rank does. A failure there -- a full disk, a
    permission -- is one the others cannot have, so they would wait at the next
    barrier forever: it takes the job down instead. On one rank there is nobody
    to wait, and the exception is the better answer."""
    try:
        yield
    except BaseException as exc:
        if MPI.COMM_WORLD.size == 1:
            raise
        sys.stderr.write("rank 0 could not %s: %s: %s\n"
                         % (what, type(exc).__name__, exc))
        sys.stderr.flush()
        MPI.COMM_WORLD.Abort(1)


def _h5py():
    import h5py
    return h5py


def read_mesh_h5(path):
    """(coordinates, topology, cell_indices) of a dolfin HDF5 mesh."""
    with _h5py().File(str(path), "r") as f:
        topo = np.array(f["mesh/topology"]).astype(np.int32)
        X = np.array(f["mesh/coordinates"])
        ci = (np.array(f["mesh/cell_indices"]).astype(np.int64)
              if "mesh/cell_indices" in f else None)
    return X, topo, ci


def read_element(path, field):
    """(degree, ncomp) from the element signature, as the loaders parse it."""
    import re
    with _h5py().File(str(path), "r") as f:
        if field not in f or "signature" not in f[field].attrs:
            return None
        sig = f[field].attrs["signature"]
    sig = sig.decode() if hasattr(sig, "decode") else str(sig)
    m = re.search(r"(triangle|tetrahedron),\s*(\d+)", sig)
    if not m:
        raise ValueError("%s: '%s' is %s, which names no degree" % (path, field, sig))
    deg = int(m.group(2))
    ncomp = 1
    if sig.startswith("VectorElement"):
        c = re.search(r"dim=(\d+)", sig)
        ncomp = int(c.group(1)) if c else (2 if m.group(1) == "triangle" else 3)
    return deg, ncomp


class DofTable:
    """The per-cell dof table of a dataset's first file, and any file's values by
    node through it.

    Only the first file of a dataset carries the table; a later one holds just
    `<field>/vector_0`, which is all the loaders open for it.
    """

    def __init__(self, first, field, topo, cell_indices):
        h5py = _h5py()
        self.field = field
        self.topo = topo
        self.degree, self.ncomp = read_element(first, field)
        nloc = topo.nv + (topo.cell_edge.shape[1] if self.degree == 2 else 0)
        with h5py.File(str(first), "r") as f:
            g = f[field]
            cd = np.array(g["cell_dofs"]).astype(np.int32)
            xc = np.array(g["x_cell_dofs"]).astype(np.int64)
            fc = np.array(g["cells"]).astype(np.int32)
        per = nloc * self.ncomp
        if not np.all(np.diff(xc) == per):
            raise ValueError("%s: '%s/x_cell_dofs' is ragged" % (first, field))
        # a vector element's cell dofs come component by component, nloc each
        rows = cd.reshape(topo.ncells, self.ncomp, nloc).transpose(0, 2, 1)
        if cell_indices is not None:
            # both files label their rows by global cell id
            row_of = np.full(int(cell_indices.max()) + 1, -1, np.int32)
            row_of[cell_indices] = np.arange(topo.ncells)
            dest = row_of[fc]
            if (dest < 0).any():
                raise ValueError("%s: '%s/cells' names a cell the mesh does not"
                                 % (first, field))
            out = np.empty_like(rows)
            out[dest] = rows
            rows = out
        self.rows = rows
        self.nodes = (np.concatenate([topo.cells, topo.nverts + topo.cell_edge], axis=1)
                      if self.degree == 2 else topo.cells)
        self.n_total = topo.nverts + (topo.nedges if self.degree == 2 else 0)

    def values(self, path):
        """The field of one stamp or component, by node.

        A cell chunk at a time: the whole gather is three times the field's own
        size in temporaries, which on a mesh of millions of cells is the largest
        allocation of the read and is paid on every rank."""
        with _h5py().File(str(path), "r") as f:
            vec = np.array(f[self.field + "/vector_0"])
        if not np.isfinite(vec).all():
            raise ValueError("%s: '%s' holds a non-finite value at dof %d"
                             % (path, self.field, int(np.flatnonzero(~np.isfinite(vec))[0])))
        values = np.full((self.n_total, self.ncomp), np.nan)
        per = max(1, READ_CHUNK_VALUES // (self.rows.shape[1] * self.ncomp))
        for c0 in range(0, len(self.rows), per):
            sl = slice(c0, c0 + per)
            values[self.nodes[sl].ravel()] = vec[self.rows[sl].reshape(-1, self.ncomp)]
        if np.isnan(values).any():
            raise ValueError("%s: '%s' leaves a node without a value" % (path, self.field))
        for c0 in range(0, len(self.rows), per):
            sl = slice(c0, c0 + per)
            if not np.array_equal(values[self.nodes[sl].ravel()],
                                  vec[self.rows[sl].reshape(-1, self.ncomp)]):
                raise ValueError("%s: cells disagree on a node of '%s'"
                                 % (path, self.field))
        return values


def read_params(path):
    """(ordered keys and values, the file's lines) of a parameter file."""
    prm, lines = {}, []
    for line in Path(path).read_text().split("\n"):
        lines.append(line)
        if "=" in line and not line.strip().startswith("#"):
            k, v = line.split("=", 1)
            prm[k.strip()] = v.strip()
    return prm, lines


def read_stamps(folder, prm):
    """(kind, entries) of a case's stamp list: `timestamps` lines are
    (["t"], file) and `freqstamps` lines (["t", "a"], file) or
    (["omega", "phi", "a"], file), the numbers kept as written."""
    kind = "freqstamps" if "freqstamps" in prm else "timestamps"
    name = prm.get(kind, kind + ".dat")
    entries = []
    for line in (folder / name).read_text().split("\n"):
        tok = line.split()
        if not tok or tok[0].startswith("#"):
            continue
        entries.append((tok[:-1], tok[-1]))
    return kind, name, entries


def read_case(params_path, periodic_tol=PERIODIC_TOL):
    """The mesh, the stamp list and the field names of a dolfin HDF5 case."""
    folder = Path(params_path).parent
    prm, lines = read_params(params_path)
    X, cells, ci = read_mesh_h5(folder / prm["mesh"])
    X = X[:, :cells.shape[1] - 1]
    kind, stampfile, stamps = read_stamps(folder, prm)
    periodic = [prm.get("periodic_%s" % a, "false") == "true" for a in "xyz"]
    return dict(folder=folder, prm=prm, lines=lines, X=X, cells=cells, cell_indices=ci,
                kind=kind, stampfile=stampfile, stamps=stamps,
                periodic=periodic[:X.shape[1]], periodic_tol=periodic_tol)


# ------------------------------------------------------------------- writing


def _signature(degree, ncomp, dim):
    cell = "triangle" if dim == 2 else "tetrahedron"
    fe = "FiniteElement('Lagrange', %s, %d)" % (cell, degree)
    return fe if ncomp == 1 else "VectorElement(%s, dim=%d)" % (fe, ncomp)


def write_mesh_h5(path, X, cells):
    """A dolfin HDF5 mesh: what `read_mesh` opens, and nothing else."""
    h5py = _h5py()
    X = np.ascontiguousarray(X, dtype=float)
    cells = np.ascontiguousarray(cells, dtype=np.int64)
    with h5py.File(str(path), "w") as f:
        g = f.create_group("mesh")
        t = g.create_dataset("topology", data=cells)
        t.attrs["celltype"] = np.bytes_("triangle" if cells.shape[1] == 3 else "tetrahedron")
        t.attrs["partition"] = np.array([0], dtype=np.uint64)
        g.create_dataset("coordinates", data=X)
        g.create_dataset("cell_indices", data=np.arange(len(cells), dtype=np.int64))


def _cell_dofs(topo, degree, ncomp):
    """The per-cell dof table of our own numbering: dof = node*ncomp + component,
    the components blocked as dolfin blocks them."""
    node = (np.concatenate([topo.cells, topo.nverts + topo.cell_edge], axis=1)
            if degree == 2 else topo.cells)
    rows = [node * ncomp + c for c in range(ncomp)]
    return np.concatenate(rows, axis=1).astype(np.int32)


def write_checkpoint(path, field, values, topo, degree, mode="w", cell_indices=None):
    """One field as a dolfin checkpoint: the dof table, the cell list and the
    values, in the numbering `write_vector` then reuses.

    The rows are the mesh's own, and `cells` labels them by global cell id, as
    dolfin does: `cell_indices` is the map the mesh beside them carries, which a
    reader composes through to find the row it wants. Without it the labels are
    `0..n-1`, which is the map a serially written mesh carries."""
    h5py = _h5py()
    values = np.ascontiguousarray(values, dtype=float)
    ncomp = values.shape[1]
    cd = _cell_dofs(topo, degree, ncomp)
    gid = (np.arange(topo.ncells) if cell_indices is None
           else np.asarray(cell_indices)).astype(np.uint64)
    if len(gid) != topo.ncells:
        raise ValueError("%s: %d cell indices for %d cells" % (path, len(gid), topo.ncells))
    with h5py.File(str(path), mode) as f:
        g = f.create_group(field)
        g.attrs["signature"] = np.bytes_(_signature(degree, ncomp, topo.dim))
        g.create_dataset("cell_dofs", data=cd.ravel())
        g.create_dataset("x_cell_dofs",
                         data=np.arange(topo.ncells + 1, dtype=np.uint64) * cd.shape[1])
        g.create_dataset("cells", data=gid)
        v = g.create_dataset("vector_0", data=values.ravel())
        v.attrs["partition"] = np.array([0], dtype=np.uint64)


def write_vector(path, field, values, degree, ncomp, dim, mode="w"):
    """A later stamp or component: `<field>/vector_0`, which is all `read_vector`
    opens, and the signature, which the frequency loader checks on every
    component against the first."""
    with _h5py().File(str(path), mode) as f:
        g = f.create_group(field)
        g.attrs["signature"] = np.bytes_(_signature(degree, ncomp, dim))
        g.create_dataset("vector_0", data=np.ascontiguousarray(values, dtype=float).ravel())


def copy_group(src, dst, name):
    """A field copied through unchanged, datasets and attributes as they are."""
    h5py = _h5py()
    with h5py.File(str(src), "r") as a:
        if name not in a:
            return False
        with h5py.File(str(dst), "a") as b:
            a.copy(a[name], b, name)
    return True


def write_params(path, lines, changes):
    """The input's parameter file with these keys set; keys it does not hold are
    appended."""
    out, seen = [], set()
    for line in lines:
        if "=" in line and not line.strip().startswith("#"):
            k = line.split("=", 1)[0].strip()
            if k in changes:
                seen.add(k)
                if changes[k] is None:
                    continue
                out.append("%s=%s" % (k, changes[k]))
                continue
        out.append(line)
    while out and not out[-1].strip():
        out.pop()
    for k, v in changes.items():
        if k not in seen and v is not None:
            out.append("%s=%s" % (k, v))
    Path(path).write_text("\n".join(out) + "\n")


def write_stamps(path, kind, entries):
    """The stamp list with the numbers as written and the new file names."""
    Path(path).write_text("\n".join(" ".join(list(n) + [f]) for n, f in entries) + "\n")


# ------------------------------------------------------------------ checking


def taylor_hood_moments(topo, U):
    """(max, rms) of |m_i| / n_i over the vertices, with m_i = int phi_i div u
    over the P1 hat of vertex i and n_i = int phi_i |div u| the scale of the
    terms that make it up; periodic images are one vertex.

    A converged Taylor-Hood velocity satisfies int q div u = 0 for every P1 q, so
    every m_i is at round-off; data that has been interpolated, projected or
    lifted from P1 has m_i of the size of n_i. The construction takes any data,
    but only compatible data has a split field whose volume mean is the
    solver's."""
    dv = cell_div(topo, U)
    M = mass_p1(topo.dim)
    w = np.einsum('k,ab,kb->ka', topo.vol, M, dv, optimize=True)
    v = topo.master[topo.cells]
    m = np.zeros(topo.nverts)
    n = np.zeros(topo.nverts)
    np.add.at(m, v.ravel(), w.ravel())
    np.add.at(n, v.ravel(), np.abs(w).ravel())
    ok = n > 1e-13 * max(n.max(), 1e-300)
    rel = np.abs(m[ok]) / n[ok]
    if not len(rel):
        return 0.0, 0.0
    return float(rel.max()), float(np.sqrt(np.mean(rel ** 2)))


def check_case(params_path, periodic_tol=PERIODIC_TOL, field=None, quiet=False):
    """The worst cell's balance in every stamp, under the criterion the loader
    applies, with the net flux through the boundary, the volume mean and the
    data's Taylor-Hood compatibility."""
    case = read_case(params_path, periodic_tol)
    topo = Topo(case["X"], case["cells"], case["periodic"], periodic_tol=periodic_tol)
    field = field or case["prm"].get("velocity_field", "u")
    first = case["folder"] / case["stamps"][0][1]
    dof = DofTable(first, field, topo, case["cell_indices"])
    worst = 0.0
    worst_th = 0.0
    for cols, name in case["stamps"]:
        U = dof.values(case["folder"] / name)
        if dof.degree == 1:
            U = p1_to_p2(topo, U)
        r, big = cell_flux(topo, U)
        ratio, facet = flux_ratios(topo, U, (r, big))
        rel = float(ratio.max())
        worst = float(np.maximum(worst, rel))
        mmax, mrms = taylor_hood_moments(topo, U)
        worst_th = float(np.maximum(worst_th, mrms))
        if not quiet and _root():
            print("  %s %s: max |net flux| %.3e, largest facet flux %.3e, worst cell %.3e of "
                  "its scale (a loader refuses above %.0e), net boundary flux %.3e of the "
                  "largest facet flux"
                  % (" ".join(cols), name, np.abs(r).max(), facet, rel, FLUX_TOL,
                     abs(float(r.sum())) / max(facet, 1e-300)))
            shift = div_moment(topo, U) / topo.vol.sum()
            print("      |int phi_i div u| / int phi_i |div u|: max %.2e, rms %.2e;"
                  " volume mean of the split field (%s), the throughput in %s, which the "
                  "macro field misses by %.2e of |u|max"
                  % (mmax, mrms,
                     " ".join("%.6g" % x for x in split_volume_mean(topo, U)),
                     _axes(topo.periodic),
                     np.abs(shift).max() / max(float(np.abs(U).max()), 1e-300)))
    if not quiet and _root():
        print("  the data %s Taylor-Hood-compatible: rms |m_i| %.2e of the scale, and a"
              " converged Taylor-Hood velocity has it at round-off"
              % ("looks" if worst_th < 1e-6 else "does not look", worst_th))
    return worst


# ------------------------------------------------------------- the whole case


def _bary_eval(topo, values, degree, cell, lam):
    """A P1 or P2 field of the macro mesh at barycentric coordinates of a cell."""
    if degree == 1:
        return np.einsum('cn,cnk->ck', lam, values[topo.cells[cell]])
    vb = np.concatenate([values[topo.cells[cell]],
                         values[topo.nverts + topo.cell_edge[cell]]], axis=1)
    B = [lam[:, i] * (2 * lam[:, i] - 1) for i in range(topo.nv)]
    B += [4 * lam[:, a] * lam[:, b] for a, b in LOC_EDGES[topo.nv]]
    return np.einsum('nc,cnk->ck', np.array(B), vb)


def _order_index(want, pos):
    """The position in pos of each wanted node; the mesh decides it, so a
    dataset builds the tree once a field rather than once a stamp."""
    from scipy.spatial import cKDTree
    d, idx = cKDTree(pos).query(want)
    if d.max() > 1e-9:
        raise ValueError("a node of the split mesh has no value: %g" % d.max())
    return idx


def split_scalar_nodes(topo, values, degree):
    """(positions, values) of a macro P1 or P2 field at the split mesh's nodes of
    the same degree; exact, since its restriction to a sub-simplex is a
    polynomial of that degree."""
    nc, nv = topo.ncells, topo.nv
    cell = np.arange(nc)
    z = topo.X[topo.cells].mean(axis=1)
    at_z = _bary_eval(topo, values, degree, cell, np.full((nc, nv), 1.0 / nv))
    if degree == 1:
        return np.vstack([topo.X, z]), np.vstack([values, at_z])
    inner_v, inner_x = [], []
    for i in range(nv):
        lam = np.full((nc, nv), 0.5 / nv)
        lam[:, i] += 0.5
        inner_v.append(_bary_eval(topo, values, degree, cell, lam))
        inner_x.append(0.5 * (topo.X[topo.cells[:, i]] + z))
    pos = np.vstack([topo.node_x, z, np.stack(inner_x, axis=1).reshape(-1, topo.dim)])
    val = np.vstack([values, at_z, np.stack(inner_v, axis=1).reshape(-1, values.shape[1])])
    return pos, val


def clean_case(params_path, out, split=False, weights="volume", write_key=True,
               periodic_tol=PERIODIC_TOL, verbose=True, penalty=None,
               rest_tol=None, boundary_weight=None, field_file=None):
    """Clean a whole dataset and write the output case. Returns a report.

    field_file names one file of the case's folder to clean in place of the
    series: every stamp then names it, so a field written beside the case's own
    becomes a steady case of its own. The report carries the topology the output
    is written on and the last file's values by its nodes, for a caller that
    writes its own view of the field.

    Every rank reads the same files, builds the same tables and holds the same
    field; only the solve is shared out, and only the first rank writes and
    prints. The report is every rank's.
    """
    case = read_case(params_path, periodic_tol)
    folder, prm = case["folder"], case["prm"]
    table = case["stamps"][0][1]      # the file that carries the dof table
    if field_file is not None:
        name = Path(field_file).name
        if not (folder / name).is_file():
            raise ValueError("%s holds no '%s'" % (folder, name))
        case["stamps"] = [(cols, name) for cols, _ in case["stamps"]]
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    u_name = prm.get("velocity_field", "u")
    p_name = prm.get("pressure_field", "p")
    phi_name = (prm.get("phase_field", "phi")
                if prm.get("include_phi", "false") == "true" else None)
    have_p = prm.get("ignore_pressure", "false") != "true"

    t0 = time.time()
    topo = Topo(case["X"], case["cells"], case["periodic"], periodic_tol=periodic_tol)
    files = []
    for _, name in case["stamps"]:
        if name not in files:
            files.append(name)
    dof = DofTable(folder / table, u_name, topo, case["cell_indices"])
    if dof.ncomp != topo.dim:
        raise ValueError("'%s' has %d components, not %d" % (u_name, dof.ncomp, topo.dim))
    # the held set is the dataset's: a node at rest in every stamp or component,
    # so every stamp is read before the first is cleaned. What fits the budget
    # is kept, so a short series is read once.
    nmax = np.zeros(topo.nnodes)
    scale = 0.0
    cached = {}
    # the budget is the node's, and every rank runs this pass
    budget = READ_CACHE_BYTES // MPI.COMM_WORLD.size
    for name in files:
        U = dof.values(folder / name)
        if dof.degree == 1:
            U = p1_to_p2(topo, U)
        nmax = np.maximum(nmax, np.linalg.norm(U, axis=1))
        scale = max(scale, float(np.abs(U).max()))
        if U.nbytes <= budget:
            cached[name] = U
            budget -= U.nbytes
    topo.set_held(at_rest_nodes(nmax, scale, rest_tol))
    read_time = time.time() - t0

    eq = Equil(topo, weights=weights, penalty=penalty,
               boundary_weight=boundary_weight)
    n_bfree = int((topo.boundary_edge & ~topo.held_edge).sum())
    root = _root()
    verbose = verbose and root
    if verbose:
        print("  %d cells, %d edges, %d held nodes, %d held midpoints; read in %.1f s"
              % (topo.ncells, topo.nedges, int(topo.held_node.sum()),
                 int(topo.held_edge.sum()), read_time))
        print("  objective %s"
              % ("s^T W s + gamma sum_K ||div(u+s)||^2, g = %g, h = %.4g, gamma = %.4g"
                 % (eq.g, eq.h, eq.gamma) if eq.gamma else "s^T W s (the smallest change)"))
        print("  %s, %d free midpoints, %d of them on the boundary at weight %g, setup "
              "%.1f s%s"
              % (eq.solver, len(eq.free), n_bfree, eq.boundary_weight, eq.setup_time,
                 ", singular (every boundary midpoint held)" if eq.singular else ""))

    if split:
        Xs, cs = split_mesh(topo)
        stopo = Topo(Xs, cs, case["periodic"], periodic_tol=periodic_tol)
        if root:
            with _root_writes("write the split mesh"):
                write_mesh_h5(out / "mesh.h5", Xs, cs)
        mesh_name = "mesh.h5"
        out_cells = None                  # the split mesh's cell_indices is the identity
    else:
        mesh_name = Path(prm["mesh"]).name
        if root:
            with _root_writes("link the mesh into the output"):
                _link_or_copy(folder / prm["mesh"], out / mesh_name)
        stopo = topo
        # the input's mesh is the output's, so its labels are the ones to write
        out_cells = case["cell_indices"]

    # the pressure and the phase field come through untouched: their datasets as
    # they are on the same mesh, resampled where --split changes the mesh
    others = []
    for name, take in ((p_name, have_p), (phi_name, phi_name is not None)):
        if not take:
            continue
        others.append((name, DofTable(folder / table, name, topo, case["cell_indices"])
                       if split else None))

    rename = {}
    values = None
    u_order, orders = None, {}
    rep = dict(before=0.0, after=0.0, facet=0.0, boundary=(0.0, 0.0, 0), bulk=(0.0, 0.0, 0),
               mean_before=None, mean_after=None, mean_shift=0.0, imbalance=0.0,
               throughput_drift=0.0, balance_after=0.0, div_rms=0.0, stepped=False,
               refine_stop="")
    t1 = time.time()
    for i, name in enumerate(files):
        U = cached.pop(name, None)
        if U is None:
            U = dof.values(folder / name)
            if dof.degree == 1:
                U = p1_to_p2(topo, U)
        Uc, info = eq.apply(U, name)
        umax = float(np.linalg.norm(Uc, axis=1).max())
        ch = change_report(topo, info["change"], umax)
        rep["before"] = max(rep["before"], info["flux_before_max"])
        rep["after"] = max(rep["after"], info["flux_after_max"])
        rep["facet"] = max(rep["facet"], info["facet_flux_max"])
        rep["div_rms"] = max(rep["div_rms"], info["div_rms"])
        # the stamp that came out worst is the one whose stop reason explains the run
        if info["balance_after"] >= rep["balance_after"]:
            rep["refine_stop"] = info["refine_stop"]
        rep["balance_after"] = max(rep["balance_after"], info["balance_after"])
        rep["imbalance"] = max(rep["imbalance"], info["imbalance"])
        rep["stepped"] = rep["stepped"] or info["stepped"]
        if rep["mean_before"] is None:
            rep["mean_before"] = info["mean_before"]
            rep["mean_after"] = info["mean_after"]
        tp = throughput_drift(topo, info["drift"])
        rep["throughput_drift"] = max(rep["throughput_drift"],
                                      float(np.abs(tp).max()) if len(tp) else 0.0)
        rep["mean_shift"] = max(rep["mean_shift"], float(np.abs(info["drift"]).max()))
        for k in ("boundary", "bulk"):
            rep[k] = (max(rep[k][0], ch[k][0]), max(rep[k][1], ch[k][1]), ch[k][2])
        if verbose:
            print("  %s: worst cell %.2e of its scale, net boundary flux %.2e of the largest "
                  "facet flux, volume mean (%s) -> (%s), throughput drift in %s %.2e of |u|max"
                  % (name, info["balance_after"], info["imbalance"],
                     " ".join("%.6g" % x for x in info["mean_before"]),
                     " ".join("%.6g" % x for x in info["mean_after"]),
                     _axes(topo.periodic), float(np.abs(tp).max()) if len(tp) else 0.0))
        if split:
            pos, val = split_nodes(topo, Uc, split_interior(topo, Uc))
            if u_order is None:
                u_order = _order_index(stopo.node_x, pos)
            values = val[u_order]
        else:
            values = Uc
        new = "u_%04d.h5" % i
        rename[name] = new
        if root:
            with _root_writes("write %s" % (out / new)):
                if i == 0:
                    write_checkpoint(out / new, u_name, values, stopo, 2,
                                     cell_indices=out_cells)
                else:
                    write_vector(out / new, u_name, values, 2, topo.dim, topo.dim)
        # a field that comes through untouched is the first rank's alone: nothing
        # of it reaches the report
        for other, odof in (others if root else []):
            with _root_writes("write '%s' of %s" % (other, out / new)):
                if odof is None:
                    if not copy_group(folder / name, out / new, other):
                        raise ValueError("%s holds no '%s' to carry through; name what is "
                                         "in it, or set ignore_pressure"
                                         % (folder / name, other))
                    continue
                pos, val = split_scalar_nodes(topo, odof.values(folder / name), odof.degree)
                if other not in orders:
                    orders[other] = _order_index(
                        stopo.node_x if odof.degree == 2 else stopo.X, pos)
                if i == 0:
                    write_checkpoint(out / new, other, val[orders[other]], stopo,
                                     odof.degree, mode="a", cell_indices=out_cells)
                else:
                    write_vector(out / new, other, val[orders[other]], odof.degree,
                                 odof.ncomp, topo.dim, mode="a")
    solve_time = time.time() - t1

    if root:
        with _root_writes("write the stamp list and the parameter file"):
            write_stamps(out / case["stampfile"], case["kind"],
                         [(cols, rename[name]) for cols, name in case["stamps"]])
            divfree = write_key and not split
            changes = {"mesh": mesh_name, "velocity_space": "P2",
                       "divfree": "true" if divfree else None}
            if divfree:
                changes["mesh_cache"] = None      # that loader refuses a cached mesh
            write_params(out / Path(params_path).name, case["lines"], changes)
    MPI.COMM_WORLD.Barrier()
    if verbose:
        print("  flux %.3e -> %.3e (largest facet flux %.3e), worst cell %.2e of its "
              "scale, %d iterations, %d refinements, %.1f s"
              % (rep["before"], rep["after"], rep["facet"], rep["balance_after"],
                 eq.iterations, eq.refinements, solve_time))
        print("  refinement stopped: %s; the worst cell is %s"
              % (rep["refine_stop"],
                 "%.3gx inside what a loader accepts" % (FLUX_TOL / rep["balance_after"])
                 if rep["balance_after"] > 0.0 else "carrying no net flux at all"))
        print("  volume mean of the reconstructed field, first stamp (%s) -> (%s); moved by "
              "at most %.2e of |u|max, the throughput in %s by %.2e; net boundary flux at "
              "most %.2e of the largest facet flux; rms ||div u||_K %.4g"
              % (" ".join("%.6g" % x for x in rep["mean_before"]),
                 " ".join("%.6g" % x for x in rep["mean_after"]), rep["mean_shift"],
                 _axes(topo.periodic), rep["throughput_drift"], rep["imbalance"],
                 rep["div_rms"]))
        print("  midpoint change / |u|max: boundary cells rms %.2e max %.2e (%d edges), "
              "elsewhere rms %.2e max %.2e (%d edges)"
              % (rep["boundary"][0], rep["boundary"][1], rep["boundary"][2],
                 rep["bulk"][0], rep["bulk"][1], rep["bulk"][2]))
        if dof.degree == 1:
            print("  the input is P1: the split field stops the trapping, but it keeps the P1"
                  " interpolant's deficit in the mean velocity, so cleaning the solver's"
                  " P2 output is preferable")
        print("  wrote %d stamps to %s" % (len(files), out))
    rep.update(cells=topo.ncells, edges=topo.nedges, held_nodes=int(topo.held_node.sum()),
               held_edges=int(topo.held_edge.sum()), boundary_free=n_bfree, solver=eq.solver,
               iterations=eq.iterations, refinements=eq.refinements, degree_in=dof.degree,
               out=out, g=eq.g, gamma=eq.gamma, h=eq.h, topo=stopo, values=values,
               boundary_weight=eq.boundary_weight,
               setup_time=eq.setup_time, solve_time=solve_time, read_time=read_time)
    eq.destroy()
    return rep


def _link_or_copy(src, dst):
    """The mesh beside the cleaned fields: a hard link where the file system
    allows one."""
    if Path(dst).exists():
        os.remove(dst)
    try:
        os.link(os.path.realpath(src), dst)
    except OSError:
        shutil.copy(src, dst)


# ----------------------------------------------------------------------- CLI


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="divfree_clean.py",
        description="Make a dolfin HDF5 velocity dataset divergence-free in every cell.",
        epilog="The cleaned case is read with divfree=true; --no-key leaves that key out, "
               "so it is read as a plain P2 field.")
    p.add_argument("params", help="the case's parameter file (dolfin_params.dat)")
    p.add_argument("--out", help="output folder; without it nothing is written")
    p.add_argument("--field-file", default=None, metavar="FILE",
                   help="clean this file of the case's folder in place of the series the "
                        "stamp list names, so a field written beside the case's own is "
                        "cleaned as a steady case of its own")
    p.add_argument("--check", action="store_true",
                   help="report each stamp's worst cell balance under the criterion a loader "
                        "applies, its net boundary flux, its volume mean and the Taylor-Hood "
                        "moments of its divergence, and stop")
    p.add_argument("--split", action="store_true",
                   help="write the full reconstruction as P2 on the barycentric split mesh, "
                        "which the plain mesh loaders read (a diagnostic; needs no dolfin)")
    p.add_argument("--weights", choices=("volume", "none"), default="volume",
                   help="weight of the midpoint change: by the volume around an edge, so a "
                        "graded mesh is not adjusted where it is fine, or unweighted")
    p.add_argument("--penalty", type=float, default=None, metavar="G",
                   help="weight of the divergence left in the cells, gamma = G h^2 with h the "
                        "median edge length of the mesh (default %g in 2D and %g in 3D); "
                        "0 is the smallest change alone" % (DEFAULT_G[2], DEFAULT_G[3]))
    p.add_argument("--boundary-weight", type=float, default=BOUNDARY_W, metavar="W",
                   help="least-change weight of a free midpoint on an exterior facet the mesh "
                        "does not pair, relative to an interior one (default %(default)g); it "
                        "keeps a moving wall's values and an open facet's flux near the data's")
    p.add_argument("--no-key", action="store_true",
                   help="do not write divfree=true into the output parameter file")
    p.add_argument("--rest-tol", type=float, default=REST_TOL, metavar="T",
                   help="a node whose |u| stays at or below this times the dataset's max |u| "
                        "in every stamp is at rest, and one at rest on an exterior facet the "
                        "mesh does not pair is held (default %(default)g)")
    p.add_argument("--periodic-tol", type=float, default=PERIODIC_TOL,
                   help="position tolerance pairing periodic images (default %(default)g, "
                        "the dolfin HDF5 loaders' own)")
    a = p.parse_args(argv)
    if a.check:
        # the whole dataset, on one rank: it solves nothing
        if _root():
            print("worst cell balance %.3e of the %.0e a loader accepts"
                  % (check_case(a.params, periodic_tol=a.periodic_tol), FLUX_TOL))
        return 0
    if not a.out:
        p.error("--out is required unless --check is given")
    clean_case(a.params, a.out, split=a.split, weights=a.weights, write_key=not a.no_key,
               periodic_tol=a.periodic_tol, penalty=a.penalty, field_file=a.field_file,
               rest_tol=a.rest_tol, boundary_weight=a.boundary_weight)
    return 0


if __name__ == "__main__":
    sys.exit(main())
