"""wall_upgrade=true on a P2 velocity (mode=triangle and tet, and the
frequency modes trianglefreq and tetfreq): the P2 field plus quartic terms
in the cells with a wall vertex at rest, and a correction in the cells with
a wall facet, a wall edge (3D) or a wall vertex (2D) at rest.

The reference here is a numpy evaluation of the rule. The quartic terms w
vanish at every P2 node and are shared by the cells around their entity:
per edge end at a wall vertex, per face over a wall edge and, in 3D, a
bubble per cell with wall facets. The correction is continuous and of
degree 4 on the cell's centroid split, zero on the cell's faces, mapped to
the cell by A from the regular simplex of edge 1, under div u = 0 and n.grad
div u = 0 on the wall facets (in their sub-cells) and div u = 0 on the wall
edges or vertices (in every sub-cell holding them): the least reference H1
field there. The loaders must reproduce it, value and gradient, at random
points (interpol): in 2D on a kinked channel periodic along it, in 3D around
a sphere in a periodic box, and on the unit square and cube walled on every
side, whose corner and edge cells have two wall facets, at a time between
two stamps. Outside the cells with a wall vertex the field is the plain P2
one bit for bit. The frequency modes read the same two files as components,
summed with their weights, against the reference summed with the same
weights. A meshio facet file (wall_facets) naming the walls reads as the
walls found at rest; naming the sides off x of a box alone, with an inlet
and an outlet on them and the first stamp at rest, it matches the reference
whose walls are those sides. Plane Poiseuille flow, which P2 holds exactly,
is carried exactly by the upgrade, periodic or through an inlet and an outlet.
"""
import os
import shutil
from itertools import permutations, product

import numpy as np
import pytest

from dumps import dump_at
from paths import app, needs
from runs import from_points, interpolation, read_checkpoint, run_app
from wall_cases import DT, dolfin_mesh, kink_wall, kinked_channel_mesh, probe, restart_from
from wall_cases import sphere_mesh  # noqa: F401  (fixture)

INTERPOL = app("interpol")
TRACERS = app("tracers")
TRACERVECTORS = app("tracervectors")

pytestmark = needs(INTERPOL, TRACERS, TRACERVECTORS)

T0 = 0.4   # the probe time, between the stamps at 0 and 1

LOCAL_EDGES = {2: [(1, 2), (0, 2), (0, 1)],
               3: [(2, 3), (1, 3), (1, 2), (0, 3), (0, 2), (0, 1)]}
# 3D: the faces; 2D: the cell
LOCAL_FACES = {2: [(0, 1, 2)], 3: [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]}


# ------------------------------------------------------------------ reference


def bary2_nodes(dim):
    """Barycentrics of the P2 nodes in dolfin's order."""
    I = np.eye(dim + 1)
    return np.array([I[i] for i in range(dim + 1)] + [(I[i] + I[j]) / 2 for i, j in LOCAL_EDGES[dim]])


def p2_basis(L):
    dim = L.shape[1] - 1
    cols = [L[:, i] * (2 * L[:, i] - 1) for i in range(dim + 1)]
    cols += [4 * L[:, i] * L[:, j] for (i, j) in LOCAL_EDGES[dim]]
    return np.array(cols).T


def p2_dbasis(L):
    """(n, nodes, nv): d basis / d lambda_m."""
    dim = L.shape[1] - 1
    nv = dim + 1
    out = np.zeros((len(L), nv + len(LOCAL_EDGES[dim]), nv))
    for i in range(nv):
        out[:, i, i] = 4 * L[:, i] - 1
    for e, (i, j) in enumerate(LOCAL_EDGES[dim]):
        out[:, nv + e, i] = 4 * L[:, j]
        out[:, nv + e, j] = 4 * L[:, i]
    return out


def bary_gradients(X, cells):
    """(ncells, nv, dim): the gradients of the barycentrics."""
    V = X[cells]
    A = np.concatenate([V.transpose(0, 2, 1), np.ones((len(cells), 1, V.shape[1]))], axis=1)
    return np.linalg.inv(A)[:, :, :-1]


def wall_slots(X, cells, periodic, tol=1e-9):
    """(cell, local vertex) of every exterior facet off the periodic seams."""
    nv = cells.shape[1]
    lo, hi = X.min(axis=0), X.max(axis=0)
    owners = {}
    for c in range(len(cells)):
        for k in range(nv):
            owners.setdefault(tuple(sorted(cells[c][j] for j in range(nv) if j != k)), []).append((c, k))
    out = []
    for f, o in owners.items():
        if len(o) != 1:
            continue
        P = X[list(f)]
        if any(periodic[d] and (np.all(np.abs(P[:, d] - lo[d]) < tol) or np.all(np.abs(P[:, d] - hi[d]) < tol))
               for d in range(X.shape[1])):
            continue
        out.append(o[0])
    return np.array(sorted(out), np.int64)


K = 4   # the correction's degree on a sub-cell

# the regular simplex of edge 1
REG = {2: np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]]),
       3: np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float) / np.sqrt(8)}


def lattice(nv, k):
    return np.array([a for a in product(range(k + 1), repeat=nv) if sum(a) == k])


def lagrange(k, A, mu, hess=False):
    """Degree-k Lagrange basis of nodes A (nloc, nv) at mu (p, nv): values (p, nloc), d/dmu (p, nloc, nv)
    and, with hess, d2/dmu2 (p, nloc, nv, nv)."""
    p, n, nv = len(mu), len(A), mu.shape[1]
    f = np.ones((p, n, nv)); df = np.zeros((p, n, nv)); ddf = np.zeros((p, n, nv))
    for l, a in enumerate(A):
        for i in range(nv):
            v, d, dd = np.ones(p), np.zeros(p), np.zeros(p)
            for m in range(a[i]):
                t = (k * mu[:, i] - m) / (m + 1); tp = k / (m + 1)
                dd = dd * t + 2 * d * tp
                d = d * t + v * tp
                v = v * t
            f[:, l, i], df[:, l, i], ddf[:, l, i] = v, d, dd
    phi = f.prod(2)
    dphi = np.empty((p, n, nv))
    for i in range(nv):
        o = [m for m in range(nv) if m != i]
        dphi[:, :, i] = df[:, :, i] * f[:, :, o].prod(2)
    if not hess:
        return phi, dphi
    H = np.empty((p, n, nv, nv))
    for i in range(nv):
        for j in range(nv):
            o = [m for m in range(nv) if m not in (i, j)]
            H[:, :, i, j] = (ddf[:, :, i] if i == j else df[:, :, i] * df[:, :, j]) * f[:, :, o].prod(2)
    return phi, dphi, H


def gauss_simplex(dim, n):
    """Collapsed Gauss rule on the unit simplex, as barycentrics (q, dim + 1) and weights summing to 1."""
    x, w = np.polynomial.legendre.leggauss(n)
    x, w = (x + 1) / 2, w / 2
    P, W = [], []
    for c in product(zip(x, w), repeat=dim):
        p, wt, rest = [], 1.0, 1.0
        for a, wa in c:
            p.append(a * rest)
            wt *= wa * rest
            rest *= 1 - a
        P.append([1 - sum(p)] + p); W.append(wt)
    W = np.array(W)
    return np.array(P), W / W.sum()


class SubSplit:
    """The centroid split of the reference cell: each sub-cell's local nodes (multi-indices over the
    macro vertices but j in order, then the centroid), the interior node of each (or -1 on the faces of
    the cell), quadrature."""

    def __init__(self, dim, k=K):
        nv = dim + 1
        self.k, self.nv, self.dim = k, nv, dim
        self.A = lattice(nv, k)
        self.verts = [[m for m in range(nv) if m != j] for j in range(nv)]
        keys, self.idx = {}, np.full((nv, len(self.A)), -1)
        for j in range(nv):
            for l, a in enumerate(self.A):
                if a[-1] == 0:
                    continue
                key = np.full(nv, a[-1])
                for s, m in enumerate(self.verts[j]):
                    key[m] += nv * a[s]
                self.idx[j, l] = keys.setdefault(tuple(key), len(keys))
        self.nint = len(keys)
        self.qp, self.qw = gauss_simplex(dim, k + 1)
        _, self.dq = lagrange(k, self.A, self.qp)

    def mu_grads(self, G):
        """(ncells, nv sub-cells, nv, dim): grad mu of sub-cell j, from bary gradients G (ncells, nv, dim)."""
        nv = self.nv
        out = np.empty((len(G), nv, nv, G.shape[2]))
        for j in range(nv):
            for s, m in enumerate(self.verts[j]):
                out[:, j, s] = G[:, m] - G[:, j]
            out[:, j, nv - 1] = nv * G[:, j]
        return out

    def mu_of(self, j, lam):
        return np.concatenate([lam[:, self.verts[j]] - lam[:, j:j + 1], self.nv * lam[:, j:j + 1]], 1)


def entities(dim):
    """The codimension-2 entities of a cell by their vertices: 3D edges 01, 02, 03, 12, 13, 23; 2D vertices."""
    return [(a, b) for a in range(4) for b in range(a + 1, 4)] if dim == 3 else [(v, v) for v in range(3)]


class WallRows:
    """The conditions: div at the degree k-1 lattice of each codimension-2 entity in every sub-cell holding
    it (a vertex itself in 2D), and on each facet f in sub-cell f, div at its degree k-1 lattice and
    n . grad div at its degree k-2 lattice. bit: the pattern bit of each row (facets, then the entities)."""

    def __init__(self, sp):
        nv, k = sp.nv, sp.k
        I = np.eye(nv)
        self.j, self.lam, self.normal, self.bit = [], [], [], []

        def add(j, lam, normal, bit):
            self.j.append(j); self.lam.append(lam); self.normal.append(normal); self.bit.append(bit)

        def face(fl, p):
            if p == 0:
                return [I[fl].mean(0)]
            return [sum(b[i] / p * I[fl[i]] for i in range(len(fl)))
                    for b in product(range(p + 1), repeat=len(fl)) if sum(b) == p]
        for q, (a, b) in enumerate(entities(sp.dim)):
            for j in range(nv):
                if j in (a, b):
                    continue
                for t in ([0.0] if a == b else np.arange(k) / (k - 1)):
                    lam = (1 - t) * I[a] + t * I[b] if a != b else I[a]
                    add(j, lam, False, nv + q)
        for f in range(nv):
            fl = [m for m in range(nv) if m != f]
            for lam in face(fl, k - 1):
                add(f, lam, False, f)
            for lam in face(fl, k - 2):
                add(f, lam, True, f)
        self.j, self.lam, self.normal, self.bit = map(np.array, (self.j, self.lam, self.normal, self.bit))
        mu = np.array([sp.mu_of(j, l[None])[0] for j, l in zip(self.j, self.lam)])
        _, self.dphi, self.H = lagrange(k, sp.A, mu, hess=True)


def reference_system(sp, wr):
    """On the regular simplex: the barycentric gradients, the split's nodes (FullNodes), the stiffness over
    every node, the rows over every node (nrows, nodes, dim) and their targets per vertex divergence."""
    nv, dim = sp.nv, sp.dim
    R = REG[dim]
    G = bary_gradients(R, np.arange(nv)[None])[0]
    Gm = sp.mu_grads(G[None])[0]
    vol = abs(np.linalg.det(R[1:] - R[:1])) / np.prod(np.arange(1, dim + 1))
    fs = FullNodes(sp)
    S = np.zeros((fs.n, fs.n))
    for j in range(nv):
        gq = np.einsum('qlm,mx->qlx', sp.dq, Gm[j])
        ii = fs.full[j]
        S[np.ix_(ii, ii)] += np.einsum('q,qlx,qkx->lk', sp.qw, gq, gq) * vol / nv
    nr = len(wr.j)
    B = np.zeros((nr, fs.n, dim)); T = np.zeros((nr, nv))
    for r in range(nr):
        j = wr.j[r]; ii = fs.full[j]
        if not wr.normal[r]:
            B[r, ii] = np.einsum('lm,mx->lx', wr.dphi[r], Gm[j])
            T[r] = -wr.lam[r]
        else:
            n = G[wr.bit[r]] / np.linalg.norm(G[wr.bit[r]])
            B[r, ii] = np.einsum('lmo,mx,oy,y->lx', wr.H[r], Gm[j], Gm[j], n)
            T[r] = -G @ n
    return G, fs, S, B, T


class FullNodes:
    """Every node of the split: sub-cell j's local node l is node full[j, l]; the interior ones first, in
    SubSplit's order; lam their barycentrics."""

    def __init__(self, sp):
        nv, k = sp.nv, sp.k
        keys, self.full, lam = {}, np.full((nv, len(sp.A)), -1), {}
        nb = sp.nint
        for j in range(nv):
            for l, a in enumerate(sp.A):
                key = np.full(nv, a[-1])
                for s, m in enumerate(sp.verts[j]):
                    key[m] += nv * a[s]
                key = tuple(key)
                if sp.idx[j, l] >= 0:
                    i = sp.idx[j, l]
                elif key in keys:
                    i = keys[key]
                else:
                    i = keys[key] = nb
                    nb += 1
                lam[i] = np.array(key) / (nv * k)
                self.full[j, l] = i
        self.n = nb
        self.lam = np.array([lam[i] for i in range(nb)])


def correction_tables(sp, wr, patterns):
    """Per wall pattern Wg (nint, dim, nrows), the least reference |grad g|^2 on the split, zero on the faces
    of the cell, under the pattern's rows with targets t: g_I = Wg t in reduced units."""
    G, fs, S, B, T = reference_system(sp, wr)
    ni = sp.nint
    Si = np.linalg.inv(S[:ni, :ni])
    BI = B[:, :ni]
    out = {}
    for p in patterns:
        act = ((p >> wr.bit) & 1).astype(float)
        BIa = BI * act[:, None, None]
        MBt = np.einsum('ab,rbx->axr', Si, BIa)
        Sch = np.einsum('rax,axs->rs', BIa, MBt)
        out[p] = np.einsum('axr,rs->axs', MBt, np.linalg.pinv(Sch, rcond=1e-11, hermitian=True)) * act
    return out, G


def vertex_classes(X, periodic, tol=1e-9):
    """Each vertex's class: the least index among its periodic images."""
    key = np.where(np.array(periodic)[None, :] & (np.abs(X - 1) < tol), 0.0, X)
    _, first, inv = np.unique(np.round(key / tol).astype(np.int64), axis=0, return_index=True,
                              return_inverse=True)
    return first[inv.ravel()]


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else None


class Upgrade:
    """The wall terms of a P2 field: quartic terms, which vanish at every P2 node. Per edge end at a wall vertex
    w, towards v: -((4 u_m - u_v).n) n lambda_w^2 lambda_v (lambda_w - lambda_v), n the facets' at w with apex
    v, else all of them at w. Per face (a, b, v) over a wall edge ab: lambda_a lambda_b lambda_v (f_a lambda_a
    + f_b lambda_b), f_a = -(2 p_a + p_b) n, f_b = -(p_a + 2 p_b) n, p = n.(4 u_m - u_v) along a -> v and b ->
    v, n the facets' at ab with apex v, else all of them there (2D: the cell over its wall edge). In 3D, per
    cell with wall facets, the bubble D lambda_0 lambda_1 lambda_2 lambda_3, D the least-norm
    solution of n_k.D = -2 n_k.(sum of the slopes towards k) over its facets k. Every wall entity at rest."""

    def __init__(self, X, cells, slots, vclass):
        self.X, self.cells, self.K = X, cells, vclass
        self.dim = X.shape[1]
        nv = self.nv = self.dim + 1
        facets = []
        for c, k in slots:
            fv = sorted(cells[c][j] for j in range(nv) if j != k)
            P = X[fv]
            if self.dim == 3:
                n = np.cross(P[1] - P[0], P[2] - P[0])
            else:
                t = P[1] - P[0]
                n = np.array([-t[1], t[0]])
            n = n / np.linalg.norm(n)
            if n @ (X[cells[c][k]] - P[0]) > 0:
                n = -n
            facets.append((tuple(fv), c, k, n))
        facets.sort(key=lambda f: f[0])
        self.facet_n = {(c, k): n for _, c, k, n in facets}
        K = vclass
        self.count, self.n_sum = {}, {}
        self.side, self.edge, self.face = {}, {}, {}
        for fv, c, k, n in facets:
            apex = K[cells[c][k]]
            cls = [K[v] for v in fv]
            for w in cls:
                if w not in self.count:
                    self.count[w], self.n_sum[w] = 0, np.zeros(self.dim)
                self.count[w] += 1
                self.n_sum[w] = self.n_sum[w] + n
                self.side[(w, apex)] = self.side.get((w, apex), 0) + n
            for i in range(len(cls)):
                for j in range(i + 1, len(cls)):
                    e = (min(cls[i], cls[j]), max(cls[i], cls[j]))
                    self.edge[e] = self.edge.get(e, 0) + n
                    self.face[e + (apex,)] = self.face.get(e + (apex,), 0) + n
        self.wall_cell = np.array([any(K[v] in self.count for v in row) for row in cells])

    def end_normal(self, w, v):
        if w not in self.count:
            return None
        return unit(self.side.get((w, v), self.n_sum[w]))

    def face_normal(self, a, b, v):
        e = (min(a, b), max(a, b))
        if e not in self.edge:
            return None
        return unit(self.face.get(e + (v,), self.edge[e]))

    def coefficients(self, U2):
        """Per cell: edge ends (nedges, 2, dim), face monomials (nfaces, 4, dim) as lambda_a lambda_b lambda_c
        times (1, lambda_a, lambda_b, lambda_c), and the bubble (dim,), zero in 2D."""
        dim, nv, K = self.dim, self.nv, self.K
        eidx = {}
        for e, (i, j) in enumerate(LOCAL_EDGES[dim]):
            eidx[(i, j)] = eidx[(j, i)] = e
        nc = len(self.cells)
        C = np.zeros((nc, len(LOCAL_EDGES[dim]), 2, dim))
        F = np.zeros((nc, len(LOCAL_FACES[dim]), 4, dim))
        D = np.zeros((nc, dim))
        for c in np.nonzero(self.wall_cell)[0]:
            row, u = self.cells[c], U2[c]
            s = lambda a, v: 4 * u[nv + eidx[(a, v)]] - u[v]               # noqa: E731
            for e, (i, j) in enumerate(LOCAL_EDGES[dim]):
                for end, (w, v) in enumerate(((i, j), (j, i))):
                    n = self.end_normal(K[row[w]], K[row[v]])
                    if n is not None:
                        C[c, e, end] = -(n @ s(w, v)) * n
            for f, tri in enumerate(LOCAL_FACES[dim]):
                for (si, sj, sv) in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
                    a, b, v = tri[si], tri[sj], tri[sv]
                    n = self.face_normal(K[row[a]], K[row[b]], K[row[v]])
                    if n is None:
                        continue
                    pa, pb = (n @ s(a, v)) * n, (n @ s(b, v)) * n
                    F[c, f, 1 + si] += -(2 * pa + pb)
                    F[c, f, 1 + sj] += -(pa + 2 * pb)
            if dim == 3:
                rows = []
                for k in range(nv):
                    if (c, k) not in self.facet_n:
                        continue
                    n = self.facet_n[(c, k)]
                    rows.append((n, -2.0 * float(n @ sum(s(i, k) for i in range(nv) if i != k))))
                if rows:
                    D[c] = np.linalg.pinv(np.array([r[0] for r in rows])) @ np.array([r[1] for r in rows])
        return C, F, D

    def weights(self, L):
        """(edge (n, nedges, 2), face (n, nfaces, 4), bubble (n,)) and their lambda derivatives, one more axis
        of nv."""
        dim, nv, n = self.dim, self.nv, len(L)
        E = np.zeros((n, len(LOCAL_EDGES[dim]), 2))
        dE = np.zeros(E.shape + (nv,))
        for e, (i, j) in enumerate(LOCAL_EDGES[dim]):
            for end, (a, b) in enumerate(((i, j), (j, i))):
                E[:, e, end] = L[:, a] ** 2 * L[:, b] * (L[:, a] - L[:, b])
                dE[:, e, end, a] = 3 * L[:, a] ** 2 * L[:, b] - 2 * L[:, a] * L[:, b] ** 2
                dE[:, e, end, b] = L[:, a] ** 3 - 2 * L[:, a] ** 2 * L[:, b]
        Fw = np.zeros((n, len(LOCAL_FACES[dim]), 4))
        dF = np.zeros(Fw.shape + (nv,))
        for f, tri in enumerate(LOCAL_FACES[dim]):
            i, j, k = tri
            B = L[:, i] * L[:, j] * L[:, k]
            dB = {i: L[:, j] * L[:, k], j: L[:, i] * L[:, k], k: L[:, i] * L[:, j]}
            Fw[:, f, 0] = B
            for m in tri:
                dF[:, f, 0, m] = dB[m]
            for slot, x in enumerate(tri):
                Fw[:, f, 1 + slot] = B * L[:, x]
                for m in tri:
                    dF[:, f, 1 + slot, m] = dB[m] * L[:, x] + (B if m == x else 0.0)
        Bw = np.prod(L, axis=1) if dim == 3 else np.zeros(n)
        dBw = np.zeros((n, nv))
        if dim == 3:
            for m in range(nv):
                dBw[:, m] = np.prod(np.delete(L, m, axis=1), axis=1)
        return (E, Fw, Bw), (dE, dF, dBw)

    def terms(self, CFD, cid, L, G=None):
        """The terms' value (n, dim) and, given the barycentric gradients G of cells cid, gradient at barycentrics
        L of cells cid."""
        C, F, D = CFD
        (E, Fw, Bw), (dE, dF, dBw) = self.weights(L)
        u = np.einsum("pea,pead->pd", E, C[cid]) + np.einsum("pfk,pfkd->pd", Fw, F[cid]) + Bw[:, None] * D[cid]
        if G is None:
            return u
        J = np.einsum("peam,pmx,pead->pdx", dE, G, C[cid]) + np.einsum("pfkm,pmx,pfkd->pdx", dF, G, F[cid]) \
            + np.einsum("pm,pmx,pd->pdx", dBw, G, D[cid])
        return u, J


class Rule:
    """The wall correction of a P2 field on a mesh, every wall entity at rest: u = P2 + w in the cells with a
    wall vertex (w the wall terms), plus in the split cells (a wall facet, a wall edge in 3D, a wall vertex in 2D)
    the field on the split that is zero on the faces of the cell, the least |grad g|^2 under the rows on the
    whole field, in reduced units on the regular simplex."""

    def __init__(self, X, cells, slots, periodic):
        self.X, self.cells = X, cells
        self.dim = X.shape[1]
        self.nv = nv = self.dim + 1
        self.G = bary_gradients(X, cells)
        self.facets = np.zeros((len(cells), nv), bool)   # facet k of a cell faces vertex k
        self.facets[slots[:, 0], slots[:, 1]] = True
        # the codimension-2 entities on a wall facet, by vertex class
        vc = vertex_classes(X, periodic)
        C = vc[cells]
        on = {}
        for c, k in slots:
            fv = [C[c, j] for j in range(nv) if j != k]
            for a, b in entities(self.dim):
                if a < len(fv) and b < len(fv):
                    on[tuple(sorted((fv[a], fv[b])))] = True
        self.pattern = self.facets @ (1 << np.arange(nv))
        for q, (a, b) in enumerate(entities(self.dim)):
            hit = np.array([tuple(sorted((C[c, a], C[c, b]))) in on for c in range(len(cells))])
            self.pattern = self.pattern | (hit.astype(np.int64) << (nv + q))
        self.split_cell = self.pattern > 0
        self.up = Upgrade(X, cells, slots, vc)
        self.wall_cell = self.up.wall_cell
        self.sp = SubSplit(self.dim)
        self.wr = WallRows(self.sp)
        patterns = np.unique(self.pattern[self.split_cell])
        self.tab, self.Gref = correction_tables(self.sp, self.wr, patterns)
        V = X[cells]
        R = REG[self.dim]
        self.A = np.einsum('cxk,ky->cxy', (V[:, 1:] - V[:, :1]).transpose(0, 2, 1),
                           np.linalg.inv((R[1:] - R[:1]).T))

    def coefficients(self, U2):
        """The wall terms' coefficients, and the split's interior node values (ncells, nint, dim) beyond P2 + w."""
        CFD = self.up.coefficients(U2)
        out = np.zeros((len(self.cells), self.sp.nint, self.dim))
        for p, tab in self.tab.items():
            c = np.nonzero(self.pattern == p)[0]
            A = self.A[c]
            out[c] = np.einsum('cxy,cay->cax', A, np.einsum('axr,cr->cax', tab, self.targets(U2, CFD, c, A)))
        return CFD, out

    def targets(self, U2, CFD, c, A):
        """The rows' targets on v = P2 + w in cells c: -div v at a value row, -(A n^).grad div v at a normal
        row, by a five-point difference in the barycentrics (exact for div v of degree <= 4)."""
        wr, n = self.wr, len(c)
        G = self.Gref
        nr = len(wr.j)
        t = np.zeros((n, nr))

        def div(L):
            m = len(L)
            cid = np.repeat(c, m)
            Lr = np.tile(L, (n, 1))
            _, J = self.evaluate(U2, (CFD, None), cid, Lr, split=False)
            return np.einsum('pii->p', J).reshape(n, m)
        val = np.nonzero(~wr.normal)[0]
        t[:, val] = -div(wr.lam[val])
        nor = np.nonzero(wr.normal)[0]
        h = 1e-2
        acc = np.zeros((n, len(nor)))
        for s, w in ((2, -1.0), (1, 8.0), (-1, -8.0), (-2, 1.0)):
            L = np.stack([wr.lam[r] + s * h * (G @ (G[wr.bit[r]] / np.linalg.norm(G[wr.bit[r]]))) for r in nor])
            acc += w * div(L)
        t[:, nor] = -acc / (12 * h)
        return t

    def evaluate(self, U2, co, cid, L, terms=True, split=True):
        """(u (n, dim), J (n, dim, dim), J[p, i, j] = du_i/dx_j) at barycentrics L of cells cid; co from
        coefficients."""
        G = self.G[cid]
        u = np.einsum("pn,pnd->pd", p2_basis(L), U2[cid])
        J = np.einsum("pnm,pmx,pnd->pdx", p2_dbasis(L), G, U2[cid])
        if not terms:
            return u, J
        CFD, coef = co
        uw, Jw = self.up.terms(CFD, cid, L, G)
        u, J = u + uw, J + Jw
        if not split:
            return u, J
        sp = self.sp
        j = np.argmin(L, axis=1)
        for jj in range(self.nv):
            m = np.nonzero(j == jj)[0]
            if len(m) == 0:
                continue
            phi, dphi = lagrange(sp.k, sp.A, sp.mu_of(jj, L[m]))
            ok = sp.idx[jj] >= 0
            C = coef[cid[m]][:, sp.idx[jj, ok]]
            u[m] += np.einsum('pl,pld->pd', phi[:, ok], C)
            gr = np.einsum('plm,pmx->plx', dphi[:, ok], sp.mu_grads(G[m])[:, jj])
            J[m] += np.einsum('pld,plx->pdx', C, gr)
        return u, J


def locate(X, cells, pts, k=64):
    """(cell, barycentrics) of each point: of the nearest cells by centroid, the one it is deepest in."""
    from scipy.spatial import cKDTree
    _, cand = cKDTree(X[cells].mean(axis=1)).query(pts, k=k)
    best_c = np.full(len(pts), -1)
    best_l = np.zeros((len(pts), cells.shape[1]))
    best = np.full(len(pts), -np.inf)
    b = np.concatenate([pts, np.ones((len(pts), 1))], axis=1)
    for j in range(k):
        c = cand[:, j]
        V = X[cells[c]]
        A = np.concatenate([V.transpose(0, 2, 1), np.ones((len(c), 1, V.shape[1]))], axis=1)
        lam = np.linalg.solve(A, b[:, :, None])[:, :, 0]
        take = lam.min(axis=1) > best
        best[take], best_c[take], best_l[take] = lam.min(axis=1)[take], c[take], lam[take]
    return best_c, best_l, best


def box_mesh(dim, n, seed=3):
    """(vertex coordinates, cells) of the unit square or cube, n boxes a side:
    a square's diagonal alternating, so each corner square has a triangle on
    two sides; a cube in the six tets around its main diagonal, so the cubes
    along the box's edges have tets on two sides. The interior vertices are
    moved at random by up to a fifth of a box, the sides stay flat."""
    g = np.linspace(0, 1, n + 1)
    X = np.stack(np.meshgrid(*[g] * dim, indexing="ij"), axis=-1).reshape(-1, dim)

    def idx(*i):
        return int(np.ravel_multi_index(i, (n + 1,) * dim))

    cells = []
    if dim == 2:
        for i in range(n):
            for j in range(n):
                a, b, c, d = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
                cells += [(a, b, d), (b, c, d)] if (i + j) % 2 == 0 else [(a, b, c), (a, c, d)]
    else:
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for p in permutations(range(3)):
                        v = [i, j, k]
                        row = [idx(*v)]
                        for d in p:
                            v[d] += 1
                            row.append(idx(*v))
                        cells.append(tuple(row))
    inner = np.all((X > 1e-12) & (X < 1 - 1e-12), axis=1)
    X[inner] += np.random.default_rng(seed).uniform(-0.2 / n, 0.2 / n, (inner.sum(), dim))
    return X, np.array(cells)


# ------------------------------------------------------------------ cases


def write_facets(folder, X, faces, markers, name="facets"):
    """faces (vertex ids into X) and their markers as meshio writes a facet
    file: name.xdmf, its data in name.h5 as /data0 to /data2."""
    import h5py
    faces = np.asarray(faces, np.int32)
    with h5py.File(folder / (name + ".h5"), "w") as f:
        f["data0"] = np.asarray(X, float)
        f["data1"] = faces
        f["data2"] = np.asarray(markers, np.int32)
    n, d = X.shape
    m, k = faces.shape
    (folder / (name + ".xdmf")).write_text(
        '<Xdmf Version="3.0"><Domain><Grid Name="Grid">'
        '<Geometry GeometryType="%s"><DataItem DataType="Float" Dimensions="%d %d" Format="HDF" Precision="8">'
        '%s.h5:/data0</DataItem></Geometry>'
        '<Topology TopologyType="%s" NumberOfElements="%d" NodesPerElement="%d"><DataItem DataType="Int" '
        'Dimensions="%d %d" Format="HDF" Precision="4">%s.h5:/data1</DataItem></Topology>'
        '<Attribute Name="markers" AttributeType="Scalar" Center="Cell"><DataItem DataType="Int" '
        'Dimensions="%d" Format="HDF" Precision="4">%s.h5:/data2</DataItem></Attribute>'
        '</Grid></Domain></Xdmf>\n'
        % ("XYZ" if d == 3 else "XY", n, d, name, "Triangle" if k == 3 else "Polyline", m, k, m, k, name, m, name))


def side_marker(P, tol=1e-12):
    """1 and 2 for a facet on x = 0 and 1, 3 and 4 on y, 5 and 6 on z, else 7."""
    for d in range(P.shape[1]):
        if np.all(np.abs(P[:, d]) < tol):
            return 1 + 2 * d
        if np.all(np.abs(P[:, d] - 1) < tol):
            return 2 + 2 * d
    return 7


class Case:
    """A P2 checkpoint case, two stamps at t = 0 and 1, read with
    wall_upgrade=mode; case.U2[k] are stamp k's P2 values per cell in dolfin's
    local order. Its walls are the exterior facets off the periodic sides that
    wall(points) holds for, all by default; facets.xdmf lists every exterior
    facet off the periodic sides, marked by side_marker."""

    def __init__(self, df, base, X, cells, periodic, exprs, mode, seed_point, wall=None):
        dim = X.shape[1]
        # ascending vertices, so dolfin's local entities are the reference's
        cells = np.sort(cells, axis=1)
        self.base, self.mode, self.seed_point = base, mode, seed_point
        self.X, self.cells, self.periodic = X, cells, periodic
        mesh = dolfin_mesh(df, X, cells)

        class Periodic(df.SubDomain):
            def inside(self, x, on):
                return bool(on and any(periodic[k] and df.near(x[k], 0) for k in range(dim))
                            and not any(periodic[k] and df.near(x[k], 1) for k in range(dim)))

            def map(self, x, y):
                for k in range(dim):
                    y[k] = x[k] - 1 if periodic[k] and df.near(x[k], 1) else x[k]

        V = df.VectorFunctionSpace(mesh, "CG", 2, constrained_domain=Periodic())
        Q = df.FunctionSpace(mesh, "CG", 1, constrained_domain=Periodic())
        dofs = np.array([[V.sub(k).dofmap().cell_dofs(c) for k in range(dim)] for c in range(len(cells))])
        coords = V.tabulate_dof_coordinates().reshape(-1, dim)
        nodes = np.einsum("nl,cld->cnd", bary2_nodes(dim), X[cells])
        # dolfin's local order is the reference's, up to the periodic images
        d = np.abs(coords[dofs[:, 0]] - nodes)
        d = np.where(np.array(periodic)[None, None, :], np.minimum(d, np.abs(d - 1)), d)
        assert d.max() < 1e-10
        self.slots = wall_slots(X, cells, periodic)
        nv = dim + 1
        faces = [[cells[c][j] for j in range(nv) if j != k] for c, k in self.slots]
        if wall is not None:
            self.slots = np.array([s for s, f in zip(self.slots, faces) if wall(X[f])], np.int64)
        on_facet = [(c, n) for c, k in self.slots for n, b in enumerate(bary2_nodes(dim)) if b[k] == 0]
        os.makedirs(base, exist_ok=True)
        with df.HDF5File(mesh.mpi_comm(), str(base / "mesh.h5"), "w") as f:
            f.write(mesh, "mesh")
        self.U2 = []
        for s, expr in enumerate(exprs):
            u = df.interpolate(df.Expression(expr, degree=4), V)
            vec = u.vector().get_local()
            # no slip: every P2 node of a wall facet
            for c, n in on_facet:
                vec[dofs[c, :, n]] = 0.0
            u.vector().set_local(vec)
            p = df.interpolate(df.Expression("x[0] + %d" % s, degree=1), Q)
            with df.HDF5File(mesh.mpi_comm(), str(base / ("up_%d.h5" % s)), "w") as f:
                f.write(u, "u")
                f.write(p, "p")
            self.U2.append(vec[dofs].transpose(0, 2, 1))
        assert len(self.U2[0][0]) == nv + len(LOCAL_EDGES[dim])
        write_facets(base, X, faces, [side_marker(X[f]) for f in faces])
        (base / "timestamps.dat").write_text("0.0 up_0.h5\n1.0 up_1.h5\n")
        self.mesh_params = ("velocity_space=P2\npressure_space=P1\nmesh=mesh.h5\n"
                            + "".join("periodic_%s=%s\n" % (a, "true" if p else "false")
                                      for a, p in zip("xyz", list(periodic) + [False] * (3 - dim))))
        self.params = "timestamps=timestamps.dat\n" + self.mesh_params
        self.coef = {}
        self.rule = Rule(X, cells, self.slots, periodic)

    def variant(self, name, extra, series="timestamps=timestamps.dat"):
        """A folder of the case's files, its parameter file naming the time
        series (series) and holding extra."""
        d = self.base / name
        if not d.exists():
            d.mkdir()
            for f in os.listdir(self.base):
                if os.path.isfile(self.base / f):
                    shutil.copy(self.base / f, d / f)
            (d / "dolfin_params.dat").write_text(series + "\n" + self.mesh_params + extra)
        return d

    def __call__(self, mode):
        """The case folder read with wall_upgrade=mode ("true" or "false")."""
        return self.variant(mode, "wall_upgrade=%s\n" % mode)

    def coefficients(self, k):
        """The reference's coefficients of stamp k, computed once."""
        if k not in self.coef:
            self.coef[k] = self.rule.coefficients(self.U2[k])
        return self.coef[k]

    def reference(self, cid, L, terms=True, weights=(1 - T0, T0)):
        """u and J from the two stamps summed with weights, by default the blend at T0."""
        out = []
        for k, U2 in enumerate(self.U2):
            coef = self.coefficients(k) if terms else None
            out.append(self.rule.evaluate(U2, coef, cid, L, terms=terms))
        return tuple(weights[0] * a + weights[1] * b for a, b in zip(*out))


@pytest.fixture(scope="module")
def channel(tmp_path_factory):
    """The kinked channel of test_wall_edge_rule.py, periodic along x, a P2 field at rest on both walls."""
    df = pytest.importorskip("dolfin", reason="writing the checkpoint needs dolfin")
    X, cells = kinked_channel_mesh(16, 8)
    g = "(x[1] - 0.1*(0.5 - fabs(x[0] - 0.5)))*(1 - x[1])"
    exprs = [["%s*(4 + 2*sin(2*pi*x[0])) + 3*%s*%s" % (g, g, g), "%s*(1.2*cos(2*pi*x[0]) + x[1])" % g],
             ["%s*(3 - cos(2*pi*x[0])) - %s*%s" % (g, g, g), "%s*(0.7*sin(2*pi*x[0]) - 2*x[1])" % g]]
    return Case(df, tmp_path_factory.mktemp("upgrade_channel"), X, cells, (True, False), exprs,
                "triangle", (0.5, 0.5))


@pytest.fixture(scope="module")
def sphere(tmp_path_factory, sphere_mesh):
    """The sphere in the periodic unit cube of test_wall_edge_rule.py, a P2 field at rest on the sphere."""
    df = pytest.importorskip("dolfin", reason="writing the checkpoint needs dolfin")
    X, cells = sphere_mesh
    d = "(sqrt(pow(x[0]-0.5, 2) + pow(x[1]-0.5, 2) + pow(x[2]-0.5, 2)) - 0.25)"
    exprs = [["%s*(1 + 0.5*sin(2*pi*x[1]))" % d, "%s*(0.6*cos(2*pi*x[2]) + 2*%s)" % (d, d),
              "%s*0.4*sin(2*pi*x[0])" % d],
             ["%s*(0.8 - 0.3*cos(2*pi*x[2]))" % d, "%s*0.5*sin(2*pi*x[0])" % d,
              "%s*(0.9*cos(2*pi*x[1]) - %s)" % (d, d)]]
    return Case(df, tmp_path_factory.mktemp("upgrade_sphere"), X, cells, (True, True, True), exprs,
                "tet", (0.5, 0.1, 0.5))


@pytest.fixture(scope="module")
def box2(tmp_path_factory):
    """box_mesh's unit square, 8 a side, walled on every side, a P2 field at rest on them."""
    df = pytest.importorskip("dolfin", reason="writing the checkpoint needs dolfin")
    X, cells = box_mesh(2, 8)
    b = "(x[0]*(1 - x[0])*x[1]*(1 - x[1]))"
    exprs = [["%s*(4 + 2*sin(2*pi*x[1])) + 30*%s*%s" % (b, b, b), "%s*(1.2*cos(2*pi*x[0]) + x[1])" % b],
             ["%s*(3 - cos(2*pi*x[0]))" % b, "%s*(0.7*sin(2*pi*x[1]) - 2*x[0] + 5*x[1]*x[1])" % b]]
    return Case(df, tmp_path_factory.mktemp("upgrade_box2"), X, cells, (False, False), exprs,
                "triangle", (0.5, 0.5))


@pytest.fixture(scope="module")
def box3(tmp_path_factory):
    """box_mesh's unit cube, 4 a side, walled on every side, a P2 field at rest on them."""
    df = pytest.importorskip("dolfin", reason="writing the checkpoint needs dolfin")
    X, cells = box_mesh(3, 4)
    b = "(x[0]*(1 - x[0])*x[1]*(1 - x[1])*x[2]*(1 - x[2]))"
    exprs = [["%s*(4 + 2*sin(2*pi*x[1]))" % b, "%s*(1.2*cos(2*pi*x[2]) + x[0])" % b,
              "%s*(0.5 + x[1]*x[0])" % b],
             ["%s*(3 - cos(2*pi*x[2]))" % b, "%s*(0.7*sin(2*pi*x[0]) - 2*x[1])" % b,
              "%s*(1 - 2*x[0]*x[2])" % b]]
    return Case(df, tmp_path_factory.mktemp("upgrade_box3"), X, cells, (False, False, False), exprs,
                "tet", (0.5, 0.5, 0.5))


def off_x(P, tol=1e-12):
    """A facet off the sides x = 0 and 1."""
    return not (np.all(np.abs(P[:, 0]) < tol) or np.all(np.abs(P[:, 0] - 1) < tol))


@pytest.fixture(scope="module")
def box2_inlet(tmp_path_factory):
    """box_mesh's unit square, 8 a side, walled on y = 0 and 1, a P2 field in through x = 0 and out through
    x = 1, at rest in the first stamp."""
    df = pytest.importorskip("dolfin", reason="writing the checkpoint needs dolfin")
    X, cells = box_mesh(2, 8)
    b = "(x[1]*(1 - x[1]))"
    exprs = [["0", "0"], ["%s*(4 + 2*sin(2*pi*x[0])) + 30*%s*%s" % (b, b, b), "%s*(1.2*cos(2*pi*x[0]) + x[1])" % b]]
    return Case(df, tmp_path_factory.mktemp("upgrade_box2_inlet"), X, cells, (False, False), exprs,
                "triangle", (0.5, 0.5), wall=off_x)


@pytest.fixture(scope="module")
def box3_inlet(tmp_path_factory):
    """box_mesh's unit cube, 4 a side, walled on the sides off x, a P2 field in through x = 0 and out through
    x = 1, at rest in the first stamp."""
    df = pytest.importorskip("dolfin", reason="writing the checkpoint needs dolfin")
    X, cells = box_mesh(3, 4)
    b = "(x[1]*(1 - x[1])*x[2]*(1 - x[2]))"
    exprs = [["0", "0", "0"],
             ["%s*(4 + 2*sin(2*pi*x[0]))" % b, "%s*(1.2*cos(2*pi*x[2]) + x[0])" % b, "%s*(0.5 + x[1]*x[0])" % b]]
    return Case(df, tmp_path_factory.mktemp("upgrade_box3_inlet"), X, cells, (False, False, False), exprs,
                "tet", (0.5, 0.5, 0.5), wall=off_x)


POISEUILLE = "6*x[1]*(1 - x[1])"


def on_y(P, tol=1e-12):
    """A facet on the sides y = 0 and 1."""
    return bool(np.all(np.abs(P[:, 1]) < tol) or np.all(np.abs(P[:, 1] - 1) < tol))


@pytest.fixture(scope="module")
def poiseuille2(tmp_path_factory):
    """Plane Poiseuille flow along x in box_mesh's unit square, 8 a side, periodic along x."""
    df = pytest.importorskip("dolfin", reason="writing the checkpoint needs dolfin")
    X, cells = box_mesh(2, 8)
    return Case(df, tmp_path_factory.mktemp("poiseuille2"), X, cells, (True, False),
                [[POISEUILLE, "0"]] * 2, "triangle", (0.5, 0.5))


@pytest.fixture(scope="module")
def poiseuille3(tmp_path_factory):
    """Plane Poiseuille flow along x in box_mesh's unit cube, 4 a side, periodic along x and z."""
    df = pytest.importorskip("dolfin", reason="writing the checkpoint needs dolfin")
    X, cells = box_mesh(3, 4)
    return Case(df, tmp_path_factory.mktemp("poiseuille3"), X, cells, (True, False, True),
                [[POISEUILLE, "0", "0"]] * 2, "tet", (0.5, 0.5, 0.5))


@pytest.fixture(scope="module")
def poiseuille2_open(tmp_path_factory):
    """Plane Poiseuille flow in through x = 0 and out through x = 1 of box_mesh's unit square, 8 a
    side, walled on y = 0 and 1 alone."""
    df = pytest.importorskip("dolfin", reason="writing the checkpoint needs dolfin")
    X, cells = box_mesh(2, 8)
    return Case(df, tmp_path_factory.mktemp("poiseuille2_open"), X, cells, (False, False),
                [[POISEUILLE, "0"]] * 2, "triangle", (0.5, 0.5), wall=on_y)


@pytest.fixture(scope="module")
def poiseuille3_open(tmp_path_factory):
    """Plane Poiseuille flow in through x = 0 and out through x = 1 of box_mesh's unit cube, 4 a
    side, walled on y = 0 and 1 alone; the sides z = 0 and 1 carry the flow along them."""
    df = pytest.importorskip("dolfin", reason="writing the checkpoint needs dolfin")
    X, cells = box_mesh(3, 4)
    return Case(df, tmp_path_factory.mktemp("poiseuille3_open"), X, cells, (False, False, False),
                [[POISEUILLE, "0", "0"]] * 2, "tet", (0.5, 0.5, 0.5), wall=on_y)


CASES = ["channel", "sphere", "box2", "box3"]


def interpolate(folder, mode, n):
    """interpol's probes of the case in folder read by mode, at T0, n random
    points: {column: array}."""
    for p in folder.glob("Interpolation*"):
        shutil.rmtree(p)
    run_app(INTERPOL, folder / "dolfin_params.dat",
            "mode=%s Nrw=%d int_order=2 t0=%r random=false seed=1" % (mode, n, T0), timeout=600)
    return interpolation(folder)


def probed(v, dim):
    """(points, u, J) of interpol's columns."""
    x = np.stack([v[a] for a in "xyz"[:dim]], axis=-1)
    u = np.stack([v["u" + a] for a in "xyz"[:dim]], axis=-1)
    J = np.stack([np.stack([v["u" + a + b] for b in "xyz"[:dim]], axis=-1) for a in "xyz"[:dim]], axis=1)
    return x, u, J


def check_against_reference(case, folder, mode, n, weights=(1 - T0, T0)):
    """The upgraded field and its gradient at interpol's points against the
    reference summed with weights; the counts of the kinds of cells
    reached."""
    dim = case.X.shape[1]
    x, u, J = probed(interpolate(folder, mode, n), dim)
    cid, L, depth = locate(case.X, case.cells, x)
    assert depth.min() > -1e-12
    ur, Jr = case.reference(cid, L, weights=weights)
    _, J2 = case.reference(cid, L, terms=False, weights=weights)
    us, Js = np.abs(ur).max(), np.abs(Jr).max()
    err_u = (np.abs(u - ur).max(axis=1) / us).max()
    err_J = (np.abs(J - Jr).max(axis=(1, 2)) / Js).max()
    assert err_u < 1e-9, err_u
    assert err_J < 1e-9, err_J
    # the correction is there to be checked: it moves the wall cells' gradient
    rule = case.rule
    wall, split = rule.wall_cell[cid], rule.split_cell[cid]
    assert np.abs(J - J2)[wall].max() > 1e-3 * Js
    assert np.abs(J - J2)[~wall].max() < 1e-12 * Js
    facets = rule.facets.sum(axis=1)[cid]
    return dict(facet=int((facets > 0).sum()), two_facets=int((facets == 2).sum()),
                no_facet=int((split & (facets == 0)).sum()), vertex=int((wall & ~split).sum()),
                err_u=err_u, err_J=err_J)


@pytest.mark.parametrize("name", CASES)
def test_the_stamped_loader_matches_the_reference(request, name):
    """mode=triangle/tet: u and grad u to 1e-9 of their scale at
    every probe, in cells with one wall facet, on the boxes with two, in
    cells with a wall edge (3D) or vertex (2D) but no wall facet, and in 3D
    in cells with a wall vertex alone."""
    case = request.getfixturevalue(name)
    reached = check_against_reference(case, case("true"), case.mode, 20000)
    assert reached["facet"] > 100 and reached["no_facet"] > 100, reached
    if case.X.shape[1] == 3:
        assert reached["vertex"] > 100, reached
    if name.startswith("box"):
        assert reached["two_facets"] > 10, reached


@pytest.mark.parametrize("name", CASES)
def test_the_frequency_loader_upgrades_the_summed_components(request, name):
    """wall_upgrade=true in mode=trianglefreq/tetfreq: the
    stamps' two files read as components, one steady and one oscillating,
    summed with their weights a cos(omega t + phi) and corrected by the rule
    a stamp is. The rule is linear in the field and both components are at
    rest on the same walls, so that is the reference's correction of each
    stamp summed with the same weights: u and grad u to 1e-9 of their
    scale."""
    case = request.getfixturevalue(name)
    lines = [(0.0, 0.0, 0.7), (2.0, 0.3, -1.3)]
    d = case.variant("freq", "t_min=0\nt_max=10\nwall_upgrade=true\n", series="freqstamps=freqstamps.dat")
    (d / "freqstamps.dat").write_text("".join("%r %r %r up_%d.h5\n" % (w, f, a, k)
                                              for k, (w, f, a) in enumerate(lines)))
    weights = [a * np.cos(w * T0 + f) for w, f, a in lines]
    reached = check_against_reference(case, d, case.mode + "freq", 20000, weights=weights)
    assert reached["facet"] > 100 and reached["no_facet"] > 100, reached


@pytest.mark.parametrize("name", CASES)
def test_off_the_wall_cells_the_field_is_plain_p2_bit_for_bit(request, name):
    """wall_upgrade=true and false read the same points to the last bit
    outside the wall cells, those with a wall vertex."""
    case = request.getfixturevalue(name)
    dim = case.X.shape[1]
    a, b = interpolate(case("true"), case.mode, 5000), interpolate(case("false"), case.mode, 5000)
    assert all((a[k] == b[k]).all() for k in "xyz")
    x, _, _ = probed(a, dim)
    cid, _, _ = locate(case.X, case.cells, x)
    off = ~case.rule.wall_cell[cid]
    nv = dim + 1
    wall_vertices = np.unique([case.cells[c][j] for c, k in case.slots for j in range(nv) if j != k])
    on_wall = np.isin(case.cells, wall_vertices).any(axis=1)
    assert off.sum() > 300
    assert (off & on_wall[cid]).sum() == 0
    for k in a:
        assert (a[k][off] == b[k][off]).all(), k
    assert any((a[k][~off] != b[k][~off]).any() for k in a)


@pytest.mark.parametrize("name", CASES)
def test_wall_facets_naming_every_wall_reads_as_the_walls_at_rest(request, name):
    """wall_facets listing every exterior facet off the periodic sides, where
    the field is at rest, reads the same points as the walls found at rest, to
    the last bit."""
    case = request.getfixturevalue(name)
    d = case.variant("facets", "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=1:7\n")
    a, b = interpolate(case("true"), case.mode, 5000), interpolate(d, case.mode, 5000)
    for k in a:
        assert (a[k] == b[k]).all(), k


@pytest.mark.parametrize("series", ["stamps", "freq"])
@pytest.mark.parametrize("name", ["box2_inlet", "box3_inlet"])
def test_wall_facets_leave_an_inlet_and_an_outlet_out(request, name, series):
    """wall_facets with wall_markers selecting the sides off x, the inlet and
    outlet on x = 0 and 1 no walls, the first stamp (component) at rest: the
    stamped and the frequency loader match the reference whose walls are
    those sides alone, u and grad u to 1e-9 of their scale, the field plain
    P2 outside the cells with a vertex on them."""
    case = request.getfixturevalue(name)
    extra = "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=3:6\n"
    if series == "stamps":
        reached = check_against_reference(case, case.variant("facets", extra), case.mode, 20000)
    else:
        lines = [(0.0, 0.0, 0.7), (2.0, 0.3, -1.3)]
        d = case.variant("facets_freq", "t_min=0\nt_max=10\n" + extra, series="freqstamps=freqstamps.dat")
        (d / "freqstamps.dat").write_text("".join("%r %r %r up_%d.h5\n" % (w, f, a, k)
                                                  for k, (w, f, a) in enumerate(lines)))
        weights = [a * np.cos(w * T0 + f) for w, f, a in lines]
        reached = check_against_reference(case, d, case.mode + "freq", 20000, weights=weights)
    assert reached["facet"] > 100 and reached["no_facet"] > 100, reached


def wall_facets(case):
    """(foot points, inward normals, |grad l_k|) at three points of every wall facet."""
    X, cells = case.X, case.cells
    dim, nv = X.shape[1], X.shape[1] + 1
    pts, normals, gk = [], [], []
    w = [(0.5, 0.5), (0.2, 0.8), (0.7, 0.3)] if dim == 2 else \
        [(1 / 3, 1 / 3, 1 / 3), (0.6, 0.2, 0.2), (0.2, 0.2, 0.6)]
    for c, k in case.slots:
        fv = [cells[c][j] for j in range(nv) if j != k]
        g = case.rule.G[c, k]
        for b in w:
            pts.append(np.array(b) @ X[fv])
            normals.append(g / np.linalg.norm(g))
            gk.append(np.linalg.norm(g))
    return np.array(pts), np.array(normals), np.array(gk)


@pytest.mark.parametrize("name", CASES)
def test_the_normal_velocity_grows_as_the_squared_wall_distance(request, name):
    """Above every wall facet, those of cells with two included, u.n /
    delta^2 is the same at l_k = 1e-7, 1e-6 and 1e-5 to a few per cent with
    wall_upgrade=true, so d(u.n)/dn is zero on the facet; P2 alone has u.n ~
    delta, the ratio growing tenfold a step. No slip: at l_k = 1e-11 the
    speed is below 1e-5 of that at 1e-5, as a field vanishing on the wall has
    it."""
    case = request.getfixturevalue(name)
    foot, n, gk = wall_facets(case)
    assert len(foot) == 3 * len(case.slots)
    levels = np.array([1e-11, 1e-7, 1e-6, 1e-5])
    deltas = levels[None, :] / gk[:, None]
    pts = (foot[:, None, :] + deltas[:, :, None] * n[:, None, :]).reshape(-1, foot.shape[1])
    normals = np.repeat(n, len(levels), axis=0)
    ratio = {}
    for mode in ("true", "false"):
        u = probe(case, mode, pts)
        un = np.sum(u * normals, axis=1).reshape(-1, len(levels))
        ratio[mode] = un[:, 1:] / deltas[:, 1:] ** 2
        speed = np.linalg.norm(u, axis=1).reshape(-1, len(levels))
        assert (speed[:, 0] < 1e-5 * speed[:, -1]).all()
    r = ratio["true"]
    assert (np.abs(r - r[:, -1:]).max(axis=1) <= 0.05 * np.abs(r).max()).all()
    r = ratio["false"]
    assert np.median(np.abs(r[:, 0]) / np.maximum(np.abs(r[:, -1]), 1e-300)) > 50


def probe_gradient(case, mode, points):
    """The velocity gradient at each point, in the order given, with
    wall_upgrade=mode: tracervectors, which evaluates it, dumping J, in a
    folder of its own."""
    d = case.variant("gradient_" + mode, "wall_upgrade=%s\n" % mode)
    folder, _ = from_points(TRACERVECTORS, d / "dolfin_params.dat", points, DT,
                            ["mode=" + case.mode, "Dm=0", "scheme=RK4"],
                            ["T=%g" % (2 * DT), "dump_intv=%g" % (2 * DT), "stat_intv=1000", "output_J=true"],
                            seed_point=case.seed_point, timeout=1200)
    ids = read_checkpoint(folder)["id"][:, 0].astype(int)
    g = dump_at(folder, 2 * DT)
    dim = points.shape[1]
    assert np.allclose(g["points"][ids][:, :dim], points, rtol=0, atol=1e-12)
    return g["J"][ids].reshape(-1, 3, 3)[:, :dim, :dim]


@pytest.mark.parametrize("name", CASES)
def test_the_divergence_grows_as_the_squared_wall_distance(request, name):
    """Above every wall facet, those of cells with two included, at l_k =
    1e-3, 2e-3 and 4e-3 (in the sub-cell of the facet), div u / l_k^2 is
    linear in l_k to 1e-8 of its scale with wall_upgrade=true: div u, cubic
    there, and its normal derivative vanish on the facet. P2 alone has div u
    = O(1) there, far above the upgraded field's at the lowest height."""
    case = request.getfixturevalue(name)
    foot, n, gk = wall_facets(case)
    levels = np.array([1e-3, 2e-3, 4e-3])
    heights = levels[None, :] / gk[:, None]
    pts = (foot[:, None, :] + heights[:, :, None] * n[:, None, :]).reshape(-1, foot.shape[1])
    div = {}
    for mode in ("true", "false"):
        J = probe_gradient(case, mode, pts)
        div[mode] = np.trace(J, axis1=1, axis2=2).reshape(-1, len(levels))
    r = div["true"] / levels ** 2
    scale = np.abs(r).max()
    assert scale > 0
    assert np.abs(r[:, 2] - 3 * r[:, 1] + 2 * r[:, 0]).max() <= 1e-8 * scale
    assert np.median(np.abs(div["false"][:, 0])) > 100 * np.abs(div["true"][:, 0]).max()


def test_wall_p2_and_divfree_are_refused_naming_wall_upgrade(channel, tmp_path):
    """wall_p2 is gone without an alias: a file holding it fails and names its
    replacement. The divergence-free split replaces the upgrade, so the two
    together are refused."""
    d = channel("false")
    for extra, word in (("wall_p2=edge\n", "wall_upgrade"), ("wall_upgrade=true\ndivfree=true\n", "wall_upgrade")):
        c = tmp_path / ("c%d" % len(extra))
        shutil.copytree(d, c)
        (c / "dolfin_params.dat").write_text(channel.params + extra)
        r = run_app(INTERPOL, c / "dolfin_params.dat", "mode=triangle Nrw=10 int_order=1 t0=0",
                    check=False)
        assert r.returncode != 0
        assert word in r.stdout + r.stderr, r.stdout + r.stderr


def near_wall_points(n, seed=5):
    """Points of the channel in the layer of cells above either wall."""
    rng = np.random.default_rng(seed)
    x = rng.uniform(0, 1, n)
    d = rng.uniform(0.01, 0.1, n)
    y = np.where(np.arange(n) % 2, kink_wall(x) + d, 1 - d)
    return np.stack([x, y], axis=1)


def ends(case, mode, points, args, scheme, dt):
    """Positions at the end of a restarted run, by id."""
    folder = restart_from(case, mode, points, args, scheme=scheme, dt=dt)
    T = float(dict(a.split("=") for a in args)["T"])
    return dump_at(folder, T)["points"][:, :2]


def test_rk4cells_ends_as_close_to_rk4_as_without_the_upgrade(channel):
    """RK4cells splits a step at the facets it meets, where the field's
    gradient jumps; in a cell the two schemes take the same RK4 step. Tracers
    above the channel's walls, stepped to T = 0.9 under both, end as close
    together with the upgrade as without it, within a factor of ten: the
    upgraded field is continuous across facets, as the P2 field is."""
    points = near_wall_points(200)
    dt, T = 0.02, 0.9
    gap = {}
    for mode in ("true", "false"):
        args = ["T=%g" % T, "dump_intv=%g" % T, "stat_intv=1000"]
        a = ends(channel, mode, points, args, "RK4", dt)
        b = ends(channel, mode, points, args, "RK4cells", dt)
        gap[mode] = np.abs(a - b).max()
    assert 0 < gap["true"] < 10 * gap["false"], gap


def test_frozen_fields_between_stamps_carry_the_correction(channel):
    """frozen_fields at a time between the stamps upgrades the blended
    velocity, the walls at rest in both: over ten steps of dt = 1e-3, frozen
    and live tracers above the walls part only by what the field's change in
    time moves them, as without the upgrade. A correction missing from either
    would part them by what it moves them, over ten times as far."""
    points = near_wall_points(200)
    dt = 1e-3
    t0 = 2 * dt   # the restart's time
    T = t0 + 10 * dt
    gap = {}
    for mode in ("true", "false"):
        args = ["T=%r" % T, "dump_intv=%r" % dt, "stat_intv=1000"]
        live = ends(channel, mode, points, args, "RK4", dt)
        frozen = ends(channel, mode, points, args + ["frozen_fields=true", "t_frozen=%r" % t0], "RK4", dt)
        gap[mode] = np.abs(live - frozen).max()
    assert gap["false"] > 0 and 0.5 < gap["true"] / gap["false"] < 2, gap


def poiseuille_errors(folder, mode, dim, n):
    """The largest error of u and grad u against plane Poiseuille flow at interpol's points,
    relative to their scales, and how many points are within a cell of a wall."""
    x, u, J = probed(interpolate(folder, mode, 5000), dim)
    y = x[:, 1]
    ue = np.zeros_like(u)
    ue[:, 0] = 6 * y * (1 - y)
    Je = np.zeros_like(J)
    Je[:, 0, 1] = 6 * (1 - 2 * y)
    near = int(((y < 1 / n) | (y > 1 - 1 / n)).sum())
    return np.abs(u - ue).max() / 1.5, np.abs(J - Je).max() / 6, near


@pytest.mark.parametrize("name,n", [("poiseuille2", 8), ("poiseuille3", 4)])
@pytest.mark.parametrize("walls", ["at rest", "wall_facets"])
def test_poiseuille_flow_is_carried_exactly(request, name, n, walls):
    """Plane Poiseuille flow is quadratic, divergence-free and zero on the walls, so P2 holds
    it exactly; the wall terms and the correction must then vanish, and the upgraded field
    and its gradient are the flow's to round-off, in the cells at the walls as everywhere,
    whether the walls are found at rest or named by wall_facets."""
    case = request.getfixturevalue(name)
    if walls == "at rest":
        d = case("true")
    else:
        d = case.variant("facets", "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=3:4\n")
    err_u, err_J, near = poiseuille_errors(d, case.mode, case.X.shape[1], n)
    assert near > 500, near
    assert err_u < 1e-12 and err_J < 1e-12, (err_u, err_J)


@pytest.mark.parametrize("name,n", [("poiseuille2_open", 8), ("poiseuille3_open", 4)])
def test_poiseuille_flow_through_an_inlet_is_carried_exactly(request, name, n):
    """The same flow entering and leaving through faces that are no walls, wall_facets naming
    y = 0 and 1 alone: the cells where a wall meets the inlet or the outlet, and in 3D the
    open sides along the flow, carry it exactly too."""
    case = request.getfixturevalue(name)
    d = case.variant("facets", "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=3:4\n")
    err_u, err_J, near = poiseuille_errors(d, case.mode, case.X.shape[1], n)
    assert near > 500, near
    assert err_u < 1e-12 and err_J < 1e-12, (err_u, err_J)
