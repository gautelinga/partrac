"""wall_upgrade=true on a P2 velocity (mode=triangle and tet): the P2 field
plus quartic terms in the cells with a wall vertex at rest.

The terms vanish at every P2 node. Per edge end at a wall vertex w, towards v:
-((4 u_m - u_v).n) n lambda_w^2 lambda_v (lambda_w - lambda_v), n the edge
end's normal (the facet's own where the edge rises from one facet, else the
mean of the facets there; none at a corner). Per face (w_i, w_j, v) over a
wall edge, lambda_i lambda_j lambda_v (f_i lambda_i + f_j lambda_j) with
f_i = -P(2 s_i + s_j), f_j = -P(s_i + 2 s_j), s the slopes 4 u_m - u_v along
w -> v and P the face's normal projector. In 3D, per cell with wall facets,
the bubble d lambda_1 lambda_2 lambda_3 lambda_4 with d the least-norm
solution of n_k.d = -2 n_k.(s_1 + s_2 + s_3) over its facets k. Together they
cancel the normal slope of u.n on a flat wall facet, so u.n grows as the
square of the wall distance, and they depend on each entity's own data, so
the field stays continuous.

The reference here is a numpy evaluation of that construction, entity by
entity, with the same normal rules. The loader must reproduce it, value and
gradient, at random points (interpol), in 2D on a kinked channel periodic
along it and in 3D around a sphere in a periodic box: cells with a wall facet,
a wall edge or a wall vertex only, on the periodic seams, at a time between
two stamps. Away from the wall cells the field is the plain P2 one bit for
bit.
"""

import os
import shutil

import numpy as np
import pytest

from dumps import dump_at
from paths import app
from runs import run_app
from wall_cases import (dolfin_mesh, kink_wall, kinked_channel_mesh, need_gmsh, periodic_sphere_mesh, probe,
                        restart_from)

INTERPOL = app("interpol")
TRACERS = app("tracers")

pytestmark = pytest.mark.skipif(not (os.path.exists(INTERPOL) and os.path.exists(TRACERS)),
                                reason="interpol and tracers are not built")

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


def periodic_classes(X, periodic, tol=1e-9):
    """One class per vertex, images on a periodic seam merged by position."""
    lo, hi = X.min(axis=0), X.max(axis=0)
    Y = X.copy()
    for d in range(X.shape[1]):
        if periodic[d]:
            Y[np.abs(X[:, d] - hi[d]) < tol, d] = lo[d]
    _, inv = np.unique(np.round((Y - lo) / tol).astype(np.int64), axis=0, return_inverse=True)
    return inv.ravel()


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


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else None


class Upgrade:
    """The wall geometry of a mesh and the quartic terms of a P2 field on it."""

    def __init__(self, X, cells, slots, vclass):
        self.X, self.cells, self.K = X, cells, vclass
        self.dim = X.shape[1]
        nv = self.nv = self.dim + 1
        # the facets in the order of their sorted vertex ids, outward normals
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
        self.count, self.corner, self.n_sum, first = {}, set(), {}, {}
        self.side, self.edge, self.face = {}, {}, {}
        for fv, c, k, n in facets:
            apex = K[cells[c][k]]
            cls = [K[v] for v in fv]
            for w in cls:
                if w not in self.count:
                    self.count[w], first[w], self.n_sum[w] = 0, n, np.zeros(self.dim)
                elif (self.dim == 2 and self.count[w] > 1) or first[w] @ n < 0.5:
                    self.corner.add(w)
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
        if w not in self.count or w in self.corner:
            return None
        return unit(self.side.get((w, v), self.n_sum[w]))

    def face_normal(self, a, b, v):
        e = (min(a, b), max(a, b))
        if e not in self.edge or a in self.corner or b in self.corner:
            return None
        return unit(self.face.get(e + (v,), self.edge[e]))

    def coefficients(self, U2):
        """Per cell: edge ends (nedges, 2, dim), face monomials (nfaces, 4, dim)
        as lambda_a lambda_b lambda_c times (1, lambda_a, lambda_b, lambda_c),
        and the 3D bubble (dim,)."""
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
                    if any(K[row[j]] in self.corner for j in range(nv) if j != k):
                        continue
                    n = self.facet_n[(c, k)]
                    rows.append((n, -2.0 * float(n @ sum(s(i, k) for i in range(nv) if i != k))))
                if rows:
                    D[c] = np.linalg.pinv(np.array([r[0] for r in rows])) @ np.array([r[1] for r in rows])
        return C, F, D

    def weights(self, L):
        """(edge (n, nedges, 2), face (n, nfaces, 4), bubble (n,)) and their
        lambda derivatives, one more axis of nv."""
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

    def evaluate(self, U2, CFD, cid, L, terms=True):
        """(u (n, dim), J (n, dim, dim), J[p, i, j] = du_i/dx_j) at barycentrics L of cells cid."""
        G = bary_gradients(self.X, self.cells[cid])
        u = np.einsum("pn,pnd->pd", p2_basis(L), U2[cid])
        J = np.einsum("pnm,pmx,pnd->pdx", p2_dbasis(L), G, U2[cid])
        if terms:
            C, F, D = CFD
            (E, Fw, Bw), (dE, dF, dBw) = self.weights(L)
            u = u + np.einsum("pea,pead->pd", E, C[cid]) + np.einsum("pfk,pfkd->pd", Fw, F[cid]) \
                + Bw[:, None] * D[cid]
            J = J + np.einsum("peam,pmx,pead->pdx", dE, G, C[cid]) \
                + np.einsum("pfkm,pmx,pfkd->pdx", dF, G, F[cid]) \
                + np.einsum("pm,pmx,pd->pdx", dBw, G, D[cid])
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


# ------------------------------------------------------------------ cases


class Case:
    """A P2 checkpoint case, two stamps at t = 0 and 1, read with
    wall_upgrade=mode; case.U2[k] are stamp k's P2 values per cell in dolfin's
    local order."""

    def __init__(self, df, base, X, cells, periodic, exprs, mode, seed_point):
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
        (base / "timestamps.dat").write_text("0.0 up_0.h5\n1.0 up_1.h5\n")
        self.params = ("velocity_space=P2\npressure_space=P1\ntimestamps=timestamps.dat\nmesh=mesh.h5\n"
                       + "".join("periodic_%s=%s\n" % (a, "true" if p else "false")
                                 for a, p in zip("xyz", list(periodic) + [False] * (3 - dim))))
        self.up = Upgrade(X, cells, self.slots, periodic_classes(X, periodic))

    def __call__(self, mode):
        """The case folder read with wall_upgrade=mode ("true" or "false")."""
        d = self.base / mode
        if not d.exists():
            d.mkdir()
            for f in os.listdir(self.base):
                if os.path.isfile(self.base / f):
                    shutil.copy(self.base / f, d / f)
            (d / "dolfin_params.dat").write_text(self.params + "wall_upgrade=%s\n" % mode)
        return d

    def reference(self, cid, L, terms=True):
        """u and J at time T0, from the two stamps."""
        out = []
        for U2 in self.U2:
            CFD = self.up.coefficients(U2) if terms else None
            out.append(self.up.evaluate(U2, CFD, cid, L, terms=terms))
        return tuple((1 - T0) * a + T0 * b for a, b in zip(*out))


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
def sphere(tmp_path_factory):
    """The sphere in the periodic unit cube of test_wall_edge_rule.py, a P2 field at rest on the sphere."""
    df = pytest.importorskip("dolfin", reason="writing the checkpoint needs dolfin")
    X, cells = periodic_sphere_mesh(need_gmsh())
    d = "(sqrt(pow(x[0]-0.5, 2) + pow(x[1]-0.5, 2) + pow(x[2]-0.5, 2)) - 0.25)"
    exprs = [["%s*(1 + 0.5*sin(2*pi*x[1]))" % d, "%s*(0.6*cos(2*pi*x[2]) + 2*%s)" % (d, d),
              "%s*0.4*sin(2*pi*x[0])" % d],
             ["%s*(0.8 - 0.3*cos(2*pi*x[2]))" % d, "%s*0.5*sin(2*pi*x[0])" % d,
              "%s*(0.9*cos(2*pi*x[1]) - %s)" % (d, d)]]
    return Case(df, tmp_path_factory.mktemp("upgrade_sphere"), X, cells, (True, True, True), exprs,
                "tet", (0.5, 0.1, 0.5))


def interpolate(case, mode, n):
    """interpol's probes of case(mode) at T0, n random points: {column: array}."""
    import h5py
    d = case(mode)
    for p in d.glob("Interpolation*"):
        shutil.rmtree(p)
    run_app(INTERPOL, d / "dolfin_params.dat",
            "mode=%s Nrw=%d int_order=2 t0=%r random=false seed=1" % (case.mode, n, T0), timeout=600)
    [out] = list(d.rglob("interpolation.h5part"))
    with h5py.File(out, "r") as h:
        return {k: np.array(h["Step#0"][k]) for k in h["Step#0"]}


def probed(v, dim):
    """(points, u, J) of interpol's columns."""
    x = np.stack([v[a] for a in "xyz"[:dim]], axis=-1)
    u = np.stack([v["u" + a] for a in "xyz"[:dim]], axis=-1)
    J = np.stack([np.stack([v["u" + a + b] for b in "xyz"[:dim]], axis=-1) for a in "xyz"[:dim]], axis=1)
    return x, u, J


def check_against_reference(case, n):
    """The upgraded field and its gradient at interpol's points against the
    reference, in the wall cells; the counts of the kinds of wall cells reached."""
    dim = case.X.shape[1]
    x, u, J = probed(interpolate(case, "true", n), dim)
    cid, L, depth = locate(case.X, case.cells, x)
    assert depth.min() > -1e-12
    wall = case.up.wall_cell[cid]
    ur, Jr = case.reference(cid, L)
    u2, J2 = case.reference(cid, L, terms=False)
    us, Js = np.abs(ur).max(), np.abs(Jr).max()
    err_u = np.abs(u - ur).max(axis=1) / us
    err_J = np.abs(J - Jr).max(axis=(1, 2)) / Js
    assert err_u.max() < 1e-12, err_u.max()
    assert err_J.max() < 1e-12, err_J.max()
    # the terms are there to be checked: they move the wall cells' gradient
    assert np.abs(J - J2)[wall].max() > 1e-2 * Js
    assert np.abs(J - J2)[~wall].max() < 1e-12 * Js
    # which wall cells: by their wall vertices, and whether one is on a seam
    K = case.up.K
    nwall = np.array([sum(K[v] in case.up.count for v in row) for row in case.cells])
    facet = np.zeros(len(case.cells), bool)
    facet[case.slots[:, 0]] = True
    images = np.bincount(K)[K] > 1
    seam_wall = np.array([any(images[v] and K[v] in case.up.count for v in row) for row in case.cells])
    return dict(points=int(wall.sum()), facet=int(facet[cid].sum()),
                vertex_only=int((nwall[cid] == 1).sum()),
                edge_only=int(((nwall[cid] == 2) & ~facet[cid]).sum()),
                seam=int(seam_wall[cid].sum()))


def test_the_channel_matches_the_reference(channel):
    """2D: u and grad u to 1e-12 of their scale at every probe in a wall cell,
    including the cells on the periodic seam and those touching the wall at a
    vertex only."""
    reached = check_against_reference(channel, 20000)
    assert reached["facet"] > 100 and reached["vertex_only"] > 100 and reached["seam"] > 10, reached


def test_the_sphere_matches_the_reference(sphere):
    """3D: the same, around a sphere meshed by gmsh in a periodic box: cells
    with a wall facet (and its bubble), a wall edge only, a wall vertex only."""
    reached = check_against_reference(sphere, 30000)
    assert reached["facet"] > 100 and reached["edge_only"] > 100 and reached["vertex_only"] > 100, reached


@pytest.mark.parametrize("name", ["channel", "sphere"])
def test_off_the_wall_cells_the_field_is_plain_p2_bit_for_bit(request, name):
    """wall_upgrade=true and false read the same points to the last bit
    outside the cells with a wall vertex."""
    case = request.getfixturevalue(name)
    dim = case.X.shape[1]
    a, b = interpolate(case, "true", 5000), interpolate(case, "false", 5000)
    assert all((a[k] == b[k]).all() for k in "xyz")
    x, _, _ = probed(a, dim)
    cid, _, _ = locate(case.X, case.cells, x)
    off = ~case.up.wall_cell[cid]
    assert off.sum() > 1000
    for k in a:
        assert (a[k][off] == b[k][off]).all(), k
    assert any((a[k][~off] != b[k][~off]).any() for k in a)


def single_facets(case):
    """(foot points, inward normals) at three points of each wall facet whose
    side edges rise from it alone, so its terms all use its own normal, and
    whose apex is off the wall, as the construction's slope cancellation
    assumes."""
    X, cells, K = case.X, case.cells, case.up.K
    dim, nv = X.shape[1], X.shape[1] + 1
    apexes = {}
    for c, k in case.slots:
        for j in range(nv):
            if j != k:
                apexes[(K[cells[c][j]], K[cells[c][k]])] = apexes.get((K[cells[c][j]], K[cells[c][k]]), 0) + 1
    pts, normals = [], []
    w = [(0.5, 0.5), (0.2, 0.8), (0.7, 0.3)] if dim == 2 else \
        [(1 / 3, 1 / 3, 1 / 3), (0.6, 0.2, 0.2), (0.2, 0.2, 0.6)]
    for c, k in case.slots:
        fv = [cells[c][j] for j in range(nv) if j != k]
        if any(apexes[(K[v], K[cells[c][k]])] > 1 for v in fv):
            continue
        if any(K[v] in case.up.corner for v in fv) or K[cells[c][k]] in case.up.count:
            continue
        for b in w:
            pts.append(np.array(b) @ X[fv])
            normals.append(-case.up.facet_n[(c, k)])
    return np.array(pts), np.array(normals)


@pytest.mark.parametrize("name", ["channel", "sphere"])
def test_the_normal_velocity_grows_as_the_squared_wall_distance(request, name):
    """Above a wall facet whose terms use its own normal, u.n / delta^2 is the
    same at delta = 1e-5, 1e-4 and 1e-3 to a few per cent with
    wall_upgrade=true, so d(u.n)/dn is zero on the facet; P2 alone has u.n ~
    delta, the ratio growing tenfold a step. No slip: at delta = 1e-9 the speed
    is below 1e-5 of that at 1e-3, as a field vanishing on the wall has it."""
    case = request.getfixturevalue(name)
    foot, n = single_facets(case)
    assert len(foot) > 30
    deltas = np.array([1e-9, 1e-5, 1e-4, 1e-3])
    pts = (foot[:, None, :] + deltas[None, :, None] * n[:, None, :]).reshape(-1, foot.shape[1])
    normals = np.repeat(n, len(deltas), axis=0)
    ratio = {}
    for mode in ("true", "false"):
        u = probe(case, mode, pts)
        un = np.sum(u * normals, axis=1).reshape(-1, len(deltas))
        ratio[mode] = un[:, 1:] / deltas[1:] ** 2
        speed = np.linalg.norm(u, axis=1).reshape(-1, len(deltas))
        assert (speed[:, 0] < 1e-5 * speed[:, -1]).all()
    r = ratio["true"]
    assert (np.abs(r - r[:, -1:]).max(axis=1) <= 0.05 * np.abs(r).max()).all()
    r = ratio["false"]
    assert np.median(np.abs(r[:, 0]) / np.maximum(np.abs(r[:, -1]), 1e-300)) > 50


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
    terms are continuous across facets, as the P2 field is, and RK4cells
    evaluates them as RK4 does."""
    points = near_wall_points(200)
    dt, T = 0.02, 0.9
    gap = {}
    for mode in ("true", "false"):
        args = ["T=%g" % T, "dump_intv=%g" % T, "stat_intv=1000"]
        a = ends(channel, mode, points, args, "RK4", dt)
        b = ends(channel, mode, points, args, "RK4cells", dt)
        gap[mode] = np.abs(a - b).max()
    assert 0 < gap["true"] < 10 * gap["false"], gap


def test_frozen_fields_between_stamps_carry_the_terms(channel):
    """frozen_fields at a time between the stamps upgrades the blended
    velocity, the walls at rest in both: over ten steps of dt = 1e-3, frozen
    and live tracers above the walls part only by what the field's change in
    time moves them, as without the upgrade. Terms missing from either would
    part them by what the terms move them, over ten times as far."""
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
