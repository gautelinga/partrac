"""python/divfree/divfree_clean.py: the flux equilibration and the per-cell
reconstruction, on arrays built here. What they check:

- the construction: every cell's net flux is zero afterwards, no held or seam
  node moves and the flux through an all-held facet is the data's, a second
  cleaning changes nothing under either objective, and a divergence-free
  quadratic survives both steps;
- the held set: a node at rest on an exterior facet the mesh does not pair, and
  nothing else, so a slow interior node stays free and a boundary at rest in one
  component only is not held; --rest-tol sets the tolerance;
- the weights, the penalty and the boundary weight, which is what keeps a moving
  wall's values and an open facet's flux near the data's, and the cases the flux
  rows find hardest: a noisy 2D lid cavity and the regular Kuhn 3D cavity;
- what is reported and never constrained: the volume mean before and after, the
  throughput drift of a periodic direction and the data's net boundary flux;
- the solve: a first solve that does not converge is refused as that, where a
  refinement pass that runs out of its own iterations keeps the step it was
  improving; the matrix handed to PETSc is sorted and within its index range;
- the criterion the tool shares with the loader: the same constants, what a
  cell's net flux is measured against, the floor under a numerically dead cell,
  the refusal when a solve misses the fluxes, and a non-finite value refused
  where it is read.

Slivers are the construction's known risk: cells of aspect ratio 20 and 100 are
cleaned here and the interior values of the reconstruction reported against the
boundary values.
"""

import os

import numpy as np
import pytest
import scipy.sparse as sp

from divfree_cases import (CASES, D, IDS, PATH_IDS, PATHS, all_held_facets,
                           blocked2d, blocked_field, channel2d, channel3d, facet_flux,
                           held_topo, pinch, random_divfree, smooth_noslip, walled,
                           write_case)
from paths import REPO


# ------------------------------------------------------------ the construction


@pytest.mark.parametrize("dim,mesh", CASES, ids=IDS)
def test_the_reference_matrix_is_unique_and_exact(dim, mesh):
    """The local divergence-free problem has as many independent constraints as
    interior dofs, so its solution is the field and not a least-squares fit, and
    it reproduces a divergence-free quadratic exactly."""
    R, info = D.reference_matrix(dim)
    assert info["rank"] == info["ndof"]
    f = random_divfree(dim, np.random.default_rng(4))
    _, B, I = D.reference_nodes(dim)
    got = (R @ f(B).ravel()).reshape(-1, dim)
    assert np.abs(got - f(I)).max() < 1e-12


@pytest.mark.parametrize("dim,mesh", CASES, ids=IDS)
def test_held_nodes_and_the_periodic_seam_are_untouched(dim, mesh):
    """No slip has to stay exact, so a held node must not move by so much as a
    bit, and the gross flux through a facet all of whose nodes are held is the
    data's; and an edge and its periodic image are one unknown, so the seam
    values stay identical and the loaders' 1e-8 seam check has nothing to
    find."""
    t, _ = walled(*mesh, dim)
    rng = np.random.default_rng(1)
    U = D.p1_to_p2(t, rng.normal(size=t.X.shape))
    U[t.held_node] = 0.0
    U = U[t.master]                                  # a periodic input
    Uc, info = D.Equil(t).apply(U)
    assert info["stepped"]
    assert t.held_edge.sum() > 0 and t.facet_boundary.sum() > 0
    assert np.abs(Uc[t.held_node] - U[t.held_node]).max() == 0.0
    assert np.abs(Uc[t.master] - Uc).max() == 0.0
    sel = all_held_facets(t)
    assert sel.sum() > 0
    a, b = facet_flux(t, U), facet_flux(t, Uc)
    assert np.array_equal(a[sel], b[sel])


@pytest.mark.parametrize("dim,mesh", CASES, ids=IDS)
@pytest.mark.parametrize("g", PATHS, ids=PATH_IDS)
def test_a_second_cleaning_changes_nothing(dim, mesh, g):
    """Cleaning is idempotent under either objective: each is refined to below
    the early exit, so a dataset that has been through the tool is a fixed point
    of it, and re-running it cannot drift the data."""
    t, _ = walled(*mesh, dim)
    rng = np.random.default_rng(2)
    U = D.p1_to_p2(t, rng.normal(size=t.X.shape))
    U[t.held_node] = 0.0
    eq = D.Equil(t, penalty=g)
    U1, i1 = eq.apply(U[t.master])
    U2, i2 = eq.apply(U1)
    assert i1["stepped"] and i1["balance_after"] <= D.FLUX_EXIT * D.FLUX_TOL
    assert not i2["stepped"]
    assert np.array_equal(U2, U1)


def cavity(n=12, noise=0.0):
    """(topology, P2 nodes) of a lid-driven cavity on channel2d's square, closed
    on all four sides: psi = x^2 (1-x)^2 (y^2 - y^3), at rest on three walls,
    the lid y = 1 moving and so free. With noise, white noise of that size times
    |u|max on every node that is not at rest, a solver's discretisation error."""
    X, cells = channel2d(n)
    t = D.Topo(X, cells, [False, False])
    x, y = t.node_x[:, 0], t.node_x[:, 1]
    f, df = x ** 2 * (1 - x) ** 2, 2 * x * (1 - x) * (1 - 2 * x)
    h, dh = y ** 2 - y ** 3, 2 * y - 3 * y ** 2
    U = np.stack([f * dh, -df * h], axis=1)
    if noise:
        scale = float(np.abs(U).max())
        moves = np.linalg.norm(U, axis=1) > D.REST_TOL * scale
        U = U + noise * scale * np.random.default_rng(3).normal(size=U.shape) * moves[:, None]
    return held_topo(X, cells, [False, False], U), U


@pytest.mark.parametrize("g", PATHS, ids=PATH_IDS)
def test_a_second_cleaning_of_a_cavity_changes_nothing(g):
    """The same fixed point where the flux rows are not all a periodic channel's:
    a closed cavity whose moving lid is its only free boundary. The
    smallest-change path must refine too, or its first balance can land above
    the early exit and a second cleaning steps again."""
    t, U = cavity()
    assert t.held_edge.any() and (t.boundary_edge & ~t.held_edge).any()
    eq = D.Equil(t, penalty=g)
    U1, i1 = eq.apply(U)
    U2, i2 = eq.apply(U1)
    assert i1["stepped"] and i1["balance_after"] <= D.FLUX_EXIT * D.FLUX_TOL
    assert not i2["stepped"]
    assert np.array_equal(U2, U1)


@pytest.mark.parametrize("dim,mesh", CASES, ids=IDS)
def test_a_divergence_free_quadratic_survives_both_steps(dim, mesh):
    """Its P2 interpolant is the field itself and already has zero fluxes, so
    the global step must leave it alone and the reconstruction must give it back
    at any point inside a cell, with zero divergence: the construction adds
    nothing of its own."""
    X, cells = mesh
    t = D.Topo(X, cells, [False] * dim)
    f = random_divfree(dim, np.random.default_rng(5))
    U = f(t.node_x)
    r, big = D.cell_flux(t, U)
    assert np.abs(r).max() < 1e-13 * big.max()
    Uc, _ = D.Equil(t).apply(U)
    assert np.abs(Uc - U).max() < 1e-11 * np.abs(U).max()
    uint = D.split_interior(t, Uc)
    rng = np.random.default_rng(6)
    k = rng.integers(0, t.ncells, 20)
    w = rng.dirichlet(np.ones(dim + 1), 20)
    pts = np.einsum('ni,nij->nj', w, t.X[t.cells[k]])
    assert np.abs(D.eval_split(t, Uc, uint, pts, k) - f(pts)).max() < 1e-11 * np.abs(U).max()
    assert np.abs(D.eval_split(t, Uc, uint, pts, k, grad=True)).max() < 1e-10 * np.abs(U).max()


# ---------------------------------------------------------------- the held set


def test_a_boundary_at_rest_in_one_component_only_is_not_held():
    """The held set is decided once a dataset. A steady inlet is at rest in the
    harmonics that vanish there, and holding those would freeze the inlet's edge
    midpoints and leave the flux unbalanced where it matters."""
    X, cells = channel2d(6)
    t = D.Topo(X, cells, [False, False])
    inflow = np.zeros_like(t.node_x)
    inflow[:, 0] = t.node_x[:, 1] * (1 - t.node_x[:, 1])            # nonzero at x = 0 and 1
    harmonic = inflow * np.sin(np.pi * t.node_x[:, 0])[:, None]     # zero there
    nmax = np.zeros(t.nnodes)
    for U in (inflow, harmonic):
        nmax = np.maximum(nmax, np.linalg.norm(U, axis=1))
    scale = max(float(np.abs(U).max()) for U in (inflow, harmonic))
    both = D.Topo(X, cells, [False, False], at_rest=D.at_rest_nodes(nmax, scale))
    one = held_topo(X, cells, [False, False], harmonic)
    # the inlet away from its corners, where the walls make it at rest anyway
    inside = np.isclose(t.node_x[:, 0], 0.0) & (t.node_x[:, 1] > 1e-12) \
        & (t.node_x[:, 1] < 1 - 1e-12)
    at_inlet = lambda topo: int((topo.held_node & inside).sum())
    assert not at_inlet(both), "the inlet was held"
    assert at_inlet(one), "the harmonic alone does look at rest there"
    assert both.held_node.sum() > 0 and both.held_node.sum() < one.held_node.sum()


def test_a_wall_at_solver_round_off_is_held_and_a_slow_interior_node_is_not():
    """An iterative Stokes solve leaves its no-slip values at round-off, not at
    zero, so the whole wall is held. The classification is the only one there
    is: a node as slow as that inside the domain -- a blocked channel's dead
    pocket -- is free, or the cells behind it would carry a net flux no step can
    take away."""
    X, cells = blocked2d(12)
    t0 = D.Topo(X, cells, [True, False])
    exact = blocked_field(t0.node_x)
    scale = float(np.abs(exact).max())
    on_wall = np.isclose(t0.node_x[:, 1], 0.0) | np.isclose(t0.node_x[:, 1], 1.0)
    rng = np.random.default_rng(9)
    # the solver's round-off, spread over four decades around 1e-12
    rough = exact.copy()
    rough[on_wall] = scale * 10.0 ** rng.uniform(-14, -10, (int(on_wall.sum()), 1)) \
        * rng.choice([-1, 1], (int(on_wall.sum()), 2))
    t = held_topo(X, cells, [True, False], rough)
    assert t.held_node[on_wall].all(), "a wall node at round-off was taken for flow"
    slow = np.linalg.norm(rough, axis=1) <= D.REST_TOL * scale
    assert (slow & ~t.boundary_node).any(), "the slit leaves no slow interior node"
    assert not (t.held_node & ~t.boundary_node).any()
    # the same held set, so it is this step the exact field is balanced by
    assert np.array_equal(t.held_node, held_topo(X, cells, [True, False], exact).held_node)
    Uc, info = D.Equil(t).apply(exact)
    assert info["stepped"] and info["balance_after"] <= D.FLUX_TOL
    assert np.abs(Uc[t.held_node] - exact[t.held_node]).max() == 0.0


def stragglers(rest_tol=D.REST_TOL):
    """(channel topology with the held set the data gives, the data): a channel
    whose walls are at round-off but for three nodes of the bottom wall at 3e-8
    of max |u|, as a solve looser there would leave them, above the default rest
    tolerance."""
    X, cells = channel2d(12)
    t0 = D.Topo(X, cells, [True, False])
    U = smooth_noslip(t0.node_x)
    scale = float(np.abs(U).max())
    y = t0.node_x[:, 1]
    on_wall = np.isclose(y, 0.0) | np.isclose(y, 1.0)
    U[on_wall] = 0.0
    few = np.nonzero(on_wall & (y < 0.5))[0][3:6]
    U[few] = [3e-8 * scale, 0.0]
    rest = D.at_rest_nodes(np.linalg.norm(U, axis=1), scale, rest_tol)
    return D.Topo(X, cells, [True, False], at_rest=rest), U


def test_rest_tol_reaches_the_held_set_from_the_command_line(tmp_path):
    """A handful of wall nodes just above the rest tolerance are free, so the
    step may move them; raising the tolerance above them holds them again, and
    the knob has to work end to end."""
    pytest.importorskip("h5py")
    t, U = stragglers()
    tw, _ = stragglers(rest_tol=1e-7)
    assert tw.held_node.sum() > t.held_node.sum()
    cfg = write_case(tmp_path / "in", t.X, t.cells, [U], [True, False])
    assert D.main([str(cfg), "--out", str(tmp_path / "out"), "--rest-tol", "1e-7"]) == 0
    case = D.read_case(tmp_path / "out" / "dolfin_params.dat")
    ct = D.Topo(case["X"], case["cells"], case["periodic"])
    Uc = D.DofTable(tmp_path / "out" / "u_0000.h5", "u", ct, case["cell_indices"]).values(
        tmp_path / "out" / "u_0000.h5")
    assert np.array_equal(Uc[tw.held_node], U[tw.held_node])
    assert D.check_case(tmp_path / "out" / "dolfin_params.dat", quiet=True) < D.FLUX_TOL


# ----------------------------------------------------------------- the weights


@pytest.mark.parametrize("dim,mesh", CASES, ids=IDS)
def test_unweighted_and_volume_weighted_both_balance_the_flux(dim, mesh):
    """The point of the global step: what leaves a cell enters it again, to the
    solver's tolerance, measured against the largest flux through one of its
    facets. --weights changes which midpoints take the change, not whether the
    fluxes end up zero."""
    t, U = walled(*mesh, dim, field=smooth_noslip)
    a, ia = D.Equil(t, weights="volume").apply(U)
    b, ib = D.Equil(t, weights="none").apply(U)
    for info in (ia, ib):
        assert info["flux_before_max"] > 1e-6 * info["facet_flux_max"], "nothing to clean"
        assert info["flux_after_max"] < 1e-10 * info["facet_flux_max"]
    assert np.abs(a - b).max() > 0.0


def test_a_non_positive_boundary_weight_is_refused():
    """The weight multiplies those rows of W and a sparse matrix drops exact
    zeros, so with the penalty off their rows of K are empty and both
    preconditioner blocks are built on a zero diagonal."""
    t, _ = walled(*channel2d(6), 2)
    for w in (0.0, -1.0):
        with pytest.raises(ValueError, match="must be positive"):
            D.Equil(t, boundary_weight=w)


# ------------------------------- the penalty, the boundary weight, the reports


def smallest_change(t, U):
    """The smallest weighted change with zero net cell fluxes, solved densely:
    the step the penalty weight 0 has to reproduce."""
    C, free, col = D.flux_matrix(t)
    d = t.dim
    w = np.bincount((t.master[t.nverts + t.cell_edge] - t.nverts).ravel(),
                    weights=np.repeat(t.vol, t.cell_edge.shape[1]), minlength=t.nedges)
    w = np.maximum(w[free], 1e-300) * np.where(t.boundary_edge[free], D.BOUNDARY_W, 1.0)
    root = np.repeat(np.sqrt(w), d)
    r, _ = D.cell_flux(t, U)
    y = np.linalg.lstsq(np.asarray((C @ sp.diags(1.0 / root)).todense()), -r, rcond=None)[0]
    s = (y / root).reshape(-1, d)
    out = U.copy()
    has = col >= 0
    out[t.nverts:][has] += s[col[has]]
    return out


def th_matrix(t):
    """int phi_i div u over the P1 hat of every vertex, as a matrix in the P2
    nodal values."""
    G = D.p2_vertex_grads(t.dim)
    coef = np.einsum('k,ab,bne,kec->kanc', t.vol, D.mass_p1(t.dim), G, t.Jinv)
    nodes = np.concatenate([t.cells, t.nverts + t.cell_edge], axis=1)
    B = np.zeros((t.nverts, t.nnodes * t.dim))
    for a in range(t.nv):
        for c in range(t.dim):
            np.add.at(B, (t.master[t.cells[:, a]][:, None], nodes * t.dim + c),
                      coef[:, a, :, c])
    return B


def split_field_mean(t, Uc):
    """The volume mean of the final field, integrated on the split mesh itself:
    the split field is an ordinary P2 field there."""
    Xs, cs = D.split_mesh(t)
    st = D.Topo(Xs, cs, [False] * t.dim)
    pos, val = D.split_nodes(t, Uc, D.split_interior(t, Uc))
    return D.volume_mean(st, val[D._order_index(st.node_x, pos)])


@pytest.mark.parametrize("dim,mesh", CASES, ids=IDS)
def test_the_penalty_weight_zero_is_the_smallest_change(dim, mesh):
    """--penalty 0 is the smallest weighted change subject to the fluxes alone,
    so it must come out of the solver as the dense least-norm solution does."""
    t, U = walled(*mesh, dim, field=smooth_noslip)
    a, _ = D.Equil(t, penalty=0.0).apply(U)
    assert np.abs(a - smallest_change(t, U)).max() < 1e-12 * np.abs(U).max()


@pytest.mark.parametrize("dim,mesh", CASES, ids=IDS)
def test_the_penalty_leaves_less_divergence_in_the_cells(dim, mesh):
    """What the reconstruction has to cancel is the divergence the global step
    leaves behind; that is what the penalty is for, so it must be smaller than the
    smallest-change step's on data that is not already divergence-free."""
    t, U = walled(*mesh, dim, field=smooth_noslip)
    plain, _ = D.Equil(t, penalty=0.0).apply(U)
    pen, _ = D.Equil(t).apply(U)
    a = np.sqrt(np.mean(D.div_norms(t, plain) ** 2))
    b = np.sqrt(np.mean(D.div_norms(t, pen) ** 2))
    assert b < (0.98 if dim == 2 else 0.5) * a
    assert np.abs(pen[t.held_node]).max() == 0.0


@pytest.mark.parametrize("dim,mesh", CASES, ids=IDS)
@pytest.mark.parametrize("g", [None, 0.0], ids=["penalised", "smallest"])
def test_the_throughput_drift_is_reported_and_small(dim, mesh, g):
    """With impermeable walls a divergence-free field keeps a uniform tracer
    density uniform, so the mean Lagrangian velocity of the tracers is the volume
    mean of the final field, and in a periodic direction that mean is the
    throughput. Nothing holds it, so the tool reports what the step did to it;
    on these smooth channel fields that is round-off."""
    t, U = walled(*mesh, dim, field=smooth_noslip)
    Uc, info = D.Equil(t, penalty=g).apply(U)
    # the mean of the split field, as the step's linear functional of the macro
    # data and as an integral over the split mesh
    assert np.abs(split_field_mean(t, Uc) - info["mean_after"]).max() < 1e-12 * np.abs(U).max()
    tp = D.throughput_drift(t, info["drift"])
    assert len(tp) == dim - 1
    assert np.abs(info["mean_before"][t.periodic]).max() > 1e-3 * np.abs(U).max()
    assert np.abs(tp).max() < 1e-12
    assert 0.0 <= info["imbalance"] < 1e-10, "a closed channel carries no net flux"


def test_a_boundary_midpoint_moves_less_than_an_interior_one_by_the_weight():
    """A free midpoint on an exterior facet the mesh does not pair costs
    BOUNDARY_W times an interior one, so the step takes its change in the
    interior instead: on an open channel, where the inlet and outlet midpoints
    are free, the largest boundary change falls by about that factor and lands
    decades below the largest interior one, while the fluxes still balance."""
    X, cells = channel2d(8)
    t = D.Topo(X, cells, [False, False])             # walls in y, open in x
    U = smooth_noslip(t.node_x)
    t = held_topo(X, cells, [False, False], U)
    free = t.boundary_edge & ~t.held_edge
    assert free.sum() > 0, "the channel has no free boundary midpoint"
    got, inner = {}, {}
    for w in (1.0, D.BOUNDARY_W):
        Uc, info = D.Equil(t, boundary_weight=w).apply(U)
        assert info["balance_after"] <= D.FLUX_TOL
        mag = np.linalg.norm(info["change"], axis=1)
        got[w], inner[w] = mag[free].max(), mag[~free & ~t.held_edge].max()
    assert got[1.0] / got[D.BOUNDARY_W] > 0.1 * D.BOUNDARY_W
    assert got[1.0] > 0.1 * inner[1.0]
    assert got[D.BOUNDARY_W] < 1e-2 * inner[D.BOUNDARY_W]


@pytest.mark.parametrize("g", PATHS, ids=PATH_IDS)
def test_a_noisy_lid_cavity_cleans_and_reports_its_imbalance(g):
    """A closed cavity whose only moving boundary is its lid, with the noise a
    solver's discretisation error leaves. Its boundary then carries a small net
    flux, which no step can remove and nothing refuses: the lid's free midpoints
    absorb it, the fluxes balance, and the imbalance is reported beside them."""
    t, U = cavity(noise=1e-6)
    Uc, info = D.Equil(t, penalty=g).apply(U)
    assert info["imbalance"] > 1e-9, "the noise leaves no imbalance to report"
    assert info["balance_after"] <= D.FLUX_TOL
    assert np.abs(Uc[t.held_node] - U[t.held_node]).max() == 0.0


def kuhn_cavity(n=4, rim_at_rest=True):
    """(topology, P2 nodes) of a cavity on channel3d's regular Kuhn lattice: at
    rest on five faces, the lid z = 1 moving in x. Its rim is at rest with the
    smooth lid 256 x^2(1-x)^2 y^2(1-y)^2, and with a uniform lid whose rim is
    overwritten to zero, which is the lattice's hardest boundary data."""
    X, cells = channel3d(n)
    t = D.Topo(X, cells, [False] * 3)
    x, y, z = t.node_x[:, 0], t.node_x[:, 1], t.node_x[:, 2]
    f, df = x ** 2 * (1 - x) ** 2, 2 * x * (1 - x) * (1 - 2 * x)
    g = y ** 2 * (1 - y) ** 2
    h, dh = z ** 2 - z ** 3, 2 * z - 3 * z ** 2
    U = 256 * np.stack([-f * g * dh, 0 * x, df * g * h], axis=1)
    if rim_at_rest:
        lid = np.isclose(z, 1.0)
        rim = lid & (np.isclose(x, 0) | np.isclose(x, 1)
                     | np.isclose(y, 0) | np.isclose(y, 1))
        U[lid] = [1.0, 0.0, 0.0]
        U[rim] = 0.0
    return held_topo(X, cells, [False] * 3, U), U


@pytest.mark.parametrize("rim_at_rest", [False, True], ids=["smooth", "rim-at-rest"])
def test_the_regular_kuhn_cavity_cleans(rim_at_rest):
    """A regular Kuhn lattice with every exterior node held has three flux rows
    a jittered mesh does not, and a lid's boundary data does not satisfy them,
    so such a domain cannot be balanced. Only the nodes at rest are held here,
    which leaves the lid free: both lids balance, the five walls at rest do not
    move, and the gross flux through them is the data's."""
    t, U = kuhn_cavity()
    sel = all_held_facets(t)
    assert sel.sum() > 0 and (t.boundary_edge & ~t.held_edge).any()
    Uc, info = D.Equil(t).apply(U)
    assert info["stepped"] and info["balance_after"] <= D.FLUX_TOL
    assert np.abs(Uc[t.held_node] - U[t.held_node]).max() == 0.0
    a, b = facet_flux(t, U), facet_flux(t, Uc)
    assert np.array_equal(a[sel], b[sel])


# ------------------------------------------------------------------- the solve


def test_a_step_that_misses_its_tolerance_refuses_and_names_the_stamp(monkeypatch):
    """A dataset must never be written from a solve that did not reach what a
    loader asks for. With no refinement passes left, a loosened MINRES misses the
    fluxes; the refusal names the stamp and the cell, and nothing is returned."""
    t, U = walled(*channel2d(), 2, field=smooth_noslip)
    monkeypatch.setattr(D, "REFINE_MAX", 0)
    with pytest.raises(RuntimeError, match=r"up_0007\.h5: the flux solve left cell"):
        D.Equil(t, minres_rtol=1e-3).apply(U, "up_0007.h5")


@pytest.mark.parametrize("g", PATHS, ids=PATH_IDS)
def test_a_solve_that_does_not_converge_is_refused_by_name(g, monkeypatch):
    """An iterative solve that stops short, on the iteration limit or on a
    breakdown, must be reported as that, naming the stamp: its output would
    otherwise reach the flux refusal and be blamed on the data."""
    t, U = walled(*CASES[0][1], 2, field=smooth_noslip)
    monkeypatch.setattr(D, "MAXITER", 1)
    eq = D.Equil(t, penalty=g)
    with pytest.raises(RuntimeError, match=r"up_0003\.h5: the KKT solve did not converge"):
        eq.apply(U, "up_0003.h5")


def dependent_rows():
    """A closed channel with a slit wall across it, one half of which is at rest
    all round its boundary and so held all round: that half's cell rows sum to
    zero, while the whole domain's do not, which is the test Equil uses to spot a
    closed domain."""
    X, cells = blocked2d(12)
    t = D.Topo(X, cells, [False, False])
    U = blocked_field(t.node_x)
    U[np.isclose(t.node_x[:, 0], 1.0), 0] = 0.1 * float(np.abs(U).max())
    return held_topo(X, cells, [False, False], U)


def test_the_schur_preconditioner_is_definite_when_the_rows_are_dependent():
    """The flux rows can be dependent in a way the closed-domain test does not
    see. The V-cycle on a singular Schur approximation is then indefinite and
    MINRES breaks down, so the approximation is always shifted."""
    t = dependent_rows()
    eq = D.Equil(t)
    assert eq.solver.startswith("MINRES + GAMG") and not eq.singular
    pc = eq.ksp.getPC().getFieldSplitSubKSP()[1].getPC()
    x, y = pc.getOperators()[0].createVecs()
    cols = []
    for i in range(eq.m):
        x.set(0.0)
        x.setValue(i, 1.0)
        x.assemble()
        pc.apply(x, y)
        cols.append(y.getArray().copy())
    P = np.column_stack(cols)
    assert np.linalg.eigvalsh(0.5 * (P + P.T)).min() > 0.0


def test_the_matrix_handed_to_petsc_is_sorted_and_within_its_index_range():
    """PETSc takes its indices in its own integer type, and a cast past its
    range wraps without a word, so a matrix too large for it is refused and the
    message names the limit. Its AIJ kernels also assume sorted column indices,
    which a sparse product does not leave, so the matrix is sorted in place
    first."""
    import types
    A = sp.random(20, 20, density=0.8, format="csr", random_state=1)
    small = types.SimpleNamespace(IntType=np.int8)
    with pytest.raises(ValueError, match="past the 127 that this PETSc's 8-bit"):
        D._petsc_csr(small, A, D.MPI.COMM_SELF)
    # every row's columns reversed, as a product can leave them
    order = np.concatenate([np.arange(A.indptr[i + 1] - 1, A.indptr[i] - 1, -1)
                            for i in range(A.shape[0])])
    U = sp.csr_matrix((A.data[order], A.indices[order], A.indptr.copy()), shape=A.shape)
    sorted_rows = lambda M: all(np.all(np.diff(M.indices[M.indptr[i]:M.indptr[i + 1]]) > 0)
                                for i in range(M.shape[0]))
    assert not sorted_rows(U)
    wide = types.SimpleNamespace(IntType=np.int64)
    indptr, indices, data = D._petsc_csr(wide, U, D.MPI.COMM_SELF)
    assert sorted_rows(U)
    assert np.array_equal(indptr, A.indptr) and np.array_equal(indices, A.indices)


def test_a_mesh_with_every_midpoint_held_is_refused_as_that():
    """With nothing to solve for, the flux rows are an empty matrix and the
    reductions over them raise a numpy message about a zero-size array; what the
    dataset needs said is that its whole boundary is at rest."""
    X = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    cells = np.array([[0, 1, 2]])
    t = D.Topo(X, cells, [False, False], at_rest=np.ones(6, bool))
    assert t.held_edge.all()
    with pytest.raises(ValueError, match="every midpoint"):
        D.Equil(t)


def test_the_nonzero_count_of_a_block_with_no_entry_at_all_is_zero():
    """A rank whose whole block of rows is empty -- every one of its cells with
    only held midpoints -- has nothing for reduceat to reduce over, and an empty
    preallocation is what PETSc has to be handed."""
    A = sp.csr_matrix((np.array([1.0]), np.array([0]), np.array([0, 0, 0, 1])),
                      shape=(3, 4))
    assert np.array_equal(D._row_counts(A, 0, 2, 0, 4), np.zeros(2, np.int32))
    assert np.array_equal(D._row_counts(A, 0, 3, 0, 4), np.array([0, 0, 1]))


def test_a_refinement_pass_that_runs_out_of_iterations_keeps_the_step():
    """A refinement pass is an improvement on a step that is already fit to
    write, so one that does not converge in its own budget is a pass that does
    not pay: the step before it stands and the run says why it stopped. Refusing
    here would throw away a dataset that already meets the criterion. The first
    solve keeps its hard refusal, having nothing to fall back to."""
    t, U = walled(*channel3d(3), 3, field=smooth_noslip)
    D.REFINE_MAXITER, budget = 1, D.REFINE_MAXITER
    try:
        Uc, info = D.Equil(t).apply(U, "up_0000.h5")
    finally:
        D.REFINE_MAXITER = budget
    assert info["stepped"] and info["refinements"] == 0
    assert "did not converge" in info["refine_stop"]
    assert info["balance_after"] <= D.FLUX_TOL


def test_a_petsc_failure_in_a_refinement_pass_is_not_taken_for_one_that_does_not_pay():
    """Only a pass that stops short of its tolerance is a pass that does not
    pay. PETSc's own errors -- memory, a broken preconditioner -- are
    RuntimeErrors too, and on several ranks one raised on a single rank and
    taken there as the end of the refinement would send that rank on while the
    others solve: it is raised as it is."""
    from petsc4py import PETSc
    t, U = walled(*channel3d(3), 3, field=smooth_noslip)
    eq = D.Equil(t, minres_rtol=1e-4)
    real, calls = eq._minres, []

    def minres(*args, **kw):
        calls.append(1)
        if len(calls) > 1:
            raise PETSc.Error(55)
        return real(*args, **kw)
    eq._minres = minres
    with pytest.raises(PETSc.Error):
        eq.apply(U, "up_0000.h5")
    assert len(calls) == 2, "no refinement pass was taken"


def test_the_refinement_reaches_its_target_and_not_the_first_right_hand_side():
    """The absolute tolerance of every solve is the criterion's own smallest net
    flux, not a fraction of the first right-hand side, so the refinement passes
    reach FLUX_REFINE of the criterion rather than stopping as many decades
    below the first residual as the relative tolerance allows. Pinned on the 3D
    channel."""
    t, U = walled(*channel3d(3), 3, field=smooth_noslip)
    eq = D.Equil(t)
    Uc, info = eq.apply(U)
    eq.destroy()
    assert info["refinements"] >= 1
    assert info["refine_stop"] == "the target was reached"
    assert info["balance_after"] < 1e-2 * D.FLUX_REFINE * D.FLUX_TOL


@pytest.mark.parametrize("dim,mesh", CASES, ids=IDS)
def test_a_loosened_solve_is_refined_until_the_fluxes_are_at_round_off(dim, mesh):
    """MINRES stops on its own estimate of the residual, which on a 3D system can
    sit decades above the true one; the refinement passes are what makes the
    fluxes trustworthy whatever the estimate said."""
    t, U = walled(*mesh, dim, field=smooth_noslip)
    loose = D.Equil(t, minres_rtol=1e-6)
    Uc, info = loose.apply(U)
    assert loose.refinements >= 1
    assert info["flux_after_max"] < 1e-11 * info["facet_flux_max"]
    tight, _ = D.Equil(t).apply(U)
    assert np.abs(Uc - tight).max() < 1e-9 * np.abs(U).max()


def test_the_constants_are_read_where_they_are_used(monkeypatch):
    """A module constant captured as a default argument is evaluated once when
    the function is defined, so setting it afterwards is silently ignored; these
    are knobs a caller does set."""
    t, _ = walled(*channel2d(6), 2)
    monkeypatch.setattr(D, "MINRES_RTOL", 1e-3)
    monkeypatch.setattr(D, "BOUNDARY_W", 7.0)
    monkeypatch.setattr(D, "REST_TOL", 1e-3)
    eq = D.Equil(t)
    assert eq.minres_rtol == 1e-3 and eq.boundary_weight == 7.0
    eq.destroy()
    assert D.at_rest_nodes(np.array([1e-4]), 1.0)[0]


def test_the_fluxes_are_computed_once_where_they_are_wanted(tmp_path, monkeypatch):
    """cell_flux is the whole cost of a stamp that takes the early exit and of
    --check, paid on every rank. The step evaluates it on the data, once a
    refinement test and once on the result, and nowhere else."""
    pytest.importorskip("h5py")
    t, U = walled(*channel2d(6), 2, field=smooth_noslip)
    calls, real = [0], D.cell_flux

    def counted(topo, V):
        calls[0] += 1
        return real(topo, V)

    monkeypatch.setattr(D, "cell_flux", counted)
    eq = D.Equil(t)
    _, info = eq.apply(U)
    eq.destroy()
    assert info["stepped"]
    assert calls[0] == 3 + info["refinements"]
    cfg = write_case(tmp_path / "in", t.X, t.cells, [U, 0.5 * U], [True, False])
    calls[0] = 0
    D.check_case(cfg, quiet=True)
    assert calls[0] == 2, "one a stamp"


def test_the_solvers_petsc_objects_do_not_outlive_the_case(tmp_path):
    """A caller that cleans several cases in one process would otherwise keep
    every KSP, nest and GAMG hierarchy alive, and every instance's entries in
    PETSc's options database, which are what its options prefix is fresh for."""
    pytest.importorskip("h5py")
    from petsc4py import PETSc
    t, _ = walled(*channel2d(6), 2)
    cfg = write_case(tmp_path / "in", t.X, t.cells, [smooth_noslip(t.node_x)],
                     [True, False])
    before = len(PETSc.Options().getAll())
    D.clean_case(cfg, tmp_path / "out", verbose=False)
    assert len(PETSc.Options().getAll()) == before
    eq = D.Equil(t)
    assert len(PETSc.Options().getAll()) > before
    eq.destroy()
    assert len(PETSc.Options().getAll()) == before


# ------------------------------- the criterion the tool shares with the loader


def test_the_tool_and_the_loader_share_their_constants():
    """What the tool refines against is what the loader refuses: the two
    constants are written in both, so they are read from the C++ source here."""
    import re
    src = open(os.path.join(REPO, "src", "interpol", "SplitInterpol.cpp")).read()
    got = {k: float(v) for k, v in
           re.findall(r"constexpr double (flux_tol|flux_floor) = ([0-9.eE+-]+);", src)}
    assert got == {"flux_tol": D.FLUX_TOL, "flux_floor": D.FLUX_FLOOR}


@pytest.mark.parametrize("dim,mesh", CASES, ids=IDS)
def test_a_dead_cell_is_measured_against_the_stamp_and_not_itself(dim, mesh):
    """A cell whose own facet fluxes are the solver's noise -- a dead-end pore
    of a bead pack -- has a net flux that is round-off over round-off against
    its own scale and nothing at all against the stamp's. The criterion the tool
    and the loader share floors a cell's scale at FLUX_FLOOR of the stamp's
    largest facet flux, so such a cell is not refused; the field here is a
    nodal function of the last coordinate alone, which every cell carries as one
    quadratic, so every net flux is exactly zero before one midpoint is moved."""
    X, cells = mesh
    t = D.Topo(X, cells, [False] * dim)
    U = np.zeros((t.nnodes, dim))
    U[:, 0] = np.where(t.node_x[:, -1] < 0.5, 1e-13, 1.0)
    r, big = D.cell_flux(t, U)
    assert np.abs(r).max() < 1e-12 * big.max(), "the field itself has to balance"
    # a midpoint no live cell touches, moved off every facet normal
    nodes = np.concatenate([t.cells, t.nverts + t.cell_edge], axis=1)
    live = np.zeros(t.nedges, bool)
    live[t.cell_edge[t.node_x[nodes, -1].max(axis=1) >= 0.5].ravel()] = True
    e = np.nonzero(~live)[0][0]
    U[t.nverts + e] += 1e-18 * 4.0 ** np.arange(dim)
    r, big = D.cell_flux(t, U)
    own = np.abs(r) / np.maximum(big, 1e-300)
    ratio, facet = D.flux_ratios(t, U)
    assert big[int(np.argmax(own))] < D.FLUX_FLOOR * facet, "the cell is not dead"
    assert own.max() > D.FLUX_TOL, "its own scale alone would refuse it"
    assert ratio.max() < 1e-2 * D.FLUX_TOL, "the stamp's scale says it is round-off"
    # so the tool takes no step on it and writes it as it stands
    Uc, info = D.Equil(t).apply(U)
    assert not info["stepped"]
    assert np.abs(Uc - U).max() == 0.0


def test_a_non_finite_value_is_refused_by_stamp(tmp_path):
    """A diverged solver run writes NaNs, and a NaN compares false with every
    tolerance, so no refusal would fire: the cleaner and --check must reject a
    non-finite value where they read it, naming the stamp, and the step must
    reject one in arrays handed to it directly."""
    pytest.importorskip("h5py")
    t, U = walled(*channel2d(6), 2, field=smooth_noslip)
    with pytest.raises(ValueError, match="up_0005.h5: the velocity holds a non-finite"):
        bad = U.copy()
        bad[t.nverts + 3, 1] = np.nan
        D.Equil(t).apply(bad, "up_0005.h5")
    bad = U.copy()
    bad[t.nverts + 3, 1] = np.inf
    cfg = write_case(tmp_path / "in", t.X, t.cells, [U, bad], [True, False])
    with pytest.raises(ValueError, match="up_1.h5: 'u' holds a non-finite value"):
        D.clean_case(cfg, tmp_path / "out", verbose=False)
    with pytest.raises(ValueError, match="up_1.h5: 'u' holds a non-finite value"):
        D.check_case(cfg, quiet=True)


@pytest.mark.parametrize("dim,mesh", CASES, ids=IDS)
def test_the_check_tells_taylor_hood_data_from_lifted_data(dim, mesh, tmp_path, capsys):
    """int phi_i div u vanishes for every vertex of a converged Taylor-Hood
    velocity and for nothing else: it is what says whether a dataset's split
    field will hold the solver's volume mean."""
    X, cells = mesh
    t = D.Topo(X, cells, [False] * dim)
    rng = np.random.default_rng(11)
    B = th_matrix(t)
    u0 = smooth_noslip(t.node_x).ravel()
    good = (u0 - B.T @ np.linalg.lstsq(B @ B.T, B @ u0, rcond=None)[0]).reshape(-1, dim)
    bad = D.p1_to_p2(t, rng.normal(size=t.X.shape))
    assert D.taylor_hood_moments(t, good)[1] < 1e-10
    assert D.taylor_hood_moments(t, bad)[1] > 1e-2
    for name, U, want in (("good", good, "looks"), ("bad", bad, "does not look")):
        cfg = write_case(tmp_path / name, X, cells, [U], [False] * dim)
        D.check_case(cfg)
        out = capsys.readouterr().out
        assert "%s Taylor-Hood-compatible" % want in out
        assert "volume mean" in out


# --------------------------------------------------------------------- slivers


def cell_conditioning(t):
    """(condition number of each cell's Jacobian, norm of its own map
    J R J^-1): how far the reconstruction can amplify any boundary data."""
    d = t.dim
    R = D.reference_matrix(d)[0]
    nb, ni = R.shape[1] // d, R.shape[0] // d
    sv = np.linalg.svd(t.J, compute_uv=False)
    nm = np.array([np.linalg.norm(np.kron(np.eye(ni), t.J[k]) @ R
                                  @ np.kron(np.eye(nb), t.Jinv[k]), 2)
                   for k in range(t.ncells)])
    return sv[:, 0] / sv[:, -1], nm


def sliver_case(dim, squash, do_pinch=False, penalty=None):
    """A channel with its bottom layer of cells flattened by `squash`, cleaned;
    returns the topology, the field, the interior values and the per-cell ratio
    of the largest interior value to the largest boundary value."""
    X, cells = channel2d(6, squash) if dim == 2 else channel3d(3, squash)
    if do_pinch:
        X = pinch(X, cells, np.random.default_rng(7))
    t, U = walled(X, cells, dim, field=smooth_noslip)
    Uc, info = D.Equil(t, penalty=penalty).apply(U)
    uint = D.split_interior(t, Uc)
    g = np.concatenate([Uc[t.cells], Uc[t.nverts + t.cell_edge]], axis=1)
    ratio = (np.linalg.norm(uint, axis=2).max(axis=1)
             / np.maximum(np.linalg.norm(g, axis=2).max(axis=1), 1e-300))
    return t, Uc, uint, ratio, info


def worst_cell_divergence(t, Uc, uint, ratio, dim):
    """The pointwise divergence in the five cells with the largest interior
    values, against the size of those values."""
    rng = np.random.default_rng(3)
    worst = np.argsort(-ratio)[:5]
    w = rng.dirichlet(np.ones(dim + 1), len(worst))
    pts = np.einsum('ni,nij->nj', w, t.X[t.cells[worst]])
    div = float(np.abs(D.eval_split(t, Uc, uint, pts, worst, grad=True)).max())
    return div / max(float(np.linalg.norm(uint[worst], axis=2).max()), 1e-300)


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("squash", [20, 100])
def test_slivers_amplify_the_interior_values_by_the_cell_conditioning(dim, squash):
    """The local problem is uniquely solvable on any non-degenerate cell, but the
    constant behind it is the cell's shape: the map J R J^-1 from the boundary
    values to the interior ones has norm about 1.2 times the condition number of
    J, so the interior values do *not* stay in proportion on a flattened cell.
    That is the law this locks in, on a layer of cells of aspect ratio 20 and
    100. In 3D a smooth no-slip field realises it with the smallest-change step
    -- the largest interior value is about the aspect ratio times the boundary
    values -- while in 2D the flattened direction is the one the field is small
    in, so the realised ratio stays near one. What the interior values
    are made of is the divergence the global step leaves behind, so the
    penalised step must not realise more of the bound than that. What survives
    either way: div u is zero to round-off however flat the cell, the fluxes
    balance and no wall node moves."""
    t, Uc, uint, ratio, info = sliver_case(dim, squash, penalty=0.0)
    cond, nm = cell_conditioning(t)
    assert cond.max() > squash, "the mesh holds no sliver"
    assert np.abs(Uc[t.held_node]).max() == 0.0
    assert info["flux_after_max"] < 1e-10 * info["facet_flux_max"]
    assert worst_cell_divergence(t, Uc, uint, ratio, dim) < 1e-8
    # the amplification is bounded by the shape alone, and the bound is tight
    assert (nm <= 1.5 * cond).all()
    assert nm.max() > 0.5 * cond.max()
    if dim == 3:
        assert 0.5 * squash < ratio.max() < 1.5 * squash
    else:
        assert ratio.max() < 2.0
    pen = sliver_case(dim, squash)[3].max()
    assert pen <= (0.3 if dim == 3 else 1.05) * ratio.max()


@pytest.mark.parametrize("dim", [2, 3])
def test_near_degenerate_cells_keep_the_field_divergence_free(dim):
    """A few vertices slid almost onto a neighbour give cells of condition number
    1e3 to 1e4, far past anything a mesher would produce. The reconstruction
    still returns a divergence-free field with no slip on the walls; its interior
    values there run into the thousands, in proportion to the conditioning, which
    is the number the loader side has to live with."""
    t, Uc, uint, ratio, info = sliver_case(dim, 100, do_pinch=True)
    cond, nm = cell_conditioning(t)
    assert cond.max() > 1e3
    assert np.abs(Uc[t.held_node]).max() == 0.0
    assert info["flux_after_max"] < 1e-10 * info["facet_flux_max"]
    assert worst_cell_divergence(t, Uc, uint, ratio, dim) < 1e-8
    assert (nm <= 3.0 * cond).all()
