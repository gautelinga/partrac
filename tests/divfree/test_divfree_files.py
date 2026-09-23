"""python/divfree/divfree_clean.py: the files it reads and writes.

The datasets and attributes the dolfin-free loaders open, laid out as dolfin
lays them out; a field read by the cell chunk rather than whole; the stamp and
frequency lists copied with the new file names; the pressure and phase field
carried through, resampled onto the split mesh under --split; the keys the
output's parameter file must and must not carry; a malformed checkpoint refused
by name. Some tests need h5py, one dolfin, to write the reference layout with.
"""

import numpy as np
import pytest

from divfree_cases import (D, channel2d, channel3d, relabel_by_global_id,
                           smooth_noslip, walled, write_case)


# ------------------------------------------------------- the checkpoint layout


def test_the_written_checkpoint_holds_what_the_loaders_open(tmp_path):
    """simplex_load reads the signature attribute, cell_dofs, x_cell_dofs, cells
    and vector_0 of the first file and vector_0 of every later one. The file the
    tool writes must hold those, with a constant dof stride, and read back as the
    values that went in."""
    h5py = pytest.importorskip("h5py")
    X, cells = channel3d(2)
    t, U = walled(X, cells, 3, field=smooth_noslip)
    D.write_checkpoint(tmp_path / "u.h5", "u", U, t, 2)
    with h5py.File(tmp_path / "u.h5", "r") as f:
        assert set(f["u"].keys()) == {"cell_dofs", "x_cell_dofs", "cells", "vector_0"}
        sig = f["u"].attrs["signature"]
        assert sig.decode() == "VectorElement(FiniteElement('Lagrange', tetrahedron, 2), dim=3)"
        xc = np.array(f["u/x_cell_dofs"])
        cd = np.array(f["u/cell_dofs"])
        assert len(xc) == t.ncells + 1 and len(set(np.diff(xc))) == 1
        assert int(np.diff(xc)[0]) == 10 * 3          # P2 nodes a tet, three components
        assert np.array_equal(np.array(f["u/cells"]), np.arange(t.ncells))
        assert cd.max() < len(np.array(f["u/vector_0"]))
    D.write_mesh_h5(tmp_path / "mesh.h5", X, cells)
    _, topo, ci = D.read_mesh_h5(tmp_path / "mesh.h5")
    back = D.DofTable(tmp_path / "u.h5", "u", t, ci).values(tmp_path / "u.h5")
    assert np.array_equal(back, U)


def test_the_written_layout_matches_a_dolfin_written_one(tmp_path):
    """The same field written by dolfin's own HDF5File: the datasets, their
    shapes and types, the group attribute and the values a loader would see must
    be the same, since nothing but this tool will write these files."""
    df = pytest.importorskip("dolfin", reason="no dolfin to write the reference with")
    h5py = pytest.importorskip("h5py")
    X, cells = channel2d(3)
    t, U = walled(X, cells, 2, field=smooth_noslip)
    mesh = df.Mesh()
    ed = df.MeshEditor()
    ed.open(mesh, "triangle", 2, 2)
    ed.init_vertices(len(X))
    ed.init_cells(len(cells))
    for i, x in enumerate(X):
        ed.add_vertex(i, x)
    for i, c in enumerate(cells):
        ed.add_cell(i, np.asarray(c, dtype=np.uintp))
    ed.close()
    V = df.VectorFunctionSpace(mesh, "CG", 2)
    u = df.Function(V)
    with df.HDF5File(mesh.mpi_comm(), str(tmp_path / "df_mesh.h5"), "w") as f:
        f.write(mesh, "mesh")
    with df.HDF5File(mesh.mpi_comm(), str(tmp_path / "df_u.h5"), "w") as f:
        f.write(u, "u")
    D.write_mesh_h5(tmp_path / "our_mesh.h5", X, cells)
    D.write_checkpoint(tmp_path / "our_u.h5", "u", U, t, 2)

    def layout(path):
        out = {}
        with h5py.File(path, "r") as f:
            f.visititems(lambda n, o: out.__setitem__(
                n, (o.shape, o.dtype.kind, sorted(o.attrs))) if isinstance(o, h5py.Dataset) else None)
            groups = {k: sorted(f[k].attrs) for k in f}
        return out, groups

    for ours, theirs in (("our_mesh.h5", "df_mesh.h5"), ("our_u.h5", "df_u.h5")):
        a, ga = layout(tmp_path / ours)
        b, gb = layout(tmp_path / theirs)
        assert a == b, "%s: %s against dolfin's %s" % (ours, a, b)
        assert ga == gb
    with h5py.File(tmp_path / "our_u.h5", "r") as f, h5py.File(tmp_path / "df_u.h5", "r") as g:
        assert f["u"].attrs["signature"] == g["u"].attrs["signature"]


def test_reading_a_field_does_not_allocate_the_whole_gather(tmp_path):
    """The dof table names ten nodes a cell where the field holds one value a
    node, so gathering the whole of it at once costs several times the field
    itself -- the largest allocation of the read, paid for every stamp and on
    every rank. A cell chunk at a time costs the chunk."""
    pytest.importorskip("h5py")
    import tracemalloc
    X, cells = channel3d(20)
    t = D.Topo(X, cells, [False] * 3)
    U = smooth_noslip(t.node_x)
    D.write_checkpoint(tmp_path / "u.h5", "u", U, t, 2)
    dof = D.DofTable(tmp_path / "u.h5", "u", t, None)
    tracemalloc.start()
    got = dof.values(tmp_path / "u.h5")
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert np.array_equal(got, U)
    assert peak < 7 * U.nbytes


# ---------------------------------------------------------------- a whole case


def test_a_time_series_is_cleaned_and_its_stamp_list_copied(tmp_path):
    """A whole case through the tool: the fluxes of every stamp are zero
    afterwards, the first file carries the dof table and the later ones only the
    values, and the stamp list keeps its times with the new file names."""
    pytest.importorskip("h5py")
    X, cells = channel2d(6)
    t = D.Topo(X, cells, [True, False])
    fields = [smooth_noslip(t.node_x), 0.5 * smooth_noslip(t.node_x)]
    cfg = write_case(tmp_path / "in", X, cells, fields, [True, False], stamps=["0", "2.5"])
    rep = D.clean_case(cfg, tmp_path / "out", write_key=False, verbose=False)
    assert rep["held_edges"] > 0 and rep["after"] < 1e-10 * rep["facet"]
    assert D.check_case(tmp_path / "out" / "dolfin_params.dat", quiet=True) < 1e-10
    lines = [l.split() for l in
             (tmp_path / "out" / "timestamps.dat").read_text().split("\n") if l.strip()]
    assert [l[0] for l in lines] == ["0", "2.5"]
    assert [l[1] for l in lines] == ["u_0000.h5", "u_0001.h5"]
    h5py = pytest.importorskip("h5py")
    with h5py.File(tmp_path / "out" / "u_0001.h5", "r") as f:
        assert set(f["u"].keys()) == {"vector_0"}, "a later stamp carries more than the values"
    assert "divfree" not in (tmp_path / "out" / "dolfin_params.dat").read_text()


def test_a_case_labelled_by_global_cell_id_cleans_to_the_same_field(tmp_path):
    """A case solved on several ranks labels its rows by global cell id, and its
    mesh's `cell_indices` is then a permutation. The output is written on that
    same mesh -- it is linked through -- so its own rows must carry those labels
    too; `0..n-1` would be refused by every reader that composes through
    `cell_indices`, which all of them do.

    The cleaned field itself cannot depend on the labelling, so it is compared
    with the same case cleaned unlabelled, node for node."""
    h5py = pytest.importorskip("h5py")
    X, cells = channel2d(6)
    t = D.Topo(X, cells, [True, False])
    fields = [smooth_noslip(t.node_x), 0.5 * smooth_noslip(t.node_x)]
    plain = write_case(tmp_path / "plain", X, cells, fields, [True, False])
    cfg = write_case(tmp_path / "in", X, cells, fields, [True, False])
    gid = relabel_by_global_id(tmp_path / "in")
    D.clean_case(plain, tmp_path / "plain_out", verbose=False)
    D.clean_case(cfg, tmp_path / "out", verbose=False)

    out = tmp_path / "out"
    with h5py.File(out / "u_0000.h5", "r") as f:
        assert np.array_equal(np.array(f["u/cells"]), gid), "the rows are not labelled"
        # the pressure comes through untouched, so the file must agree with itself
        assert np.array_equal(np.array(f["p/cells"]), gid)
    case = D.read_case(out / "dolfin_params.dat")
    assert np.array_equal(case["cell_indices"], gid), "the input's mesh is the output's"
    ct = D.Topo(case["X"], case["cells"], case["periodic"])
    assert D.check_case(out / "dolfin_params.dat", quiet=True) < 1e-10
    for name in ("u_0000.h5", "u_0001.h5"):
        U = D.DofTable(out / "u_0000.h5", "u", ct, case["cell_indices"]).values(out / name)
        ref = D.DofTable(tmp_path / "plain_out" / "u_0000.h5", "u", ct, None).values(
            tmp_path / "plain_out" / name)
        assert np.array_equal(U, ref), name


@pytest.mark.parametrize("form", ["ta", "omega"])
def test_frequency_components_round_trip_with_new_file_names(tmp_path, form):
    """A frequency file is a dataset too: every component is cleaned against the
    same wall set, and the list is copied with its numbers as written -- both
    line forms, `t a file` and `omega phi a file` -- and only the names changed.
    Each component keeps its signature, which the frequency loader checks on
    every file against the first."""
    h5py = pytest.importorskip("h5py")
    X, cells = channel2d(6)
    t = D.Topo(X, cells, [True, False])
    fields = [smooth_noslip(t.node_x), 0.3 * smooth_noslip(t.node_x)]
    cols = [["0.0", "1.0"], ["0.125", "0.7"]] if form == "ta" \
        else [["0.0", "0.0", "1.0"], ["6.2831853", "-1.5707963", "0.7"]]
    cfg = write_case(tmp_path / "in", X, cells, fields, [True, False], freq=cols)
    D.clean_case(cfg, tmp_path / "out", write_key=False, verbose=False)
    lines = [l.split() for l in
             (tmp_path / "out" / "freqstamps.dat").read_text().split("\n") if l.strip()]
    assert [l[:-1] for l in lines] == cols
    assert [l[-1] for l in lines] == ["u_0000.h5", "u_0001.h5"]
    assert D.check_case(tmp_path / "out" / "dolfin_params.dat", quiet=True) < 1e-10
    with h5py.File(tmp_path / "out" / "u_0001.h5", "r") as f:
        assert "signature" in f["u"].attrs


def test_the_report_carries_what_is_no_longer_constrained(tmp_path):
    """The volume mean, the throughput drift and the data's net boundary flux
    are reported and never held, so a whole series has to bring them out: the
    means of the first stamp, and the largest drift and imbalance of any of
    them. A stamp that is balanced already takes no step and drifts nothing."""
    pytest.importorskip("h5py")
    t, _ = walled(*channel2d(6), 2)
    U = D.p1_to_p2(t, np.random.default_rng(2).normal(size=t.X.shape))
    U[t.held_node] = 0.0
    U = U[t.master]
    eq = D.Equil(t)
    Uc, info = eq.apply(U)
    _, again = eq.apply(Uc)
    assert info["stepped"] and not again["stepped"]
    assert np.abs(again["drift"]).max() == 0.0
    cfg = write_case(tmp_path / "in", t.X, t.cells, [U, Uc], [True, False])
    rep = D.clean_case(cfg, tmp_path / "out", verbose=False)
    assert np.array_equal(rep["mean_before"], info["mean_before"])
    assert np.array_equal(rep["mean_after"], info["mean_after"])
    assert rep["throughput_drift"] == abs(float(info["drift"][0])) > 0.0
    assert rep["imbalance"] == max(info["imbalance"], again["imbalance"])


def test_a_stamp_without_the_field_to_carry_through_is_refused(tmp_path):
    """The pressure and the phase field are copied out of the stamp's own file,
    which with field_file need not be the one that carries them, and a later
    stamp of some datasets holds only its velocity. Saying nothing leaves the
    output declaring a pressure field it does not have, which a loader finds out
    much later."""
    h5py = pytest.importorskip("h5py")
    t, _ = walled(*channel2d(6), 2)
    U = smooth_noslip(t.node_x)
    cfg = write_case(tmp_path / "in", t.X, t.cells, [U, 0.5 * U], [True, False])
    with h5py.File(tmp_path / "in" / "up_1.h5", "r+") as f:
        del f["p"]
    with pytest.raises(ValueError, match=r"up_1\.h5 holds no 'p'"):
        D.clean_case(cfg, tmp_path / "out", verbose=False)


def test_a_cached_mesh_key_does_not_survive_into_a_cleaned_case(tmp_path):
    """The divergence-free loader refuses mesh_cache by design, so a case solved
    with it would clean into one no loader reads. Without the key the output is
    an ordinary P2 case and the cache is the caller's business again."""
    pytest.importorskip("h5py")
    t, _ = walled(*channel2d(6), 2)
    cfg = write_case(tmp_path / "in", t.X, t.cells, [smooth_noslip(t.node_x)],
                     [True, False])
    cfg.write_text(cfg.read_text() + "mesh_cache=true\n")
    D.clean_case(cfg, tmp_path / "out", verbose=False)
    assert "divfree=true" in (tmp_path / "out" / "dolfin_params.dat").read_text()
    assert "mesh_cache" not in (tmp_path / "out" / "dolfin_params.dat").read_text()
    D.clean_case(cfg, tmp_path / "out2", write_key=False, verbose=False)
    assert "mesh_cache=true" in (tmp_path / "out2" / "dolfin_params.dat").read_text()


def test_a_p1_case_comes_out_as_p2(tmp_path):
    """P1 data has no edge midpoints; the tool takes the mean of the ends and the
    output is a P2 case, so its parameter file must say so."""
    pytest.importorskip("h5py")
    X, cells = channel2d(6)
    t = D.Topo(X, cells, [True, False])
    cfg = write_case(tmp_path / "in", X, cells, [smooth_noslip(t.node_x)],
                     [True, False], degree=1)
    rep = D.clean_case(cfg, tmp_path / "out", write_key=False, verbose=False)
    assert rep["degree_in"] == 1
    assert "velocity_space=P2" in (tmp_path / "out" / "dolfin_params.dat").read_text()
    assert D.check_case(tmp_path / "out" / "dolfin_params.dat", quiet=True) < 1e-10


# ------------------------------------------------------------------ field_file


def test_one_field_file_is_cleaned_as_a_steady_case_of_its_own(tmp_path):
    """`field_file` cleans one file of the folder in place of the series the
    parameter file names, which is what a directory of several fields on one
    mesh needs. Every stamp of the output names the one cleaned file; the file cleaned
    is the one asked for, which its untouched vertex values say; and the report
    carries that field by node, on the topology the output is written on."""
    pytest.importorskip("h5py")
    X, cells = channel2d(6)
    t = D.Topo(X, cells, [True, False])
    fields = [smooth_noslip(t.node_x), 0.5 * smooth_noslip(t.node_x)]
    cfg = write_case(tmp_path / "in", X, cells, fields, [True, False], stamps=["0", "2.5"])
    rep = D.clean_case(cfg, tmp_path / "out", field_file="up_1.h5", verbose=False)
    lines = [l.split() for l in
             (tmp_path / "out" / "timestamps.dat").read_text().split("\n") if l.strip()]
    assert [l[0] for l in lines] == ["0", "2.5"]
    assert [l[1] for l in lines] == ["u_0000.h5", "u_0000.h5"]
    assert D.check_case(tmp_path / "out" / "dolfin_params.dat", quiet=True) < 1e-10
    # only the midpoints move, so the vertex values tell the two inputs apart
    assert rep["topo"].nnodes == t.nnodes and rep["values"].shape == (t.nnodes, 2)
    assert np.abs(rep["values"][:t.nverts] - fields[1][:t.nverts]).max() == 0.0
    with pytest.raises(ValueError, match="up_7.h5"):
        D.clean_case(cfg, tmp_path / "out", field_file="up_7.h5", verbose=False)


def test_field_file_reaches_the_tool_from_the_command_line(tmp_path):
    """field_file is a supported mode with users of its own, so --help has to
    list it and the flag has to arrive where the library argument does. Only the
    midpoints move, so the vertex values tell the two inputs apart."""
    pytest.importorskip("h5py")
    t, _ = walled(*channel2d(6), 2)
    fields = [smooth_noslip(t.node_x), 0.5 * smooth_noslip(t.node_x)]
    cfg = write_case(tmp_path / "in", t.X, t.cells, fields, [True, False])
    assert D.main([str(cfg), "--out", str(tmp_path / "out"), "--field-file",
                   str(tmp_path / "in" / "up_1.h5")]) == 0
    names = (tmp_path / "out" / "timestamps.dat").read_text().split()[1::2]
    assert names == ["u_0000.h5", "u_0000.h5"]
    case = D.read_case(tmp_path / "out" / "dolfin_params.dat")
    ct = D.Topo(case["X"], case["cells"], case["periodic"])
    out = tmp_path / "out" / "u_0000.h5"
    Uc = D.DofTable(out, "u", ct, case["cell_indices"]).values(out)
    assert np.abs(Uc[:ct.nverts] - fields[1][:ct.nverts]).max() == 0.0


# --------------------------------------------------------------------- --split


def test_the_split_option_writes_a_divergence_free_field(tmp_path):
    """--split writes the full reconstruction as P2 on the barycentric split,
    which the plain mesh loaders read. Read back through that mesh, its
    pointwise divergence is zero, which the cleaned P2 on the original mesh is
    not."""
    pytest.importorskip("h5py")
    X, cells = channel2d(5)
    t = D.Topo(X, cells, [True, False])
    cfg = write_case(tmp_path / "in", X, cells, [smooth_noslip(t.node_x)], [True, False])
    D.clean_case(cfg, tmp_path / "out", split=True, write_key=False, verbose=False)
    case = D.read_case(tmp_path / "out" / "dolfin_params.dat")
    st = D.Topo(case["X"], case["cells"], case["periodic"])
    assert st.ncells == t.ncells * 3 and st.nverts == t.nverts + t.ncells
    U = D.DofTable(tmp_path / "out" / "u_0000.h5", "u", st, case["cell_indices"]).values(
        tmp_path / "out" / "u_0000.h5")
    rng = np.random.default_rng(8)
    k = rng.integers(0, st.ncells, 40)
    w = rng.dirichlet(np.ones(3), 40)
    pts = np.einsum('ni,nij->nj', w, st.X[st.cells[k]])
    assert np.abs(D.p2_eval(st, U, pts, k, grad=True)).max() < 1e-9 * np.abs(U).max()


def test_the_split_option_writes_its_own_mesh_and_its_own_labels(tmp_path):
    """--split writes a mesh of its own, whose cells are new and whose
    `cell_indices` is the identity, so a labelled input's labels must not follow
    the field onto it."""
    h5py = pytest.importorskip("h5py")
    X, cells = channel2d(5)
    t = D.Topo(X, cells, [True, False])
    cfg = write_case(tmp_path / "in", X, cells, [smooth_noslip(t.node_x)], [True, False])
    relabel_by_global_id(tmp_path / "in")
    D.clean_case(cfg, tmp_path / "out", split=True, write_key=False, verbose=False)
    case = D.read_case(tmp_path / "out" / "dolfin_params.dat")
    assert np.array_equal(case["cell_indices"], np.arange(t.ncells * 3))
    with h5py.File(tmp_path / "out" / "u_0000.h5", "r") as f:
        assert np.array_equal(np.array(f["u/cells"]), np.arange(t.ncells * 3))
    st = D.Topo(case["X"], case["cells"], case["periodic"])
    U = D.DofTable(tmp_path / "out" / "u_0000.h5", "u", st, case["cell_indices"]).values(
        tmp_path / "out" / "u_0000.h5")
    assert np.isfinite(U).all()


def quadratic_phi(x):
    """A quadratic scalar, which a P2 field holds exactly on either mesh."""
    return (x[:, 0] ** 2 - x[:, 0] * x[:, 1] + 0.5 * x[:, 1] ** 2 + 0.3)[:, None]


def test_the_split_option_resamples_a_p2_phase_field_exactly(tmp_path):
    """A field carried through keeps its degree on the split mesh, whose P2 space
    holds the macro field exactly: the phase field of every stamp reads back
    there as the quadratic it was, the later stamps' values in the first
    stamp's numbering."""
    pytest.importorskip("h5py")
    X, cells = channel2d(5)
    t = D.Topo(X, cells, [False, False])
    U = smooth_noslip(t.node_x)
    phi = [quadratic_phi(t.node_x), 2.0 * quadratic_phi(t.node_x)]
    cfg = write_case(tmp_path / "in", X, cells, [U, 0.5 * U], [False, False], phi=phi)
    D.clean_case(cfg, tmp_path / "out", split=True, write_key=False, verbose=False)
    out = tmp_path / "out"
    case = D.read_case(out / "dolfin_params.dat")
    st = D.Topo(case["X"], case["cells"], case["periodic"])
    dof = D.DofTable(out / "u_0000.h5", "phi", st, case["cell_indices"])
    assert dof.degree == 2 and st.ncells == 3 * t.ncells
    for k, name in enumerate(("u_0000.h5", "u_0001.h5")):
        want = (1.0 + k) * quadratic_phi(st.node_x)
        assert np.abs(dof.values(out / name) - want).max() < 1e-12, name


def test_a_p1_series_past_the_read_cache_cleans_as_one_read_into_it(tmp_path, monkeypatch):
    """The stamps that fit the read budget are kept from the pass that finds the
    held set; the others are read again for the cleaning, and a P1 one lifted
    to P2 again. Either way the stamp cleaned is the same."""
    pytest.importorskip("h5py")
    X, cells = channel2d(6)
    t = D.Topo(X, cells, [True, False])
    fields = [smooth_noslip(t.node_x), 0.5 * smooth_noslip(t.node_x)]
    cfg = write_case(tmp_path / "in", X, cells, fields, [True, False], degree=1)
    a = D.clean_case(cfg, tmp_path / "cached", verbose=False)
    monkeypatch.setattr(D, "READ_CACHE_BYTES", 0)
    b = D.clean_case(cfg, tmp_path / "reread", verbose=False)
    assert a["degree_in"] == b["degree_in"] == 1
    case = D.read_case(tmp_path / "cached" / "dolfin_params.dat")
    ct = D.Topo(case["X"], case["cells"], case["periodic"])
    dof = D.DofTable(tmp_path / "cached" / "u_0000.h5", "u", ct, case["cell_indices"])
    for name in ("u_0000.h5", "u_0001.h5"):
        assert np.array_equal(dof.values(tmp_path / "cached" / name),
                              dof.values(tmp_path / "reread" / name)), name


def poke_checkpoint(folder, kind, ncells):
    """Break the first stamp's velocity, or the mesh, the way `kind` names."""
    import h5py
    with h5py.File(folder / "up_0.h5", "r+") as f:
        g = f["u"]
        if kind == "ragged":
            g["x_cell_dofs"][1] = g["x_cell_dofs"][1] - 1
        elif kind == "unknown cell":
            # global ids with gaps, as a partitioned run can leave them
            with h5py.File(folder / "mesh.h5", "r+") as m:
                m["mesh/cell_indices"][...] = 2 * np.arange(ncells)
            g["cells"][...] = 2 * np.arange(ncells)
            g["cells"][0] = 1
        elif kind == "cells disagree":
            cd = g["cell_dofs"]
            cd[0], cd[1] = cd[1], cd[0]
        elif kind == "no degree":
            g.attrs["signature"] = np.bytes_("FiniteElement('Real', interval)")
    if kind == "node without a value":
        with h5py.File(folder / "mesh.h5", "r+") as m:
            X = np.array(m["mesh/coordinates"])
            del m["mesh/coordinates"]
            m["mesh/coordinates"] = np.vstack([X, [[0.5, 2.0]]])


MALFORMED = {
    "ragged": r"'u/x_cell_dofs' is ragged",
    "unknown cell": r"'u/cells' names a cell the mesh does not",
    "node without a value": r"up_0\.h5: 'u' leaves a node without a value",
    "cells disagree": r"up_0\.h5: cells disagree on a node of 'u'",
    "no degree": r"'u' is FiniteElement\('Real', interval\), which names no degree",
    "components": r"'u' has 1 components, not 2",
}


@pytest.mark.parametrize("kind", sorted(MALFORMED))
def test_a_malformed_checkpoint_is_refused_by_name(tmp_path, kind):
    """The dof table is read as the loaders read it, so a file they would refuse
    or misread -- a ragged or unknown row, a node no cell gives a value or two
    cells give two, an element with no degree or the wrong number of components
    -- is refused before anything is cleaned, naming the file and what is wrong
    in it."""
    pytest.importorskip("h5py")
    X, cells = channel2d(3)
    t = D.Topo(X, cells, [False, False])
    U = np.random.default_rng(12).normal(size=(t.nnodes, 2))
    cfg = write_case(tmp_path / "in", X, cells, [U[:, :1] if kind == "components" else U],
                     [False, False])
    if kind != "components":
        poke_checkpoint(tmp_path / "in", kind, t.ncells)
    with pytest.raises(ValueError, match=MALFORMED[kind]):
        D.clean_case(cfg, tmp_path / "out", verbose=False)
