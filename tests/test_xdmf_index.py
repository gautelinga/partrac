"""The XDMF index, mesh.xdmf, that a run writes beside its dumps.

A run that dumps keeps mesh.xdmf in its folder: one temporal collection with
a uniform grid per dump, in time order, so ParaView shows every dump so far
with no post-processing. Each grid is checked against the dump it names: the
time is the group's, every DataItem's file, dataset and dimensions exist as
stated, the topology follows from faces or edges, and every dumped dataset
other than points, faces and edges is an attribute, per node or per face or
edge. Where partractools is importable, the index agrees with what its
make_xdmf.py writes for the same dumps; where meshio is, it reads the index.
A restart continues the folder's one history: grids past the restart time are
dropped, and the index is complete and valid after every run.
"""

import os
import pathlib
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from paths import REPO, app, needs
from runs import checkpoint_folder, copy_example, run_app

PARTRAC = app("partrac")
TRACERS = app("tracers")
SPATIAL = app("tracervectors_spatial")
SINE = os.path.join(REPO, "data_example", "sine_flow", "expr_params.dat")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

MESH = ("mode=analytic x0=0.5 y0=0.5 z0=0.5 Nrw=50 Nrw_max=200000 Dm=0 int_order=1 dt=0.0625 "
        "stat_intv=1e9 random=false seed=1 dump_intv=0.125 T=0.5 ds_min=0.002 refine=true "
        "refine_intv=0.0625 integrate_tau=true tau_intv=0.0625 tau_max=0")
STRIP = "init_mode=strip_x La=0.5 ds_max=0.01"
SHEET = "init_mode=sheet_xy La=0.5 Lb=0.5 ds_init=0.05 ds_max=0.05 ds_min=0.01"
POINTS = ("mode=analytic init_mode=points_xyz x0=0 y0=0 z0=0 Nrw=20 Nrw_max=20 Dm=0 int_order=1 "
          "scheme=RK4 dt=0.1 stat_intv=1e9 random=false seed=1")
NOT_ATTRIBUTES = {"points", "faces", "edges"}


def run_folder(d):
    """The one folder under d that holds a mesh.xdmf."""
    [index] = list(pathlib.Path(d).rglob("mesh.xdmf"))
    return index.parent


def grids(folder):
    """The uniform grids of folder/mesh.xdmf, which must parse and hold one temporal collection."""
    root = ET.parse(folder / "mesh.xdmf").getroot()
    [collection] = root.find("Domain").findall("Grid")
    assert collection.get("GridType") == "Collection"
    assert collection.get("CollectionType") == "Temporal"
    return collection.findall("Grid")


def item_location(item):
    """The file and dataset a DataItem names."""
    assert item.get("Format") == "HDF"
    name, path = item.text.strip().split(":")
    assert path.startswith("/"), path
    return name, path[1:]


def check_index(folder):
    """Every grid of folder/mesh.xdmf against the dump it names; [(time, file)] in order.

    Every dataset of the group but points, faces and edges is an attribute,
    per face or edge where it has a row per cell and the cells and nodes
    differ in number, else per node.
    """
    h5py = pytest.importorskip("h5py")
    assert not [f for f in os.listdir(folder) if f.startswith("mesh.xdmf") and f != "mesh.xdmf"]
    out = []
    for g in grids(folder):
        assert g.get("GridType") == "Uniform"
        time = g.find("Time").get("Value")
        geometry = g.find("Geometry")
        assert geometry.get("GeometryType") == "XYZ"
        name, path = item_location(geometry.find("DataItem"))
        group = path.rsplit("/", 1)[0]
        assert float(group) == float(time) and path == group + "/points"
        with h5py.File(folder / name, "r") as h:
            assert group in h, (name, group)
            # every DataItem names a dataset of the shape and type it states
            for item in g.iter("DataItem"):
                n, p = item_location(item)
                assert n == name and p in h[group].parent, p
                d = h[p]
                assert [int(x) for x in item.get("Dimensions").split()] == list(d.shape), p
                kind = {"Float": "f", "Int": "i", "UInt": "u"}[item.get("NumberType")]
                assert (d.dtype.kind, d.dtype.itemsize) == (kind, int(item.get("Precision"))), p
            dsets = h[group]
            nodes = len(dsets["points"])
            topology = g.find("Topology")
            cells = "faces" if "faces" in dsets else "edges" if "edges" in dsets else None
            ncells = len(dsets[cells]) if cells else None
            if cells == "faces":
                assert (topology.get("TopologyType"), topology.get("NodesPerElement")) == ("Triangle", "3")
            elif cells == "edges":
                assert (topology.get("TopologyType"), topology.get("NodesPerElement")) == ("Polyline", "2")
            else:
                assert topology.get("TopologyType") == "Polyvertex"
                assert topology.get("NodesPerElement") == str(nodes)
            if cells:
                assert ncells != nodes   # the centre is told apart by rows
                assert topology.get("NumberOfElements") == str(ncells)
                assert item_location(topology.find("DataItem")) == (name, group + "/" + cells)
            attributes = {a.get("Name"): a for a in g.findall("Attribute")}
            assert set(attributes) == set(dsets) - NOT_ATTRIBUTES
            for k, a in attributes.items():
                rows, cols = dsets[k].shape
                assert item_location(a.find("DataItem")) == (name, group + "/" + k)
                assert a.get("Center") == ("Cell" if cells and rows == ncells else "Node"), k
                assert rows == (ncells if a.get("Center") == "Cell" else nodes), k
                assert a.get("AttributeType") == {1: "Scalar", 3: "Vector", 9: "Tensor"}[cols], k
        out.append((float(time), name))
    times = [t for t, _ in out]
    assert times == sorted(times) and len(set(times)) == len(times), times
    return out


def dumped(folder):
    """[(time, file)] of every dump group in the folder's dump files, in time order."""
    h5py = pytest.importorskip("h5py")
    out = []
    for f in sorted(folder.glob("data_from_t*.h5")):
        with h5py.File(f, "r") as h:
            out += [(float(k), f.name) for k in h]
    return sorted(out)


def same_as_make_xdmf(folder, tmp):
    """Each grid has make_xdmf.py's time, topology, geometry and attributes for the same dumps."""
    try:
        import partractools.make_xdmf   # noqa: F401
    except ImportError:
        return
    mx = tmp / "make_xdmf"
    mx.mkdir()
    for f in folder.iterdir():
        if f.name.startswith(("data_from_t", "params_from_t")):
            shutil.copy(f, mx / f.name)
    subprocess.run([sys.executable, "-m", "partractools.make_xdmf", str(mx)], check=True,
                   capture_output=True)
    ours, theirs = grids(folder), grids(mx)
    assert len(ours) == len(theirs)

    def item(e):
        d = e.find("DataItem")
        return d.get("Dimensions"), d.text.strip().replace(":/", ":")

    for a, b in zip(ours, theirs):
        assert float(a.find("Time").get("Value")) == float(b.find("Time").get("Value"))
        ta, tb = a.find("Topology"), b.find("Topology")
        for k in ("TopologyType", "NumberOfElements", "NodesPerElement"):
            assert (ta.get(k) or "").lower() == (tb.get(k) or "").lower(), k
        if tb.find("DataItem") is not None:
            assert item(ta) == item(tb)
        assert item(a.find("Geometry")) == item(b.find("Geometry"))
        attributes = {x.get("Name"): x for x in a.findall("Attribute")}
        for x in b.findall("Attribute"):
            y = attributes[x.get("Name")]
            assert (y.get("AttributeType"), y.get("Center")) == (x.get("AttributeType"), x.get("Center"))
            assert item(y) == item(x)


def meshio_reads(folder):
    """meshio's time series reader takes the index: every step's time and node data."""
    try:
        import meshio
    except ImportError:
        return
    reader = meshio.xdmf.TimeSeriesReader(str(folder / "mesh.xdmf"))
    points, _ = reader.read_points_cells()
    assert points.shape[1] == 3
    times = []
    for k in range(reader.num_steps):
        t, point_data, _ = reader.read_data(k)
        assert "c" in point_data
        times.append(t)
    assert times == [t for t, _ in check_index(folder)]


@needs(PARTRAC)
@pytest.mark.parametrize("init", [STRIP, SHEET], ids=["strip", "sheet"])
def test_a_strip_or_sheet_run_indexes_every_dump(tmp_path, init):
    """A refining strip or sheet: a grid per dump, in order, Polyline or
    Triangle topology with the per-edge or per-face fields as cell data."""
    run_app(PARTRAC, copy_example(SINE, tmp_path / "c"), MESH, init)
    folder = run_folder(tmp_path / "c")
    index = check_index(folder)
    assert index == dumped(folder)
    assert [t for t, _ in index] == pytest.approx([0.125 * k for k in range(5)])
    cells = {"Polyline": "dl", "Triangle": "dA"}
    [kind] = {g.find("Topology").get("TopologyType") for g in grids(folder)}
    assert {a.get("Name"): a.get("Center") for a in grids(folder)[-1].findall("Attribute")}[cells[kind]] == "Cell"
    same_as_make_xdmf(folder, tmp_path)
    meshio_reads(folder)


@needs(TRACERS)
def test_a_point_cloud_is_indexed_as_vertices(tmp_path):
    """Tracers with no mesh: a Polyvertex grid per dump, every field per node."""
    run_app(TRACERS, copy_example(POISEUILLE, tmp_path / "c"), POINTS, "T=0.4 dump_intv=0.1")
    folder = run_folder(tmp_path / "c")
    index = check_index(folder)
    assert index == dumped(folder)
    assert len(index) == 5
    same_as_make_xdmf(folder, tmp_path)


@needs(SPATIAL)
def test_a_march_indexes_its_per_node_tau(tmp_path):
    """tracervectors_spatial dumps each particle's local time as tau, the name
    a strip or sheet gives its per-edge or per-face tau: indexed by where the
    run wrote it, it is a node attribute, as every other field of the point
    cloud."""
    run_app(SPATIAL, copy_example(POISEUILLE, tmp_path / "c"),
            "mode=analytic init_mode=points_x x0=0 y0=0 z0=0 Nrw=50 Nrw_max=50 stat_intv=1e9 "
            "checkpoint_intv=1e9 random=false seed=1 dxn=0.01 dump_intv=0.05 Ln=0.1 T=1e9")
    folder = run_folder(tmp_path / "c")
    index = check_index(folder)
    assert index == dumped(folder) and len(index) == 3
    for g in grids(folder):
        assert {a.get("Name"): a.get("Center") for a in g.findall("Attribute")}["tau"] == "Node"


@needs(TRACERS)
def test_chunked_dumps_are_indexed_in_their_files(tmp_path):
    """With two dumps a file, each grid names the file its group is in."""
    run_app(TRACERS, copy_example(POISEUILLE, tmp_path / "c"), POINTS,
            "T=0.5 dump_intv=0.1 dump_chunk_size=2")
    folder = run_folder(tmp_path / "c")
    index = check_index(folder)
    assert index == dumped(folder)
    assert [f for _, f in index] == ["data_from_t%f.h5" % t for t in (0., 0., .2, .2, .4, .4)]


@needs(TRACERS)
def test_restarts_continue_one_history(tmp_path):
    """Runs stopped and resumed, each time to a later T, then a resume from an
    earlier checkpoint: after every run the index is complete and valid, and
    after the last it holds the one history, with no grid of the abandoned
    run past the restart and none twice."""
    params = copy_example(POISEUILLE, tmp_path / "c")
    args = [POINTS, "dump_intv=0.1 dump_chunk_size=4"]
    run_app(TRACERS, params, args, "T=0.2")
    folder = run_folder(tmp_path / "c")
    assert check_index(folder) == dumped(folder)
    # the final checkpoint is one step past T
    saved = tmp_path / "saved"
    shutil.copytree(folder / "Checkpoints", saved)
    run_app(TRACERS, params, args, "T=0.5", "restart_folder=%s" % checkpoint_folder(tmp_path / "c"))
    index = check_index(folder)
    assert [t for t, _ in index] == pytest.approx([0.1 * k for k in range(6)])
    run_app(TRACERS, params, args, "T=0.8", "restart_folder=%s" % folder)
    assert [t for t, _ in check_index(folder)] == pytest.approx([0.1 * k for k in range(9)])
    # back to t = 0.3, to a shorter T
    shutil.rmtree(folder / "Checkpoints")
    shutil.copytree(saved, folder / "Checkpoints")
    run_app(TRACERS, params, args, "T=0.6", "restart_folder=%s" % folder)
    index = check_index(folder)
    assert [t for t, _ in index] == pytest.approx([0.1 * k for k in range(7)])
    assert [f for t, f in index if t > 0.25] == ["data_from_t0.300000.h5"] + ["data_from_t0.400000.h5"] * 3


@needs(TRACERS)
def test_a_restart_on_a_chunk_boundary_dumps_into_its_start_file(tmp_path):
    """A restart whose first step starts a new chunk writes that dump into the
    file it opened at the start, not a second file of the same name, and the
    chunks after it follow as before."""
    params = copy_example(POISEUILLE, tmp_path / "c")
    args = [POINTS, "dump_intv=0.1 dump_chunk_size=2"]
    # the final checkpoint, at t = 0.4, is on a chunk boundary
    run_app(TRACERS, params, args, "T=0.3")
    folder = run_folder(tmp_path / "c")
    run_app(TRACERS, params, args, "T=0.7", "restart_folder=%s" % folder)
    index = check_index(folder)
    assert index == dumped(folder)
    assert [t for t, _ in index] == pytest.approx([0.1 * k for k in range(8)])
    assert [f for _, f in index] == ["data_from_t%f.h5" % t for t in (0., 0., .2, .2, .4, .4, .6, .6)]


@needs(TRACERS)
def test_a_restart_that_does_not_dump_at_its_start_drops_the_grids_of_the_file_it_empties(tmp_path):
    """A restart on a chunk boundary recreates, empty, the chunk file that held
    the abandoned run's dump at the restart time. With a longer dump interval
    the restart does not dump there again, and the index must keep no grid in
    that file: every grid names a dump that exists."""
    params = copy_example(POISEUILLE, tmp_path / "c")
    args = [POINTS, "dump_intv=0.1 dump_chunk_size=2"]
    run_app(TRACERS, params, args, "T=0.5")
    folder = run_folder(tmp_path / "c")
    # the same run's checkpoint at t = 0.4, from a shorter copy
    other = copy_example(POISEUILLE, tmp_path / "x")
    run_app(TRACERS, other, args, "T=0.3")
    shutil.rmtree(folder / "Checkpoints")
    shutil.copytree(checkpoint_folder(tmp_path / "x") / "Checkpoints", folder / "Checkpoints")
    run_app(TRACERS, params, args, "T=0.5 dump_intv=0.5", "restart_folder=%s" % folder)
    index = check_index(folder)
    assert index == dumped(folder)
    assert index == [(0., "data_from_t0.000000.h5"), (.1, "data_from_t0.000000.h5"), (.2, "data_from_t0.200000.h5"),
                     (.3, "data_from_t0.200000.h5"), (.5, "data_from_t0.400000.h5")]


@needs(TRACERS)
def test_a_restart_that_never_steps_keeps_the_grid_at_its_start(tmp_path):
    """A restart whose T is before its checkpoint's time runs no step, so it
    dumps nothing; the grid at the restart time, in an older file, is no
    dump it replaces and stays: the index is the files' dumps up to then."""
    params = copy_example(POISEUILLE, tmp_path / "c")
    args = [POINTS, "dump_intv=0.1"]
    run_app(TRACERS, params, args, "T=0.6")
    folder = run_folder(tmp_path / "c")
    # the same run's checkpoint at t = 0.3, from a shorter copy
    other = copy_example(POISEUILLE, tmp_path / "x")
    run_app(TRACERS, other, args, "T=0.2")
    shutil.rmtree(folder / "Checkpoints")
    shutil.copytree(checkpoint_folder(tmp_path / "x") / "Checkpoints", folder / "Checkpoints")
    run_app(TRACERS, params, args, "T=0.1", "restart_folder=%s" % folder)
    index = check_index(folder)
    assert index == [d for d in dumped(folder) if d[0] < 0.35]
    assert [t for t, _ in index] == pytest.approx([0., .1, .2, .3])


@needs(TRACERS)
@pytest.mark.parametrize("damage", ["missing", "cut"])
def test_a_restart_rebuilds_a_missing_or_damaged_index(tmp_path, damage):
    """With no mesh.xdmf (a run from before the index) or one cut short
    mid-grid, the resume rebuilds it from the dump files. Where two files hold
    a time, an abandoned run's and the one resumed over it, the history is
    the file started last but not after that time, whatever the files'
    modification times say: they are scrambled here."""
    params = copy_example(POISEUILLE, tmp_path / "c")
    args = [POINTS, "dump_intv=0.1"]
    run_app(TRACERS, params, args, "T=0.6")
    folder = run_folder(tmp_path / "c")
    # the same run's checkpoint at t = 0.3, from a shorter copy
    other = copy_example(POISEUILLE, tmp_path / "x")
    run_app(TRACERS, other, args, "T=0.2")
    shutil.rmtree(folder / "Checkpoints")
    shutil.copytree(checkpoint_folder(tmp_path / "x") / "Checkpoints", folder / "Checkpoints")
    run_app(TRACERS, params, args, "T=0.5", "restart_folder=%s" % folder)
    # the abandoned run's file the newest, the resumed one's the oldest
    os.utime(folder / "data_from_t0.000000.h5", (2e9, 2e9))
    os.utime(folder / "data_from_t0.300000.h5", (1e9, 1e9))
    index = folder / "mesh.xdmf"
    if damage == "missing":
        index.unlink()
    else:
        text = index.read_text()
        index.write_text(text[:text.index("<Attribute", text.index('<Time Value="0.200000"'))])
    run_app(TRACERS, params, args, "T=0.8", "restart_folder=%s" % folder)
    index = check_index(folder)
    assert [t for t, _ in index] == pytest.approx([0.1 * k for k in range(9)])
    assert [f for _, f in index] == (["data_from_t0.000000.h5"] * 3 + ["data_from_t0.300000.h5"] * 3
                                     + ["data_from_t0.600000.h5"] * 3)
