"""The checkpoint file: complete, written whole, and refused when it does not fit.

A checkpoint is Checkpoints/checkpoint.h5 beside the text params.dat. Both are
written in full beside the last ones as checkpoint.h5.tmp and params.dat.tmp,
then moved over them, checkpoint.h5 first and params.dat last. A run killed or
failing while writing keeps its previous checkpoint; only a kill between the
two moves leaves a pair from different checkpoints, which a resume refuses,
naming params.dat.tmp as the params.dat that belongs to the new checkpoint.h5.
A resume that finds such a pair, no checkpoint at all, a dataset missing or of
the wrong length, or an index past what it indexes, stops with a message
instead of running on from a state nobody wrote. The runs are a sheet with
injection in partrac, which writes every dataset but the carried elements', on
the plane Poiseuille example; restarts matching an uninterrupted run are
tested with each app in test_restart, test_transport, test_intervals,
test_march and test_walkers.
"""

import os
import shutil

import numpy as np
import pytest

from paths import REPO, app
from runs import (checkpoint_file, checkpoint_folder, copy_example, read_checkpoint, run_app,
                  write_checkpoint, write_text_checkpoint)

PARTRAC = app("partrac")
TRACERS = app("tracers")
HAGEN = os.path.join(REPO, "data_example", "hagen_poiseuille", "expr_params.dat")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

BASE = ("mode=analytic init_mode=uniform_x x0=0 y0=0 z0=0 Nrw=21 Nrw_max=200000 "
        "inject=true inject_edges=true T_inject=1e10 inject_intv=0.05 "
        "ds_max=0.2 ds_min=0.02 refine=true coarsen=true Dm=0 int_order=1 "
        "dt=0.01 stat_intv=1e9 dump_intv=1e9 checkpoint_intv=0.03 T=0.2 "
        "random=false seed=1").split()

# every dataset of this run's checkpoint, and the column count of each
DATASETS = {"points": 3, "id": 1, "c": 1,
            "edges": 2, "dl0": 1, "edge_tau": 1, "edge_rho_prev": 1,
            "face_edges": 3, "dA0": 1, "face_tau": 1, "face_rho_prev": 1,
            "positions_inj": 3, "edges_inj": 2, "dl0_inj": 1, "edge_tau_inj": 1,
            "edge_rho_prev_inj": 1, "edges_inlet": 1, "nodes_inlet": 1}

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
needs_tracers = pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")


def params_t(folder, name="params.dat"):
    """t in the params.dat (or `name`) of the Checkpoints folder."""
    for line in (folder / name).read_text().splitlines():
        if line.startswith("t="):
            return float(line[2:])
    raise AssertionError("no t in " + str(folder / name))


def resume(case, *args, check=True, restart=None):
    """Resume the run in case from its checkpoint (or the folder restart) to
    T=0.3, with args overriding; the process."""
    return run_app(PARTRAC, case / "expr_params.dat", BASE, "T=0.3", *args,
                   "restart_folder=%s" % (restart or checkpoint_folder(case)), check=check)


def refused(r, message):
    """Assert the run r stopped with exit code 2 and message in its stderr."""
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr[-1000:]
    assert message in r.stderr, r.stderr


@pytest.fixture
def case(tmp_path):
    """A sheet with injection run to T, checkpointing every 3 steps; its case folder."""
    run_app(PARTRAC, copy_example(HAGEN, tmp_path / "case"), BASE)
    return tmp_path / "case"


@needs_partrac
def test_a_run_leaves_one_whole_checkpoint_and_nothing_half_written(case):
    """After a run that checkpointed several times the Checkpoints folder holds
    checkpoint.h5 and params.dat only: no temporary file, no text files. The
    checkpoint has every dataset, one row per particle, edge, face or inlet
    entry, every edge on existing nodes and every face on existing edges, and
    its time is the one params.dat resumes at."""
    h5py = pytest.importorskip("h5py")
    f = checkpoint_file(case)
    assert sorted(p.name for p in f.parent.iterdir()) == ["checkpoint.h5", "params.dat"]
    assert not list(case.rglob("*.tmp"))
    ck = read_checkpoint(case)
    assert set(ck) == set(DATASETS)
    for k, cols in DATASETS.items():
        assert ck[k].ndim == 2 and ck[k].shape[1] == cols, k
    n, ne, nf = len(ck["points"]), len(ck["edges"]), len(ck["face_edges"])
    assert nf > 0, "the sheet has no faces: the test would not see them"
    for k in ("id", "c"):
        assert len(ck[k]) == n, k
    for k in ("dl0", "edge_tau", "edge_rho_prev"):
        assert len(ck[k]) == ne, k
    for k in ("dA0", "face_tau", "face_rho_prev"):
        assert len(ck[k]) == nf, k
    for k in ("dl0_inj", "edge_tau_inj", "edge_rho_prev_inj"):
        assert len(ck[k]) == len(ck["edges_inj"]), k
    assert ck["edges"].max() < n and ck["face_edges"].max() < ne
    assert ck["edges_inj"].max() < len(ck["positions_inj"])
    assert ck["nodes_inlet"].max() < n and ck["edges_inlet"].max() < ne
    assert sorted(ck["id"][:, 0]) == sorted(set(ck["id"][:, 0]))
    with h5py.File(f, "r") as h:
        assert h.attrs["t"] == params_t(f.parent)


@needs_partrac
def test_a_failed_checkpoint_write_keeps_the_previous_checkpoint(case, tmp_path):
    """A checkpoint whose checkpoint.h5.tmp cannot be written stops the run and
    leaves checkpoint.h5 and params.dat as they were, params.dat.tmp beside
    them; so does a kill before the two files are moved in place. Resumed from
    that folder, with a torn checkpoint.h5.tmp in it as a kill leaves, the run
    ends where an uninterrupted one does. Were params.dat moved first, the
    folder would pair the new params.dat with the old checkpoint.h5, and the
    run could not be resumed at all."""
    folder = checkpoint_file(case).parent
    h5, params = (folder / "checkpoint.h5").read_bytes(), (folder / "params.dat").read_text()
    (folder / "checkpoint.h5.tmp").mkdir()
    refused(resume(case, check=False), "cannot write the checkpoint")
    assert (folder / "checkpoint.h5").read_bytes() == h5
    assert (folder / "params.dat").read_text() == params
    assert (folder / "params.dat.tmp").exists()

    (folder / "checkpoint.h5.tmp").rmdir()
    (folder / "checkpoint.h5.tmp").write_bytes(b"not a checkpoint")
    resume(case)
    assert not list(case.rglob("*.tmp"))
    cont = copy_example(HAGEN, tmp_path / "cont")
    run_app(PARTRAC, cont, BASE, "T=0.3")
    a, b = read_checkpoint(cont.parent), read_checkpoint(case)
    assert set(a) == set(b)
    for k in a:
        assert np.array_equal(a[k], b[k]), k
    assert params_t(folder) == params_t(checkpoint_file(cont.parent).parent)


@needs_partrac
@pytest.mark.skipif(not os.path.exists("/dev/full"), reason="no /dev/full to fail a write")
def test_a_params_file_that_cannot_be_written_replaces_nothing(case):
    """A params.dat.tmp whose write fails (here: a link to /dev/full, which
    refuses every byte) stops the run with a message, and params.dat and
    checkpoint.h5 are left as they were. Moved in place unchecked, the empty
    file would be the params.dat a resume reads."""
    folder = checkpoint_file(case).parent
    h5, params = (folder / "checkpoint.h5").read_bytes(), (folder / "params.dat").read_text()
    os.symlink("/dev/full", folder / "params.dat.tmp")
    r = resume(case, check=False)
    assert r.returncode != 0, r.stdout[-1000:] + r.stderr[-1000:]
    assert "could not write" in r.stderr, r.stderr
    assert not (folder / "params.dat").is_symlink()
    assert (folder / "params.dat").read_text() == params
    assert (folder / "checkpoint.h5").read_bytes() == h5


@needs_partrac
@pytest.mark.parametrize("newer", ["checkpoint.h5", "params.dat"])
def test_a_pair_from_two_checkpoints_is_refused(case, newer):
    """checkpoint.h5 one checkpoint newer than params.dat, as a kill between
    moving the one and the other leaves, or params.dat newer than
    checkpoint.h5: the resume stops with exit code 2, naming both times. In
    the first case params.dat.tmp holds the params.dat written with the
    checkpoint, and the message says so, so the pair can be put together by
    hand."""
    folder = checkpoint_file(case).parent
    h5, params = (folder / "checkpoint.h5").read_bytes(), (folder / "params.dat").read_text()
    resume(case)
    if newer == "checkpoint.h5":
        shutil.copy(folder / "params.dat", folder / "params.dat.tmp")
        (folder / "params.dat").write_text(params)
    else:
        (folder / "checkpoint.h5").write_bytes(h5)
    r = resume(case, "T=0.4", check=False)
    refused(r, "its params.dat at t =")
    hint = "params.dat.tmp beside it is the params.dat written with it"
    assert (hint in r.stderr) == (newer == "checkpoint.h5"), r.stderr


@needs_partrac
def test_a_checkpoints_folder_without_a_checkpoint_is_refused(case):
    """A Checkpoints folder with params.dat but neither checkpoint.h5 nor the
    text files of an older run stops the resume with exit code 2. Read as a
    text checkpoint with its files missing, it would resume with no particles."""
    restart = checkpoint_folder(case)
    checkpoint_file(case).unlink()
    refused(resume(case, check=False, restart=restart), "no checkpoint.h5 or positions.pos")


@needs_tracers
def test_a_text_checkpoint_larger_than_the_particle_budget_is_refused(tmp_path):
    """A text checkpoint of 100 particles resumed with Nrw_max=20 stops with
    exit code 2, as an HDF5 one does, instead of writing past the particle
    arrays."""
    args = ("mode=analytic init_mode=points_z x0=0 y0=0 z0=0 Nrw=100 Nrw_max=100 Dm=0 "
            "int_order=1 stat_intv=1e9 dump_intv=1e9 random=false seed=1 dt=0.1 T=0.1 "
            "checkpoint_intv=0.1")
    d = copy_example(POISEUILLE, tmp_path / "text").parent
    run_app(TRACERS, d / "expr_params.dat", args)
    assert len(read_checkpoint(d)["points"]) == 100
    restart = checkpoint_folder(d)
    write_text_checkpoint(restart / "Checkpoints")
    r = run_app(TRACERS, d / "expr_params.dat", args, "T=0.2 Nrw=20 Nrw_max=20",
                "restart_folder=%s" % restart, check=False)
    refused(r, "100 particles, more than Nrw_max = 20")


@needs_partrac
@pytest.mark.parametrize("damage,message", [
    ("t", "its params.dat at t ="),
    ("rows", "'dl0' is"),
    ("missing", "no dataset 'c'"),
    ("edge", "joins nodes"),
    ("nodes_inlet", "of 'nodes_inlet' names node"),
    ("edges_inlet", "of 'edges_inlet' names edge"),
    ("garbage", "cannot read the checkpoint"),
])
def test_a_checkpoint_that_does_not_fit_is_refused(case, damage, message):
    """A resume stops with exit code 2 and a message naming the fault when
    checkpoint.h5 and params.dat are from different checkpoints, a dataset is
    one row short or missing, an edge names a node the checkpoint does not
    have, an inlet entry names a node or edge past the end, or the file is not
    HDF5 at all, which is reported without the HDF5 library's error stack.
    Read as it is, each would run on from a state no run wrote."""
    h5py = pytest.importorskip("h5py")
    f = checkpoint_file(case)
    ck = read_checkpoint(case)
    if damage == "t":
        with h5py.File(f, "r+") as h:
            h.attrs["t"] = h.attrs["t"] - 0.01
    elif damage == "rows":
        with h5py.File(f, "r+") as h:
            del h["dl0"]
            h.create_dataset("dl0", data=ck["dl0"][:-1])
    elif damage == "missing":
        with h5py.File(f, "r+") as h:
            del h["c"]
    elif damage == "edge":
        e = ck["edges"]
        e[0, 1] = len(ck["points"])
        write_checkpoint(case, edges=e)
    elif damage == "nodes_inlet":
        v = ck["nodes_inlet"]
        v[0] = len(ck["points"])
        write_checkpoint(case, nodes_inlet=v)
    elif damage == "edges_inlet":
        v = ck["edges_inlet"]
        v[0] = len(ck["edges"])
        write_checkpoint(case, edges_inlet=v)
    else:
        f.write_bytes(b"not a checkpoint")
    r = resume(case, check=False)
    refused(r, message)
    assert "HDF5-DIAG" not in r.stderr, r.stderr
