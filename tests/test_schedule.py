"""dump_intv, stat_intv and checkpoint_intv as schedules, through a run.

A schedule is comma-separated segments in absolute time: X every X to the end,
X:U every X until U, logN:S:U and logN:U N a decade (src/base/schedule.hpp).
Each time is snapped to the nearest step and the times on one step are one
event, and the run's start t0 is one too. Here dt = 0.003 divides neither 0.1 nor the log times, so every dump,
statistics row and index entry must sit on the step nearest its time; a
restart, mid-uniform or mid-log, must give the dumps of the run it splits; the
schedule string goes through the parameter files unchanged; a schedule whose
first segment starts later dumps the start positions at t0; and dump_chunk_size
counts a schedule's dumps.
"""

import math
import os
import pathlib
import xml.etree.ElementTree as ET

import pytest

from dumps import all_dumps, read_stats
from paths import REPO, app, needs
from runs import checkpoint_folder, continuous_and_resumed, copy_example, run_app

TRACERS = app("tracers")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
DT = 0.003
POINTS = ("mode=analytic init_mode=points_xyz x0=0 y0=0 z0=0 Nrw=20 Nrw_max=20 Dm=0 int_order=1 "
          "scheme=RK4 dt=%r random=false seed=1 checkpoint_intv=1e9" % DT)
DUMPS = "0.1:1,log5:10,1"
STATS = "log10:0.1:10,2"
needs_tracers = needs(TRACERS)


def uniform(a, X, U):
    """a + n X below U, then U."""
    out = []
    while a + len(out) * X < U - 1e-9:
        out.append(a + len(out) * X)
    return out + [U]


def log(N, S, U):
    """S 10^(k/N) up to U."""
    out = []
    k = 0
    while S * 10 ** (k / N) <= U * (1 + 1e-12):
        out.append(S * 10 ** (k / N))
        k += 1
    return out


def snapped(times, T, dt=DT, t0=0.):
    """The times t0 + k dt of the steps nearest the times from t0 to T, once
    each; half way, the later."""
    steps = sorted({math.floor((t - t0) / dt + 0.5 + 1e-6) for t in times if t >= t0})
    return [t0 + k * dt for k in steps if t0 + k * dt <= T + dt / 2]


def dump_schedule(T):
    """The dump times of DUMPS to T."""
    return snapped(uniform(0., 0.1, 1.) + log(5, 1., 10.) + uniform(10., 1., T + 1)[:-1], T)


def stat_schedule(T):
    """The statistics times of STATS to T: t0 = 0, then from 0.1."""
    return snapped([0.] + log(10, 0.1, 10.) + uniform(10., 2., T + 2)[:-1], T)


def index_times(d):
    """The times of mesh.xdmf's grids under d."""
    [index] = list(pathlib.Path(d).rglob("mesh.xdmf"))
    [collection] = ET.parse(index).getroot().find("Domain").findall("Grid")
    return [float(g.find("Time").get("Value")) for g in collection.findall("Grid")]


@needs_tracers
def test_dumps_rows_and_index_follow_their_schedules(tmp_path):
    """Every 0.1 to 1, 5 a decade to 10, then every 1; statistics 10 a decade
    from 0.1 to 10, then every 2: the h5 groups, mesh.xdmf and the tdata rows
    are on the snapped times, the switches once."""
    T = 12.
    params = copy_example(POISEUILLE, tmp_path)
    run_app(TRACERS, params, POINTS, "T=%r dump_intv=%s stat_intv=%s" % (T, DUMPS, STATS))
    want = dump_schedule(T)
    assert len(want) == 11 + 5 + 2
    assert sorted(all_dumps(tmp_path)) == pytest.approx(want, abs=1e-9)
    assert index_times(tmp_path) == pytest.approx(want, abs=1e-9)
    rows = read_stats(tmp_path)
    assert list(rows["t"]) == pytest.approx(stat_schedule(T), abs=1e-9)
    # the start, then 0.1 and 0.126 at 33 and 42 steps
    assert list(rows["t"])[:3] == pytest.approx([0., 0.099, 0.126], abs=1e-9)


@needs_tracers
@pytest.mark.parametrize("stop", [0.45, 3.0])
def test_a_resume_mid_segment_dumps_as_the_run_it_splits(tmp_path, stop):
    """Stopped mid-uniform (0.45) or mid-log (3, between 2.512 and 3.981) and
    resumed: the same dumps at the same times, the schedule read back from the
    checkpoint's parameters."""
    T = 12.
    args = [POINTS, "dump_intv=%s stat_intv=%s" % (DUMPS, STATS)]
    cont, split = continuous_and_resumed(TRACERS, POISEUILLE, tmp_path, args,
                                         ["T=%r" % stop, "checkpoint_intv=%r" % stop],
                                         ["T=%r" % T])
    ck = checkpoint_folder(split)
    text = (ck / "Checkpoints" / "params.dat").read_text()
    assert "\ndump_intv=%s\n" % DUMPS in text and "\nstat_intv=%s\n" % STATS in text
    a, b = all_dumps(cont), all_dumps(split)
    assert sorted(b) == pytest.approx(sorted(a), abs=1e-9)
    assert sorted(a) == pytest.approx(dump_schedule(T), abs=1e-9)
    for t in a:
        assert (a[t]["points"] == b[t]["points"]).all(), t
    assert index_times(split) == pytest.approx(dump_schedule(T), abs=1e-9)


@needs_tracers
def test_a_schedule_ends_at_its_last_until(tmp_path):
    """With an until on the last segment nothing comes after it: dumps every
    0.1 to 0.5 in a run to 2."""
    params = copy_example(POISEUILLE, tmp_path)
    run_app(TRACERS, params, POINTS, "T=2 dump_intv=0.1:0.5 stat_intv=0.5:1")
    assert sorted(all_dumps(tmp_path)) == pytest.approx(snapped(uniform(0., 0.1, 0.5), 2.), abs=1e-9)
    assert list(read_stats(tmp_path)["t"]) == pytest.approx(snapped([0., 0.5, 1.], 2.), abs=1e-9)


@needs_tracers
def test_a_checkpoint_schedule_resumes(tmp_path):
    """checkpoint_intv takes a schedule too, written to and read back from the checkpoint."""
    params = copy_example(POISEUILLE, tmp_path)
    run_app(TRACERS, params, POINTS, "T=1 dump_intv=0.5:2 stat_intv=0 checkpoint_intv=0.1:0.3,log5:10")
    ck = checkpoint_folder(tmp_path)
    assert "\ncheckpoint_intv=0.1:0.3,log5:10\n" in (ck / "Checkpoints" / "params.dat").read_text()
    run_app(TRACERS, params, POINTS, "T=2 dump_intv=0.5:2 stat_intv=0", "restart_folder=%s" % ck)
    assert sorted(all_dumps(tmp_path)) == pytest.approx(snapped([0., 0.5, 1.0, 1.5, 2.0], 2.), abs=1e-9)


@needs_tracers
@pytest.mark.parametrize("dumps, first", [("log20:0.1:0.8,0.5", 0.1), ("log20:1:1.6,0.5", 1.)])
def test_a_schedule_that_starts_later_dumps_the_start(tmp_path, dumps, first):
    """The first segment starts at 0.1 or 1, after t0 = 0: the start positions
    are dumped at t0 all the same, then nothing until the first segment."""
    main = tmp_path / "main"
    params = copy_example(POISEUILLE, main)
    run_app(TRACERS, params, POINTS, "T=2 t0=0 dump_intv=%s stat_intv=%s" % (dumps, dumps))
    times = sorted(all_dumps(main))
    assert times[:2] == pytest.approx([0., snapped([first], 2.)[0]], abs=1e-9)
    assert list(read_stats(main)["t"])[:2] == pytest.approx(times[:2], abs=1e-9)
    # the start positions, as a bare interval dumps them at t0
    other = copy_example(POISEUILLE, tmp_path / "bare")
    run_app(TRACERS, other, POINTS, "T=0.003 dump_intv=1e9 stat_intv=0")
    assert (all_dumps(tmp_path / "bare")[0.]["points"] == all_dumps(main)[0.]["points"]).all()


@needs_tracers
def test_a_resume_dumps_no_extra_start(tmp_path):
    """On a restart t0 is the run's start: the resumed run, whose start is no
    event, dumps only the scheduled times after it."""
    params = copy_example(POISEUILLE, tmp_path)
    args = [POINTS, "dump_intv=log10:0.1:1,0.5"]
    run_app(TRACERS, params, args, "T=0.45 checkpoint_intv=0.45")
    run_app(TRACERS, params, args, "T=2", "restart_folder=%s" % checkpoint_folder(tmp_path))
    want = snapped([0.] + log(10, 0.1, 1.) + [1.5, 2.], 2.)
    assert sorted(all_dumps(tmp_path)) == pytest.approx(want, abs=1e-9)


@needs_tracers
def test_chunks_count_a_schedules_dumps(tmp_path):
    """dump_chunk_size counts a schedule's dumps, t0's first: the start file
    takes the first four."""
    h5py = pytest.importorskip("h5py")
    params = copy_example(POISEUILLE, tmp_path)
    run_app(TRACERS, params, POINTS, "T=2 dump_intv=log10:0.1:1,0.5 dump_chunk_size=4 stat_intv=0")
    want = snapped([0.] + log(10, 0.1, 1.) + [1.5, 2.], 2.)
    assert len(want) == 14
    files = sorted(tmp_path.rglob("data_from_t*.h5"))
    sizes = []
    for f in files:
        with h5py.File(f, "r") as h:
            sizes.append(len(h))
    assert sizes == [4, 4, 4, 2]
    assert files[0].name == "data_from_t0.000000.h5"
    assert [f.name for f in files[1:]] == ["data_from_t%f.h5" % want[k] for k in (4, 8, 12)]


@needs_tracers
def test_a_bad_schedule_is_refused_with_its_segment_and_rule(tmp_path):
    """The parameters name the key, the segment and the rule (the rules: test_schedule.cpp)."""
    params = copy_example(POISEUILLE, tmp_path)
    r = run_app(TRACERS, params, POINTS, "T=1 dump_intv=0.1:1,abc", check=False)
    assert r.returncode == 2, r.stdout + r.stderr
    assert "dump_intv: 'abc': not a schedule" in r.stderr, r.stderr


@needs_tracers
def test_chunks_count_from_the_run_start(tmp_path):
    """t0 = 0.5 and 10 a decade from 0.1: the times below t0 are no dumps and
    not counted, so four dumps a file give files of 4, 4 and 2."""
    h5py = pytest.importorskip("h5py")
    params = copy_example(POISEUILLE, tmp_path)
    run_app(TRACERS, params, POINTS, "t0=0.5 T=3 dump_intv=log10:0.1:2,0.5 dump_chunk_size=4 stat_intv=0")
    want = snapped([0.5] + log(10, 0.1, 2.) + [2., 2.5, 3.], 3., t0=0.5)
    assert len(want) == 10
    assert sorted(all_dumps(tmp_path)) == pytest.approx(want, abs=1e-9)
    sizes = []
    for f in sorted(tmp_path.rglob("data_from_t*.h5")):
        with h5py.File(f, "r") as h:
            sizes.append(len(h))
    assert sizes == [4, 4, 2]


@needs_tracers
def test_a_resume_with_another_step_passes_the_times_before_it(tmp_path):
    """Stopped at 0.453 with dt = 0.003 and resumed with dt = 0.2: the dump at
    0.4, nearest the resumed run's first step, was written before; the resumed
    run dumps 0.6 and 0.8 on its own steps, and 1 too."""
    params = copy_example(POISEUILLE, tmp_path)
    args = [POINTS, "dump_intv=0.2:1,10 stat_intv=0"]
    run_app(TRACERS, params, args, "T=0.45 checkpoint_intv=0.45")
    run_app(TRACERS, params, args, "dt=0.2 T=1.2", "restart_folder=%s" % checkpoint_folder(tmp_path))
    assert sorted(all_dumps(tmp_path)) == pytest.approx([0., 0.201, 0.399, 0.653, 0.853, 1.053], abs=1e-9)
