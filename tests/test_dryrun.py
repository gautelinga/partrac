"""Parameter validation and initializer construction in partrac, mostly via --check.

--check parses the parameters, builds the interpolator, particle set, topology
and initializer, then stops before the time loop without writing anything. That
exercises the parameter reads inside each initializer, which are otherwise only
reached when someone runs that particular init_mode, and lets invalid input be
rejected before any compute time is spent. The input is the plane Poiseuille
example; a few tests need a short full run because the quantity they check is
only used inside the time loop.
"""

import os
import shutil
import subprocess

import pytest

from paths import REPO, app

PARTRAC = app("partrac")
EXAMPLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

# domain is [-1, 1]^3
BASE = ("mode=analytic Nrw=50 Nrw_max=5000 ds_max=0.4 ds_min=0.1 "
        "Dm=0 dt=0.01 T=0.1 int_order=1").split()

# every init_mode, with the parameters it needs to construct
INIT_MODES = [
    ("point", []),
    ("uniform_x", []),
    ("uniform_y", []),
    ("uniform_z", []),
    ("strip_x", ["La=0.5"]),
    ("sheet_xy", ["La=0.5", "Lb=0.5", "ds_init=0.2"]),
    ("ellipsoid_xy", ["La=0.5", "Lb=0.5"]),
    ("pair_xyz", ["ds_init=0.1"]),
    ("pairs_xyz", ["ds_init=0.1"]),
    ("points_xy", ["init_weight=none", "ds_init=0.1"]),
    ("randomgaussianstrip_x_y", ["La=0.5", "Lb=0.1"]),
    ("randomgaussiancircle_xy", ["La=0.5", "Lb=0.1"]),
]

# invalid input that must be rejected before anything runs, with the key the
# error message must name
BAD_CASES = [
    (["init_mode=ellipsoid_xy"], "La"),                         # conditionally required
    (["init_mode=strip_x"], "La"),                              # conditionally required
    (["init_mode=uniform_x", "exit_plane=x"], "filter_intv"),   # conditionally required
    (["init_mode=uniform_x", "inject=true"], "inject_intv"),    # conditionally required
    (["init_mode=uniform_x", "Nrw=1"], "Nrw"),                  # no interval to step
    (["init_mode=uniform_x", "ds_ini=0.1"], "ds_init"),         # typo
    (["init_mode=uniform_x", "nx=0"], "nx"),                    # unknown key
    (["init_mode=uniform_x", "Nrw_max=-1"], "Nrw_max"),         # negative size
    (["init_mode=uniform_x", "verbose=ture"], "verbose"),       # not a boolean
    (["init_mode=uniform_x", "Dm=10meters"], "Dm"),             # not a number
    (["init_mode=uniform"], "init_mode"),                       # key[1] out of bounds
    (["init_mode=randomgaussianstrip_x", "La=1", "Lb=0.1"], "init_mode"),  # needs key[2]
    (["init_mode=from_nowhere:x"], "init_mode"),                # left init_state null
    (["init_mode=from_file"], "init_mode"),                     # no path
]


@pytest.fixture
def case(tmp_path):
    """A scratch copy of the example, since output lands next to the input file."""
    shutil.copy(EXAMPLE, tmp_path / "expr_params.dat")
    return tmp_path


def run_check(case_dir, args):
    """Run partrac --check on the case; return the process result."""
    return subprocess.run([PARTRAC, str(case_dir / "expr_params.dat"), "--check"] + args,
                          capture_output=True, text=True, timeout=300)


def run_full(case_dir, args):
    """Run partrac on the case through its time loop; return the process result."""
    return subprocess.run([PARTRAC, str(case_dir / "expr_params.dat")] + args,
                          capture_output=True, text=True, timeout=300)


def dumped_params(case_dir):
    """The parameter file the run wrote, as key -> string value."""
    files = list(case_dir.glob("RandomWalkers/*/0/params_from_t*.dat"))
    assert len(files) == 1, files
    return dict(l.split("=", 1) for l in files[0].read_text().splitlines() if "=" in l)


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
@pytest.mark.parametrize("init_mode,extra", INIT_MODES, ids=[m for m, _ in INIT_MODES])
def test_init_mode_constructs(case, init_mode, extra):
    """Every init_mode builds its initializer from valid parameters and --check
    reports success, so a user choosing any of them does not hit a parse error
    only after the rest of the setup has run."""
    r = run_check(case, BASE + ["init_mode=" + init_mode] + extra)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "Check OK" in r.stdout


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
@pytest.mark.parametrize("extra,expected", BAD_CASES,
                         ids=[e[0].split("=")[0] + ":" + k for e, k in BAD_CASES])
def test_bad_parameters_are_rejected(case, extra, expected):
    """Missing, misspelled, unknown, malformed or out-of-range parameters make
    --check fail with an error naming the offending key. Otherwise a typo would
    be silently ignored or a run would start with a meaningless value."""
    r = run_check(case, BASE + extra)
    assert r.returncode != 0
    assert expected in r.stderr


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_dry_run_writes_nothing(case):
    """--check leaves the case directory untouched, so validating a case never
    creates output folders or overwrites results from an earlier run."""
    before = set(os.listdir(case))
    r = run_check(case, BASE + ["init_mode=uniform_x"])
    assert r.returncode == 0
    assert set(os.listdir(case)) == before


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_dump_interval_is_clamped_to_the_timestep(case):
    """An output interval shorter than dt is raised to dt, and the raised value is
    recorded. Intervals become integer step counts, and rounding down to zero
    steps would never write or divide by zero."""
    # a full run: the step count is only used inside the time loop, which
    # --check skips
    args = [a for a in BASE if not a.startswith(("dt=", "T="))]
    r = run_full(case, args + ["init_mode=uniform_x", "dt=0.4", "T=1.2",
                               "dump_intv=0.1", "stat_intv=0.1",
                               "random=false", "seed=1"])
    assert r.returncode == 0, r.stdout + r.stderr
    prm = dumped_params(case)
    assert float(prm["dump_intv"]) == 0.4
    assert float(prm["stat_intv"]) == 0.4


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_nrw_is_an_input_and_the_counts_are_recorded_separately(case):
    """Nrw stays the requested count in the dumped parameters, while Nrw_init and
    Nrw_current record how many particles were placed and remain. A dumped
    parameter file must replay as the same run, not a smaller one."""
    args = [a for a in BASE if not a.startswith(("ds_max=", "ds_min=", "Nrw="))]
    # a strip longer than the domain, so part of it is dropped
    r = run_full(case, args + ["init_mode=strip_x", "Nrw=1000", "La=3.0",
                               "ds_max=1e9", "ds_min=1e-9", "refine=false",
                               "coarsen=false", "random=false", "seed=1",
                               "dump_intv=1e9", "stat_intv=1e9",
                               "checkpoint_intv=1e9"])
    assert r.returncode == 0, r.stdout + r.stderr
    prm = dumped_params(case)

    assert int(prm["Nrw"]) == 1000                  # the request survives
    assert int(prm["Nrw_init"]) < 1000              # some of the strip was dropped
    assert int(prm["Nrw_current"]) == int(prm["Nrw_init"])

    written = list(case.glob("**/Checkpoints/positions.pos"))
    assert len(written) == 1
    n = len([l for l in written[0].read_text().splitlines() if l.strip()])
    assert n == int(prm["Nrw_current"])


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_an_interval_of_zero_turns_that_output_off(case):
    """An interval of 0 turns that output off entirely, at t = 0 too, and writes
    no file; it is not clamped to every step and does not divide by zero. With
    positive intervals the same run does write the files."""
    args = [a for a in BASE if not a.startswith("T=")]
    r = run_full(case, args + ["T=0.05", "init_mode=uniform_x", "random=false",
                               "seed=1", "dump_intv=0", "stat_intv=0",
                               "checkpoint_intv=0"])
    assert r.returncode == 0, r.stdout + r.stderr
    written = {p.name for p in case.rglob("*") if p.is_file()}
    assert not [f for f in written if f.startswith(("data_from_t", "tdata_from_t"))]

    # control: with the intervals on, those files are there
    other = case / "on"
    other.mkdir()
    shutil.copy(EXAMPLE, other / "expr_params.dat")
    r = run_full(other, args + ["T=0.05", "init_mode=uniform_x", "random=false",
                                "seed=1", "dump_intv=0.05", "stat_intv=0.05"])
    assert r.returncode == 0, r.stdout + r.stderr
    on = {p.name for p in other.rglob("*") if p.is_file()}
    assert [f for f in on if f.startswith("data_from_t")]
    assert [f for f in on if f.startswith("tdata_from_t")]


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_a_negative_interval_is_rejected(case):
    """A negative output interval is rejected by --check with a message saying so,
    since it has no meaning as a step count."""
    r = run_check(case, BASE + ["init_mode=uniform_x", "dump_intv=-1"])
    assert r.returncode != 0
    assert "negative" in r.stderr


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_a_huge_interval_does_not_overflow_the_step_count(case):
    """A very large interval, the usual way to make an output effectively never
    happen, runs normally. Its step count must not overflow into a negative
    value that breaks the time loop."""
    # 1e9/0.005 = 2e11 steps does not fit a 32-bit int
    args = [a for a in BASE if not a.startswith(("dt=", "T="))]
    r = run_full(case, args + ["dt=0.005", "T=0.05", "init_mode=uniform_x",
                               "random=false", "seed=1", "dump_intv=1e9",
                               "stat_intv=1e9", "checkpoint_intv=1e9"])
    assert r.returncode == 0, r.stdout + r.stderr


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_help_lists_required_parameters(case):
    """--help succeeds and lists the required parameters, so a user can find out
    what a case file must contain without reading the source."""
    r = subprocess.run([PARTRAC, "--help"], capture_output=True, text=True, timeout=60)
    assert r.returncode == 0
    assert "Required:" in r.stdout
    for key in ("init_mode", "Nrw_max", "ds_max", "mode"):
        assert key in r.stdout
