"""Parameter validation and initializer construction, via partrac --check.

--check parses the parameters, builds the interpolator, particle set, topology
and initializer, then stops before the time loop without writing anything. That
exercises the parameter reads inside each initializer, which are otherwise only
reached when someone runs that particular init_mode.
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

# invalid input that must be rejected before anything runs
BAD_CASES = [
    (["init_mode=ellipsoid_xy"], "La"),                         # conditionally required
    (["init_mode=strip_x"], "La"),                              # conditionally required
    (["init_mode=uniform_x", "exit_plane=x"], "filter_intv"),   # conditionally required
    (["init_mode=uniform_x", "inject=true"], "inject_intv"),    # conditionally required
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
    """A scratch copy of the example: output would land next to the input file."""
    shutil.copy(EXAMPLE, tmp_path / "expr_params.dat")
    return tmp_path


def run_check(case_dir, args):
    return subprocess.run([PARTRAC, str(case_dir / "expr_params.dat"), "--check"] + args,
                          capture_output=True, text=True, timeout=300)


def run_full(case_dir, args):
    return subprocess.run([PARTRAC, str(case_dir / "expr_params.dat")] + args,
                          capture_output=True, text=True, timeout=300)


def dumped_params(case_dir):
    files = list(case_dir.glob("RandomWalkers/*/0/params_from_t*.dat"))
    assert len(files) == 1, files
    return dict(l.split("=", 1) for l in files[0].read_text().splitlines() if "=" in l)


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
@pytest.mark.parametrize("init_mode,extra", INIT_MODES, ids=[m for m, _ in INIT_MODES])
def test_init_mode_constructs(case, init_mode, extra):
    r = run_check(case, BASE + ["init_mode=" + init_mode] + extra)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "Check OK" in r.stdout


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
@pytest.mark.parametrize("extra,expected", BAD_CASES,
                         ids=[e[0].split("=")[0] + ":" + k for e, k in BAD_CASES])
def test_bad_parameters_are_rejected(case, extra, expected):
    r = run_check(case, BASE + extra)
    assert r.returncode != 0
    assert expected in r.stderr


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_dry_run_writes_nothing(case):
    before = set(os.listdir(case))
    r = run_check(case, BASE + ["init_mode=uniform_x"])
    assert r.returncode == 0
    assert set(os.listdir(case)) == before


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_dump_interval_is_clamped_to_the_timestep(case):
    # dump_intv/dt becomes an integer step count, so a dump_intv below dt has to
    # be raised to dt rather than rounding down to zero. Needs a real run: the
    # step count is only used inside the time loop, which --check skips.
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
    # Nrw used to be overwritten with the number the initializer managed to
    # place, so a dumped parameter file replayed as a smaller run every time.
    # It stays the request; Nrw_init and Nrw_current carry what happened.
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
def test_help_lists_required_parameters(case):
    r = subprocess.run([PARTRAC, "--help"], capture_output=True, text=True, timeout=60)
    assert r.returncode == 0
    assert "Required:" in r.stdout
    for key in ("init_mode", "Nrw_max", "ds_max", "mode"):
        assert key in r.stdout
