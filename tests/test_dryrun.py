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

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARTRAC = os.path.join(REPO, "bin", "partrac")
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

# each of these used to be accepted silently and then misbehave at run time
BAD_CASES = [
    (["init_mode=ellipsoid_xy"], "La"),            # divided by zero at La = 0
    (["init_mode=strip_x"], "La"),                 # collapsed onto a single point
    (["init_mode=uniform_x", "exit_plane=x"], "filter_intv"),   # modulo by zero
    (["init_mode=uniform_x", "inject=true"], "inject_intv"),    # modulo by zero
    (["init_mode=uniform_x", "ds_ini=0.1"], "ds_init"),         # typo, was ignored
    (["init_mode=uniform_x", "nx=0"], "nx"),                    # removed key
    (["init_mode=uniform_x", "Nrw_max=-1"], "Nrw_max"),         # wrapped to 2^64-1
    (["init_mode=uniform_x", "verbose=ture"], "verbose"),       # was read as false
    (["init_mode=uniform_x", "Dm=10meters"], "Dm"),             # was read as 10
]


@pytest.fixture
def case(tmp_path):
    """A scratch copy of the example: output would land next to the input file."""
    shutil.copy(EXAMPLE, tmp_path / "expr_params.dat")
    return tmp_path


def run_check(case_dir, args):
    return subprocess.run([PARTRAC, str(case_dir / "expr_params.dat"), "--check"] + args,
                          capture_output=True, text=True, timeout=300)


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
def test_help_lists_required_parameters(case):
    r = subprocess.run([PARTRAC, "--help"], capture_output=True, text=True, timeout=60)
    assert r.returncode == 0
    assert "Required:" in r.stdout
    for key in ("init_mode", "Nrw_max", "ds_max", "mode"):
        assert key in r.stdout
