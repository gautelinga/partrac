"""A user's mistake that the parameter schema cannot see ends the run the way
a parameter error does: exit code 2 and a message on stderr, from the one
place each app reports (partrac::report_errors), not an exit() where it was
found. Any placement of the initial state that misses the domain is such a
mistake; the schema only knows the keys, not the domain."""

import os
import shutil
import subprocess

import pytest

from paths import REPO, app

PARTRAC = app("partrac")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
ARGS = ("mode=analytic dt=0.01 T=0.02 Nrw=10 Nrw_max=100 ds_max=0.1 ds_min=1e-9 "
        "ds_init=0.01 Dm=0 int_order=1 refine=false coarsen=false dump_intv=1e9 "
        "stat_intv=1e9 checkpoint_intv=1e9").split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")


def run(tmp_path, extra):
    shutil.copy(POISEUILLE, tmp_path / "expr_params.dat")
    return subprocess.run([PARTRAC, str(tmp_path / "expr_params.dat")] + ARGS + extra,
                          capture_output=True, text=True, timeout=120)


@needs_partrac
@pytest.mark.parametrize("extra,message", [
    (["init_mode=strip_x", "x0=100", "y0=100", "z0=100", "La=0.1"], "strip not inside domain"),
    (["init_mode=ellipsoid_z", "x0=100", "y0=100", "z0=100", "La=0.1", "Lb=0.1"],
     "ellipsoid not inside domain"),
])
def test_an_initial_state_outside_the_domain_is_reported(tmp_path, extra, message):
    r = run(tmp_path, extra)
    assert r.returncode == 2, r.stdout[-400:] + r.stderr
    assert "Error: " + message in r.stderr
    assert message not in r.stdout
