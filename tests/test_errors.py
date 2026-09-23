"""A user's mistake that the parameter schema cannot see ends the run the way
a parameter error does: exit code 2 and a message on stderr, from the one
place each app reports (partrac::report_errors), not an exit() where it was
found. The schema knows the keys, not what they name: an initial state that
misses the domain, a file that is not there, an expression or an element the
loader does not know, an XDMF grid without a time. And a misspelt key in an
expression file, which that expression's schema refuses by name."""

import os
import re

import pytest

from paths import REPO, app
from runs import copy_case, copy_example, run_app

PARTRAC = app("partrac")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
ARGS = ("mode=analytic dt=0.01 T=0.02 Nrw=10 Nrw_max=100 ds_max=0.1 ds_min=1e-9 "
        "ds_init=0.01 Dm=0 int_order=1 refine=false coarsen=false dump_intv=1e9 "
        "stat_intv=1e9 checkpoint_intv=1e9").split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")


SEED = ["init_mode=strip_x", "La=0.1", "x0=0", "y0=0", "z0=0"]


def run(tmp_path, extra, params=None):
    """partrac on a copy of the Poiseuille example, or on the parameter file `params`."""
    if params is None:
        params = copy_example(POISEUILLE, tmp_path)
    return run_app(PARTRAC, params, ARGS, extra, check=False, timeout=120)


def reported(r, message):
    """Exit code 2 and `message` after "Error: " on stderr, not on stdout."""
    assert r.returncode == 2, r.stdout[-400:] + r.stderr
    assert "Error: " in r.stderr and message in r.stderr, r.stderr
    assert message not in r.stdout


@needs_partrac
@pytest.mark.parametrize("extra,message", [
    (["init_mode=strip_x", "x0=100", "y0=100", "z0=100", "La=0.1"], "strip not inside domain"),
    (["init_mode=ellipsoid_z", "x0=100", "y0=100", "z0=100", "La=0.1", "Lb=0.1"],
     "ellipsoid not inside domain"),
    (["init_mode=pairs_xy", "x0=100", "y0=100", "z0=100"], "pair centre is not inside the domain"),
    (["init_mode=randomgaussiancircle_x", "La=0.1", "Lb=0.01", "x0=100", "y0=100", "z0=100"],
     "no points inside the domain"),
])
def test_an_initial_state_outside_the_domain_is_reported(tmp_path, extra, message):
    reported(run(tmp_path, extra), message)


@needs_partrac
def test_an_init_mode_the_schema_lets_through_is_reported(tmp_path):
    """The schema checks the number of directions, not the mode's name."""
    reported(run(tmp_path, ["init_mode=bogus_x", "La=0.1", "x0=0", "y0=0", "z0=0"]),
             "unknown init_mode: bogus_x")


@needs_partrac
def test_a_missing_parameter_file_is_reported(tmp_path):
    reported(run(tmp_path, SEED, params=tmp_path / "missing.dat"), "no such file")


@needs_partrac
def test_a_missing_positions_file_is_reported(tmp_path):
    reported(run(tmp_path, ["init_mode=from_file:%s" % (tmp_path / "none.h5")]), "no such file")


@needs_partrac
@pytest.mark.parametrize("edit,message", [
    (lambda text: "".join(l for l in text.splitlines(True) if not l.startswith("expression")),
     "no expression= in"),
    (lambda text: "".join("expression=nonsense\n" if l.startswith("expression") else l
                          for l in text.splitlines(True)),
     "unknown expression nonsense"),
])
def test_an_analytic_file_without_a_known_expression_is_reported(tmp_path, edit, message):
    with open(POISEUILLE) as f:
        (tmp_path / "expr_params.dat").write_text(edit(f.read()))
    reported(run(tmp_path, SEED, params=tmp_path / "expr_params.dat"), message)


@needs_partrac
def test_a_misspelt_key_in_the_expression_file_stops_the_run(tmp_path):
    """A key the expression does not read stops partrac before it runs, naming
    the key and the one it resembles.

    Read loosely, u_innf would be ignored and the flow would run with whatever
    u_inf the file also sets, or stop later on a missing key without saying
    which line was wrong.
    """
    with open(POISEUILLE) as f:
        (tmp_path / "expr_params.dat").write_text(f.read().replace("u_inf=", "u_innf="))
    r = run(tmp_path, SEED, params=tmp_path / "expr_params.dat")
    assert r.returncode != 0
    out = r.stdout + r.stderr
    assert "unknown parameter 'u_innf'" in out, out
    assert "did you mean 'u_inf'" in out, out
    assert "missing required parameter 'u_inf'" in out, out
    assert not list(tmp_path.rglob("tdata_from_t*.dat"))


def mesh_case(src, tmp_path, edit):
    """A copy of the mesh case in src with dolfin_params.dat passed through edit."""
    f = copy_case(src, tmp_path / "case") / "dolfin_params.dat"
    f.write_text(edit(f.read_text()))
    return f


@needs_partrac
def test_an_element_the_loader_does_not_know_is_reported(mesh_dir, tmp_path):
    f = mesh_case(mesh_dir("tet"), tmp_path,
                  lambda t: "".join("velocity_space=P7\n" if l.startswith("velocity_space") else l
                                    for l in t.splitlines(True)))
    reported(run(tmp_path, ["mode=tet", "init_mode=points_xyz", "init_weight=uniform",
                            "x0=0.5", "y0=0.5", "z0=0.5"], params=f),
             "unrecognized velocity element: P7")


@needs_partrac
def test_an_xdmf_grid_without_a_time_is_reported(xdmf_dir, tmp_path):
    d = copy_case(xdmf_dir, tmp_path / "case")
    u = d / "u.xdmf"
    u.write_text(re.sub(r"<Time [^>]*/>", "", u.read_text(), count=1))
    reported(run(tmp_path, ["mode=xdmftriangle", "init_mode=points_xy", "init_weight=uniform",
                            "x0=0.5", "y0=0.5", "z0=0"], params=d / "dolfin_params.dat"),
             "XDMF: a grid without a time")
