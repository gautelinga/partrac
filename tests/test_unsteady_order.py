"""The second-order explicit step on an unsteady field converges at second order.

scheme=explicit int_order=2 steps by x += u dt + (du/dt + (grad u) u) dt^2/2,
so a field that changes in time must give du/dt (PointValues::A) as well as
u; without it the step is first order wherever the field is unsteady. The
field is the ABC flow with A oscillating, A_amp and A_omega made large so the
unsteady part dominates; the reference is RK4 at a step eight times finer
than the finest explicit one, whose error is far below theirs.
"""

import os
import re

import numpy as np
import pytest

from dumps import dump_at
from paths import REPO, app, needs
from runs import copy_example, run_app

TRACERS = app("tracers")
ABC = os.path.join(REPO, "data_example", "abc_flow_unsteady", "expr_params.dat")
T = 1.0
ARGS = ("mode=analytic init_mode=points_xyz init_weight=uniform Nrw=20 Nrw_max=20 Dm=0 "
        "T=%g dump_intv=%g stat_intv=1e9 checkpoint_intv=1e9 random=false seed=1" % (T, T)).split()


def final_positions(tmp_path, name, *args):
    """Positions at T, in id order, of a run on the strongly forced ABC flow."""
    params = copy_example(ABC, tmp_path / name)
    text = re.sub(r"^A_amp=.*$", "A_amp=1.0", params.read_text(), flags=re.M)
    params.write_text(re.sub(r"^A_omega=.*$", "A_omega=3.0", text, flags=re.M))
    run_app(TRACERS, params, ARGS, *args)
    return dump_at(tmp_path / name, T)["points"]


@needs(TRACERS)
def test_the_second_order_explicit_step_converges_at_second_order_on_an_unsteady_flow(tmp_path):
    """Halving dt divides the error of the positions at T by about 4."""
    ref = final_positions(tmp_path, "ref", "scheme=RK4 int_order=1 dt=0.00125")
    err = []
    for dt in (0.04, 0.02, 0.01):
        x = final_positions(tmp_path, "dt%g" % dt, "scheme=explicit int_order=2 dt=%g" % dt)
        err.append(np.abs(x - ref).max())
    ratios = [err[0]/err[1], err[1]/err[2]]
    assert all(3.4 < r < 4.6 for r in ratios), (err, ratios)
