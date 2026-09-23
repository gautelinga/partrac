"""The march in path length: static_space_stepper and tracervectors_spatial.

Both apps advance particles by a fixed path-length step dxn rather than a time
step. The flow is plane Poiseuille, u = (0, 0, 1.5 (1 - x^2)): particles placed
on the x axis march along z and never move in x, so each particle's speed is
the same at every step and which of them are slower than u_eps is known from
the first dump.
"""

import os

import numpy as np
import pytest

from dumps import all_dumps, by_id
from paths import REPO, app
from runs import continuous_and_resumed, copy_example, run_app

SPATIAL = app("tracervectors_spatial")
STEPPER = app("static_space_stepper")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

BASE = {
    SPATIAL: ("mode=analytic init_mode=points_x x0=0 y0=0 z0=0 Nrw=400 Nrw_max=400 "
              "stat_intv=1e9 checkpoint_intv=1e9 random=false seed=1"),
    # the stepper floors its intervals to multiples of dt, so dt is set to a
    # value that divides every interval used here
    STEPPER: ("mode=analytic init_mode=points_x x0=0 y0=0 z0=0 Nrw=400 Nrw_max=400 "
              "ds_init=0 init_weight=uniform int_order=1 dx_max=1e9 T=1e9 dt=0.005 "
              "stat_intv=1e9 checkpoint_intv=1e9 random=false seed=1"),
}


def run(binary, d, extra):
    """Run binary on plane Poiseuille in d, with BASE overridden by extra; returns d."""
    run_app(binary, copy_example(POISEUILLE, d), BASE[binary], extra, timeout=600)
    return d


def groups(d):
    """Path length -> the dump's datasets in id order, with the sorted ids."""
    return {xn: dict(by_id(g), id=np.sort(g["id"][:, 0]))
            for xn, g in all_dumps(d, raw=True).items()}


@pytest.mark.skipif(not os.path.exists(SPATIAL), reason="tracervectors_spatial is not built")
@pytest.mark.parametrize("outside", ["remove", "mark", "ignore"])
def test_a_particle_slower_than_u_eps_follows_outside(tmp_path, outside):
    """A particle with |u| <= u_eps cannot cover a path-length step in any
    finite time, so the march treats it as unable to move and outside= decides
    its fate: removed, kept in place and marked c = 2, or kept in place unmarked.
    Every other particle must still take its step along z."""
    # u_eps = 0.1 stops the particles near the walls, where 1.5 (1 - x^2) <= 0.1;
    # dxn = Ln = dump_intv gives exactly one step, dumped before and after
    d = run(SPATIAL, tmp_path / outside,
            "u_eps=0.1 dxn=0.005 Ln=0.005 dump_intv=0.005 outside=" + outside)
    g = groups(d)
    first, second = sorted(g)[:2]
    a, b = g[first], g[second]
    slow = np.linalg.norm(a["u"], axis=1) <= 0.1
    assert slow.any() and not slow.all()
    if outside == "remove":
        assert np.array_equal(b["id"], a["id"][~slow])
        return
    assert np.array_equal(b["id"], a["id"])
    assert np.array_equal(b["points"][slow], a["points"][slow])
    assert np.all(b["points"][~slow][:, 2] > a["points"][~slow][:, 2])
    c = b["c"][:, 0]
    if outside == "mark":
        assert np.all(c[slow] == 2.0) and not np.any(c[~slow] == 2.0)
    else:
        assert not np.any(c == 2.0)


@pytest.mark.parametrize("binary", [SPATIAL, STEPPER],
                         ids=["tracervectors_spatial", "static_space_stepper"])
def test_a_resumed_march_is_identical_to_one_never_stopped(tmp_path, binary):
    """A march stopped at a checkpoint and resumed must give bit-identical dumps
    to one run straight through. The checkpoint carries each particle's
    integration time and the step count resumes; otherwise a long march split
    over several jobs would restart its integration times from zero or dump at
    the wrong path lengths."""
    if not os.path.exists(binary):
        pytest.skip(os.path.basename(binary) + " is not built")
    # the final checkpoint is written one step past Ln: xn = 0.2, step 20
    cont, split = continuous_and_resumed(binary, POISEUILLE, tmp_path,
                                         [BASE[binary], "dxn=0.01 dump_intv=0.1"],
                                         "Ln=0.19", "Ln=0.4")
    a, b = groups(cont), groups(split)
    for xn in (0.3, 0.4):
        assert xn in a, "the march never reached xn = %g: the loop end dropped its last steps" % xn
        assert xn in b, ("the resumed march has no dump at xn = %g (it dumped at %s): the step "
                         "count or the path length did not resume" % (xn, sorted(b)))
        assert set(a[xn]) == set(b[xn])
        for k in a[xn]:
            assert np.array_equal(a[xn][k], b[xn][k]), "%s differs at xn = %g after a restart" % (k, xn)
