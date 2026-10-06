"""frozen_fields on the analytic input: the expression at t_frozen at every time.

The sine flow switches direction every tau = 0.5, so a field that followed the
run's time would turn within the runs below, and one frozen at another time
than t_frozen would point elsewhere. The loaders that read files have their
frozen tests beside their other tests; the default and the refusal, shared by
every loader, and the clamp of t_frozen to the flow's times, are here.
"""

import os

import numpy as np
import pytest

from dumps import dump_at
from paths import REPO, app
from runs import checkpoint_folder, continuous_and_resumed, copy_example, run_app

TRACERS = app("tracers")
SINE = os.path.join(REPO, "data_example", "sine_flow", "expr_params.dat")
ARGS = ("mode=analytic init_mode=points_xy x0=0.5 y0=0.5 z0=0.5 Nrw=50 Nrw_max=50 Dm=0 int_order=1 "
        "scheme=RK4 dt=0.05 stat_intv=1e9 checkpoint_intv=1e9 random=false seed=1 frozen_fields=true")


@pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")
def test_a_frozen_expression_does_not_depend_on_the_time(tmp_path):
    """Frozen at t = 0.7, in another half period than either start, a run
    from t0 = 0 and the same run from t0 = 1.1, each across a switch of the
    sine flow, end at the same points bit for bit."""
    ends = []
    for t0 in (0.0, 1.1):
        params = copy_example(SINE, tmp_path / ("t0_%g" % t0))
        run_app(TRACERS, params, ARGS, "t_frozen=0.7", "t0=%r T=%r dump_intv=0.4" % (t0, t0 + 0.8))
        start, end = dump_at(params.parent, t0)["points"], dump_at(params.parent, t0 + 0.8)["points"]
        assert np.abs(end - start).max() > 0.1
        ends.append(end)
    assert np.array_equal(ends[0], ends[1])


@pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")
def test_t_frozen_is_clamped_to_the_flows_times(tmp_path):
    """A t_frozen before the flow's t_min = 0 freezes it at 0, and one past
    its t_max = 1e7 at 1e7, the run going on from t0 = 0 either way; the time
    frozen at is the one recorded."""
    ends = {}
    for name, t_frozen in (("below", -1), ("start", 0), ("above", 1e30), ("end", 1e7)):
        params = copy_example(SINE, tmp_path / name)
        r = run_app(TRACERS, params, ARGS, "t0=0 T=0.4 dump_intv=0.4 t_frozen=%r" % t_frozen)
        ends[name] = dump_at(params.parent, 0.4)["points"]
        assert "Fields frozen at t = %g" % max(0, min(t_frozen, 1e7)) in r.stdout
    assert np.array_equal(ends["below"], ends["start"])
    assert np.array_equal(ends["above"], ends["end"])
    assert not np.array_equal(ends["start"], ends["end"])


@pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")
def test_frozen_fields_default_to_the_start_time(tmp_path):
    """Without t_frozen the fields are frozen at t0: a run from t0 = 0.7 ends
    where the run from t0 = 0 frozen at 0.7 does, after as many steps."""
    ends = []
    for name, args in (("default", "t0=0.7 T=1.5"), ("given", "t0=0 T=0.8 t_frozen=0.7")):
        params = copy_example(SINE, tmp_path / name)
        run_app(TRACERS, params, ARGS, args, "dump_intv=0.8")
        ends.append(dump_at(params.parent, 1.5 if name == "default" else 0.8)["points"])
    assert np.array_equal(ends[0], ends[1])


@pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")
def test_a_resumed_run_stays_frozen_at_the_first_start(tmp_path):
    """The default t_frozen is recorded as the start time, so a run resumed
    later in time keeps the field of the first start, bit for bit the run
    never stopped."""
    cont, split = continuous_and_resumed(TRACERS, SINE, tmp_path, [ARGS, "t0=0.7 dump_intv=0.4"],
                                         ["T=1.1", "checkpoint_intv=0.4"], ["T=2.3", "checkpoint_intv=1e9"])
    a, b = dump_at(cont, 2.3)["points"], dump_at(split, 2.3)["points"]
    assert np.array_equal(a, b)
