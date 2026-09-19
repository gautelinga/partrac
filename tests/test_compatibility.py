"""The partrac binary starts on this system and rejects a missing input file."""

import os
import subprocess

import pytest

"""
These tests just tests that partrac runs on the current system.
"""


from paths import REPO, app

PARTRAC = app("partrac")


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_run():
    """partrac without arguments loads, prints its usage hint and exits nonzero.

    This catches a binary that cannot start at all (missing shared libraries,
    wrong architecture), and makes sure scripts see a usage error as a failure.
    """
    r = subprocess.run(PARTRAC, capture_output=True, text=True, timeout=60)
    assert "Specify an input file." in r.stdout
    # a missing input file is a usage error, not a successful run
    assert r.returncode != 0


if __name__ == "__main__":
    test_run()
