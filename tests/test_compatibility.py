import os
import re
import subprocess

import numpy as np
import pytest

"""
These tests just tests that partrac runs on the current system.
"""


REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARTRAC = os.path.join(REPO, "bin", "partrac")


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_run():
    out = subprocess.check_output(PARTRAC, shell=True).decode("utf-8")
    assert "Specify an input file." in out


if __name__ == "__main__":
    test_run()
