import os
import re
import subprocess

import numpy as np
import pytest

"""
These tests just tests that partrac runs on the current system.
"""


def test_run():
    cmd = "partrac"
    out = subprocess.check_output(cmd, shell=True)
    ref = b'Specify an input timestamps file.\n'
    assert out == ref


if __name__ == "__main__":
    test_run()
