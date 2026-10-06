#!/bin/bash
# Remakes this fixture with OpenFOAM (Foundation, 14): the lid-driven cavity
# of the tutorials (20x20 cells, one thick, kEpsilon), blockMesh and foamRun
# from 0/ to t = 10, where U and p are kept; then cellpoint_U.txt, OpenFOAM's
# own cellPoint U at the 200 points interpol probes; split_expected.h5, the
# reference split, by tests/openfoam_expected.py. Seconds.
set -e
cd "$(dirname "$0")"
if [ -z "$WM_PROJECT_DIR" ]; then source /opt/openfoam14/etc/bashrc > /dev/null 2>&1 || true; fi
rm -rf constant/polyMesh 10 postProcessing cellpoint_U.txt
blockMesh > /dev/null
foamRun > /dev/null
find 10 -mindepth 1 ! -name U ! -name p -exec rm -rf {} +
python3 probe_points.py write 200 2d
foamPostProcess -func probesCP -time 10 > /dev/null
python3 probe_points.py collect
rm -rf postProcessing system/probesCP
python3 ../../tests/openfoam_expected.py .
