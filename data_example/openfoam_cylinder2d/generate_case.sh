#!/bin/bash
# Remakes this fixture with OpenFOAM (Foundation, 14): the box [0,1]^2,
# cyclic in x and y, one cell thick (front and back empty), less one
# circular cylinder (no-slip), blockMesh's O-grid of 2048 hexes;
# ../stokes_case.py's solve, foamRun's steady Stokes-like flow (nu = 1)
# under a body force of 1 at 20 degrees to x, converged to residuals 1e-9;
# the mesh and the last time's U kept, converted to binary and gzipped.
# cylinder_case.py writes the dictionaries and prints the at-rest
# quantities. Seconds.
set -e
cd "$(dirname "$0")"
if [ -z "$WM_PROJECT_DIR" ]; then source /opt/openfoam14/etc/bashrc > /dev/null 2>&1 || true; fi
rm -rf constant system 0 [1-9]* postProcessing
python3 cylinder_case.py dicts
blockMesh > /dev/null
foamRun > /dev/null
T=$(foamListTimes -latestTime | tail -1)
foamPostProcess -func writeCellVolumes -time "$T" > /dev/null
python3 cylinder_case.py measure "$T"
find "$T" -mindepth 1 ! -name U -exec rm -rf {} +
rm -rf 0 postProcessing constant/fvModels constant/physicalProperties constant/momentumTransport \
    system/*MeshDict
sed -i 's/^writeFormat .*/writeFormat     binary;/' system/controlDict
foamFormatConvert > /dev/null
gzip -9n constant/polyMesh/* "$T"/*
