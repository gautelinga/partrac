#!/bin/bash
# Remakes this fixture with OpenFOAM (Foundation, 14): blockMesh's 12^3
# channel, its points jittered and the linear field written by
# channel_fields.py, then cellpoint_U.txt, OpenFOAM's own cellPoint U at the
# 200 points interpol probes; the mesh and the fields converted to binary and
# gzipped, which OpenFOAM reads as they are; split_expected.h5, the reference
# split, by tests/openfoam_expected.py. Seconds.
set -e
cd "$(dirname "$0")"
if [ -z "$WM_PROJECT_DIR" ]; then source /opt/openfoam14/etc/bashrc > /dev/null 2>&1 || true; fi
rm -rf constant 0 postProcessing cellpoint_U.txt
sed -i 's/^writeFormat .*/writeFormat     ascii;/' system/controlDict
blockMesh > /dev/null
python3 channel_fields.py jitter
mkdir 0
foamPostProcess -func writeCellCentres -time 0 > /dev/null
python3 channel_fields.py fields
rm -f 0/C 0/Ccx 0/Ccy 0/Ccz
python3 probe_points.py write 200
foamPostProcess -func probesCP -time 0 > /dev/null
python3 probe_points.py collect
rm -rf postProcessing system/probesCP
sed -i 's/^writeFormat .*/writeFormat     binary;/' system/controlDict
foamFormatConvert > /dev/null
gzip -9n constant/polyMesh/* 0/*
python3 ../../tests/openfoam_expected.py .
