#!/bin/bash
# Removes from this fixture everything the tests do not read: run output,
# anything OpenFOAM or Python left, the solver's dictionaries.
# Keeps the mesh, time 548's U, controlDict, the parameter file and the
# generator. -n lists what it would remove.
set -e
cd "$(dirname "$0")"
keep=" generate_case.sh clean_case.sh cylinder_case.py partrac_params.dat system/controlDict 
       constant/polyMesh/boundary.gz constant/polyMesh/faces.gz constant/polyMesh/neighbour.gz
       constant/polyMesh/owner.gz constant/polyMesh/points.gz 548/U.gz "
while IFS= read -r -d '' f; do
  f=${f#./}
  case "$keep" in *" $f "*|*" $f"$'\n'*) continue ;; esac
  echo "remove $f"
  [ "$1" = -n ] || rm -f -- "$f"
done < <(find . -type f -print0)
[ "$1" = -n ] || find . -mindepth 1 -type d -empty -delete
