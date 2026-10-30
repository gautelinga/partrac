#!/bin/bash
# Cleans the solved OpenFOAM fixtures (their clean_case.sh) and Python's
# bytecode here. -n lists what it would remove.
set -e
cd "$(dirname "$0")"
for c in openfoam_cylinder2d openfoam_sphere3d; do ./$c/clean_case.sh "$@"; done
while IFS= read -r -d '' d; do
  echo "remove ${d#./}/"
  [ "$1" = -n ] || rm -rf -- "$d"
done < <(find . -name __pycache__ -type d -prune -print0)
