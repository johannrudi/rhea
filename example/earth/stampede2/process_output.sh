#!/bin/bash

BASEDIR="${1:-.}"

# merge and pack text output
merged_file="earth_solution_surf.txt"
packed_file="txt_surf.tgz"
cat "$BASEDIR/txt/earth_solution_surf*.txt" > "$BASEDIR/$merged_file"
tar -czf "$BASEDIR/$packed_file" "$BASEDIR/$merged_file"

# pack text output of Aleutian subduction
packed_file="txt_aleutian"
tar -czf "$BASEDIR/$packed_file" "aleutian"

# pack vtk output
packed_file="vtk_surf.tgz"
tar -czf "$BASEDIR/$packed_file" \
  "$BASEDIR/vtk/earth_solution?face1*" \
  "$BASEDIR/vtk/earth_solution_[0-9]*face1*"
