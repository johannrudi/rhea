#!/bin/bash

in_files=$@
out_dir="data_cartesian"

mkdir -p "$out_dir"

for in_file in $in_files; do
  out_file="$out_dir/${in_file##*/}"
  echo "Transform: $in_file > $out_file"
  cat $in_file | perl transform_spherical.pl > "$out_file"
done
