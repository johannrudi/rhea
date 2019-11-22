#!/bin/bash

# get environment
script_path=$0
script_dir=$(dirname $script_path)

# get input files from arguments
in_files=$@

for in_file in $in_files; do
  # set output file path
  out_file="${in_file%.*}.png"

  # create plot
  echo "Plot $in_file > $out_file"
  python $script_dir/plot_velocities.py $in_file $out_file
done
