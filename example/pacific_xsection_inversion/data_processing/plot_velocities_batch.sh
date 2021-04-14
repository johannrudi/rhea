#!/bin/bash

# get environment
script_path=$0
script_dir=$(dirname $script_path)

# get rhea file from arguments
rhea_out_path='None'
while getopts r: flag
do
  case "${flag}" in
    r) rhea_out_path=${OPTARG};;
  esac
done
echo "Rhea file $rhea_out_path"

# get all arguments
in_args=$@

for in_file in $in_args; do
  if [ -f $in_file ] && [ $in_file != $rhea_out_path ]; then
    # set output file path
    out_file="${in_file%.*}.png"

    # create plot
    echo "Plot $in_file > $out_file"
    python3 $script_dir/plot_velocities.py $in_file $rhea_out_path $out_file
  fi
done
