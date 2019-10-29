#!/bin/bash

in_files_list=$1
out_file=$2

# purge out file
> $out_file

# loop over all files from input list
while read -r -a strarr; do
  # get label
  label="${strarr[0]}"

  # get number of lines (points) from coordinates file
  count=$(cat "${strarr[1]}" | wc -l | tr -d '[:space:]')

  # output status
  echo "process label: $label, count: $count"

  # write the current label for `count` times into the output file
  printf "${label}\n%.0s" $(seq $count) >> $out_file
done < "$in_files_list"
