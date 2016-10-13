#!/bin/bash -e

###############################################################################
# Batch processing of files for workload plotting.
#
# Author:             Johann Rudi <johann@ices.utexas.edu>
###############################################################################

# set constant parameters
EXTENSION='txt'
EXTENSION_PETSC='petsc.txt'

########################################
# Process arguments
########################################

function print_usage()
{
  echo "Usage: `basename $0` -r PLOT_SCRIPT -f FILEPATH" \
       "[-p] [-d OUT_DPI] [-o OUT_DIR]"
}

# initialize parameters
plot_script=''
filepath=''
plot_petsc_file=0
out_dpi=0
out_directory=''

# set parameters from arguments
while getopts ":hr:f:pd:o:" opt; do
  case $opt in
    h) # help
      print_usage; exit 1
      ;;
    r) # path to Python plotting script
      plot_script=$OPTARG
      ;;
    f) # path to files that should be rendered
      filepath=$OPTARG
      ;;
    p) # plot memory usage recorded by petsc
      plot_petsc_file=1
      ;;
    d) # DPI for output pictures
      out_dpi=$OPTARG
      ;;
    o) # path to output directory
      out_directory=$OPTARG
      ;;
    \?) # invalid option
      echo "Invalid option: -$OPTARG" >&2; exit 1
      ;;
    :) # option requires argument
      echo "Option -$OPTARG requires an argument" >&2; exit 1
      ;;
  esac
done

# extract input directory and file name pattern
if [ -d "$filepath" ]; then  # if path is directory
  echo "Error: '$filepath' cannot be a directory"
  exit 1
elif echo "$filepath" | egrep -iq "\/"; then  # if path (to file) has dir's
  in_directory="${filepath%/*}/"
  filename_pattern=$(basename "$filepath")
else  # if path is just a filename
  in_directory="$(pwd)/"
  filename_pattern="$filepath"
fi

# check parameters generated from arguments
if [ "$filepath" = '' ] || [ "$plot_script" = '' ]; then
  print_usage
  exit 1
fi

if [ ! -d "$in_directory" ]; then
  echo "Error: input directory '$in_directory' does not exist" >&2
  exit 1
fi

if [ ! -f "$plot_script" ]; then
  echo "Error: rendering script '$plot_script' does not exist" >&2
  exit 1
fi

if [[ "$out_directory" != '' ]]; then
  if [ ! -d "$out_directory" ]; then
    echo "Error: output directory '$out_directory' does not exist" >&2
    exit 1
  else
    out_directory="${out_directory%/}/"
  fi
fi

########################################
# Start rendering
########################################

# set extension for file search
if echo "$filename_pattern" | egrep -iq "\.${EXTENSION}$"; then
  extension_pattern=''
else
  extension_pattern="*.${EXTENSION}"
fi

# status output
echo "Start plotting of data in directory '$in_directory' with pattern" \
     "'${filename_pattern}${extension_pattern}'"

# plot all files matching the name pattern
for file in $(find "$in_directory" -maxdepth 1 -type f \
              -name "${filename_pattern}${extension_pattern}");
do
  # skip if file matches petsc file
  if echo "$file" | egrep -iq "\.${EXTENSION_PETSC}$"; then
    continue
  fi

  # status output
  echo "========================================"
  echo "Plot dataset '$file'"

  # extract basename
  file_basename=$(basename "$file")
  file_basename=${file_basename%.*}

  # set execution path and arguments
  exec_path="$plot_script"
  exec_args="$file"
  if [ 1 -le "$plot_petsc_file" ] && \
     [ -f "${in_directory}${file_basename}.${EXTENSION_PETSC}" ]; then
    # set petsc file
    exec_args+=" -p ${in_directory}${file_basename}.${EXTENSION_PETSC}"
  fi
  if [ 0 -lt "$out_dpi" ]; then
    # set output DPI
    exec_args+=" -d ${out_dpi}"
  fi
  if [ "$out_directory" != '' ]; then
    # set output file path
    exec_args+=" -o ${out_directory}${file_basename}"
  fi

  # create plot
  eval "$exec_path $exec_args"
done

# status output
echo "========================================"
echo "Done plotting of data in directory '$in_directory'"

