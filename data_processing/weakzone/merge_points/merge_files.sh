#!/bin/bash

in_files_list=$1
out_file=$2

xargs < "$in_files_list" cat > "$out_file"

# Note: `xargs` sends the file names from the list to `cat` in batches (so they
# won't overflow the command line).
