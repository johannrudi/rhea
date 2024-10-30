#!/bin/bash

function first_word()
{
  echo $1 | awk '{print $1;}'
}

function filter()
{
  local pattern=$1
  local input_string=$2
  local output_string=()

  for s in $input_string; do
    if [[ "$s" == $pattern ]] ; then
      output_string+=("$s")
    fi
  done
  echo "${output_string[@]}"
}

function filter_out()
{
  local pattern=$1
  local input_string=$2
  local output_string=()

  for s in $input_string; do
    if [[ "$s" == $pattern ]] ; then
      continue
    fi
    output_string+=("$s")
  done
  echo "${output_string[@]}"
}

function get_python_ldflags()
{
  echo $(filter -L* "$(python3-config --ldflags --embed)")
}

function get_python_libs()
{
  echo $(filter_out -L* "$(python3-config --ldflags --embed)")
}

function print_sep_line()
{
  echo "========================================"
}

function print_modules()
{
  print_sep_line
  echo "Modules:"
  module -t list 2>&1 | sort
}

function print_setup()
{
  print_sep_line
  echo "Configure setup for Rhea"
  echo "- Site:             $1"
  echo "- Optimized build:  $2"
  echo "- Code directory:   $3"
  echo "- Build directory:  $4"

  echo "- FET dir:          $5"
  echo "- PETSc dir:        $6"
  echo "- Rhea-kit dir:     $7"

  echo "- MPI dir:          $8"
  echo "- OpenMP flag:      $9"
  echo "- Valgrind include: ${10}"
  echo "- Python ldflags:   ${11}"
  echo "- Python libs:      ${12}"
  echo "- BLAS libs:        ${13}"
}
