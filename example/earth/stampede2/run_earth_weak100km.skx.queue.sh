#!/bin/bash

#SBATCH -J earth_weak100km
#SBATCH -o earth_weak100km.o%j
#SBATCH -e earth_weak100km.e%j
#SBATCH -p skx-normal
#SBATCH -N 64             # Total number of nodes (now required)
#SBATCH -n 3072           # Total number of mpi tasks
#SBATCH -t 10:00:00
#SBATCH -A TG-DPP130002
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johann@ices.utexas.edu

###############################################################################
# Author:             Johann Rudi <johann@ices.utexas.edu>
###############################################################################

########################################
# Constants
########################################

declare -r BUILD_DIR="$HOME/build/perf/rhea"
declare -r EXEC_RELPATH="example/earth/rhea_earth"
declare -r CODE_DIR="$HOME/code/rhea"

declare -r JOB_DIR="$SCRATCH/runs/rhea/earth_weak100km_YYYY-MM-DD"

declare -r BIN_DIR="$JOB_DIR/bin"
declare -r TXT_DIR="$JOB_DIR/txt"
declare -r VTK_DIR="$JOB_DIR/vtk"
declare -r ALEUTIAN_DIR="$JOB_DIR/aleutian"

declare -r OMPSIZE=2

########################################
# Functions
########################################

function print_path_environment()
{
  echo "PWD:    `pwd`"
  echo "Script: $0"
}

function print_job_environment()
{
  # job status output
  echo "Job id:                   $SLURM_JOB_ID"
  echo "Job name:                 $SLURM_JOB_NAME"
  echo "Number of nodes:          $SLURM_NNODES"
  echo "Number of MPI tasks:      $SLURM_NTASKS"
  echo "Number of tasks per node: $SLURM_TASKS_PER_NODE"
  echo "Number of OpenMP threads: $OMP_NUM_THREADS"
  echo "Queue:                    $SLURM_QUEUE"
  echo "Allocation:               $SLURM_TACC_ACCOUNT"
  echo "Directory of submission:  $SLURM_SUBMIT_DIR"
  echo "========================================"
  echo "List of nodes: $SLURM_NODELIST"
  echo "========================================"

  # get number of tasks per node
  length=`expr match "$SLURM_TASKS_PER_NODE" '^[0-9]*'`
  n_tasks_per_node=${SLURM_TASKS_PER_NODE:0:$length}

  # get total memory per node in kilobytes
  node_memory=`free -k | grep ^Mem: | awk '{ print $2; }'`
  echo "Total memory per node: $node_memory kilobytes"
}

function print_program_environment()
{
  # output of code version/commit
  cd $CODE_DIR; echo "$(git log -3)"; cd - &> /dev/null;
  echo "========================================"

  # print list of modules
  modules=$(module list 2>&1)
  modules=${modules##+++ }
  echo "${modules%%++*}"
}

########################################
# Set Environment
########################################

# echo commands in stderr
set -x

# load modules
#module restore default

# set number of OpenMP threads per MPI rank
export OMP_NUM_THREADS=$OMPSIZE

########################################
# Print Environment
########################################

echo "========================================"
print_path_environment
echo "========================================"
print_job_environment
echo "========================================"
print_program_environment
echo "========================================"

########################################
# Run executable
########################################

# create output paths
mkdir -p "$BIN_DIR"
mkdir -p "$TXT_DIR"
mkdir -p "$VTK_DIR"
mkdir -p "$ALEUTIAN_DIR"

# set execution options
options_path="$JOB_DIR/input/options.ini"
exec_args="-f $options_path"
exec_path="$BUILD_DIR/$EXEC_RELPATH"

echo "Launch main executable"

# launch executable
out_base="${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
out_log="$JOB_DIR/${out_base}.out"
err_log="$JOB_DIR/${out_base}.err"
#ibrun tacc_affinity $exec_path $exec_args 1>> $out_log 2>$err_log
ibrun $exec_path $exec_args 1>> $out_log 2>$err_log

echo "========================================"
