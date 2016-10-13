#!/bin/bash -e
#SBATCH -J __JOB_NAME__
#SBATCH -o __JOB_NAME__.o%j
#SBATCH -e __JOB_NAME__.e%j
#SBATCH -p __QUEUE__
#SBATCH -N __NUM_NODES__
#SBATCH -n __MPISIZE__
#SBATCH -t __RUNTIME__
#SBATCH -A __PROJECT__
#SBATCH --mail-type=ALL
#SBATCH --mail-user=__EMAIL__

###############################################################################
# Template script for running ymir on Stampede.
#
# Replace the following strings:
#   __JOB_NAME__
#   __QUEUE__
#   __NUM_NODES__
#   __MPISIZE__
#   __OMPSIZE__
#   __RUNTIME__
#   __PROJECT__
#   __EMAIL__
#   __HOME_DIR__
#   __SOURCE_PATH__
#   __EXECUTABLE_PATH__
#   __YMIR_OPTIONS__
#
#   __VIS_JOB_PATH__
#   __VIS_RENDER_SCRIPT__
#   __CONCURRENT_VIS__
#   __CONCURRENT_VIS_FILEPATH__
#   __CONCURRENT_VIS_OUT_DIR__
#   __INPUT_VIS__
#   __INPUT_VIS_FILEPATH__
#   __INPUT_VIS_OUT_DIR__
#   __SOLUTION_VIS__
#   __SOLUTION_VIS_FILEPATH__
#   __SOLUTION_VIS_OUT_DIR__
#
#   __PROFILE_HPCTOOLKIT__
#   __PROFILE_TAU_LIB__
#   __PROFILE_TAU_SOURCE__
#   __PROFILE_IPM__
#   __PROFILE_MASSIF__
#
# Author:             Johann Rudi <johann@ices.utexas.edu>
###############################################################################

# set constant parameters
declare -r REL_MEMORY_LIMIT=0.98

########################################
# General settings & information
########################################

# load standard modules
module restore ymir

# echo commands in stderr
set -x

# path output
echo `pwd`
echo $0

echo "========================================"

# output of code version/commit
cd "__SOURCE_PATH__"; echo `git log -3`; cd -;

echo "========================================"

########################################
# Launch concurrent visualization of nonlinear iterations
########################################

# set common parameters for all stages of visualization
vis_job_path="__VIS_JOB_PATH__"
vis_render_script="__VIS_RENDER_SCRIPT__"
declare -i vis_concurrent='__CONCURRENT_VIS__'
declare -i vis_input='__INPUT_VIS__'
declare -i vis_solution='__SOLUTION_VIS__'

# launch concurrent visualization
if [ 1 -le "$vis_concurrent" ]; then
  echo "Launch concurrent visualization"
  module restore default

  vis_filepath="__CONCURRENT_VIS_FILEPATH__"
  vis_out_dir="__CONCURRENT_VIS_OUT_DIR__"

  vis_job_options="-r '$vis_render_script' -f '$vis_filepath' "
  vis_job_options+="-d -j $SLURM_JOB_ID -n __MPISIZE__"
  if [ -d "$vis_out_dir" ]; then
    vis_job_options+=" -o $vis_out_dir"
  fi

  # submit visualization job from login node
  ssh login2 "cd '${vis_job_path%/*}/'; " \
             "sbatch '$vis_job_path' $vis_job_options"

  echo "========================================"
fi

########################################
# Set environment
########################################

# set number of OpenMP threads per MPI rank
export OMP_NUM_THREADS="__OMPSIZE__"

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

# set memory limits to 95% of total memory to prevent node death;
# this overwrites .bashrc in home directory!
#memory_limit=`echo "$REL_MEMORY_LIMIT * $node_memory / $n_tasks_per_node" | bc`
#memory_limit=`echo "$REL_MEMORY_LIMIT * $node_memory / 1" | bc`
#ulimit -v $memory_limit -m $memory_limit
#echo "ulimit -v $memory_limit -m $memory_limit" > "__HOME_DIR__/.bashrc"
#echo "Memory limit set to $memory_limit kilobytes per task"

echo "========================================"

# set threshold for enabling on-demand connection management scheme
# (can help avoiding an MPI malloc bug if set larger than jobsize)
export MV2_ON_DEMAND_THRESHOLD=256
echo "Threshold for on-demand connection management set to" \
     "$MV2_ON_DEMAND_THRESHOLD"

# enable only UD transport; do not use any RC/XRC connections
export MV2_USE_UD_HYBRID=1
export MV2_USE_ONLY_UD=1
echo "Enable native InfiniBand Unreliable Datagram (UD) based asynchronous" \
     "connection management"

# set up profiling
declare -i profile_hpctoolkit='__PROFILE_HPCTOOLKIT__'
declare -i profile_tau_lib='__PROFILE_TAU_LIB__'
declare -i profile_tau_source='__PROFILE_TAU_SOURCE__'
declare -i profile_ipm='__PROFILE_IPM__'
declare -i profile_massif='__PROFILE_MASSIF__'
if [ 1 -le "$profile_hpctoolkit" ]; then
  echo "Enable profiling with HPCToolkit"
  module load hpctoolkit

elif [ 1 -le "$profile_tau_lib" ] || [ 1 -le "$profile_tau_source" ]; then
  echo "Enable profiling with TAU"
  module load tau
  export TAU_PROFILE=1
  export TAU_METRICS=LINUX_TIMERS  # more :PAPI_FP_OPS:PAPI_L2_DCM
  export TAU_VERBOSE=1         # print TAU measuring steps to stderr
  export TAU_TRACK_HEADROOM=0  # (`1` causes crashes on Stampede)
  export TAU_TRACK_HEAP=1
 #export TAU_TRACK_IO_PARAMS=1
  export TAU_TRACK_MEMORY_LEAKS=1
  export TAU_THROTTLE=1                # reduce performance monitoring overhead
  export TAU_THROTTLE_NUMCALLS=1000000 # throttle if executes more than X
  export TAU_THROTTLE_PERCALL=5        # time per function call less than X ms
  export TAU_CALLPATH=1         # enable callpath profiling
  export TAU_CALLPATH_DEPTH=100 # depth of callpath

 #export TAU_TRACE=1
 #export TAU_SYNCHRONIZE_CLOCKS

elif [ 1 -le "$profile_ipm" ]; then
  echo "Enable profiling with IPM"
  module load ipm
  export LD_PRELOAD="$TACC_IPM_LIB/libipm.so"
  export IPM_REPORT='full'      # none, terse (default), full
  export IPM_MPI_THRESHOLD=0.1  # report routines with >10% of MPI total time

elif [ 1 -le "$profile_massif" ]; then
  echo "Enable profiling with Valgrind Massif"
  module load valgrind
fi

# print list of modules
modules=$(module list 2>&1)
modules=${modules##+++ }
echo ${modules%%++*}

echo "========================================"

########################################
# Run executable
########################################

echo "Launch main executable"

exec_path="__EXECUTABLE_PATH__"
ymir_options_path="__YMIR_OPTIONS__"
exec_args="-o $ymir_options_path"

if [ 1 -le "$profile_hpctoolkit" ]; then
  ibrun tacc_affinity hpcrun $exec_path $exec_args
elif [ 1 -le "$profile_tau_lib" ]; then
  ibrun tacc_affinity tauex $exec_path $exec_args
elif [ 1 -le "$profile_massif" ]; then
  ibrun tacc_affinity valgrind --tool=massif \
      --massif-out-file="massif.o${SLURM_JOB_ID}.p%q{MV2_COMM_WORLD_RANK}" \
      $exec_path $exec_args
# ibrun tacc_affinity \
#     valgrind --tool=massif --massif-out-file="massif.o${SLURM_JOB_ID}.p%p" \
#     --time-unit=B \
#     $exec_path $exec_args
else
  ibrun tacc_affinity $exec_path $exec_args
fi

echo "========================================"

########################################
# Launch visualization of input
########################################
if [ 1 -le "$vis_input" ]; then
  echo "Launch input visualization"
  module restore default

  vis_filepath="__INPUT_VIS_FILEPATH__"
  vis_out_dir="__INPUT_VIS_OUT_DIR__"

  vis_job_options="-r '$vis_render_script' -f '$vis_filepath'"
  if [ -d "$vis_out_dir" ]; then
    vis_job_options+=" -o $vis_out_dir"
  fi

  # submit visualization job from login node
  ssh login2 "cd '${vis_job_path%/*}/'; " \
             "sbatch '$vis_job_path' $vis_job_options"

  echo "========================================"
fi

########################################
# Launch visualization of solution
########################################
if [ 1 -le "$vis_solution" ]; then
  echo "Launch solution visualization"
  module restore default

  vis_filepath="__SOLUTION_VIS_FILEPATH__"
  vis_out_dir="__SOLUTION_VIS_OUT_DIR__"

  vis_job_options="-r '$vis_render_script' -f '$vis_filepath'"
  if [ -d "$vis_out_dir" ]; then
    vis_job_options+=" -o $vis_out_dir"
  fi

  # submit visualization job from login node
  ssh login2 "cd '${vis_job_path%/*}/'; " \
             "sbatch '$vis_job_path' $vis_job_options"

  echo "========================================"
fi

########################################
# Postprocessing of profiling files
########################################
if [ 1 -le "$profile_tau_lib" ] || [ 1 -le "$profile_tau_source" ]; then
  echo "Postprocessing of TAU profiling files"
  module load tau

  # pack profiling files
  paraprof --pack tau.ppk

  # delete single profile output
  if [ -f "profile.0.0.0" ]; then
    rm profile.*.*.*
  fi

  # delete multi-profile output directories
  #TODO
fi
