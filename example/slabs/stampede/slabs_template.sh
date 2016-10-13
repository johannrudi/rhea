#!/bin/bash -e

#######################################
# Author:             Johann Rudi <johann@ices.utexas.edu>
#######################################

PROBLEM_TYPE='nl'       # lin, nl
PROBLEM_DOMAIN='earth'  # earth, slice, brick

if [ "$PROBLEM_TYPE" = 'lin' ]; then
  script_name="Convergence of linear solver for "
elif [ "$PROBLEM_TYPE" = 'nl' ]; then
  script_name="Convergence of nonlinear solver for "
else
  echo "Unknown problem type \"${PROBLEM_TYPE}\""
  exit 1
fi

if [ "$PROBLEM_DOMAIN" = 'earth' ]; then
  script_name+="global mantle flow"
elif [ "$PROBLEM_DOMAIN" = 'slice' ]; then
  script_name+="2-plates problem on slice"
elif [ "$PROBLEM_DOMAIN" = 'brick' ]; then
  script_name+="2-plates problem on brick"
else
  echo "Unknown problem domain \"${PROBLEM_DOMAIN}\""
  exit 1
fi

#######################################
# Set script environment
#######################################

script_base_name=$(basename "$0")
script_base_name="${script_base_name%.*}"
script_path="${HOME}/jobs/ymir/${PROBLEM_TYPE}_conv/"
called_from_path="$(pwd)/"
logfile="${called_from_path}${script_base_name}.log"
script_launch_all_filename="${script_base_name}_launch_all_jobs.sh"

#######################################
# Set system environment
#######################################

SYSTEM='stampede'  # stampede, bgq

# set available memory per process (in Bytes)
if [ "$SYSTEM" = 'stampede' ]; then # if Stampede
  arch_mem_max='2000000000.0'  # 2 GB
elif [ "$SYSTEM" = 'bgq' ]; then # if BG/Q
  arch_mem_max='900000000.0'  # 900 MB
else
  echo "Unknown system \"${SYSTEM}\""
  exit 1
fi

#######################################
# Set varying parameters for convergence tests
#######################################

# physics parameters
viscosity_decay=(17.5)
init_amr_n_elements_max=('0')
init_amr_visc_tol_max=('1e-4')
init_amr_post_uniform_n_steps=('0')
amr_n_elements_max=('0')
amr_visc_tol_max=('1e-4')
amr_visc_dr_tol_max=('1e+6')
amr_strain_rate_tol_max=('1e-4')

weak_factor='1e-5'                   # fine: 1e-5, coarse: 1e-4, scoarse: 1e-4
viscosity_for_init_nl_stokes='4'

viscosity_model=('UWYL_SHIFT_LREG')
viscosity_min=('1e-2')               # fine: 1e-2, coarse: 1e-2, scoarse: 1e-1
viscosity_max=('1e+4')               # fine: 1e+4, coarse: 1e+4, scoarse: 1e+3
viscosity_temp_max=('-1')            # fine: 1e+7, coarse: 1e+7, scoarse: 1e+7

# discretization parameters
order=(2)
pressure_space=(0)  # 0: P_{k-1}, 1: Q_{k-2}, 2: stabilized (k=1 only)
maxlevel=(15)

# nonlinear solver parameters
if [ "$PROBLEM_TYPE" = 'lin' ]; then
  nonlinear_type='-1'
  nonlinear_forcing_max='1.0e-6'
else
  nonlinear_type='1'  # 0: Picard, 1: Newton, 2: Picard-Newton
  nonlinear_forcing_max='1.0e-4'
fi
nonlinear_primaldual=(0)
nonlinear_primaldual_scal=(0)
nonlinear_step_length_min='-1'
nonlinear_norm_type='4'

# linear solver parameters
stress_gmg='1'
stress_gmg_smoothing_iter='3'
stress_gmg_cycle_type='0'  # 0: V-cycle, 1: W-cycle
stress_gmg_num_cycles='1'
stress_amg_drop_tol='0.001'
stress_amg_agg_nsmooths='1'
stress_amg_smoother='pbJacobi'
stress_amg_smoothing_iter='3'
stress_amg_vcycles='1'

schur_type=(1)
bfbt_type=(3)
nl_scaling=(0)

bbt_stiff_pc='1'
bbt_gmg='0'
bbt_gmg_smoothing_iter='3'
bbt_gmg_cycle_type='0'  # 0: V-cycle, 1: W-cycle
bbt_gmg_num_cycles='1'
bbt_amg_drop_tol='0.001'
bbt_amg_agg_nsmooths='1'
bbt_amg_smoother='Jacobi'
bbt_amg_smoothing_iter='3'
bbt_amg_vcycles='1'

solve_stress_block_only=(0)
solve_press_bbt_only=(0)  # Required: schur_type = 1
solve_cnode_bbt_only=(0)  # Required: schur_type = 1

#######################################
# Set constant parameters during tests
#######################################

# physics parameters
if [ "$PROBLEM_DOMAIN" = 'earth' ]; then
  domain_shape='shell'
elif [ "$PROBLEM_DOMAIN" = 'slice' ]; then
  domain_shape='shell_slice'
elif [ "$PROBLEM_DOMAIN" = 'brick' ]; then
  domain_shape='brick'
fi
boundary_condition='1'
boundary_default_dir_scale='0'
if [ "$PROBLEM_DOMAIN" = 'earth' ]; then # if earth domain
  read_temp_weak_bin_only='1'
  path_temp_bin="${HOME}/research/mantle_data/rhea1/shell_temp_k2_ll456.bin"
  path_temp_txt="${HOME}/research/mantle_data/rhea1/shell_temp_k2_ll456.txt"
  path_weak_bin="${HOME}/research/mantle_data/rhea1/shell_weakzone0324.bin"
  path_weak_txt="${HOME}/research/mantle_data/rhea1/shell_weakzone0324.txt"
  temp_import_background_plate_age='50.0e6'  # 50 My
  temp_import_plate_age_min='1.0e6'          #  1 My
  weakzone_size='2541799'
  weak_thickness='10.0e3'            # fine: 10e3, coarser: 20e3, 50e3, 100e3
  weak_thickness_const='4.0e3'       # fine:  4e3, coarser:  5e3, 10e3,  15e3
else # if slice or brick domain
  weak_thickness='10.0e3'            # fine: 10e3, coarse: 50e3
  weak_thickness_const='4.0e3'       # fine:  4e5, coarse: 15e3
  weak_ridge_smoothwidth='10.0e3'
fi # end if domain
if [ "$PROBLEM_TYPE" = 'lin' ]; then
  viscosity='1'
  viscosity_scaling='4.0e3'
  lower_mantle='1'
  lower_mantle_scaling='2.0e6'
  stress_exp='1'
  rhs_scaling='8.0e7'
  stress_yield='-1'
  stress_yielding_reg='0'
else
  if [ "$PROBLEM_DOMAIN" = 'earth' ]; then # if earth domain
    viscosity='2'
    viscosity_scaling='8.0e6'
    lower_mantle='1'
    lower_mantle_scaling='5.0e5'
    stress_exp='3'
    rhs_scaling='2.34e9'
    stress_yield='8.93e7'
    stress_yielding_reg='0.01'
  else # if slice or brick domain
    viscosity='2'
    viscosity_scaling='3.0e7'
    lower_mantle='1'
    lower_mantle_scaling='2.0e6'
    stress_exp='3'
    rhs_scaling='2.34e9'
    stress_yield='8.12e7'
    stress_yielding_reg='0.01'
  fi # end if domain
fi

# discretization parameters
if [ "$PROBLEM_DOMAIN" = 'earth' ]; then # if earth domain
  minlevel='4'
  refine_radius='0.76, 0.96'
  import_mesh_order='2'
  import_mesh_minlevel='4'
  init_amr_weak_indicator='15'  # 0 or 15
  init_amr_weak_elem_res="$weak_thickness"
else # if slice or brick domain
  minlevel='1'
  refine_radius='0.76, 0.96'
  import_mesh_order='0'
  import_mesh_minlevel='0'
  init_amr_weak_indicator='13'  # 0 or 13
  init_amr_weak_elem_res="$weak_thickness"
fi # end if domain
refine_surface_maxdist='0.0'
refine_surface_elem_res='16.0e3'
refine_lm_um_maxdist='32.0e3'
refine_lm_um_elem_res='16.0e3'
init_amr_override_order='2'
init_amr_visc_indicator='5'
init_amr_visc_tol_min='0'
init_amr_rhs_indicator='0'  # 0, 16, or 17
init_amr_rhs_tol_min='0'
init_amr_rhs_tol_max='1e-6'
init_amr_rhs_norm_shift='1e3'
amr_rel_threshold='0.05'
amr_visc_indicator='5'
amr_visc_tol_min='0'
amr_visc_dr_indicator='3'
amr_visc_dr_tol_min='0'
amr_strain_rate_indicator='10'
amr_strain_rate_tol_min='0'
mesh_partitioning_type='0'

# nonlinear solver parameters
nl_maxiter='40'
nl_rtol='1.0e-7'
nl_forcing_exponent='1.618'  # superlinear: 1.618, quadratic: 2
nl_forcing_max_progressive_iter='20'
nl_forcing_total_min='1.0e-10'
nl_forcing_saveguard='0'
nl_step_reduction_type='0'
nl_switch_picard_step_length_min='0.01'
nl_switch_picard_after_amr='0'
nl_switch_picard_init='1'
nl_switch_picard_maxiter='3'
nl_switch_picard_rtol='1.0e-6'
nl_norm_Hminus1_mass_scaling='0'
nl_initial_guess_type='0'
nl_schur_diag_type='0'
nl_project_out_nullspace='0'  # 0: off, 1: proj on, 2: proj res, 3: proj symm
nl_enforce_unscaled_reduction='0'
nl_grid_cont_init_threshold='0.005'
nl_grid_cont_init_forcing='0'
nl_grid_cont_init_steps='4'
nl_grid_cont_skipsteps='0'
nl_grid_cont_maxsteps='8'
nl_visc_bounds_cont_min='1.0e-0'
nl_visc_bounds_cont_max='1.0e+2'
nl_visc_bounds_cont_steps='0'

# linear solver parameters
krylov_maxiter='600'
gmres_restart='150'

stress_dir_scale='1'
stress_gmg_coarse_coeff_pw_const='0'
stress_gmg_default_dir_scale='0'
stress_gmg_smoother_matrix_type='2'  # 1: diag, 2: point-block diag
stress_gmg_smoother_pos_lumping='0'
stress_gmg_smoother_submesh='0'
stress_gmg_smooth_coarse_solve='0'
stress_eig_est_max_it='10'
stress_nl_pc_without_tensor='0'

bfbt_inner_uscale_force_lumping='0'
bfbt_inner_uscale_force_submesh='0'
bbt_inner_uscale_low_order='0'
bbt_stiff_pc_mass_scaling='0'
bbt_stiff_pc_coeff_dr_max='1.0e9'
bbt_stiff_pc_modal_space_smoother_iter_multiply='1'
bbt_eig_est_max_it='10'

stiff_gmg='1'
stiff_gmg_coarse_coeff_pw_const='0'
stiff_gmg_smoother_matrix_type='1'  # 1: diag
stiff_gmg_smoother_pos_lumping='0'
stiff_gmg_smoother_submesh='0'
stiff_gmg_smooth_coarse_solve='0'
stiff_gmg_coarse_force_mass_scaling='0.0'  # 1.0e-4
stiff_gmg_vcycles='1'

Hminus1_norm_smoothing_iter='1'

chebyshev_est_eig_monitor='1'  # 0: off, 1: on

# null space projections
if [ "$PROBLEM_DOMAIN" = 'earth' ]; then # if earth domain
  stress_project_out_mean_rotation='2'  # 0: off, 1: on, 2: res, 3: symm
  stokes_project_out_mean_rotation='2'
  stress_gmg_cycle_project_out_rot='3'
  stress_gmg_level_project_out_rot='0'
  stress_setup_project_out_rot='1'  # 0: off, 1: on
else # if slice or brick domain
  stress_project_out_mean_rotation='0'
  stokes_project_out_mean_rotation='0'
  stress_gmg_cycle_project_out_rot='0'
  stress_gmg_level_project_out_rot='0'
  stress_setup_project_out_rot='0'
fi # end if domain
bbt_invert_project_out_mean_pressure='0'
bbt_project_out_mean_pressure='0'
bbt_gmg_project_out_mean_pressure='1'
bfbt_project_out_stokes_nullspace='0'
stokes_project_out_mean_pressure='0'  # 0: off, 1: on, 2: res, 3: symm

# GMG
gmg_type='1'  # 0: h-multrigrid, 1: ph-multigrid
gmg_max_levels='20'
gmg_coarsening_strategy='0'
gmg_coarsening_reduction_min='0.75'
gmg_mesh_minlevel='0';

# GMG partitioning and coarse problem size
gmg_n_nodes_per_proc_min='1800'
gmg_n_nodes_per_proc_max='2000'
gmg_reserve_procs_coarse='0'
if [ "$PROBLEM_DOMAIN" = 'earth' ]; then # if earth domain
  gmg_n_procs_coarse='512'
  gmg_n_cnodes_coarse='200000'
  gmg_n_elements_coarse='400'
else # if slice or brick domain
  gmg_n_procs_coarse='32'
  gmg_n_cnodes_coarse='20000'
  gmg_n_elements_coarse='400'
fi # end if domain

# AMG
amg_type='GAMG'
amg_max_levels='5'

# AMG partitioning and coarse problem size
amg_stress_mem_fraction='0.1'   # use this fraction of available memory for AMG
amg_stress_dof_reduction='0.04' # approx reduction of DOF at coarser levels
amg_stress_nzz_per_dof='70.0'   # avg #nonzeros per DOF at level 0
amg_stress_nzz_increase='4.0'   # appox increase of #nonzeros at coarser levels
amg_stress_mem_per_dof='96.0'   # appox memory per DOF (in Bytes)

amg_stress_dof_per_proc=$(printf '%.0f' $(echo "3.0*$gmg_n_cnodes_coarse*$amg_stress_dof_reduction/$gmg_n_procs_coarse" | bc -l))
amg_stress_dof_coarse=$(printf '%.0f' $(echo "sqrt( $amg_stress_mem_fraction*$arch_mem_max/$amg_stress_mem_per_dof - $amg_stress_dof_per_proc*$amg_stress_nzz_per_dof*$amg_max_levels*$amg_stress_nzz_increase )" | bc -l))

# AMG mesh
if [ 1 -le "$gmg_type" ]; then
  submesh='0'
else
  submesh='1'
fi
submesh_interp='0'
submesh_gll='1'
submesh_reweight='0'

# implementation parameters
optimized_matvecs='1'

# output
bin_out='1'
bin_out_input='0'
bin_out_nl_iter='1'
vtk_out='1'
vtk_out_input='0'
vtk_out_nl_iter='0'

# performance monitoriing and profiling
monitor_performance='1' # measure runtimes and flops
workload_out='1'        # write workload to file
profile_hpctoolkit='0'  # HPCToolkit
profile_tau_lib='0'     # TAU with instrumentation: library interposition
profile_tau_source='0'  # TAU with instrumentation: source transformation (PDT)
profile_ipm='0'         # IPM
profile_massif='0'      # Valgrind Massif (not working)

# performance test of matvecs
run_tensor_test_num='0'
run_matvec_test_stress_num='0'
run_matvec_test_stress_nl='0'
run_matvec_test_stress_dirty='1'
run_matvec_test_stress_optimized='1'
run_matvec_test_stokes_num='0'
run_matvec_test_stokes_nl='0'
run_matvec_test_stokes_dirty='1'
run_matvec_test_stokes_optimized='1'
run_matvec_test_bbt_num='0'
run_matvec_test_bbt_optimized='1'
run_matvec_test_stiff_num='0'
run_matvec_test_stiff_dirty='1'
run_matvec_test_stiff_optimized='1'

# job parameters
dry_run='1'
launch_multiple_jobs='0'
launch_vis_jobs='0'
queue='normal'
mpisize='4096'          # fine: >4096, coarse: 3072, scoarse: 512
ompsize='1'             # 1, 2, 4, 8, 16
vis_mpisize='256'       # fine: ??,    coarse: 160,  scoarse: 32
max_runtime='01:00:00'  # fine: ??h,   coarse: ??h,  scoarse: 12h
project_name='TG-DPP130002'
email='johann@ices.utexas.edu'

# paths for input
path_code="${HOME}/code/ymir/"
if [ 1 -le "$profile_tau_source" ]; then
  path_build="${HOME}/build/prof/ymir_tau/"
else
  #path_build="${HOME}/build/dev/ymir/"
  path_build="${HOME}/build/perf/ymir/"
fi
if [ 0 -lt "$run_tensor_test_num" ] || \
   [ 0 -lt "$run_matvec_test_stress_num" ] || \
   [ 0 -lt "$run_matvec_test_stokes_num" ] || \
   [ 0 -lt "$run_matvec_test_bbt_num" ]    || \
   [ 0 -lt "$run_matvec_test_stiff_num" ]; then
  executable="example/slabs/ymir_slabs_perf_matvec"
else
  executable="example/slabs/ymir_slabs"
fi
path_job_template="${path_code}example/slabs/stampede/job_template.sh"
path_vis_template="${HOME}/code/paraview/stampede_render_pv.sh"
if [ "$PROBLEM_DOMAIN" = 'earth' ]; then # if earth domain
  path_render_script="${HOME}/code/paraview/earth_render_pv.py"
else # if slice or brick domain
  path_render_script="${HOME}/code/paraview/slice_render_pv.py"
fi # end if domain
path_render_other="${HOME}/code/paraview/rhea_pv.py"
dir_input='input'

dir_input_override=''    # override directory to input files
dir_options_override=''  # override directory to option files

# paths for output
path_output="${SCRATCH}/"  # or "${WORK}/"
path_output+="runs/ymir/${PROBLEM_TYPE}_conv/${script_base_name}_$(date +%F)/"
dir_output_bin='bin'
dir_output_vis='vtk'
dir_output_render='png'

#######################################
# Functions for setting Petsc options
#######################################

function petsc_options_set_amg()
{
  prefix=$1
  amg_type=$2
  amg_max_num_levels=$3
  amg_drop_tol=$4
  amg_agg_n_smooth=$5
  amg_dof_per_proc=$6
  amg_dof_coarse=$7
  amg_num_vcycles=$8

  opt="
  ### AMG Galerkin coarse grid generation
  -${prefix}_pc_mg_levels $amg_max_num_levels
"
  if [ "$amg_type" = 'GAMG' ]; then
    opt+="
    # GAMG
    -${prefix}_pc_type gamg
    -${prefix}_pc_gamg_verbose 10
    -${prefix}_pc_gamg_threshold $amg_drop_tol
    -${prefix}_pc_gamg_agg_nsmooths $amg_agg_n_smooth
    -${prefix}_pc_gamg_process_eq_limit $amg_dof_per_proc
    -${prefix}_pc_gamg_coarse_eq_limit $amg_dof_coarse
   #-${prefix}_pc_gamg_repartition true"
  else
    opt+="
    # ML
    -${prefix}_pc_type ml
    -${prefix}_pc_ml_PrintLevel 10
    -${prefix}_pc_ml_CoarsenScheme Uncoupled
    -${prefix}_pc_ml_Threshold $amg_drop_tol
    -${prefix}_pc_ml_maxCoarseSize $amg_dof_coarse
   #-${prefix}_pc_ml_Reusable
   #-${prefix}_pc_ml_KeepAggInfo
    -${prefix}_pc_ml_repartition
    -${prefix}_pc_ml_repartitionMinPerProc $amg_dof_per_proc
    -${prefix}_pc_ml_repartitionPutOnSingleProc $amg_dof_coarse
   #-${prefix}_pc_ml_repartitionType ParMETIS  # Zoltan or ParMETIS"
  fi
  opt+="\n"

  opt+="
  ### AMG algorithmic options
  -${prefix}_pc_mg_type multiplicative
  -${prefix}_pc_mg_multiplicative_cycles $amg_num_vcycles"

  printf "$opt"
}

function petsc_options_set_mg_smoother()
{
  prefix=$1
  n_iter=$2
  smooth_accel=$3
  cheby_est_eig=$4
  cheby_est_eig_random=$5
  cheby_est_eig_n_iter=$6
  cheby_est_eig_monitor=$7
  smooth_type=$8
  smooth_icc_shift=$9
  smooth_use_amat=${10}

  # set KSP options

  opt="
    ### KSP type
"
  if [ "$smooth_accel" = 'CG' ]; then # if smooth accel with CG
    opt+="
      # CG
      -${prefix}_ksp_type cg
      -${prefix}_ksp_cg_single_reduction  # combine 2 i.prod. into 1 allreduce

    ### KSP tolerances
    -${prefix}_ksp_rtol 1.0e-16
    -${prefix}_ksp_atol 0.0
    -${prefix}_ksp_max_it $n_iter
    -${prefix}_ksp_norm_type none"
  elif [ "$smooth_accel" = 'Chebyshev' ]; then # if smooth accel Chebyshev
    opt+="
      # Chebyshev smoother
      -${prefix}_ksp_type chebyshev"
    if [ -n "$cheby_est_eig" ]; then # if estimate eigenvalues
      opt+="
      -${prefix}_ksp_chebyshev_estimate_eigenvalues"
      if [ -n "$cheby_est_eig_random" ]; then # if estimate with random RHS
        opt+="
      -${prefix}_ksp_chebyshev_estimate_eigenvalues_random"
      fi
      opt+="
      -${prefix}_est_ksp_type gmres
      -${prefix}_est_ksp_pc_side left  # same as chebyshev
      -${prefix}_est_ksp_gmres_restart $cheby_est_eig_n_iter
      -${prefix}_est_ksp_max_it $cheby_est_eig_n_iter"
      if [ 1 -le "$cheby_est_eig_monitor" ]; then # if monitor singular values
        opt+="
      -${prefix}_est_ksp_monitor_singular_value"
      else
        opt+="
      -${prefix}_est_ksp_norm_type none"
      fi
    fi
    opt+="\n"
    opt+="
    ### KSP tolerances
    -${prefix}_ksp_rtol 1.0e-16
    -${prefix}_ksp_atol 0.0
    -${prefix}_ksp_max_it $n_iter
    -${prefix}_ksp_norm_type none"
  else # if no smooth accel
    opt+="
      # none
      -${prefix}_ksp_type preonly

    ### KSP tolerances
    -${prefix}_ksp_max_it 1"
  fi # end if smoothing acceleration
  opt+="\n"

  # set PC options

  opt+="
    ### PC type
"
  if [ "$smooth_type" = 'Jacobi' ]; then # if point Jacobi smoother
    opt+="
      # Point Jacobi
      -${prefix}_pc_type jacobi"
  elif [ "$smooth_type" = 'pbJacobi' ]; then # if point block Jacobi smoother
    opt+="
      # Point block Jacobi
      -${prefix}_pc_type pbjacobi"
  elif [ "$smooth_type" = 'SOR' ]; then # if SOR smoother
    opt+="
      # SOR
      -${prefix}_pc_type sor"
  elif [ "$smooth_type" = 'BlockSSOR' ]; then # if block Jac + SSOR smoother
    opt+="
      # Block Jacobi with SSOR as subsolver
      -${prefix}_pc_type bjacobi
      -${prefix}_sub_ksp_type preonly
      -${prefix}_sub_pc_type sor
      -${prefix}_sub_pc_sor_symmetric"
  elif [ "$smooth_type" = 'BlockSSORe' ]; then # if block Jac + SSOR Eisenstat
    opt+="
      # Block Jacobi with SSOR (with Eisenstat's trick to reduce computation)
      -${prefix}_pc_type bjacobi
      -${prefix}_sub_ksp_type preonly
      -${prefix}_sub_pc_type eisenstat"
  elif [ "$smooth_type" = 'BlockICC' ]; then # if block ICC smoother
    opt+="
      # Block Jacobi with ICC(0) as subsolver
      -${prefix}_pc_type bjacobi
      -${prefix}_sub_ksp_type preonly
      -${prefix}_sub_pc_type icc
      -${prefix}_sub_pc_factor_levels 0
      -${prefix}_sub_pc_factor_in_place"
    if [ -n "$smooth_icc_shift" ]; then # if shift not NULL
      opt+="
      -${prefix}_sub_pc_factor_shift_type $smooth_icc_shift"
    fi
  else # otherwise do not set pc_type
    opt+="
      #-${prefix}_pc_type NOT_SET"
  fi # end if smoother type

  if [ -n "$smooth_use_amat" ]; then # if use matrix that defines the lin sys
    opt+="\n"
    opt+="
    -${prefix}_pc_use_amat $smooth_use_amat"
  fi

  printf "$opt"
}

function petsc_options_set_mg_coarse()
{
  prefix=$1
  n_iter=$2
  coarse_ksp=$3
  coarse_type=$4
  coarse_shift=$5

  opt="
    ### KSP type
"
  if [ "$coarse_ksp" = 'CG' ]; then # if coarse solve with CG
    opt+="
      # CG
      -${prefix}_mg_coarse_ksp_type cg
      -${prefix}_mg_coarse_ksp_cg_single_reduction

    ### KSP tolerances
    -${prefix}_mg_coarse_ksp_rtol 1.0e-16
    -${prefix}_mg_coarse_ksp_atol 0.0
    -${prefix}_mg_coarse_ksp_max_it $n_iter
    -${prefix}_mg_coarse_ksp_norm_type none"
  else # if direct solve (only this seems to work here!)
    opt+="
      # none
      -${prefix}_mg_coarse_ksp_type preonly

    ### KSP tolerances
    -${prefix}_mg_coarse_ksp_max_it 1"
  fi # end if coarse KSP
  opt+="\n"

  opt+="
    ### PC type
"
  if [ "$coarse_type" = 'Cholesky' ]; then # if Cholesky coarse solver
    opt+="
      # Block Jacobi with Cholesky factorization as subsolver
      -${prefix}_mg_coarse_pc_type bjacobi
      -${prefix}_mg_coarse_sub_ksp_type preonly
      -${prefix}_mg_coarse_sub_pc_type cholesky"
    if [ -n "$coarse_shift" ]; then # if shift not NULL
      opt+="
      -${prefix}_mg_coarse_sub_pc_factor_shift_type $coarse_shift"
    fi
  else # otherwise LU coarse solver
    opt+="
      # Block Jacobi with LU factorization as subsolver
      -${prefix}_mg_coarse_pc_type bjacobi
      -${prefix}_mg_coarse_sub_ksp_type preonly
      -${prefix}_mg_coarse_sub_pc_type lu"
    if [ -n "$coarse_shift" ]; then # if shift not NULL
      opt+="
      -${prefix}_mg_coarse_sub_pc_factor_shift_type $coarse_shift"
    fi
  fi # end if smoother type

  printf "$opt"
}

#######################################
# Run all tests
#######################################

echo "Start: $script_name"

# collect all input and output in this directory
mkdir "$path_output"

# set number of nodes
num_nodes="$(echo "$mpisize * $ompsize / 16" | bc)"

# create script that launches all simulations
script_launch_all="${path_output}${script_launch_all_filename}"
touch "$script_launch_all"
chmod u+x "$script_launch_all"
echo "#!/bin/bash -e
#SBATCH -J ${script_base_name}
#SBATCH -o ${script_base_name}.o%j
#SBATCH -e ${script_base_name}.e%j
#SBATCH -p ${queue}
#SBATCH -N ${num_nodes}
#SBATCH -n ${mpisize}
#SBATCH -t ${max_runtime}
#SBATCH -A ${project_name}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${email}
" 1>>"$script_launch_all"

# create directory for all input data and copy input files
mkdir "${path_output}${dir_input}"
path_executable="${path_output}${dir_input}/$(basename $executable)"
cp "${path_build}${executable}" "$path_executable"
if [ -n "$dir_input_override" ]; then # if override input directory
  path_executable="${dir_input_override}$(basename $executable)"
fi
if [ "$PROBLEM_DOMAIN" = 'earth' ]; then # if earth domain
  temperature_bin="${path_output}${dir_input}/$(basename $path_temp_bin)"
  temperature_txt="${path_output}${dir_input}/$(basename $path_temp_txt)"
  weakzone_bin="${path_output}${dir_input}/$(basename $path_weak_bin)"
  weakzone_txt="${path_output}${dir_input}/$(basename $path_weak_txt)"
  cp "$path_temp_bin" "$temperature_bin"
  cp "$path_weak_bin" "$weakzone_bin"
  if [ "$read_temp_weak_bin_only" -le 0 ]; then # if import from txt file
    cp "$path_temp_txt" "$temperature_txt"
    cp "$path_weak_txt" "$weakzone_txt"
  fi
  if [ -n "$dir_input_override" ]; then # if override input directory
    temperature_bin="${dir_input_override}$(basename $path_temp_bin)"
    temperature_txt="${dir_input_override}$(basename $path_temp_txt)"
    weakzone_bin="${dir_input_override}$(basename $path_weak_bin)"
    weakzone_txt="${dir_input_override}$(basename $path_weak_txt)"
  fi
fi

# output to logfile
touch "$logfile"
echo "# $script_name
#
# physics parameters:
#   domain-shape = $domain_shape
#   boundary-cond = $boundary_condition
#   boundary-set-default-dirichlet-scale = $boundary_default_dir_scale
#   viscosity = $viscosity
#   viscosity-model = ${viscosity_model[*]}
#   viscosity-min = ${viscosity_min[*]}
#   viscosity-max = ${viscosity_max[*]}
#   viscosity-temp-max = ${viscosity_temp_max[*]}
#   viscosity-scaling = $viscosity_scaling
#   viscosity-temp-decay = ${viscosity_decay[*]}
#   viscosity-lower-mantle = $lower_mantle
#   viscosity-lower-mantle-scaling = $lower_mantle_scaling
#   viscosity-stress-exponent = $stress_exp
#   viscosity-stress-yield = $stress_yield
#   viscosity-yielding-reg = $stress_yielding_reg
#   right-hand-side-scaling = $rhs_scaling
#
# discretiztion parameters:
#   order = ${order[*]}
#   pressure space = ${pressure_space[*]}
#   maxlevel = ${maxlevel[*]}
#   init-amr-num-elements-max = ${init_amr_n_elements_max[*]}
#   init-amr-visc-tol-max = ${init_amr_visc_tol_max[*]}
#   init-amr-post-uniform-num-steps = ${init_amr_post_uniform_n_steps[*]}
#   amr-num-elements-max = ${amr_n_elements_max[*]}
#   amr-visc-tol-max = ${amr_visc_tol_max[*]}
#   amr-visc-dr-tol-max = ${amr_visc_dr_tol_max[*]}
#   amr-strain-rate-tol-max = ${amr_strain_rate_tol_max[*]}
#
# solver parameters:
#   Nonlinear primal-dual type = ${nonlinear_primaldual[*]}
#   Nonlinear primal-dual scaling = ${nonlinear_primaldual_scal[*]}
#   Nonlinear maxiter = $nl_maxiter
#   Nonlinear rtol = $nl_rtol
#   Krylov maxiter = $krylov_maxiter
#   GMRES restart = $gmres_restart
#   GMG type = $gmg_type
#   GMG max h-levels = $gmg_max_levels
#   AMG type = $amg_type
#   AMG levels = $amg_max_levels
#   schur-type = ${schur_type[*]}
#   BFBT type = ${bfbt_type[*]}
#   solve-stress-block-only = ${solve_stress_block_only[*]}
#   solve-press-bbt-only = ${solve_press_bbt_only[*]}
#   solve-cnode-bbt-only = ${solve_cnode_bbt_only[*]}
#
# job parameters:
#   queue = $queue
#   mpisize = $mpisize
#   ompsize = $ompsize
#   #nodes = $num_nodes
#   max runtime = $max_runtime
#" 1>>"$logfile"

# count elements
num_primaldual=${#nonlinear_primaldual[*]}
num_orders=${#order[*]}
num_visc_decays=${#viscosity_decay[*]}
num_visc_models=${#viscosity_model[*]}
num_schur_types=${#schur_type[*]}
num_solve_only=${#solve_stress_block_only[*]}

for (( p=0; p<$num_primaldual; p++ )); do
for (( k=0; k<$num_orders; k++ )); do
for (( d=0; d<$num_visc_decays; d++ )); do
for (( m=0; m<$num_visc_models; m++ )); do
for (( s=0; s<$num_schur_types; s++ )); do

for weak in $weak_factor ; do
for visc_init in $viscosity_for_init_nl_stokes ; do
for nl_type in $nonlinear_type ; do
for nl_forcing_max in $nonlinear_forcing_max ; do
for nl_step_length_min in $nonlinear_step_length_min ; do

for nl_norm in $nonlinear_norm_type ; do
for stress_use_gmg in $stress_gmg ; do
for stress_gmg_n_smoothing in $stress_gmg_smoothing_iter ; do
for stress_gmg_cycle_t in $stress_gmg_cycle_type ; do
for stress_gmg_n_cycles in $stress_gmg_num_cycles ; do

for stress_drop_tol in $stress_amg_drop_tol ; do
for stress_agg_nsmooths in $stress_amg_agg_nsmooths ; do
for stress_smoother in $stress_amg_smoother ; do
for stress_amg_n_smoothing in $stress_amg_smoothing_iter ; do
for stress_vcycles in $stress_amg_vcycles ; do

for bbt_use_stiff_pc in $bbt_stiff_pc ; do
for bbt_use_gmg in $bbt_gmg ; do
for bbt_gmg_n_smoothing in $bbt_gmg_smoothing_iter ; do
for bbt_gmg_cycle_t in $bbt_gmg_cycle_type ; do
for bbt_gmg_n_cycles in $bbt_gmg_num_cycles ; do

for bbt_drop_tol in $bbt_amg_drop_tol ; do
for bbt_agg_nsmooths in $bbt_amg_agg_nsmooths ; do
for bbt_smoother in $bbt_amg_smoother ; do
for bbt_smoothing_iter in $bbt_amg_smoothing_iter ; do
for bbt_vcycles in $bbt_amg_vcycles ; do

for (( o=0; o<$num_solve_only; o++ )); do

  # skip order 1 and Schur type BFBT, because it's not implemented
  if [ "${order[k]}" -eq 1 ] && [ "${schur_type[s]}" -eq 1 ]; then
    continue
  fi

  # set basename
  base_name="${PROBLEM_DOMAIN}"
  base_name+="_bc${boundary_condition}"
 #base_name+="_weak${weak}_${weak_thickness}m"
 #base_name+="_visc${viscosity_model[m]}"
 #base_name+="_iAMRtol${init_amr_visc_tol_max[d]}"
 #base_name+="_postRefine${init_amr_post_uniform_n_steps[d]}"
  base_name+="_L${maxlevel[k]}_k${order[k]}"
 #base_name+="_nl${nl_type}_pd${nonlinear_primaldual[p]}"
 #base_name+="_Agmg${stress_use_gmg}"
 #base_name+="_BBTstiff${bbt_use_stiff_pc}"
 #base_name+="_BBTgmg${bbt_use_gmg}"
 #base_name+="_dropTol${stress_drop_tol}"
 #base_name+="_amgVcycl${stress_vcycles}_amgSmooth${stress_amg_n_smoothing}"
  base_name+="_schur${schur_type[s]}"
  base_name+="_bfbt${bfbt_type[s]}"
  if [ 1 -le "${solve_stress_block_only[o]}" ]; then
    base_name+="_Aonly"
  fi
  if [ 1 -le "${solve_press_bbt_only[o]}" ]; then
    base_name+="_BBTonly"
  fi
  if [ 1 -le "${solve_cnode_bbt_only[o]}" ]; then
    base_name+="_cnodeBBTonly"
  fi

  # output
  echo "Process: $base_name"
  if [ 1 -le "${solve_stress_block_only[o]}" ]; then
    echo "Warning: Solving viscous stress block only!"
  fi
  if [ 1 -le "${solve_press_bbt_only[o]}" ]; then
    echo "Warning: Solving pressure BB^T only!"
  fi
  if [ 1 -le "${solve_cnode_bbt_only[o]}" ]; then
    echo "Warning: Solving continuous approximation of BB^T only!"
  fi

  # set directory name
  dir_name="$base_name"

  # create directories for all output data
  mkdir "${path_output}${dir_name}"
  mkdir "${path_output}${dir_name}/${dir_output_bin}"
  mkdir "${path_output}${dir_name}/${dir_output_vis}"
  mkdir "${path_output}${dir_name}/${dir_output_render}"

  # copy relevant files
  cp "${path_build}config.log" "${path_output}${dir_name}/config.log"

  # convert numbers from scientific notation
  nl_forcing_max_dec=$(echo "$nl_forcing_max" | sed -e 's/[eE]+*/\*10\^/')
  stiff_gmg_coarse_force_mass_scaling_dec=$(echo "$stiff_gmg_coarse_force_mass_scaling" | sed -e 's/[eE]+*/\*10\^/')

  # set dependent parameters
  nl_forcing_saveguard_threshold=$(echo "0.1 * $nl_forcing_max_dec" | bc -l)
  krylov_rtol=$nl_forcing_max
  if [ 0 -lt "${nl_scaling[s]}" ]; then
    stokes_scaling='1'
  else
    stokes_scaling='0'
  fi
  bbt_stiff_pc_modal_space_smoother_iter=$(echo "$bbt_gmg_n_smoothing * $bbt_stiff_pc_modal_space_smoother_iter_multiply" | bc -l)

  # set option filepaths
  ymir_options_filepath="${path_output}${dir_name}/${base_name}.ini"
  petsc_options_filepath="${path_output}${dir_name}/${base_name}.petsc"

  #
  # Create ymir options file
  #

  ymir_options=""

  ymir_options+="###### Global Options ######\n"
  ymir_options+="
[Options]
# Override physics options
  viscosity-override-const = 0
  weakzone-override-factor = 0
# Performance and load balance statistics
  monitor-performance = $monitor_performance"
  if [ 1 -le "$workload_out" ]; then
    ymir_options+="\n  workload-out-file = ${path_output}${dir_name}/workload"
  fi

  if [ 1 -le "$bin_out" ]; then
    ymir_options+="
# Binary output
  bin-out = \"${path_output}${dir_name}/${dir_output_bin}/${base_name}\"
  bin-out-input = $bin_out_input
  bin-out-nonlinear-iter = $bin_out_nl_iter
  bin-out-solution = 1"
  fi

  if [ 1 -le "$vtk_out" ]; then
    ymir_options+="
# VTK output
  vtk-out = \"${path_output}${dir_name}/${dir_output_vis}/${base_name}\"
# vtk-out-input-state = 1
  vtk-out-input = $vtk_out_input
  vtk-out-nonlinear-iter = $vtk_out_nl_iter"
  else
    ymir_options+="
# VTK output
# vtk-out = \"${path_output}${dir_name}/${dir_output_vis}/${base_name}\"
# vtk-out-input-state = 1"
  fi

  ymir_options+="
# exit-after-mesh-setup = 1"

  if [ 0 -lt "$run_tensor_test_num" ] || \
     [ 0 -lt "$run_matvec_test_stress_num" ] || \
     [ 0 -lt "$run_matvec_test_stokes_num" ] || \
     [ 0 -lt "$run_matvec_test_bbt_num" ]    || \
     [ 0 -lt "$run_matvec_test_stiff_num" ]; then
    ymir_options+="
# Performance test of matvec
  tensor-test-num = $run_tensor_test_num
  matvec-test-stress-num = $run_matvec_test_stress_num
  matvec-test-stress-nl = $run_matvec_test_stress_nl
  matvec-test-stress-dirty = $run_matvec_test_stress_dirty
  matvec-test-stress-optimized = $run_matvec_test_stress_optimized
  matvec-test-stokes-num = $run_matvec_test_stokes_num
  matvec-test-stokes-nl = $run_matvec_test_stokes_nl
  matvec-test-stokes-dirty = $run_matvec_test_stokes_dirty
  matvec-test-stokes-optimized = $run_matvec_test_stokes_optimized
  matvec-test-bbt-num = $run_matvec_test_bbt_num
  matvec-test-bbt-optimized = $run_matvec_test_bbt_optimized
  matvec-test-stiff-num = $run_matvec_test_stiff_num
  matvec-test-stiff-dirty = $run_matvec_test_stiff_dirty
  matvec-test-stiff-optimized = $run_matvec_test_stiff_optimized"
  fi

  ymir_options+="\n\n###### Mantle Options ######\n"
  ymir_options+="
[Mantle]
# Physics
  domain-shape = $domain_shape
  boundary-cond = $boundary_condition
  boundary-set-default-dirichlet-scale = $boundary_default_dir_scale
# p4est-import = "
  if [ "$PROBLEM_DOMAIN" = 'earth' ]; then # if earth domain
    ymir_options+="
  temperature = $temperature_bin
  temp-background-plate-age = $temp_import_background_plate_age
  temp-import-plate-age-min = $temp_import_plate_age_min
  weakzone = $weakzone_bin
  weakzone-import-pointcloud-size = $weakzone_size
  weakzone-import-thickness = $weak_thickness
  weakzone-import-thickness-const = $weak_thickness_const
  weakzone-import-weak-factor = $weak

  # only required once for creating temperature and weak zone .bin files"
    if [ "$read_temp_weak_bin_only" -le 0 ]; then # if import from txt file
      ymir_options+="
  temp-import-read-textfile = $temperature_txt
  weakzone-import-read-textfile = $weakzone_txt\n"
    else
      ymir_options+="
# temp-import-read-textfile = $temperature_txt
# weakzone-import-read-textfile = $weakzone_txt\n"
    fi
    ymir_options+="
# velocity-import =
# pressure-import = "
  else # if slice or brick domain
    ymir_options+="
  domain-brick-dx = 1
  domain-brick-dy = 16
  domain-brick-dz = 16
  temperature = 2plates_poly2
  temp-2plates-trench-longitude = 0.13
  temp-2plates-dip-angle = 5.0
  temp-2plates-subd-depth = 400.0e3
  temp-2plates-subd-width = 500.0e3
  temp-2plates-subd-edge-width = 1.0e3
  temp-2plates-subd-edge-smoothwidth = 40.0e3
  temp-2plates-subd-plate-velocity = 4.0e-2
  temp-2plates-subd-plate-initial-age = 1.0e6
  temp-2plates-over-plate-age = 40.0e6
  weakzone = 2plates_poly2
  weakzone-2plates-subdu-depth = 60.0e3
  weakzone-2plates-subdu-thickness = $weak_thickness
  weakzone-2plates-subdu-thickness-const = $weak_thickness_const
  weakzone-2plates-subdu-weak-factor = $weak
  weakzone-2plates-ridge-depth = 15.0e3
  weakzone-2plates-ridge-width = 15.0e3
  weakzone-2plates-ridge-smoothwidth = $weak_ridge_smoothwidth
  weakzone-2plates-ridge-weak-factor = $weak"
  fi # end if domain
  ymir_options+="
  viscosity = $viscosity
  viscosity-for-init-nl-stokes = $visc_init
  viscosity-model = ${viscosity_model[m]}
  viscosity-min = ${viscosity_min[m]}
  viscosity-max = ${viscosity_max[m]}
  viscosity-temp-max = ${viscosity_temp_max[m]}
  viscosity-scaling = $viscosity_scaling
  viscosity-temp-decay = ${viscosity_decay[d]}
  viscosity-lower-mantle = $lower_mantle
  viscosity-lower-mantle-scaling = $lower_mantle_scaling
  viscosity-stress-exponent = $stress_exp
  viscosity-stress-yield = $stress_yield
  viscosity-yielding-reg = $stress_yielding_reg
  viscosity-coarsen-eval = 0
  viscosity-p-coarsen-eval = 0
  viscosity-coarsen-type = 0
  weakzone-coarsen-type = 0
  right-hand-side-scaling = $rhs_scaling
  right-hand-side-random = 0
  right-hand-side-multiply-in-weakzone = 1
  plume-type = 0
# Discretization
  order = ${order[k]}
  minlevel = $minlevel
  maxlevel = ${maxlevel[k]}
  refine = uniform
  refine-radius = $refine_radius
  import-mesh-order = $import_mesh_order
  import-mesh-minlevel = $import_mesh_minlevel
  refine-surface-max-dist = $refine_surface_maxdist
  refine-surface-element-resolution = $refine_surface_elem_res
  refine-lower-upper-interface-max-dist = $refine_lm_um_maxdist
  refine-lower-upper-interface-element-resolution = $refine_lm_um_elem_res
  init-amr-rel-threshold = 0
  init-amr-num-elements-max = ${init_amr_n_elements_max[d]}
  init-amr-override-order = $init_amr_override_order
  init-amr-lower-mantle = 0
  init-amr-visc-indicator = $init_amr_visc_indicator
  init-amr-visc-tol-min = $init_amr_visc_tol_min
  init-amr-visc-tol-max = ${init_amr_visc_tol_max[d]}"
  if [ "$PROBLEM_DOMAIN" = 'earth' ]; then # if earth domain
    ymir_options+="
  init-amr-weak-import-indicator = $init_amr_weak_indicator
  init-amr-weak-import-element-resolution = $init_amr_weak_elem_res"
  else # if slice or brick domain
    ymir_options+="
  init-amr-weak-subdu-indicator = $init_amr_weak_indicator
  init-amr-weak-subdu-element-resolution = $init_amr_weak_elem_res"
  fi # end if domain
  ymir_options+="
  init-amr-rhs-indicator = $init_amr_rhs_indicator
  init-amr-rhs-tol-min = $init_amr_rhs_tol_min
  init-amr-rhs-tol-max = $init_amr_rhs_tol_max
  init-amr-rhs-norm-shift = $init_amr_rhs_norm_shift
  init-amr-post-uniform-num-steps = ${init_amr_post_uniform_n_steps[d]}
  amr-rel-threshold = $amr_rel_threshold
  amr-num-elements-max = ${amr_n_elements_max[d]}
  amr-lower-mantle = 0
  amr-visc-indicator = $amr_visc_indicator
  amr-visc-tol-min = $amr_visc_tol_min
  amr-visc-tol-max = ${amr_visc_tol_max[d]}
  amr-visc-dr-indicator = $amr_visc_dr_indicator
  amr-visc-dr-tol-min = $amr_visc_dr_tol_min
  amr-visc-dr-tol-max = ${amr_visc_dr_tol_max[d]}
  amr-strain-rate-indicator = $amr_strain_rate_indicator
  amr-strain-rate-tol-min = $amr_strain_rate_tol_min
  amr-strain-rate-tol-max = ${amr_strain_rate_tol_max[d]}
  amr-log-maxlevel = 1
  mesh-partitioning-type = $mesh_partitioning_type
# Nonlinear solver
  nonlinear-type = $nl_type
  nonlinear-primaldual-type = ${nonlinear_primaldual[p]}
  nonlinear-primaldual-scaling-type = ${nonlinear_primaldual_scal[p]}
  nonlinear-maxiter = $nl_maxiter
  nonlinear-rtol = $nl_rtol
  nonlinear-forcing-exponent = $nl_forcing_exponent
  nonlinear-forcing-max = $nl_forcing_max
  nonlinear-forcing-max-progressive-iter = $nl_forcing_max_progressive_iter
  nonlinear-forcing-total-min = $nl_forcing_total_min
  nonlinear-forcing-saveguard = $nl_forcing_saveguard
  nonlinear-forcing-saveguard-threshold = $nl_forcing_saveguard_threshold
  nonlinear-step-length-reduction-type = $nl_step_reduction_type
  nonlinear-step-length-reduction-min = 0.2
  nonlinear-step-length-reduction-max = 0.8
  nonlinear-step-length-reduction-reg = 1.0e-12
  nonlinear-step-length-min = $nl_step_length_min
  nonlinear-switch-picard-step-length-min = $nl_switch_picard_step_length_min
  nonlinear-switch-picard-after-amr = $nl_switch_picard_after_amr
  nonlinear-switch-picard-init = $nl_switch_picard_init
  nonlinear-switch-picard-maxiter = $nl_switch_picard_maxiter
  nonlinear-switch-picard-rtol = $nl_switch_picard_rtol
  nonlinear-norm-type = $nl_norm
  nonlinear-norm-Hminus1-mass-scaling = $nl_norm_Hminus1_mass_scaling
  nonlinear-initial-guess-type = $nl_initial_guess_type
  nonlinear-schur-diag-type = $nl_schur_diag_type
  nonlinear-project-out-nullspace = $nl_project_out_nullspace
  nonlinear-scaling-type = ${nl_scaling[s]}
  nonlinear-enforce-unscaled-reduction = $nl_enforce_unscaled_reduction
  nonlinear-grid-continuation-init-threshold = $nl_grid_cont_init_threshold
  nonlinear-grid-continuation-init-forcing = $nl_grid_cont_init_forcing
  nonlinear-grid-continuation-init-steps = $nl_grid_cont_init_steps
  nonlinear-grid-continuation-skipsteps = $nl_grid_cont_skipsteps
  nonlinear-grid-continuation-maxsteps = $nl_grid_cont_maxsteps
  nonlinear-viscosity-bounds-continuation-min = $nl_visc_bounds_cont_min
  nonlinear-viscosity-bounds-continuation-max = $nl_visc_bounds_cont_max
  nonlinear-viscosity-bounds-continuation-steps = $nl_visc_bounds_cont_steps
  nonlinear-check-derivative = 0
  log-physics-stats = 1
# Krylov solver
  krylov-maxiter = $krylov_maxiter
  krylov-atol = 0.0
  krylov-rtol = $krylov_rtol
  solve-stress-block-only = ${solve_stress_block_only[o]}
  solve-press-bbt-only = ${solve_press_bbt_only[o]}
  solve-cnode-bbt-only = ${solve_cnode_bbt_only[o]}"

  ymir_options+="\n\n###### Petsc Options ######\n"
  if [ -n "$dir_options_override" ]; then # if override options directory
    ymir_options+="
[PETSc]
  file = ${dir_options_override}$(basename $petsc_options_filepath)"
  else
    ymir_options+="
[PETSc]
  file = $petsc_options_filepath"
  fi

  ymir_options+="\n\n###### Viscous Stress ######\n"
  ymir_options+="
[Stress:Op]
  set-dirichlet-scale = $stress_dir_scale
  log-coeff-range = 1
  nsp-test = 0
  project-out-transl = 0
  project-out-rot = $stress_project_out_mean_rotation
  log-nsp = 0
  monitor-performance = $monitor_performance
[Stress:PC]
  picard = $stress_nl_pc_without_tensor
  gmg = $stress_use_gmg
  amg-rbm = 1
  amg-reblock = 0
  ksp-setup-immediately = 1
  setup-project-out-rot = $stress_setup_project_out_rot
  optimized-matvec = $optimized_matvecs
  monitor-performance = $monitor_performance"

  ymir_options+="\n\n###### Stiffness ######\n"
  ymir_options+="
[Stiff:Op]
  log-coeff-range = 1
  log-coeff-anisotropy = 1
  project-out-mean = 0
  log-nsp = 0
  monitor-performance = $monitor_performance
[Stiff:PC]
  gmg = $stiff_gmg
  ksp-setup-immediately = 1
  optimized-matvec = $optimized_matvecs
  monitor-performance = $monitor_performance"

  ymir_options+="\n\n###### Stokes ######\n"
  ymir_options+="
[Pressure:Elem]
  default-space = ${pressure_space[k]}

[Stokes:Op]
  project-out-rot = $stokes_project_out_mean_rotation
  project-out-mean = $stokes_project_out_mean_pressure
  log-nsp = 0
  monitor-performance = $monitor_performance
[Stokes:PC]
  block-type = 1
  schur-type = ${schur_type[s]}
  project-out-rot = 0
  project-out-mean = 0
  use-scaling = $stokes_scaling
  log-residual-comp = 0
  optimized-matvec = $optimized_matvecs
  monitor-performance = $monitor_performance
[Stokes:PC:BFBT]
  bfbt-type = ${bfbt_type[s]}
  inner-uscale-force-lumping = $bfbt_inner_uscale_force_lumping
  inner-uscale-force-submesh = $bfbt_inner_uscale_force_submesh
  project-out-stokes-nullspace = $bfbt_project_out_stokes_nullspace
  log-mean = 0
  monitor-performance = $monitor_performance
[Stokes:PC:BBT]
  stiffness-pc = ${bbt_use_stiff_pc}
  stiffness-pc-mass-scaling = $bbt_stiff_pc_mass_scaling
  stiffness-pc-project-out-mean = 1
  stiffness-pc-log-mean-dbg = 0
  stiffness-pc-coeff-dr-max = $bbt_stiff_pc_coeff_dr_max
  stiffness-pc-coarsen-coeff-type = 0
  stiffness-pc-override-coeff-with-viscosity = 0
  stiffness-pc-num-cycles = 1
  inner-uscale-low-order = $bbt_inner_uscale_low_order
  gmg = $bbt_use_gmg
  amg-reblock = 0
  ksp-setup-immediately = 1
  invert-project-out-mean = $bbt_invert_project_out_mean_pressure
  project-out-mean = $bbt_project_out_mean_pressure
  log-mean = 0
  optimized-matvec = $optimized_matvecs
  monitor-performance = $monitor_performance"

  ymir_options+="\n\n###### GMG ######\n"
  ymir_options+="
[GMG:HierarchyMesh]
  gmg-type = $gmg_type
  max-h-levels = $gmg_max_levels
  coarsening-strategy = $gmg_coarsening_strategy
  coarsening-strategy-coeff-jump = 1.0e3
  coarsening-reduction-min = $gmg_coarsening_reduction_min
  mesh-minlevel = $gmg_mesh_minlevel"
  if [ 1 -le "$gmg_reserve_procs_coarse" ]; then
    ymir_options+="
  num-procs-coarse = $gmg_n_procs_coarse"
  else
    ymir_options+="
  num-procs-coarse = 0"
  fi
  ymir_options+="
  num-cnodes-coarse = $gmg_n_cnodes_coarse
  num-elements-coarse = $gmg_n_elements_coarse
  partition-nodes-per-proc-min = $gmg_n_nodes_per_proc_min
  partition-nodes-per-proc-max = $gmg_n_nodes_per_proc_max
  partition-subset-type = 1"
  if [ "$PROBLEM_TYPE" = 'lin' ]; then
    ymir_options+="
  partition-skip-num-levels = ${init_amr_post_uniform_n_steps[d]}"
  fi
  ymir_options+="
  reduce-mpi-comm = 1
  monitor-performance = $monitor_performance
[GMG:HierarchyStress]
  coefficient-interp-linear = 0
  coefficient-min = 2.0e-2
  coefficient-max = 2.0e+4
  coefficient-bound-min-linear-interp = 0
  coefficient-piecewise-constant = $stress_gmg_coarse_coeff_pw_const
  set-default-dirichlet-scale = $stress_gmg_default_dir_scale
  smoother-matrix-type = $stress_gmg_smoother_matrix_type
  smoother-matrix-positive-lumping = $stress_gmg_smoother_pos_lumping
  smoother-matrix-submesh = $stress_gmg_smoother_submesh
  smooth-coarse-solve = $stress_gmg_smooth_coarse_solve
  cycle-project-out-rot = $stress_gmg_cycle_project_out_rot
  level-project-out-rot = $stress_gmg_level_project_out_rot
  cycle-type = $stress_gmg_cycle_t
  num-cycles = $stress_gmg_n_cycles
  monitor-performance = $monitor_performance"
  if [ 1 -le "$workload_out" ]; then
    ymir_options+="
  workload-out-file = ${path_output}${dir_name}/workload_gmg_stress"
  fi
  ymir_options+="
[GMG:HierarchyStiff]
  coefficient-interp-linear = 0
# coefficient-min = 1.0e-12
# coefficient-max = 1.0e+12
  coefficient-bound-min-linear-interp = 0
  coefficient-piecewise-constant = $stiff_gmg_coarse_coeff_pw_const
  smoother-matrix-type = $stiff_gmg_smoother_matrix_type
  smoother-matrix-positive-lumping = $stiff_gmg_smoother_pos_lumping
  smoother-matrix-submesh = $stiff_gmg_smoother_submesh
  smooth-coarse-solve = $stiff_gmg_smooth_coarse_solve
  coarse-force-mass-scaling = $stiff_gmg_coarse_force_mass_scaling
  cycle-project-out-mean = 0
  level-project-out-mean = 0
  cycle-type = 0
  num-cycles = $stiff_gmg_vcycles
  monitor-performance = $monitor_performance
[GMG:HierarchyBBT]
  project-out-mean = $bbt_gmg_project_out_mean_pressure
  restrict-inner-scaling = 0
  cycle-type = $bbt_gmg_cycle_t
  num-cycles = $bbt_gmg_n_cycles
  monitor-performance = $monitor_performance"

  ymir_options+="\n\n###### Submesh (AMG mesh) ######\n"
  ymir_options+="
[Submesh]
  use-submesh = $submesh
  interp = $submesh_interp
  use-gll = $submesh_gll
  reweight = $submesh_reweight
  simple = 1"

  # write ymir options to file
  printf "${ymir_options}\n" 1>$ymir_options_filepath

  #
  # Create Petsc options file
  #

  petsc_options=""

  # set global options
  petsc_options+="\n\n###### Global ######\n"
  petsc_options+="
 #-log_summary  # summarize the program’s performance
 #-malloc_dump  # dumps unfreed memory during call to PetscFinalize
"
  if [ 1 -le "$workload_out" ]; then # if track Petsc memory usage
    petsc_options+="  -malloc  # enable Petsc malloc logging\n"
  fi

  #
  # Petsc options: mass_
  #

  petsc_options+="\n\n###### Mass Matrix Inversion KSP ######\n"
  petsc_options+="
  ### KSP type

    # CG
    -mass_ksp_type cg
    -mass_ksp_cg_single_reduction  # combine 2 inner prods into 1 allreduce
    -mass_ksp_pc_side left

  ### KSP tolerances
  -mass_ksp_rtol 1.0e-6
  -mass_ksp_atol 0
  -mass_ksp_max_it 25

  ### KSP monitoring
 #-mass_ksp_view
 #-mass_ksp_monitor
 #-mass_ksp_monitor_true_residual
 #-mass_ksp_monitor_singular_value
"

  petsc_options+="\n\n###### Mass Matrix Inversion PC ######\n"
  petsc_options+="
  # inner preconditioner uses the matrix that defines the linear system
  -mass_pc_use_amat true
"

  #
  # Petsc options: Hminus1_norm_
  #

  if [ 0 -le "$nl_type" ] || \
     [ 4 -eq "${bfbt_type[s]}" ] || \
     [ 7 -eq "${bfbt_type[s]}" ]; then # if H^-1 norm used

    petsc_options+="\n\n###### H^-1 norm KSP ######\n"

    if [ 0 -le "$nl_type" ]; then # if nonlinear solver
      petsc_options+="
  ### KSP type

    # CG
    -Hminus1_norm_ksp_type cg
    -Hminus1_norm_ksp_cg_single_reduction
    -Hminus1_norm_ksp_pc_side left

  ### KSP tolerances
  -Hminus1_norm_ksp_rtol 1.0e-6
  -Hminus1_norm_ksp_atol 0.0
  -Hminus1_norm_ksp_max_it 50
"
    else
      petsc_options+="
  ### KSP type
  -Hminus1_norm_ksp_type preonly

  ### KSP tolerances
  -Hminus1_norm_ksp_max_it 1
"
    fi # end if nonlinear solver

    petsc_options+="\n\n###### H^-1 norm PC ######\n"

    if [ 1 -le "$stiff_gmg" ]; then # if use GMG for stiffness
      Hminus1_gmg_prefix="Hminus1_norm_gmg"
      Hminus1_amg_prefix="${Hminus1_gmg_prefix}_coarse"
      petsc_options+="
  # inner preconditioner uses the matrix that defines the linear system
  -Hminus1_norm_pc_use_amat true

  ### GMG smoother
"
      petsc_options+="$(petsc_options_set_mg_smoother \
                            "${Hminus1_gmg_prefix}_levels" \
                            "$Hminus1_norm_smoothing_iter" \
                            "Chebyshev" "est_eig" "" "10" "0" \
                            "PCMAT" "" \
                            "true")\n"

      petsc_options+="
  ### GMG coarse KSP

    ### KSP type (use AMG preconditioner only)
    -${Hminus1_amg_prefix}_ksp_type preonly

    ### KSP tolerances
    -${Hminus1_amg_prefix}_ksp_max_it 1

  ### GMG coarse PC
"
    else # if use AMG for stiffness
      Hminus1_amg_prefix='Hminus1_norm'
    fi # end if use GMG for stiffness

    petsc_options+="
  # inner preconditioner uses the matrix that defines the linear system
  -${Hminus1_amg_prefix}_pc_use_amat false
"
    petsc_options+="$(petsc_options_set_amg \
                          $Hminus1_amg_prefix \
                          $amg_type \
                          $amg_max_levels \
                          "0.0" \
                          "1" \
                          $amg_stress_dof_per_proc \
                          $amg_stress_dof_coarse \
                          "1")\n"

    petsc_options+="
  ### AMG smoother
"
    petsc_options+="$(petsc_options_set_mg_smoother \
                          "${Hminus1_amg_prefix}_mg_levels" \
                          "$Hminus1_norm_smoothing_iter" \
                          "Chebyshev" "" "" "10" "" \
                          "Jacobi" "" \
                          "false")\n"

    petsc_options+="
  ### AMG coarse solver
"
    petsc_options+="$(petsc_options_set_mg_coarse \
                          "${Hminus1_amg_prefix}_mg_coarse" \
                          "1" "preonly" "Cholesky")\n"
  fi # end if H^-1 norm used

  #
  # Petsc options: stress_
  #

  petsc_options+="\n\n###### Stress KSP ######\n"

  if [ "$PROBLEM_TYPE" = 'lin' ] && \
     [ 1 -le "${solve_stress_block_only[o]}" ]; then # if solve stress only
    petsc_options+="
  ### KSP type

    # CG
   #-stress_ksp_type cg
   #-stress_ksp_cg_single_reduction
   #-stress_ksp_pc_side left

    # GMRES
    -stress_ksp_type gmres
    -stress_ksp_gmres_restart $gmres_restart
    -stress_ksp_pc_side right

  ### KSP tolerances
  -stress_ksp_rtol $krylov_rtol
  -stress_ksp_atol 0.0
  -stress_ksp_max_it $krylov_maxiter

  ### KSP monitoring
  -stress_ksp_view
  -stress_ksp_monitor
 #-stress_ksp_monitor_range
 #-stress_ksp_monitor_singular_value
"
  else
    petsc_options+="
  ### KSP type
  -stress_ksp_type preonly

  ### KSP tolerances
  -stress_ksp_max_it 1
"
  fi # end if solve stress block only

  petsc_options+="\n\n###### Stress PC ######\n"

  if [ 1 -le "$stress_use_gmg" ]; then # if use GMG for viscous stress
    stress_gmg_prefix="stress_gmg"
    stress_amg_prefix="${stress_gmg_prefix}_coarse"
    petsc_options+="
  # inner preconditioner uses the matrix that defines the linear system
  -stress_pc_use_amat true

  ### GMG smoother
"
    petsc_options+="$(petsc_options_set_mg_smoother \
                          "${stress_gmg_prefix}_levels" \
                          "$stress_gmg_n_smoothing" \
                          "Chebyshev" "est_eig" "" "$stress_eig_est_max_it" \
                          "$chebyshev_est_eig_monitor" \
                          "PCMAT" "" \
                          "true")\n"

    petsc_options+="
  ### GMG coarse KSP

    ### KSP type (use AMG preconditioner only)
    -${stress_amg_prefix}_ksp_type preonly

    ### KSP tolerances
    -${stress_amg_prefix}_ksp_max_it 1

  ### GMG coarse PC
"
  else # if use AMG only for viscous stress
    stress_amg_prefix='stress'
  fi # end if use GMG for viscous stress

  petsc_options+="
  # inner preconditioner uses the matrix that defines the linear system
  -${stress_amg_prefix}_pc_use_amat false
"
  petsc_options+="$(petsc_options_set_amg \
                        $stress_amg_prefix \
                        $amg_type \
                        $amg_max_levels \
                        $stress_drop_tol \
                        $stress_agg_nsmooths \
                        $amg_stress_dof_per_proc \
                        $amg_stress_dof_coarse \
                        $stress_vcycles)\n"

  petsc_options+="
  ### AMG smoother
"
  petsc_options+="$(petsc_options_set_mg_smoother \
                        "${stress_amg_prefix}_mg_levels" \
                        "$stress_amg_n_smoothing" \
                        "Chebyshev" "" "" "$stress_eig_est_max_it" \
                        "$chebyshev_est_eig_monitor" \
                        "$stress_smoother" "" \
                        "false")\n"

  petsc_options+="
  ### AMG coarse solver
"
  petsc_options+="$(petsc_options_set_mg_coarse \
                        "${stress_amg_prefix}_mg_coarse" \
                        "1" "preonly" "Cholesky")\n"

  #
  # Petsc options: stokes_ and schur_
  #

  if [ "$PROBLEM_TYPE" = 'nl' ] || \
     [ 0 -ge "${solve_stress_block_only[o]}" ]; then # if Stokes solve

    petsc_options+="\n\n###### Stokes KSP ######\n"
    petsc_options+="
  ### KSP type

    # GMRES
    -stokes_ksp_type gmres
    -stokes_ksp_gmres_restart $gmres_restart
   #-stokes_ksp_gmres_haptol 1e-30
    -stokes_ksp_gmres_classicalgramschmidt  # fast
   #-stokes_ksp_gmres_modifiedgramschmidt   # more stable, but slower
    -stokes_ksp_gmres_cgs_refinement_type refine_always  # more accurate G-S
    -stokes_ksp_pc_side right

    # FGMRES (use if inner Krylov cycles are used anywhere)
   #-stokes_ksp_type fgmres
   #-stokes_ksp_pc_side right

    # MINRES
   #-stokes_ksp_type minres

  ### KSP tolerances (if commented out, then ymir options are passed)
 #-stokes_ksp_rtol 1.0e-16
 #-stokes_ksp_atol 0.0
 #-stokes_ksp_max_it 1000

  ### KSP monitoring
  -stokes_ksp_view
  -stokes_ksp_monitor
 #-stokes_ksp_monitor_true_residual
 #-stokes_ksp_monitor_range
 #-stokes_ksp_monitor_singular_value
 #-stokes_ksp_monitor_eigenvalues
"

    petsc_options+="\n\n###### Stokes PC ######\n"
    petsc_options+="  # handled by ymir\n"

    petsc_options+="\n\n###### Schur KSP ######\n"
    petsc_options+="
  ### KSP type
  -schur_ksp_type preonly

  ### KSP tolerances
  -schur_ksp_max_it 1
  -schur_ksp_norm_type none
"

    petsc_options+="\n\n###### Schur PC ######\n"
    petsc_options+="  # handled by ymir\n"

    petsc_options+="\n\n###### Stress KSP inside Schur ######\n"
    petsc_options+="
  ### KSP type
  -schur_stress_ksp_type preonly

  ### KSP tolerances
  -schur_stress_ksp_max_it 1
  -schur_stress_ksp_norm_type none
"

    petsc_options+="\n\n###### Stress PC inside Schur ######\n"
    petsc_options+="  # handled by ymir\n"
  fi # end if Stokes solve

  #
  # Petsc options: bbt_
  #

  if [ 1 -le "${schur_type[s]}" ]; then
    # if using BFBT pressure Schur preconditioner
    petsc_options+="\n\n###### BB^T KSP ######\n"

    if [ "$PROBLEM_TYPE" = 'lin' ] && \
       [ 1 -le "${solve_press_bbt_only[o]}" ]; then # if solve BB^T only
      petsc_options+="
  ### KSP type

    # CG
   #-bbt_ksp_type cg
   #-bbt_ksp_pc_side left

    # GMRES
    -bbt_ksp_type gmres
    -bbt_ksp_gmres_restart $gmres_restart
    -bbt_ksp_pc_side right

  ### KSP tolerances
  -bbt_ksp_rtol $krylov_rtol
  -bbt_ksp_atol 0.0
  -bbt_ksp_max_it $krylov_maxiter

  ### KSP monitoring
  -bbt_ksp_view
  -bbt_ksp_monitor
 #-bbt_ksp_monitor_range
 #-bbt_ksp_monitor_singular_value
"
    else
      petsc_options+="
  ### KSP type
  -bbt_ksp_type preonly

  ### KSP tolerances
  -bbt_ksp_max_it 1
"
    fi # end if solve BB^T only

    petsc_options+="\n\n###### BB^T PC ######\n"

    if [ 1 -le "$bbt_use_stiff_pc" ]; then # if use Stiffness PC for BB^T
      bbt_smooth_prefix='bbt_smooth'
      bbt_gmg_prefix="bbt_stiff_gmg"
      bbt_amg_prefix="${bbt_gmg_prefix}_coarse"

      petsc_options+="
  # inner preconditioner uses the matrix that defines the linear system
  -bbt_pc_use_amat true

  ### BB^T smoother
"
      petsc_options+="$(petsc_options_set_mg_smoother \
                            "$bbt_smooth_prefix" \
                            "$bbt_stiff_pc_modal_space_smoother_iter" \
                            "Chebyshev" "est_eig" "" "$bbt_eig_est_max_it" \
                            "$chebyshev_est_eig_monitor" \
                            "PCMAT" "" \
                            "true")\n"

      petsc_options+="\n\n###### BB^T stiffness KSP ######\n"

      if [ "$PROBLEM_TYPE" = 'lin' ] && \
         [ 1 -le "${solve_cnode_bbt_only[o]}" ]; then # if solve cnode BB^T only
        petsc_options+="
  ### KSP type
  -bbt_stiff_ksp_type gmres
  -bbt_stiff_ksp_gmres_restart $gmres_restart
  -bbt_stiff_ksp_pc_side right

  ### KSP tolerances
  -bbt_stiff_ksp_rtol $krylov_rtol
  -bbt_stiff_ksp_atol 0.0
  -bbt_stiff_ksp_max_it $krylov_maxiter

  ### KSP monitoring
  -bbt_stiff_ksp_view
  -bbt_stiff_ksp_monitor
 #-bbt_stiff_ksp_monitor_range
 #-bbt_stiff_ksp_monitor_singular_value
"
      else
        petsc_options+="
  ### KSP type
  -bbt_stiff_ksp_type preonly

  ### KSP tolerances
  -bbt_stiff_ksp_max_it 1
"
      fi # end if solve cnode BB^T only

      petsc_options+="\n\n###### BB^T stiffness PC ######\n"
      petsc_options+="
  # inner preconditioner uses the matrix that defines the linear system
  -bbt_stiff_pc_use_amat true

  ### GMG smoother
"
      petsc_options+="$(petsc_options_set_mg_smoother \
                            "${bbt_gmg_prefix}_levels" \
                            "$bbt_gmg_n_smoothing" \
                            "Chebyshev" "est_eig" "" "$bbt_eig_est_max_it" \
                            "$chebyshev_est_eig_monitor" \
                            "PCMAT" "" \
                            "true")\n"

      petsc_options+="
  ### GMG coarse KSP

    # AMG preconditioner
    -${bbt_amg_prefix}_ksp_type preonly

    ### KSP tolerances
    -${bbt_amg_prefix}_ksp_max_it 1

  ### GMG coarse PC
"
    elif [ 1 -le "$bbt_use_gmg" ]; then # if use GMG for BB^T
      bbt_gmg_prefix="bbt_gmg"
      bbt_amg_prefix="${bbt_gmg_prefix}_coarse"

      petsc_options+="
  # inner preconditioner uses the matrix that defines the linear system
  -bbt_pc_use_amat true

  ### GMG smoother
"
      petsc_options+="$(petsc_options_set_mg_smoother \
                            "${bbt_gmg_prefix}_levels" \
                            "$bbt_gmg_n_smoothing" \
                            "Chebyshev" "est_eig" "" "$bbt_eig_est_max_it" \
                            "$chebyshev_est_eig_monitor" \
                            "PCMAT" "" \
                            "true")\n"

      petsc_options+="
  ### GMG coarse KSP

    ### KSP type (use AMG preconditioner only)
    -${bbt_amg_prefix}_ksp_type preonly

    ### KSP tolerances
    -${bbt_amg_prefix}_ksp_max_it 1

  ### GMG coarse PC
"
    else # if use AMG for BB^T
      bbt_amg_prefix='bbt'
    fi # end if use Stiffness/GMG/AMG for BB^T

    petsc_options+="
  # inner preconditioner uses the matrix that defines the linear system
  -${bbt_amg_prefix}_pc_use_amat false
"
    petsc_options+="$(petsc_options_set_amg \
                          $bbt_amg_prefix \
                          $amg_type \
                          $amg_max_levels \
                          $bbt_drop_tol \
                          $bbt_agg_nsmooths \
                          $amg_stress_dof_per_proc \
                          $amg_stress_dof_coarse \
                          $bbt_vcycles)\n"

    petsc_options+="
  ### AMG smoother
"
    petsc_options+="$(petsc_options_set_mg_smoother \
                          "${bbt_amg_prefix}_mg_levels" \
                          "$bbt_smoothing_iter" \
                          "Chebyshev" "" "" "$bbt_eig_est_max_it" \
                          "$chebyshev_est_eig_monitor" \
                          "$bbt_smoother" "" \
                          "false")\n"

    petsc_options+="
  ### AMG coarse solver
"
   #if [ 1 -le "$bbt_use_stiff_pc" ] && \
   #   (( $(echo "0 < ($stiff_gmg_coarse_force_mass_scaling_dec)" | bc -l) ));
   #then # if BB^T matrix invertible
   #  bbt_amg_coarse_shift=''
   #else # if BB^T matrix (near) singular
   #  bbt_amg_coarse_shift='positive_definite'
   #fi
    bbt_amg_coarse_shift='positive_definite'
    petsc_options+="$(petsc_options_set_mg_coarse \
                          "${bbt_amg_prefix}_mg_coarse" \
                          "1" "preonly" "Cholesky" $bbt_amg_coarse_shift)\n"
  fi # end if using BFBT pressure Schur preconditioner

  # write Petsc options to file
  printf "${petsc_options}\n" 1>$petsc_options_filepath

  #
  # Create job scripts
  #

  job_script="${path_output}${dir_name}/${base_name}.sh"

  # copy template job script
  cp "$path_job_template" "$job_script"

  # replace placeholder strings in template
  sed -i "s/__JOB_NAME__/${base_name}/g" "$job_script"
  sed -i "s/__QUEUE__/${queue}/g" "$job_script"
  sed -i "s/__NUM_NODES__/${num_nodes}/g" "$job_script"
  sed -i "s/__MPISIZE__/${mpisize}/g" "$job_script"
  sed -i "s/__OMPSIZE__/${ompsize}/g" "$job_script"
  sed -i "s/__RUNTIME__/${max_runtime}/g" "$job_script"
  sed -i "s/__PROJECT__/${project_name}/g" "$job_script"
  sed -i "s/__EMAIL__/${email}/g" "$job_script"
  sed -i "s/__HOME_DIR__/${HOME//\//\\/}/g" "$job_script"
  sed -i "s/__SOURCE_PATH__/${path_code//\//\\/}/g" "$job_script"
  sed -i "s/__EXECUTABLE_PATH__/${path_executable//\//\\/}/g" "$job_script"
  if [ -n "$dir_options_override" ]; then # if override options directory
    ymir_options_filepath="$(basename $ymir_options_filepath)"
    ymir_options_filepath="${dir_options_override}${ymir_options_filepath}"
  fi
  sed -i "s/__YMIR_OPTIONS__/${ymir_options_filepath//\//\\/}/g" "$job_script"

  # create visualization script
  if [ 1 -le "$vtk_out" ]; then
    vis_script="${path_output}${dir_name}/vis_${base_name}.sh"
    render_script="${path_output}${dir_name}/$(basename $path_render_script)"
    render_other="${path_output}${dir_name}/$(basename $path_render_other)"
    vis_in_dir="${path_output}${dir_name}/${dir_output_vis}/"
    vis_out_dir="${path_output}${dir_name}/${dir_output_render}/"

    sed -i "s/__VIS_JOB_PATH__/${vis_script//\//\\/}/g" "$job_script"
    sed -i "s/__VIS_RENDER_SCRIPT__/${render_script//\//\\/}/g" "$job_script"

    nl_vis_filepath="${vis_in_dir}${base_name}_nl_itn"
    if [ 1 -le "$launch_vis_jobs" ] && [ "$PROBLEM_TYPE" = 'nl' ] && \
       [ 1 -le "$vtk_out_nl_iter" ]; then
      sed -i "s/__CONCURRENT_VIS__/1/g" "$job_script"
    else
      sed -i "s/__CONCURRENT_VIS__/0/g" "$job_script"
    fi
    sed -i "s/__CONCURRENT_VIS_FILEPATH__/${nl_vis_filepath//\//\\/}/g" \
           "$job_script"
    sed -i "s/__CONCURRENT_VIS_OUT_DIR__/${vis_out_dir//\//\\/}/g" \
           "$job_script"

    in_vis_filepath="${vis_in_dir}${base_name}_input.pvtu"
    if [ 1 -le "$launch_vis_jobs" ] && [ 1 -le "$vtk_out_input" ]; then
      sed -i "s/__INPUT_VIS__/1/g" "$job_script"
    else
      sed -i "s/__INPUT_VIS__/0/g" "$job_script"
    fi
    sed -i "s/__INPUT_VIS_FILEPATH__/${in_vis_filepath//\//\\/}/g" \
           "$job_script"
    sed -i "s/__INPUT_VIS_OUT_DIR__/${vis_out_dir//\//\\/}/g" \
           "$job_script"

    sol_vis_filepath="${vis_in_dir}${base_name}_solution.pvtu"
    if [ 1 -le "$launch_vis_jobs" ]; then
      sed -i "s/__SOLUTION_VIS__/1/g" "$job_script"
    else
      sed -i "s/__SOLUTION_VIS__/0/g" "$job_script"
    fi
    sed -i "s/__SOLUTION_VIS_FILEPATH__/${sol_vis_filepath//\//\\/}/g" \
           "$job_script"
    sed -i "s/__SOLUTION_VIS_OUT_DIR__/${vis_out_dir//\//\\/}/g" \
           "$job_script"

    cp "$path_vis_template" "$vis_script"
    cp "$path_render_script" "$render_script"
    cp "$path_render_other" "$render_other"

    sed -i "s/__JOB_NAME__/vis_${base_name}/g" "$vis_script"
    sed -i "s/__NUM_NODES__/$(echo "$vis_mpisize / 16" | bc)/g" "$vis_script"
    sed -i "s/__MPISIZE__/${vis_mpisize}/g" "$vis_script"
    sed -i "s/__PROJECT__/${project_name}/g" "$vis_script"
    sed -i "s/__EMAIL__/${email}/g" "$vis_script"
    sed -i "s/__VIS_JOB_PATH__/${vis_script//\//\\/}/g" "$vis_script"
  fi

  # set profiling
  if [ 1 -le "$profile_hpctoolkit" ]; then # if profile with HPCToolkit
    sed -i "s/__PROFILE_HPCTOOLKIT__/1/g" "$job_script"
  elif [ 1 -le "$profile_tau_lib" ]; then # if profile with TAU
    sed -i "s/__PROFILE_TAU_LIB__/1/g" "$job_script"
  elif [ 1 -le "$profile_tau_source" ]; then # if profile with TAU (PDT)
    sed -i "s/__PROFILE_TAU_SOURCE__/1/g" "$job_script"
  elif [ 1 -le "$profile_ipm" ]; then # if profile with IPM
    sed -i "s/__PROFILE_IPM__/1/g" "$job_script"
  elif [ 1 -le "$profile_massif" ]; then # if profile with Valgrind Massif
    sed -i "s/__PROFILE_MASSIF__/1/g" "$job_script"
  fi

  # add this simulation to script that launches all simulations
  printf "echo 'Launch simulation ${base_name}'\n\n" 1>>"$script_launch_all"
  printf "cd ${path_output}${dir_name}\n\n" 1>>"$script_launch_all"
  printf "$job_script > ${path_output}${dir_name}/${base_name}.out\n\n" \
  1>>"$script_launch_all"
  printf "echo '========================================'\n\n" \
  1>>"$script_launch_all"

  # add this simulation to logfile
  echo "${dir_name}/${base_name}" 1>>"$logfile"

  # submit job launching only this simulation
  if [ "$dry_run" -le 0 ] && [ 1 -le "$launch_multiple_jobs" ]; then
    echo "Submit job: $base_name"
    cd "${path_output}${dir_name}"
    sbatch "$job_script"
  fi

done
done
done
done
done

done
done
done
done
done

done
done
done
done
done

done
done
done
done
done

done
done
done
done
done

done
done
done
done
done

done

# submit script that launches all jobs
if [ "$dry_run" -le 0 ] && [ "$launch_multiple_jobs" -le 0 ]; then
  echo "Submit one job for all simulations"
  cd "${path_output}"
  sbatch "$script_launch_all"
fi

# move logfile of this script
mv "$logfile" "$path_output"

# copy this script
cp "${script_path}$(basename "$0")" "$path_output"

echo "Finished: $script_name"

