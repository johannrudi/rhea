/*
  This file is part of the ymir Library.
  ymir is a C library for modeling ice sheets

  Copyright (C) 2015 Carsten Burstedde, Toby Isaac, Johann Rudi, Georg Stadler,
                     Lucas Wilcox.

  The ymir Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The ymir Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the ymir Library.  If not, see <http://www.gnu.org/licenses/>.

  ---

  This example runs perfomace tests.

*/

#include <slabs_setup.h>
#include <slabs_discretization_extended.h>
#include <ymir_gmg.h>
#include <ymir_perf_counter.h>
#include <ymir_monitor.h>
#include <p8est_extended.h>
#ifdef YMIR_DEBUG
# include <ymir_vtk.h>
#endif

#if defined(__bgq__)
#include <ymir_bgq.h>
#endif

/* default physics options TODO move down */
#define SL_DOMAIN_SHAPE "shell"
#define SL_BRICK_DX (1)
#define SL_BRICK_DY (16)
#define SL_BRICK_DZ (16)

#define SL_WEAKZONE_TYPE_NAME "NONE"
#define SL_WEAKZONE_IMPORT_THICKNESS (20.0e3)
#define SL_WEAKZONE_IMPORT_THICKNESS_CONST (5.0e3)
#define SL_WEAKZONE_IMPORT_WEAK_FACTOR (1.0e-5)
#define SL_WEAKZONE_2PLATES_SUBDU_LONGITUDE (-100.0)      /* rhea1:  0.13 */
#define SL_WEAKZONE_2PLATES_SUBDU_DIP_ANGLE (-1.0)        /* rhea1:   N/A */
#define SL_WEAKZONE_2PLATES_SUBDU_DEPTH (80.0e3)          /* rhea1: 50 km */
#define SL_WEAKZONE_2PLATES_SUBDU_WIDTH (-1.0)            /* rhea1:   N/A */
#define SL_WEAKZONE_2PLATES_SUBDU_THICKNESS (20.0e3)      /* rhea1: 20 km */
#define SL_WEAKZONE_2PLATES_SUBDU_THICKNESS_CONST (5.0e3) /* rhea1: 10 km */
#define SL_WEAKZONE_2PLATES_SUBDU_WEAK_FACTOR (1.0e-5)    /* rhea1:  1e-5 */
#define SL_WEAKZONE_2PLATES_RIDGE_DEPTH (30.0e3)          /* rhea1: 30 km */
#define SL_WEAKZONE_2PLATES_RIDGE_WIDTH (30.0e3)          /* rhea1: 30 km */
#define SL_WEAKZONE_2PLATES_RIDGE_SMOOTHWIDTH (10.0e3)    /* rhea1:  5 km */
#define SL_WEAKZONE_2PLATES_RIDGE_WEAK_FACTOR (1.0e-5)    /* rhea1:  1e-5 */

#define SL_VISCOSITY_MODEL_TYPE "UWYL"
#define SL_VISCOSITY_IIE_REG (1.0e-10)
#define SL_VISCOSITY_MIN (-1.0)
#define SL_VISCOSITY_MAX (-1.0)
#define SL_VISCOSITY_SCALING (1.0)
#define SL_VISCOSITY_TEMP_DECAY (7.0)
#define SL_VISCOSITY_LOWER_MANTLE (0)
#define SL_VISCOSITY_LOWER_MANTLE_SCALING (1.0)
#define SL_VISCOSITY_LOWER_MANTLE_TEMP_DECAY (-1.0)
#define SL_VISCOSITY_STRESS_EXPONENT (3.0)
#define SL_VISCOSITY_STRESS_YIELD (-1.0)

#define SL_RHS_SCALING (1.0)

#define SL_CUBE_INNER_PLUME_CENTER_X (0.5)
#define SL_CUBE_INNER_PLUME_CENTER_Y (0.5)
#define SL_CUBE_INNER_PLUME_CENTER_Z (0.5)
#define SL_CUBE_INNER_PLUME_DECAY (10.0)
#define SL_CUBE_INNER_PLUME_SCALING (-1.0)

#define SL_BRICK_INNER_PLUME_CENTER_X (0.5)
#define SL_BRICK_INNER_PLUME_CENTER_Y (0.5)
#define SL_BRICK_INNER_PLUME_CENTER_Z (0.4)
#define SL_BRICK_INNER_PLUME_DECAY (0.2)
#define SL_BRICK_INNER_PLUME_SCALING (-0.5)

#define SL_SHELL_INNER_PLUME_CENTER_X (0.0)
#define SL_SHELL_INNER_PLUME_CENTER_Y (0.0)
#define SL_SHELL_INNER_PLUME_CENTER_Z (0.5 * (SL_SHELL_RADIUS_BOTTOM + \
                                              SL_SHELL_RADIUS_TOP))
#define SL_SHELL_INNER_PLUME_DECAY (5.0)
#define SL_SHELL_INNER_PLUME_SCALING (1.0)

#define SL_SHELL_RISING_PLUME_CENTER_X (0.0)
#define SL_SHELL_RISING_PLUME_CENTER_Y (0.0)
#define SL_SHELL_RISING_PLUME_CENTER_Z (0.9 * SL_SHELL_RADIUS_BOTTOM - \
                                        0.1 * SL_SHELL_RADIUS_TOP)
#define SL_SHELL_RISING_PLUME_DECAY (5.0)
#define SL_SHELL_RISING_PLUME_SCALING (7.6818e4)

#define SL_SHELL_SLICE_INNER_PLUME_CENTER_X (0.0)
#define SL_SHELL_SLICE_INNER_PLUME_CENTER_Y (0.0)
#define SL_SHELL_SLICE_INNER_PLUME_CENTER_Z (0.8 * SL_SHELL_RADIUS_BOTTOM + \
                                             0.2 * SL_SHELL_RADIUS_TOP)
#define SL_SHELL_SLICE_INNER_PLUME_DECAY (30.0)
#define SL_SHELL_SLICE_INNER_PLUME_SCALING (-0.5)

/* default Krylov solver options TODO move down */
#define SL_KRYLOV_MAXITER (300)
#define SL_KRYLOV_ATOL (0.0)
#define SL_KRYLOV_RTOL (1.0e-1)
#define SL_KRYLOV_GMRES_NUM_VECS (50)

/* initialize physics options: general */

#define SLABS_DEFAULT_P4EST_IMPORT_FILENAME NULL

char               *domain_shape_name = SL_DOMAIN_SHAPE;
int                 brick_dx = SL_BRICK_DX;
int                 brick_dy = SL_BRICK_DY;
int                 brick_dz = SL_BRICK_DZ;
int                 bc_type = SL_VEL_BC_DIRICHLET_ALL;
int                 bc_default_dir_scale = 0; //TODO

char               *slabs_p4est_import_filename =
  SLABS_DEFAULT_P4EST_IMPORT_FILENAME;

/* initialize physics options: temperature */

#define SLABS_DEFAULT_TEMP_TYPE_NAME "NONE"
#define SLABS_DEFAULT_TEMP_BACKGROUND_PLATE_AGE (50.0e6)    /* rhea1: 60 Myr */
#define SLABS_DEFAULT_TEMP_IMPORT_PLATE_AGE_MIN (1.0e6)
#define SLABS_DEFAULT_TEMP_IMPORT_FILENAME_TXT NULL
#define SLABS_DEFAULT_TEMP_IMPORT_VERIFICATION_OUT NULL
#define SLABS_DEFAULT_TEMP_IMPORT_WRITE_COORD_PATH NULL
#define SLABS_DEFAULT_TEMP_2PL_TRENCH_LONGITUDE (0.13)      /* rhea1: 0.13   */
#define SLABS_DEFAULT_TEMP_2PL_DIP_ANGLE (5.0)              /* rhea1: --     */
#define SLABS_DEFAULT_TEMP_2PL_SUBD_DEPTH (400.0e3)         /* rhea1: 400 km */
#define SLABS_DEFAULT_TEMP_2PL_SUBD_WIDTH (300.0e3)         /* rhea1: --     */
#define SLABS_DEFAULT_TEMP_2PL_SUBD_EDGE_WIDTH (1.0e3)      /* rhea1: --     */
#define SLABS_DEFAULT_TEMP_2PL_SUBD_EDGE_SMOOTHWIDTH (40.0e3)/* rhea1: --     */
#define SLABS_DEFAULT_TEMP_2PL_SUBD_PLATE_VELOCITY (4.0e-2)  /* rhea1: 4 cm/y */
#define SLABS_DEFAULT_TEMP_2PL_SUBD_PLATE_INITIAL_AGE (1.0e6)/* rhea1: --     */
#define SLABS_DEFAULT_TEMP_2PL_OVER_PLATE_AGE (40.0e6)       /* rhea1: 40 Myr */

char               *temp_type_name = SLABS_DEFAULT_TEMP_TYPE_NAME;
double              temp_back_plate_age =
  SLABS_DEFAULT_TEMP_BACKGROUND_PLATE_AGE;
double              temp_import_plate_age_min =
  SLABS_DEFAULT_TEMP_IMPORT_PLATE_AGE_MIN;
char               *temp_import_filename_txt =
  SLABS_DEFAULT_TEMP_IMPORT_FILENAME_TXT;
char               *temp_import_verification_out =
  SLABS_DEFAULT_TEMP_IMPORT_VERIFICATION_OUT;
char               *temp_import_write_coord_path =
  SLABS_DEFAULT_TEMP_IMPORT_WRITE_COORD_PATH;
double              temp_2pl_trench_lon =
  SLABS_DEFAULT_TEMP_2PL_TRENCH_LONGITUDE;
double              temp_2pl_dip_angle = SLABS_DEFAULT_TEMP_2PL_DIP_ANGLE;
double              temp_2pl_subd_depth = SLABS_DEFAULT_TEMP_2PL_SUBD_DEPTH;
double              temp_2pl_subd_width = SLABS_DEFAULT_TEMP_2PL_SUBD_WIDTH;
double              temp_2pl_subd_edge_width =
  SLABS_DEFAULT_TEMP_2PL_SUBD_EDGE_WIDTH;
double              temp_2pl_subd_edge_smoothwidth =
  SLABS_DEFAULT_TEMP_2PL_SUBD_EDGE_SMOOTHWIDTH;
double              temp_2pl_subd_plate_vel =
  SLABS_DEFAULT_TEMP_2PL_SUBD_PLATE_VELOCITY;
double              temp_2pl_subd_plate_init_age =
  SLABS_DEFAULT_TEMP_2PL_SUBD_PLATE_INITIAL_AGE;
double              temp_2pl_over_plate_age =
  SLABS_DEFAULT_TEMP_2PL_OVER_PLATE_AGE;

/* initialize physics options: weak zone */

char               *weakzone_type_name = SL_WEAKZONE_TYPE_NAME;
char               *weakzone_import_filename_txt = NULL; //TODO
char               *weakzone_import_verification_out = NULL; //TODO
int                 weakzone_import_pointcloud_size = 0; //TODO
double              weakzone_import_thickness = SL_WEAKZONE_IMPORT_THICKNESS;
double              weakzone_import_thickness_const =
  SL_WEAKZONE_IMPORT_THICKNESS_CONST;
double              weakzone_import_weak_factor =
  SL_WEAKZONE_IMPORT_WEAK_FACTOR;
double              weakzone_2pl_subdu_lon =
  SL_WEAKZONE_2PLATES_SUBDU_LONGITUDE;
double              weakzone_2pl_subdu_dip_angle =
  SL_WEAKZONE_2PLATES_SUBDU_DIP_ANGLE;
double              weakzone_2pl_subdu_depth = SL_WEAKZONE_2PLATES_SUBDU_DEPTH;
double              weakzone_2pl_subdu_width = SL_WEAKZONE_2PLATES_SUBDU_WIDTH;
double              weakzone_2pl_subdu_thickness =
  SL_WEAKZONE_2PLATES_SUBDU_THICKNESS;
double              weakzone_2pl_subdu_thickness_const =
  SL_WEAKZONE_2PLATES_SUBDU_THICKNESS_CONST;
double              weakzone_2pl_subdu_weak_factor =
  SL_WEAKZONE_2PLATES_SUBDU_WEAK_FACTOR;
double              weakzone_2pl_ridge_depth = SL_WEAKZONE_2PLATES_RIDGE_DEPTH;
double              weakzone_2pl_ridge_width = SL_WEAKZONE_2PLATES_RIDGE_WIDTH;
double              weakzone_2pl_ridge_smoothwidth =
  SL_WEAKZONE_2PLATES_RIDGE_SMOOTHWIDTH;
double              weakzone_2pl_ridge_weak_factor =
  SL_WEAKZONE_2PLATES_RIDGE_WEAK_FACTOR;
int                 weakzone_2pl_ridge_mesh_align = 0; //TODO

/* initialize physics options: viscosity */

#define SLABS_DEFAULT_VELOCITY_IMPORT_FILENAME NULL
#define SLABS_DEFAULT_PRESSURE_IMPORT_FILENAME NULL

char               *slabs_velocity_import_filename =
  SLABS_DEFAULT_VELOCITY_IMPORT_FILENAME;
char               *slabs_pressure_import_filename =
  SLABS_DEFAULT_PRESSURE_IMPORT_FILENAME;

int                 viscosity_type = SL_VISCOSITY_CONST;
int                 viscosity_type_for_init_nl_stokes =
  SL_VISCOSITY_INIT_NL_STOKES_DEFAULT;
char               *visc_model_type_name = SL_VISCOSITY_MODEL_TYPE;
double              visc_IIe_reg = SL_VISCOSITY_IIE_REG;
double              visc_min = SL_VISCOSITY_MIN;
double              visc_max = SL_VISCOSITY_MAX;
double              visc_temp_max = SL_VISCOSITY_MAX; //TODO
double              visc_scaling = SL_VISCOSITY_SCALING;
double              visc_temp_decay = SL_VISCOSITY_TEMP_DECAY;
int                 visc_lower_mantle = SL_VISCOSITY_LOWER_MANTLE;
double              visc_lower_mantle_scaling =
  SL_VISCOSITY_LOWER_MANTLE_SCALING;
double              visc_lower_mantle_temp_decay =
  SL_VISCOSITY_LOWER_MANTLE_TEMP_DECAY;
double              visc_stress_exp = SL_VISCOSITY_STRESS_EXPONENT;
double              visc_stress_yield = SL_VISCOSITY_STRESS_YIELD;
double              visc_yield_reg = 0.0;

int                 visc_coarsen_eval = 0; //TODO
int                 visc_p_coarsen_eval = 0; //TODO
int                 visc_coarsen_type = SL_VISCOSITY_COARSEN_EVAL;
int                 weakzone_coarsen_type = SL_WEAKZONE_COARSEN_EVAL;
double              visc_lower_upper_transition_zone_incr = 0.0; //TODO

/* initialize physics options: right-hand side */

double              rhs_scaling = SL_RHS_SCALING;
int                 rhs_random = 0; //TODO
int                 rhs_multiply_in_weak_zone = 0; //TODO
int                 plume_type = SL_PLUME_NONE; //TODO
double              plume_center_x = 0.0; //TODO
double              plume_center_y = 0.0; //TODO
double              plume_center_z = 0.0; //TODO
double              plume_decay = -1.0; //TODO
double              plume_scaling = 0.0; //TODO

/* initialize discretization options */

#define SLABS_DEFAULT_DISCR_ORDER (2)
#define SLABS_DEFAULT_DISCR_REFINEMENT_MINLEVEL (0)
#define SLABS_DEFAULT_DISCR_REFINEMENT_MAXLEVEL (15)
#define SLABS_DEFAULT_DISCR_REFINEMENT_TYPE "uniform"
#define SLABS_DEFAULT_DISCR_IMPORT_ORDER (0)
#define SLABS_DEFAULT_DISCR_IMPORT_MINLEVEL (0)

int                 slabs_order = ymir_n (SLABS_DEFAULT_DISCR_ORDER);
int                 slabs_minlevel = SLABS_DEFAULT_DISCR_REFINEMENT_MINLEVEL;
int                 slabs_maxlevel = SLABS_DEFAULT_DISCR_REFINEMENT_MAXLEVEL;
char               *slabs_refine = SLABS_DEFAULT_DISCR_REFINEMENT_TYPE;
int                 slabs_import_order = SLABS_DEFAULT_DISCR_IMPORT_ORDER;
int                 slabs_import_minlevel = SLABS_DEFAULT_DISCR_IMPORT_MINLEVEL;

                    //TODO all below
char               *refine_radius_str = NULL;
double              refine_surface_maxdist = 0.0;
double              refine_surface_elem_res = 0.0;
double              refine_lm_um_interface_maxdist = 0.0;
double              refine_lm_um_interface_elem_res = 0.0;
int                 enforce_refinement_lm_um_interface = 0;

int                 init_amr_max_steps = 0;

double              init_amr_rel_threshold = 0.0;
double              init_amr_n_elements_max = 0.0;
int                 init_amr_override_order = 0;
int                 init_amr_lower_mantle = 0;

int                 init_amr_visc_indicator_type = SL_AMR_INDICATOR_NONE;
double              init_amr_visc_tol_min = 0.0;
double              init_amr_visc_tol_max = 0.0;
double              init_amr_visc_in_plates_tol_min = 0.0;
double              init_amr_visc_in_plates_tol_max = 0.0;

int                 init_amr_weak_subdu_indicator_type = SL_AMR_INDICATOR_NONE;
double              init_amr_weak_subdu_tol_min = 0.0;
double              init_amr_weak_subdu_tol_max = 0.0;
double              init_amr_weak_subdu_elem_res = 0.0;
int                 init_amr_weak_ridge_indicator_type = SL_AMR_INDICATOR_NONE;
double              init_amr_weak_ridge_tol_min = 0.0;
double              init_amr_weak_ridge_tol_max = 0.0;
double              init_amr_weak_ridge_elem_res = 0.0;

int                 init_amr_weak_import_indicator_type = SL_AMR_INDICATOR_NONE;
double              init_amr_weak_import_tol_min = 0.0;
double              init_amr_weak_import_tol_max = 0.0;
double              init_amr_weak_import_elem_res = 0.0;

int                 init_amr_rhs_indicator_type = SL_AMR_INDICATOR_NONE;
double              init_amr_rhs_tol_min = 0.0;
double              init_amr_rhs_tol_max = 0.0;
double              init_amr_rhs_norm_shift = 1.0;

int                 init_amr_post_uniform_n_steps = 0;

int                 amr_max_steps = 0;
double              amr_rel_threshold = 0.0;
double              amr_n_elements_max = 0.0;
int                 amr_lower_mantle = 0;

int                 amr_visc_indicator_type = SL_AMR_INDICATOR_NONE;
double              amr_visc_tol_min = 0.0;
double              amr_visc_tol_max = 0.0;

int                 amr_visc_dr_indicator_type = SL_AMR_INDICATOR_NONE;
double              amr_visc_dr_tol_min = 0.0;
double              amr_visc_dr_tol_max = 0.0;

int                 amr_strain_rate_indicator_type = SL_AMR_INDICATOR_NONE;
double              amr_strain_rate_tol_min = 0.0;
double              amr_strain_rate_tol_max = 0.0;

int                 amr_log_maxlevel = 0;
int                 mesh_partitioning_type = SL_MESH_PARTITIONING_ELEM;

/* initialize nonlinear solver options */

#define SLABS_DEFAULT_NL_SOLVER_TYPE (SL_NL_SOLVER_NONE)
#define SLABS_DEFAULT_NL_SOLVER_PRIMALDUAL_TYPE (SL_NL_SOLVER_PRIMALDUAL_NONE)
#define SLABS_DEFAULT_NL_SOLVER_PRIMALDUAL_SCAL_TYPE \
  (SL_NL_SOLVER_PRIMALDUAL_SCAL_NONE)
#define SLABS_DEFAULT_NL_SOLVER_MAXITER (40)
#define SLABS_DEFAULT_NL_SOLVER_RTOL (1.0e-6)
#define SLABS_DEFAULT_NL_SOLVER_INIT_GUESS_TYPE \
  (SL_NL_SOLVER_INITIAL_GUESS_ZERO)

#define SLABS_DEFAULT_NL_SOLVER_RESUME_AT_ITER (0)
#define SLABS_DEFAULT_NL_SOLVER_RESUME_AT_TIME (0.0)
#define SLABS_DEFAULT_NL_SOLVER_RESUME_PREV_RES (0.0)
#define SLABS_DEFAULT_NL_SOLVER_RESUME_INIT_RES (0.0)

#define SLABS_DEFAULT_NL_SOLVER_FORCING_EXPONENT (1.618) /* (1 + sqrt(5)) / 2 */
#define SLABS_DEFAULT_NL_SOLVER_FORCING_MAX (0.5)
#define SLABS_DEFAULT_NL_SOLVER_FORCING_MAX_PROGRESSIVE_ITER (0)
#define SLABS_DEFAULT_NL_SOLVER_FORCING_TOTAL_MIN (0.0)
#define SLABS_DEFAULT_NL_SOLVER_FORCING_SAVEGUARD (1)
#define SLABS_DEFAULT_NL_SOLVER_FORCING_SAVEGUARD_THRESHOLD (0.1)

#define SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_REDUCTION_TYPE \
  (SL_NL_SOLVER_STEP_REDUCTION_CONST)
#define SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_REDUCTION_MIN (0.2)
#define SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_REDUCTION_MAX (0.8)
#define SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_REDUCTION_REG (1.0e-10)
#define SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_DESCEND_COND_RELAX (1.0e-4)
#define SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_MIN (-1.0)

#define SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_STEP_LENGTH_MIN (-1.0)
#define SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_AFTER_AMR (0)
#define SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_INIT (0)
#define SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_MAXITER (3)
#define SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_RTOL (0.1)

#define SLABS_DEFAULT_NL_SOLVER_NORM_TYPE (SL_NORM_VEC_L2)
#define SLABS_DEFAULT_NL_SOLVER_NORM_HMINUS1_MASS_SCALING (0.0)

#define SLABS_DEFAULT_NL_SOLVER_SCHUR_DIAG_TYPE \
  (SL_NL_STOKES_PROB_SCHUR_DIAG_INV_VISC_PMASS)
#define SLABS_DEFAULT_NL_SOLVER_SCALING_TYPE \
  (SL_NL_STOKES_PROB_SCALING_NONE)
#define SLABS_DEFAULT_NL_SOLVER_PROJECT_NULLSPACE \
  (SL_NL_SOLVER_PROJECT_NULLSPACE_NONE)
#define SLABS_DEFAULT_NL_SOLVER_ENFORCE_UNSCALED_REDUCTION (0)

int                 nl_solver_type = SLABS_DEFAULT_NL_SOLVER_TYPE;
int                 nl_solver_primaldual_type =
  SLABS_DEFAULT_NL_SOLVER_PRIMALDUAL_TYPE;
int                 nl_solver_primaldual_scal_type =
  SLABS_DEFAULT_NL_SOLVER_PRIMALDUAL_SCAL_TYPE;
int                 nl_maxiter = SLABS_DEFAULT_NL_SOLVER_MAXITER;
double              nl_rtol = SLABS_DEFAULT_NL_SOLVER_RTOL;
int                 nl_initial_guess_type =
  SLABS_DEFAULT_NL_SOLVER_INIT_GUESS_TYPE;

int                 nl_resume_at_iter =
  SLABS_DEFAULT_NL_SOLVER_RESUME_AT_ITER;
double              nl_resume_at_time =
  SLABS_DEFAULT_NL_SOLVER_RESUME_AT_TIME;
double              nl_resume_prev_res =
  SLABS_DEFAULT_NL_SOLVER_RESUME_PREV_RES;
double              nl_resume_init_res =
  SLABS_DEFAULT_NL_SOLVER_RESUME_INIT_RES;

double              nl_forcing_exponent =
  SLABS_DEFAULT_NL_SOLVER_FORCING_EXPONENT;
double              nl_forcing_max =
  SLABS_DEFAULT_NL_SOLVER_FORCING_MAX;
int                 nl_forcing_max_progressive_iter =
  SLABS_DEFAULT_NL_SOLVER_FORCING_MAX_PROGRESSIVE_ITER;
double              nl_forcing_total_min =
  SLABS_DEFAULT_NL_SOLVER_FORCING_TOTAL_MIN;
int                 nl_forcing_saveguard =
  SLABS_DEFAULT_NL_SOLVER_FORCING_SAVEGUARD;
double              nl_forcing_threshold =
  SLABS_DEFAULT_NL_SOLVER_FORCING_SAVEGUARD_THRESHOLD;

int                 nl_step_length_reduction_type =
  SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_REDUCTION_TYPE;
double              nl_step_length_reduction_min =
  SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_REDUCTION_MIN;
double              nl_step_length_reduction_max =
  SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_REDUCTION_MAX;
double              nl_step_length_reduction_reg =
  SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_REDUCTION_REG;
double              nl_step_length_descend_relax =
  SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_DESCEND_COND_RELAX;
double              nl_step_length_min =
  SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_MIN;

double              nl_switch_picard_step_length_min =
  SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_STEP_LENGTH_MIN;
int                 nl_switch_picard_after_amr =
  SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_AFTER_AMR;
int                 nl_switch_picard_init =
  SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_INIT;
int                 nl_switch_picard_maxiter =
  SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_MAXITER;
double              nl_switch_picard_rtol =
  SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_RTOL;

int                 nl_norm_type =
  SLABS_DEFAULT_NL_SOLVER_NORM_TYPE;
double              nl_norm_Hminus1_mass_scaling =
  SLABS_DEFAULT_NL_SOLVER_NORM_HMINUS1_MASS_SCALING;

int                 nl_schur_diag_type =
  SLABS_DEFAULT_NL_SOLVER_SCHUR_DIAG_TYPE;
int                 nl_scaling_type =
  SLABS_DEFAULT_NL_SOLVER_SCALING_TYPE;
int                 nl_project_out_nullspace =
  SLABS_DEFAULT_NL_SOLVER_PROJECT_NULLSPACE;
int                 nl_enforce_unscaled_reduction =
  SLABS_DEFAULT_NL_SOLVER_ENFORCE_UNSCALED_REDUCTION;

                    //TODO all below
double              nl_grid_cont_init_amr_threshold = 0.0;
double              nl_grid_cont_init_forcing = 0.0;
double              nl_grid_cont_init_forcing_exp = 0.0;
int                 nl_grid_cont_init_steps = 0;
int                 nl_grid_cont_skipsteps = 0;
int                 nl_grid_cont_maxsteps = 0;
double              nl_visc_bounds_cont_min = -1.0;
double              nl_visc_bounds_cont_max = -1.0;
int                 nl_visc_bounds_cont_steps = 0;

int                 nl_check_derivative = 0;
int                 log_physics_stats = 0;

/* initialize linear solver options */

int                 krylov_type = SL_KRYLOV_GMRES; //TODO
int                 gmres_num_vecs = SL_KRYLOV_GMRES_NUM_VECS;
int                 krylov_maxiter = SL_KRYLOV_MAXITER;
double              krylov_rtol = SL_KRYLOV_RTOL;
double              krylov_atol = SL_KRYLOV_ATOL;

int                 lin_solve_stress_block_only = 0; //TODO
int                 lin_solve_press_bbt_only = 0; //TODO
int                 lin_solve_cnode_bbt_only = 0; //TODO

/**
 * Defines options and adds them as sub-options.
 */
void
slabs_setup_add_suboptions (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Mantle";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /* physics options */
  YMIR_OPTIONS_S, "domain-shape", '\0', &domain_shape_name, SL_DOMAIN_SHAPE,
    "Shape of domain: unitcube, brick, shell, shell_chunk, shell_slice",
  YMIR_OPTIONS_I, "domain-brick-dx", '\0', &brick_dx, SL_BRICK_DX,
    "Dimension of the brick domain in x",
  YMIR_OPTIONS_I, "domain-brick-dy", '\0', &brick_dy, SL_BRICK_DY,
    "Dimension of the brick domain in y",
  YMIR_OPTIONS_I, "domain-brick-dz", '\0', &brick_dz, SL_BRICK_DZ,
    "Dimension of the brick domain in z",
  YMIR_OPTIONS_I, "boundary-cond", '\0', &bc_type, SL_VEL_BC_DIRICHLET_ALL,
    "Boundary condition: 0: Dirichlet all; "
    "1: Dirichlet in normal direction; "
    "2: Dirichlet norm (if shell: two inner points fixed for rotation inv); "
    "3: Dirichlet all on inner sphere & Dirichlet norm on outer sphere; "
    "4: Dirichlet all on side faces & Neumann on top, bottom faces; "
    "5: Dirichlet norm on side faces & Dirichlet all on top, bottom faces",
  YMIR_OPTIONS_I, "boundary-set-default-dirichlet-scale", '\0', &bc_default_dir_scale, 0,
    "Scale the heterogeneous Dirichlet BC component",

  YMIR_OPTIONS_S, "p4est-import", '\0',
    &slabs_p4est_import_filename, NULL,
    "Path to file for importing p4est mesh",

  YMIR_OPTIONS_S, "temperature", '\0',
    &temp_type_name, SLABS_DEFAULT_TEMP_TYPE_NAME,
    "Type of temperature: path to binary file to import from, or: "
    "hot_core, cold_plate, hot_core_cold_plate, 2plates",
  YMIR_OPTIONS_D, "temp-background-plate-age", '\0',
    &temp_back_plate_age, SLABS_DEFAULT_TEMP_BACKGROUND_PLATE_AGE,
    "Bachground temperature descibed by plate age [yr]",
  YMIR_OPTIONS_D, "temp-import-plate-age-min", '\0',
    &temp_import_plate_age_min, SLABS_DEFAULT_TEMP_IMPORT_PLATE_AGE_MIN,
    "Min age of plates to avoid overly sharp gradients at ridges [y]",
  YMIR_OPTIONS_S, "temp-import-read-textfile", '\0',
    &temp_import_filename_txt, SLABS_DEFAULT_TEMP_IMPORT_FILENAME_TXT,
    "Filename of a text file that contains temperature field "
    "(binary filename needs to be provided for option 'temperature')",
  YMIR_OPTIONS_S, "temp-import-verification-out", '\0',
    &temp_import_verification_out, SLABS_DEFAULT_TEMP_IMPORT_VERIFICATION_OUT,
    "Filename of output for verifying imported temperature data",
  YMIR_OPTIONS_S, "temp-import-write-coordinates", '\0',
    &temp_import_write_coord_path, SLABS_DEFAULT_TEMP_IMPORT_WRITE_COORD_PATH,
    "Path to a text file for output of node coordinates",
  YMIR_OPTIONS_D, "temp-2plates-trench-longitude", '\0',
    &temp_2pl_trench_lon, SLABS_DEFAULT_TEMP_2PL_TRENCH_LONGITUDE,
    "2plates temp: Longitude of trench in interval (-pi/8, pi/8)",
  YMIR_OPTIONS_D, "temp-2plates-dip-angle", '\0',
    &temp_2pl_dip_angle, SLABS_DEFAULT_TEMP_2PL_DIP_ANGLE,
    "2plates temp: Dip angle of subducting plate (in degrees < 0)",
  YMIR_OPTIONS_D, "temp-2plates-subd-depth", '\0',
    &temp_2pl_subd_depth, SLABS_DEFAULT_TEMP_2PL_SUBD_DEPTH,
    "2plates temp: Maximal depth of subducting plate inside of mantle [m]",
  YMIR_OPTIONS_D, "temp-2plates-subd-width", '\0',
    &temp_2pl_subd_width, SLABS_DEFAULT_TEMP_2PL_SUBD_WIDTH,
    "2plates temp: Maximal width of subducting zone inside of mantle [m]",
  YMIR_OPTIONS_D, "temp-2plates-subd-edge-width", '\0',
    &temp_2pl_subd_edge_width, SLABS_DEFAULT_TEMP_2PL_SUBD_EDGE_WIDTH,
    "2plates temp: Width for subducting plate's top edge [m]",
  YMIR_OPTIONS_D, "temp-2plates-subd-edge-smoothwidth", '\0',
    &temp_2pl_subd_edge_smoothwidth,
    SLABS_DEFAULT_TEMP_2PL_SUBD_EDGE_SMOOTHWIDTH,
    "2plates weak zone: Width of smoothing of subd. plate's top edge [m]",
  YMIR_OPTIONS_D, "temp-2plates-subd-plate-velocity", '\0',
    &temp_2pl_subd_plate_vel, SLABS_DEFAULT_TEMP_2PL_SUBD_PLATE_VELOCITY,
    "2plates temp: Velocity of subducting plate [m/y]",
  YMIR_OPTIONS_D, "temp-2plates-subd-plate-initial-age", '\0',
    &temp_2pl_subd_plate_init_age,
    SLABS_DEFAULT_TEMP_2PL_SUBD_PLATE_INITIAL_AGE,
    "2plates temp: Age of subducting plate at left boundary [y]",
  YMIR_OPTIONS_D, "temp-2plates-over-plate-age", '\0',
    &temp_2pl_over_plate_age, SLABS_DEFAULT_TEMP_2PL_OVER_PLATE_AGE,
    "2plates temp: Age of overriding plate [y]",

  YMIR_OPTIONS_S, "weakzone", '\0',
    &weakzone_type_name, SL_WEAKZONE_TYPE_NAME,
    "Type of weak fault zone: path to binary file to import from, or: "
    "2plates_poly2, 2plates_rhea1",
  YMIR_OPTIONS_S, "weakzone-import-read-textfile", '\0',
    &weakzone_import_filename_txt, NULL,
    "Filename of a text file that contains weak zone point cloud data "
    "(binary filename needs to be provided for option 'weakzone')",
  YMIR_OPTIONS_S, "weakzone-import-verification-out", '\0',
    &weakzone_import_verification_out, NULL,
    "Filename of output for verifying imported weak zone data",
  YMIR_OPTIONS_I, "weakzone-import-pointcloud-size", '\0',
    &weakzone_import_pointcloud_size, 0,
    "Import weak zone: Number of points of an imported weak zone",
  YMIR_OPTIONS_D, "weakzone-import-thickness", '\0',
    &weakzone_import_thickness, SL_WEAKZONE_IMPORT_THICKNESS,
    "Import weak zone: Width about center of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-import-thickness-const", '\0',
    &weakzone_import_thickness_const, SL_WEAKZONE_IMPORT_THICKNESS_CONST,
    "Import weak zone: Width of smoothing of edges of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-import-weak-factor", '\0',
    &weakzone_import_weak_factor, SL_WEAKZONE_IMPORT_WEAK_FACTOR,
    "Import weak zone: Min value of weak zone factor",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-longitude", '\0',
    &weakzone_2pl_subdu_lon, SL_WEAKZONE_2PLATES_SUBDU_LONGITUDE,
    "2plates weak zone: Longitude in interval (-pi/8, pi/8), "
    "where weak zone begins",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-dip-angle", '\0',
    &weakzone_2pl_subdu_dip_angle, SL_WEAKZONE_2PLATES_SUBDU_DIP_ANGLE,
    "2plates weak zone: Dip angle of weak zone (in degrees > 0)",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-depth", '\0',
    &weakzone_2pl_subdu_depth, SL_WEAKZONE_2PLATES_SUBDU_DEPTH,
    "2plates weak zone: Depth of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-width", '\0',
    &weakzone_2pl_subdu_width, SL_WEAKZONE_2PLATES_SUBDU_WIDTH,
    "2plates weak zone: Width of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-thickness", '\0',
    &weakzone_2pl_subdu_thickness, SL_WEAKZONE_2PLATES_SUBDU_THICKNESS,
    "2plates weak zone: Width at center of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-thickness-const", '\0',
    &weakzone_2pl_subdu_thickness_const,
    SL_WEAKZONE_2PLATES_SUBDU_THICKNESS_CONST,
    "2plates weak zone: Width of smoothing of edges of weak zone [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-subdu-weak-factor", '\0',
    &weakzone_2pl_subdu_weak_factor, SL_WEAKZONE_2PLATES_SUBDU_WEAK_FACTOR,
    "2plates weak zone: Value of weak zone factor",
  YMIR_OPTIONS_D, "weakzone-2plates-ridge-depth", '\0',
    &weakzone_2pl_ridge_depth, SL_WEAKZONE_2PLATES_RIDGE_DEPTH,
    "2plates weak zone: Depth of weak zone in left corner of domain [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-ridge-width", '\0',
    &weakzone_2pl_ridge_width, SL_WEAKZONE_2PLATES_RIDGE_WIDTH,
    "2plates weak zone: Width of weak zone in left corner of domain [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-ridge-smoothwidth", '\0',
    &weakzone_2pl_ridge_smoothwidth, SL_WEAKZONE_2PLATES_RIDGE_SMOOTHWIDTH,
    "2plates weak zone: Smoothing width of edges of weak zone in corner [m]",
  YMIR_OPTIONS_D, "weakzone-2plates-ridge-weak-factor", '\0',
    &weakzone_2pl_ridge_weak_factor, SL_WEAKZONE_2PLATES_RIDGE_WEAK_FACTOR,
    "2plates weak zone: Value of weak zone factor for weak zone in corner",
  YMIR_OPTIONS_I, "weakzone-2plates-ridge-mesh-align", '\0',
    &weakzone_2pl_ridge_mesh_align, 0,
    "2plates weak zone: Align ridge weak zone with mesh elements",

  YMIR_OPTIONS_S, "velocity-import", '\0',
    &slabs_velocity_import_filename, SLABS_DEFAULT_VELOCITY_IMPORT_FILENAME,
    "Path to file for importing velocity field",
  YMIR_OPTIONS_S, "pressure-import", '\0',
    &slabs_pressure_import_filename, SLABS_DEFAULT_PRESSURE_IMPORT_FILENAME,
    "Path to file for importing pressure field",

  YMIR_OPTIONS_I, "viscosity", '\0',
    &viscosity_type, SL_VISCOSITY_CONST,
    "Viscosity type: 0: constant=1; 1: linear, temperature dependent; "
    "2: nonlinear, temperature & strain rate dependent",
  YMIR_OPTIONS_I, "viscosity-for-init-nl-stokes", '\0',
    &viscosity_type_for_init_nl_stokes, SL_VISCOSITY_INIT_NL_STOKES_DEFAULT,
    "Viscosity type for initial viscosity of nonlinear Stokes problem",
  YMIR_OPTIONS_S, "viscosity-model", '\0',
    &visc_model_type_name, SL_VISCOSITY_MODEL_TYPE,
    "Viscosity model: WYUL, UWYL, UYWL, UYWL_RSHIFT, UWL_IIE_REG",
  YMIR_OPTIONS_D, "viscosity-IIe-regularization", '\0',
    &visc_IIe_reg, SL_VISCOSITY_IIE_REG,
    "Regularization for 2nd invariant by adding const to aviod division by 0",
  YMIR_OPTIONS_D, "viscosity-min", '\0',
    &visc_min, SL_VISCOSITY_MIN,
    "Lower bound for viscosity",
  YMIR_OPTIONS_D, "viscosity-max", '\0',
    &visc_max, SL_VISCOSITY_MAX,
    "Upper bound for viscosity",
  YMIR_OPTIONS_D, "viscosity-temp-max", '\0',
    &visc_temp_max, SL_VISCOSITY_MAX,
    "Upper bound for temperature dependent viscosity for model UWYUL",
  YMIR_OPTIONS_D, "viscosity-scaling", '\0',
    &visc_scaling, SL_VISCOSITY_SCALING,
    "Scaling factor for viscosity",
  YMIR_OPTIONS_D, "viscosity-temp-decay", '\0',
    &visc_temp_decay, SL_VISCOSITY_TEMP_DECAY,
    "Decay of temperature dependent viscosity",
  YMIR_OPTIONS_I, "viscosity-lower-mantle", '\0',
    &visc_lower_mantle, SL_VISCOSITY_LOWER_MANTLE,
    "Include lower mantle with linear rheology",
  YMIR_OPTIONS_D, "viscosity-lower-mantle-scaling", '\0',
    &visc_lower_mantle_scaling, SL_VISCOSITY_LOWER_MANTLE_SCALING,
    "Scaling factor for viscosity in lower mantle",
  YMIR_OPTIONS_D, "viscosity-lower-mantle-temp-decay", '\0',
    &visc_lower_mantle_temp_decay, SL_VISCOSITY_LOWER_MANTLE_TEMP_DECAY,
    "Decay of temperature dependent viscosity in lower mantle",
  YMIR_OPTIONS_D, "viscosity-stress-exponent", '\0',
    &visc_stress_exp, SL_VISCOSITY_STRESS_EXPONENT,
    "Stress exponent for the second invariant of the strain rate",
  YMIR_OPTIONS_D, "viscosity-stress-yield", '\0',
    &visc_stress_yield, SL_VISCOSITY_STRESS_YIELD,
    "Yielding stress for nonlinear rheology",
  YMIR_OPTIONS_D, "viscosity-yielding-reg", '\0',
    &visc_yield_reg, 0.0,
    "Regularization for yielding regions",

  YMIR_OPTIONS_I, "viscosity-coarsen-eval", '\0',
    &visc_coarsen_eval, 0,
    "Evaluate viscosity on coarse mesh in case of p/h-coarsening",
  YMIR_OPTIONS_I, "viscosity-p-coarsen-eval", '\0',
    &visc_p_coarsen_eval, 0,
    "Evaluate viscosity on coarse mesh in case of p-coarsening",
  YMIR_OPTIONS_I, "viscosity-coarsen-type", '\0',
    &visc_coarsen_type, SL_VISCOSITY_COARSEN_EVAL,
    "Type of coarsening for viscosity",
  YMIR_OPTIONS_I, "weakzone-coarsen-type", '\0',
    &weakzone_coarsen_type, SL_WEAKZONE_COARSEN_EVAL,
    "Type of coarsening for weak zones",
  YMIR_OPTIONS_D, "viscosity-lower-upper-transition-zone-incr", '\0',
    &visc_lower_upper_transition_zone_incr, 0.0,
    "Create cont transition zone instead of discont LM/UM interface",

  YMIR_OPTIONS_D, "right-hand-side-scaling", '\0', &rhs_scaling, SL_RHS_SCALING,
    "Scaling factor for right-hand side",
  YMIR_OPTIONS_I, "right-hand-side-random", '\0', &rhs_random, 0,
    "Set random right-hand side",
  YMIR_OPTIONS_I, "right-hand-side-multiply-in-weakzone", '\0', &rhs_multiply_in_weak_zone, 0,
    "Multiply in weak zone to right-hand side",
  YMIR_OPTIONS_I, "plume-type", '\0', &plume_type, SL_PLUME_NONE,
    "Add plume to temperature field. 0: none 1: inner plume; "
    "2: rising plume. (Plumes overridden by plume settings below)",
  YMIR_OPTIONS_D, "plume-center-x", '\0', &plume_center_x, 0.0,
    "Right-hand side plume: x-coordinate of center",
  YMIR_OPTIONS_D, "plume-center-y", '\0', &plume_center_y, 0.0,
    "Right-hand side plume: y-coordinate of center",
  YMIR_OPTIONS_D, "plume-center-z", '\0', &plume_center_z, 0.0,
    "Right-hand side plume: z-coordinate of center",
  YMIR_OPTIONS_D, "plume-decay", '\0', &plume_decay, -1.0,
    "Right-hand side plume decay",
  YMIR_OPTIONS_D, "plume-scaling", '\0', &plume_scaling, 0.0,
    "Right-hand side plume scaling",

  /* discretization options */
#ifndef YMIR_N
  YMIR_OPTIONS_I, "order", 'N', &slabs_order, SLABS_DEFAULT_DISCR_ORDER,
    "Approximation order",
#endif
  YMIR_OPTIONS_I, "minlevel", 'l', &slabs_minlevel, SLABS_DEFAULT_DISCR_REFINEMENT_MINLEVEL,
    "Minumum mesh refinement level",
  YMIR_OPTIONS_I, "maxlevel", 'L', &slabs_maxlevel, SLABS_DEFAULT_DISCR_REFINEMENT_MAXLEVEL,
    "Maximum mesh refinement level",

  YMIR_OPTIONS_S, "refine", 'r', &slabs_refine, SLABS_DEFAULT_DISCR_REFINEMENT_TYPE,
    "Mesh refinement type",
  YMIR_OPTIONS_S, "refine-radius", '\0',
    &refine_radius_str, NULL,
    "Array of (comma separated) radii for increasing refinement of mesh, "
    "assuming the shell to have thickness = 1",

  YMIR_OPTIONS_I, "import-mesh-order", '\0',
    &slabs_import_order, SLABS_DEFAULT_DISCR_IMPORT_ORDER,
    "Mesh for data inport: polynomial order of discretization",
  YMIR_OPTIONS_I, "import-mesh-minlevel", '\0',
    &slabs_import_minlevel, SLABS_DEFAULT_DISCR_IMPORT_MINLEVEL,
    "Mesh for data inport: minimum mesh refinement level",

  YMIR_OPTIONS_D, "refine-surface-max-dist", '\0',
    &refine_surface_maxdist, 0.0,
    "Refine at surface: distance for refinement",
  YMIR_OPTIONS_D, "refine-surface-element-resolution", '\0',
    &refine_surface_elem_res, 0.0,
    "Refine at surface: side length of elem [m]",
  YMIR_OPTIONS_D, "refine-lower-upper-interface-max-dist", '\0',
    &refine_lm_um_interface_maxdist, 0.0,
    "Refine about the lower/upper mantle interface: distance for refinement",
  YMIR_OPTIONS_D, "refine-lower-upper-interface-element-resolution", '\0',
    &refine_lm_um_interface_elem_res, 0.0,
    "Refine about the lower/upper mantle interface: side length of elem [m]",
  YMIR_OPTIONS_I, "enforce-refinement-lower-upper-interface", '\0',
    &enforce_refinement_lm_um_interface, 0,
    "Enforce refinement about the LM/UM interface during GMG coarsening",

  YMIR_OPTIONS_I, "init-amr-max-steps", '\0',
    &init_amr_max_steps, 0,
    "Init AMR: max number of AMR steps",
  YMIR_OPTIONS_D, "init-amr-rel-threshold", '\0',
    &init_amr_rel_threshold, 0.0,
    "Init AMR: Keep same mesh below this rel. threshold for marked quadrants",
  YMIR_OPTIONS_D, "init-amr-num-elements-max", '\0',
    &init_amr_n_elements_max, 0.0,
    "Init AMR: Max number of elements in mesh",
  YMIR_OPTIONS_I, "init-amr-override-order", '\0',
    &init_amr_override_order, 0,
    "Init AMR: overrider polynomial order of discretization during AMR",
  YMIR_OPTIONS_I, "init-amr-lower-mantle", '\0',
    &init_amr_lower_mantle, 0,
    "Initial AMR: activate AMR in lower mantle",

  YMIR_OPTIONS_I, "init-amr-visc-indicator", '\0',
    &init_amr_visc_indicator_type, SL_AMR_INDICATOR_NONE,
    "Initial AMR w.r.t. viscosity: type of refinement indicator",
  YMIR_OPTIONS_D, "init-amr-visc-tol-min", '\0',
    &init_amr_visc_tol_min, 0.0,
    "Initial AMR w.r.t. viscosity: lower tol under which mesh is coarsened",
  YMIR_OPTIONS_D, "init-amr-visc-tol-max", '\0',
    &init_amr_visc_tol_max, 0.0,
    "Initial AMR w.r.t. viscosity: upper tol above which mesh is refined",
  YMIR_OPTIONS_D, "init-amr-visc-inside-plates-tol-min", '\0',
    &init_amr_visc_in_plates_tol_min, 0.0,
    "Initial AMR w.r.t. viscosity inside plate: lower tol for coarsening",
  YMIR_OPTIONS_D, "init-amr-visc-inside-plates-tol-max", '\0',
    &init_amr_visc_in_plates_tol_max, 0.0,
    "Initial AMR w.r.t. viscosity inside plate: upper tol for refinement",

  YMIR_OPTIONS_I, "init-amr-weak-subdu-indicator", '\0',
    &init_amr_weak_subdu_indicator_type, SL_AMR_INDICATOR_NONE,
    "Initial AMR w.r.t. subduction weak zone: type of refinement indicator",
  YMIR_OPTIONS_D, "init-amr-weak-subdu-tol-min", '\0',
    &init_amr_weak_subdu_tol_min, 0.0,
    "Initial AMR w.r.t. subduction weak zone: mesh coarsened below this tol",
  YMIR_OPTIONS_D, "init-amr-weak-subdu-tol-max", '\0',
    &init_amr_weak_subdu_tol_max, 0.0,
    "Initial AMR w.r.t. subduction weak zone: mesh refined above this tol",
  YMIR_OPTIONS_D, "init-amr-weak-subdu-element-resolution", '\0',
    &init_amr_weak_subdu_elem_res, 0.0,
    "Initial AMR w.r.t. subduction weak zone: side length of element [m]",

  YMIR_OPTIONS_I, "init-amr-weak-ridge-indicator", '\0',
    &init_amr_weak_ridge_indicator_type, SL_AMR_INDICATOR_NONE,
    "Initial AMR w.r.t. ridge weak zone: type of refinement indicator",
  YMIR_OPTIONS_D, "init-amr-weak-ridge-tol-min", '\0',
    &init_amr_weak_ridge_tol_min, 0.0,
    "Initial AMR w.r.t. ridge weak zone: mesh coarsened below this tol",
  YMIR_OPTIONS_D, "init-amr-weak-ridge-tol-max", '\0',
    &init_amr_weak_ridge_tol_max, 0.0,
    "Initial AMR w.r.t. ridge weak zone: mesh refined above this tol",
  YMIR_OPTIONS_D, "init-amr-weak-ridge-element-resolution", '\0',
    &init_amr_weak_ridge_elem_res, 0.0,
    "Initial AMR w.r.t. ridge weak zone: side length of element [m]",

  YMIR_OPTIONS_I, "init-amr-weak-import-indicator", '\0',
    &init_amr_weak_import_indicator_type, SL_AMR_INDICATOR_NONE,
    "Initial AMR w.r.t. imported weak zone: type of refinement indicator",
  YMIR_OPTIONS_D, "init-amr-weak-import-tol-min", '\0',
    &init_amr_weak_import_tol_min, 0.0,
    "Initial AMR w.r.t. imported weak zone: mesh coarsened below this tol",
  YMIR_OPTIONS_D, "init-amr-weak-import-tol-max", '\0',
    &init_amr_weak_import_tol_max, 0.5,
    "Initial AMR w.r.t. imported weak zone: mesh refined above this tol",
  YMIR_OPTIONS_D, "init-amr-weak-import-element-resolution", '\0',
    &init_amr_weak_import_elem_res, 0.0,
    "Initial AMR w.r.t. imported weak zone: side length of element [m]",

  YMIR_OPTIONS_I, "init-amr-rhs-indicator", '\0',
    &init_amr_rhs_indicator_type, SL_AMR_INDICATOR_NONE,
    "Initial AMR w.r.t. rhs: type of refinement indicator",
  YMIR_OPTIONS_D, "init-amr-rhs-tol-min", '\0',
    &init_amr_rhs_tol_min, 0.0,
    "Initial AMR w.r.t. rhs: mesh coarsened below this tol",
  YMIR_OPTIONS_D, "init-amr-rhs-tol-max", '\0',
    &init_amr_rhs_tol_max, 0.0,
    "Initial AMR w.r.t. rhs: mesh refined above this tol",
  YMIR_OPTIONS_D, "init-amr-rhs-norm-shift", '\0',
    &init_amr_rhs_norm_shift, 1.0,
    "Initial AMR w.r.t. rhs: shift norm of rhs by this value (avoid DIV 0)",

  YMIR_OPTIONS_I, "init-amr-post-uniform-num-steps", '\0',
    &init_amr_post_uniform_n_steps, 0,
    "Post init AMR: additional number of steps of uniform refinement",

  YMIR_OPTIONS_I, "amr-max-steps", '\0',
    &amr_max_steps, 0,
    "AMR: max number of AMR steps",
  YMIR_OPTIONS_D, "amr-rel-threshold", '\0',
    &amr_rel_threshold, 0.0,
    "AMR: Keep same mesh below this rel. threshold for marked quadrants",
  YMIR_OPTIONS_D, "amr-num-elements-max", '\0',
    &amr_n_elements_max, 0.0,
    "AMR: Max number of elements in mesh",
  YMIR_OPTIONS_I, "amr-lower-mantle", '\0',
    &amr_lower_mantle, 0,
    "AMR: activate AMR in lower mantle",

  YMIR_OPTIONS_I, "amr-visc-indicator", '\0',
    &amr_visc_indicator_type, SL_AMR_INDICATOR_NONE,
    "AMR w.r.t. viscosity: type of refinement indicator",
  YMIR_OPTIONS_D, "amr-visc-tol-min", '\0',
    &amr_visc_tol_min, 0.0,
    "AMR w.r.t. viscosity: lower tol below which mesh is coarsened",
  YMIR_OPTIONS_D, "amr-visc-tol-max", '\0',
    &amr_visc_tol_max, 0.0,
    "AMR w.r.t. viscosity: upper tol above which mesh is refined",

  YMIR_OPTIONS_I, "amr-visc-dr-indicator", '\0',
    &amr_visc_dr_indicator_type, SL_AMR_INDICATOR_NONE,
    "AMR w.r.t. viscosity dynamic range: type of refinement indicator",
  YMIR_OPTIONS_D, "amr-visc-dr-tol-min", '\0',
    &amr_visc_dr_tol_min, 0.0,
    "AMR w.r.t. viscosity dynamic range: lower tol, coarsen mesh if below",
  YMIR_OPTIONS_D, "amr-visc-dr-tol-max", '\0',
    &amr_visc_dr_tol_max, 0.0,
    "AMR w.r.t. viscosity dynamic range: upper tol, refine mesh if above",

  YMIR_OPTIONS_I, "amr-strain-rate-indicator", '\0',
    &amr_strain_rate_indicator_type, SL_AMR_INDICATOR_NONE,
    "AMR w.r.t. strain rate: type of refinement indicator",
  YMIR_OPTIONS_D, "amr-strain-rate-tol-min", '\0',
    &amr_strain_rate_tol_min, 0.0,
    "AMR w.r.t. strain rate: lower tol below which mesh is coarsened",
  YMIR_OPTIONS_D, "amr-strain-rate-tol-max", '\0',
    &amr_strain_rate_tol_max, 0.0,
    "AMR w.r.t. strain rate: upper tol above which mesh is refined",

  YMIR_OPTIONS_I, "amr-log-maxlevel", '\0',
    &amr_log_maxlevel, 0,
    "AMR: print max refinement level of p4est mesh",
  YMIR_OPTIONS_I, "mesh-partitioning-type", '\0',
    &mesh_partitioning_type, SL_MESH_PARTITIONING_ELEM,
    "type of partitioning of the mesh among processors",

  /* nonlinear solver options */
  YMIR_OPTIONS_I, "nonlinear-type", '\0',
    &nl_solver_type, SLABS_DEFAULT_NL_SOLVER_TYPE,
    "Nonlinear solver: 0: Picard/fixed point; 1: Newton, 2: Picard-Newton",
  YMIR_OPTIONS_I, "nonlinear-primaldual-type", '\0',
    &nl_solver_primaldual_type, SLABS_DEFAULT_NL_SOLVER_PRIMALDUAL_TYPE,
    "Nonlinear solver primal-dual type: 0: standard, "
    "1: stress, 2: normalized strain",
  YMIR_OPTIONS_I, "nonlinear-primaldual-scaling-type", '\0',
    &nl_solver_primaldual_scal_type,
    SLABS_DEFAULT_NL_SOLVER_PRIMALDUAL_SCAL_TYPE,
    "Primal-dual scaling type: 0: none, 1: inv. max norm, 2: normalized",
  YMIR_OPTIONS_I, "nonlinear-maxiter", '\0',
    &nl_maxiter, SLABS_DEFAULT_NL_SOLVER_MAXITER,
    "Nonlinear solver maximum number of iterations",
  YMIR_OPTIONS_D, "nonlinear-rtol", '\0',
    &nl_rtol, SLABS_DEFAULT_NL_SOLVER_RTOL,
    "Nonlinear solver relative tolerance",
  YMIR_OPTIONS_I, "nonlinear-initial-guess-type", '\0',
    &nl_initial_guess_type, SLABS_DEFAULT_NL_SOLVER_INIT_GUESS_TYPE,
    "Choose initial guess for nonlinear solver",

  YMIR_OPTIONS_I, "nonlinear-resume-at-iter", '\0',
    &nl_resume_at_iter, SLABS_DEFAULT_NL_SOLVER_RESUME_AT_ITER,
    "Resume nonlinear solver at this iteration",
  YMIR_OPTIONS_D, "nonlinear-resume-at-time", '\0',
    &nl_resume_at_time, SLABS_DEFAULT_NL_SOLVER_RESUME_AT_TIME,
    "Resume nonlinear solver at this elapsed time",
  YMIR_OPTIONS_D, "nonlinear-resume-prev-res", '\0',
    &nl_resume_prev_res, SLABS_DEFAULT_NL_SOLVER_RESUME_PREV_RES,
    "Resume nonlinear solver with this previous residual",
  YMIR_OPTIONS_D, "nonlinear-resume-init-res", '\0',
    &nl_resume_init_res, SLABS_DEFAULT_NL_SOLVER_RESUME_INIT_RES,
    "Resume nonlinear solver with this initial residual",

  YMIR_OPTIONS_D, "nonlinear-forcing-exponent", '\0',
    &nl_forcing_exponent, SLABS_DEFAULT_NL_SOLVER_FORCING_EXPONENT,
    "Exponent for forcing for the linear solver during nonlinear iterations",
  YMIR_OPTIONS_D, "nonlinear-forcing-max", '\0',
    &nl_forcing_max, SLABS_DEFAULT_NL_SOLVER_FORCING_MAX,
    "Max forcing, i.e., max Krylov relative tolerance",
  YMIR_OPTIONS_I, "nonlinear-forcing-max-progressive-iter", '\0',
    &nl_forcing_max_progressive_iter,
    SLABS_DEFAULT_NL_SOLVER_FORCING_MAX_PROGRESSIVE_ITER,
    "#iterations (~ variance) for progressive max forcing",
  YMIR_OPTIONS_D, "nonlinear-forcing-total-min", '\0',
    &nl_forcing_total_min, SLABS_DEFAULT_NL_SOLVER_FORCING_TOTAL_MIN,
    "Min forcing, i.e., min Krylov relative tolerance over all nonlinear steps",
  YMIR_OPTIONS_I, "nonlinear-forcing-saveguard", '\0',
    &nl_forcing_saveguard, SLABS_DEFAULT_NL_SOLVER_FORCING_SAVEGUARD,
    "Saveguarding of Krylov rel tol during nl iter to avoid oversolving",
  YMIR_OPTIONS_D, "nonlinear-forcing-saveguard-threshold", '\0',
    &nl_forcing_threshold, SLABS_DEFAULT_NL_SOLVER_FORCING_SAVEGUARD_THRESHOLD,
    "Threshold for saveguarding of Krylov rel tol during nonlinear iter",

  YMIR_OPTIONS_I, "nonlinear-step-length-reduction-type", '\0',
    &nl_step_length_reduction_type,
    SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_REDUCTION_TYPE,
    "Type of step reduction algorithm for step search during nonlinear iter",
  YMIR_OPTIONS_D, "nonlinear-step-length-reduction-min", '\0',
    &nl_step_length_reduction_min,
    SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_REDUCTION_MIN,
    "Minimal allowed step reduction for step search during nonlinear iter",
  YMIR_OPTIONS_D, "nonlinear-step-length-reduction-max", '\0',
    &nl_step_length_reduction_max,
    SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_REDUCTION_MAX,
    "Maximal allowed step reduction for step search during nonlinear iter",
  YMIR_OPTIONS_D, "nonlinear-step-length-reduction-reg", '\0',
    &nl_step_length_reduction_reg,
    SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_REDUCTION_REG,
    "Step reduction regularization for step search during nonlinear iter",
  YMIR_OPTIONS_D, "nonlinear-step-length-descend-condition-relaxation", '\0',
    &nl_step_length_descend_relax,
    SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_DESCEND_COND_RELAX,
    "Relaxation factor for the decrease condition during step length search",
  YMIR_OPTIONS_D, "nonlinear-step-length-min", '\0',
    &nl_step_length_min, SLABS_DEFAULT_NL_SOLVER_STEPLENGTH_MIN,
    "Minimal allowed step length for backtracking at nonlinear iterations",

  YMIR_OPTIONS_D, "nonlinear-switch-picard-step-length-min", '\0',
    &nl_switch_picard_step_length_min,
    SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_STEP_LENGTH_MIN,
    "Condition to switch from Newton to Picard: min Newton step length",
  YMIR_OPTIONS_I, "nonlinear-switch-picard-after-amr", '\0',
    &nl_switch_picard_after_amr,
    SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_AFTER_AMR,
    "Switch from Newton to Picard: switch after AMR",
  YMIR_OPTIONS_I, "nonlinear-switch-picard-init", '\0',
    &nl_switch_picard_init,
    SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_INIT,
    "Switch from Newton to Picard: switch at the beginning of nl solve",
  YMIR_OPTIONS_I, "nonlinear-switch-picard-maxiter", '\0',
    &nl_switch_picard_maxiter, SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_MAXITER,
    "Switch from Newton to Picard: max number of Picard steps",
  YMIR_OPTIONS_D, "nonlinear-switch-picard-rtol", '\0',
    &nl_switch_picard_rtol, SLABS_DEFAULT_NL_SOLVER_SWITCH_PICARD_RTOL,
    "Switch from Newton to Picard: relative tolerance for set of "
    "consecutive Picard steps",

  YMIR_OPTIONS_I, "nonlinear-norm-type", '\0',
    &nl_norm_type, SLABS_DEFAULT_NL_SOLVER_NORM_TYPE,
    "Norm for convergence checks within nonlinear solver",
  YMIR_OPTIONS_D, "nonlinear-norm-Hminus1-mass-scaling", '\0',
    &nl_norm_Hminus1_mass_scaling,
    SLABS_DEFAULT_NL_SOLVER_NORM_HMINUS1_MASS_SCALING,
    "Weight factor of mass matrix of H^-1 norm operator",

  YMIR_OPTIONS_I, "nonlinear-schur-diag-type", '\0',
    &nl_schur_diag_type, SLABS_DEFAULT_NL_SOLVER_SCHUR_DIAG_TYPE,
    "Type of diagonal pressure Schur complement approximation.",
  YMIR_OPTIONS_I, "nonlinear-scaling-type", '\0',
    &nl_scaling_type, SLABS_DEFAULT_NL_SOLVER_SCALING_TYPE,
    "Blockwise diagonal scaling of Stokes matrix",
  YMIR_OPTIONS_I, "nonlinear-project-out-nullspace", '\0',
    &nl_project_out_nullspace, SLABS_DEFAULT_NL_SOLVER_PROJECT_NULLSPACE,
    "Type of projection to remove nullspaces during nonlinear Stokes solve",
  YMIR_OPTIONS_I, "nonlinear-enforce-unscaled-reduction", '\0',
    &nl_enforce_unscaled_reduction,
    SLABS_DEFAULT_NL_SOLVER_ENFORCE_UNSCALED_REDUCTION,
    "Enforce the unscaled residual reduction during the linear solves",

  YMIR_OPTIONS_D, "nonlinear-grid-continuation-init-threshold", '\0',
    &nl_grid_cont_init_amr_threshold, 0.0,
    "Grid continuation at beginning of nl solve: relative AMR tolerance",
  YMIR_OPTIONS_D, "nonlinear-grid-continuation-init-forcing", '\0',
    &nl_grid_cont_init_forcing, 0.0,
    "Grid continuation at beginning of nl solve: relative Krylov tolerance",
  YMIR_OPTIONS_D, "nonlinear-grid-continuation-init-forcing-exp", '\0',
    &nl_grid_cont_init_forcing_exp, 0.0,
    "Grid continuation at beginning of nl solve: forcing exponent",
  YMIR_OPTIONS_I, "nonlinear-grid-continuation-init-steps", '\0',
    &nl_grid_cont_init_steps, 0,
    "Grid continuation at beginning of nl solve: number of nl. steps",
  YMIR_OPTIONS_I, "nonlinear-grid-continuation-skipsteps", '\0',
    &nl_grid_cont_skipsteps, 0,
    "Grid continuation: skip AMR for this number of nonlinear steps",
  YMIR_OPTIONS_I, "nonlinear-grid-continuation-maxsteps", '\0',
    &nl_grid_cont_maxsteps, 0,
    "Grid continuation: max number of nonlinear steps",
  YMIR_OPTIONS_D, "nonlinear-viscosity-bounds-continuation-min", '\0',
    &nl_visc_bounds_cont_min, -1.0,
    "Viscosity bounds continuation: min viscosity bound",
  YMIR_OPTIONS_D, "nonlinear-viscosity-bounds-continuation-max", '\0',
    &nl_visc_bounds_cont_max, -1.0,
    "Viscosity bounds continuation: max viscosity bound",
  YMIR_OPTIONS_I, "nonlinear-viscosity-bounds-continuation-steps", '\0',
    &nl_visc_bounds_cont_steps, 0,
    "Viscosity bounds continuation: number of nonlinear steps",

  YMIR_OPTIONS_I, "nonlinear-check-derivative", '\0',
    &nl_check_derivative, 0,
    "Check Newton derivative after construction",
  YMIR_OPTIONS_I, "log-physics-stats", '\0',
    &log_physics_stats, 0,
    "Monitor physics statistics",

  /* Linear solver options */
  YMIR_OPTIONS_I, "krylov-type", '\0', &krylov_type, SL_KRYLOV_GMRES,
    "0: MINRES, 1: GMRES",
  YMIR_OPTIONS_I, "krylov-gmres-num-vecs", '\0', &gmres_num_vecs, SL_KRYLOV_GMRES_NUM_VECS,
    "Number of GMRES vectors",
  YMIR_OPTIONS_I, "krylov-maxiter", '\0', &krylov_maxiter, SL_KRYLOV_MAXITER,
    "Krylov maximum number of iterations",
  YMIR_OPTIONS_D, "krylov-atol", '\0', &krylov_atol, SL_KRYLOV_ATOL,
    "Krylov absolute tolerance",
  YMIR_OPTIONS_D, "krylov-rtol", '\0', &krylov_rtol, SL_KRYLOV_RTOL,
    "Krylov relative tolerance",

  YMIR_OPTIONS_I, "solve-stress-block-only", '\0', &lin_solve_stress_block_only, 0,
    "0: Solve linear Stokes system, 1: Solve linear viscous stress block",
  YMIR_OPTIONS_I, "solve-press-bbt-only", '\0', &lin_solve_press_bbt_only, 0,
    "Solve pressure Laplacian, BB^T",
  YMIR_OPTIONS_I, "solve-cnode-bbt-only", '\0', &lin_solve_cnode_bbt_only, 0,
    "Solve approximation of pressure Laplacian, BB^T, on continuous space",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}

/**
 * Converts a string into a double array.
 */
static double *
slabs_convert_str_to_double (char *in, int *num_values)
{
  char               *str_prev, *str_next;
  int                 j = 0;
  double             *out;

  /* count number of floating point numbers */
  *num_values = 0;
  str_next = in;
  while ( (str_next = strpbrk (str_next, "+-0123456789.eE")) != NULL ) {
    (*num_values)++;
    /* look for separator */
    if ( (str_next = strpbrk (str_next, ", ")) == NULL ) {
      break;
    }
  }

  /* allocate memory for array */
  out = YMIR_ALLOC (double, *num_values);

  /* fill in array */
  str_next = in;
  str_prev = in;
  while (j < *num_values) {
    out[j] = strtod (str_prev, &str_next);
    if (str_next == str_prev) {
      str_next++;
      str_prev++;
    }
    else {
      j++;
      str_prev = str_next;
    }
  }

  YMIR_ASSERT (j == *num_values);

  /* return array with double values */
  return out;
}

/**
 * Processes options.
 */
void
slabs_setup_process_options (slabs_physics_options_t *physics_options,
                             slabs_discr_options_t *discr_options,
                             slabs_nl_solver_options_t *solver_options)
{
  const char         *this_fn_name = "slabs_setup_process_options";
  slabs_domain_shape_t  domain_shape;
  slabs_temp_t        temp_type;
  char               *temp_import_filename_bin = NULL;
  slabs_weakzone_t    weakzone_type;
  char               *weakzone_import_filename_bin = NULL;
  slabs_viscosity_model_t  visc_model_type;

  /*
   * Physics Options
   */

  /* set shape of domain & mangll forest from input */
  if (strcmp (domain_shape_name, "unitcube") == 0) {
    domain_shape = SL_DOMAIN_CUBE;
  }
  else if (strcmp (domain_shape_name, "brick") == 0) {
    domain_shape = SL_DOMAIN_BRICK;
  }
  else if (strcmp (domain_shape_name, "shell") == 0) {
    domain_shape = SL_DOMAIN_SHELL;
  }
  else if (strcmp (domain_shape_name, "shell_chunk") == 0) {
    domain_shape = SL_DOMAIN_SHELL_CHUNK;
    domain_shape_name = "shell";
  }
  else if (strcmp (domain_shape_name, "shell_slice") == 0) {
    domain_shape = SL_DOMAIN_SHELL_SLICE;
    domain_shape_name = "shell";
  }
  else { /* if unknown domain name */
    YMIR_ABORTF ("Invalid domain name `%s`", domain_shape_name);
  }

  /* store domain and BC options */
  physics_options->domain_shape = domain_shape;
  switch (domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_SHELL:
  case SL_DOMAIN_SHELL_CHUNK:
    physics_options->domain_brick_dx = 1;
    physics_options->domain_brick_dy = 1;
    physics_options->domain_brick_dz = 1;
    break;

  case SL_DOMAIN_BRICK:
  case SL_DOMAIN_SHELL_SLICE:
    physics_options->domain_brick_dx = brick_dx;
    physics_options->domain_brick_dy = brick_dy;
    physics_options->domain_brick_dz = brick_dz;
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
  slabs_physics_compute_domain_setup (physics_options);

  physics_options->bc_type = (slabs_vel_bc_t) bc_type;
  physics_options->bc_default_dirichlet_scale = bc_default_dir_scale;
  physics_options->p4est_import_filename = slabs_p4est_import_filename;

  /* set type of temperature */
  if (strcmp (temp_type_name, "NONE") == 0) {
    temp_type = SL_TEMP_NONE;
  }
  else if (strcmp (temp_type_name, "hot_core") == 0) {
    temp_type = SL_TEMP_HOT_CORE;
  }
  else if (strcmp (temp_type_name, "cold_plate") == 0) {
    temp_type = SL_TEMP_COLD_PLATE;
  }
  else if (strcmp (temp_type_name, "hot_core_cold_plate") == 0) {
    temp_type = SL_TEMP_HOT_CORE_COLD_PLATE;
  }
  else if (strcmp (temp_type_name, "2plates_poly2") == 0) {
    temp_type = SL_TEMP_2PLATES_POLY2;
  }
  else if (strcmp (temp_type_name, "2plates_rhea1") == 0) {
    temp_type = SL_TEMP_2PLATES_RHEA1;
  }
  else { /* if import from file */
    temp_type = SL_TEMP_IMPORT_FILE;
    temp_import_filename_bin = temp_type_name;
  }

  /* store temperature options */
  physics_options->temperature_type = temp_type;
  physics_options->temp_background_plate_age = temp_back_plate_age;

  physics_options->temp_import_plate_age_min = temp_import_plate_age_min;
  physics_options->temp_import_filename_txt = temp_import_filename_txt;
  physics_options->temp_import_filename_bin = temp_import_filename_bin;
  physics_options->temp_import_verification_out = temp_import_verification_out;

  physics_options->temp_2plates_trench_longitude = temp_2pl_trench_lon;
  physics_options->temp_2plates_dip_angle = temp_2pl_dip_angle;
  physics_options->temp_2plates_subd_depth = temp_2pl_subd_depth;
  physics_options->temp_2plates_subd_width = temp_2pl_subd_width;
  physics_options->temp_2plates_subd_edge_width = temp_2pl_subd_edge_width;
  physics_options->temp_2plates_subd_edge_smoothwidth =
    temp_2pl_subd_edge_smoothwidth;
  physics_options->temp_2plates_subd_plate_velocity = temp_2pl_subd_plate_vel;
  physics_options->temp_2plates_subd_plate_initial_age =
    temp_2pl_subd_plate_init_age;
  physics_options->temp_2plates_over_plate_age = temp_2pl_over_plate_age;

  /* set type of weak zone */
  if (strcmp (weakzone_type_name, "NONE") == 0) {
    weakzone_type = SL_WEAKZONE_NONE;
  }
  else if (strcmp (weakzone_type_name, "2plates_poly2") == 0) {
    weakzone_type = SL_WEAKZONE_2PLATES_POLY2;
  }
  else if (strcmp (weakzone_type_name, "2plates_rhea1") == 0) {
    weakzone_type = SL_WEAKZONE_2PLATES_RHEA1;
  }
  else { /* if import from file */
    weakzone_type = SL_WEAKZONE_IMPORT_FILE;
    weakzone_import_filename_bin = weakzone_type_name;
  }

  /* store weak zone options */
  physics_options->weakzone_type = weakzone_type;
  physics_options->weakzone_import_filename_bin = weakzone_import_filename_bin;
  physics_options->weakzone_import_filename_txt = weakzone_import_filename_txt;
  physics_options->weakzone_import_verification_out =
    weakzone_import_verification_out;
  physics_options->weakzone_import_n_points =
    weakzone_import_pointcloud_size;
  physics_options->weakzone_import_thickness = weakzone_import_thickness;
  physics_options->weakzone_import_thickness_const =
    weakzone_import_thickness_const;
  physics_options->weakzone_import_weak_factor = weakzone_import_weak_factor;

  physics_options->weakzone_2plates_subdu_longitude =
    weakzone_2pl_subdu_lon;
  physics_options->weakzone_2plates_subdu_dip_angle =
    weakzone_2pl_subdu_dip_angle;
  physics_options->weakzone_2plates_subdu_depth =
    weakzone_2pl_subdu_depth;
  physics_options->weakzone_2plates_subdu_width =
    weakzone_2pl_subdu_width;
  physics_options->weakzone_2plates_subdu_thickness =
    weakzone_2pl_subdu_thickness;
  physics_options->weakzone_2plates_subdu_thickness_const =
    weakzone_2pl_subdu_thickness_const;
  physics_options->weakzone_2plates_subdu_weak_factor =
    weakzone_2pl_subdu_weak_factor;
  physics_options->weakzone_2plates_ridge_depth =
    weakzone_2pl_ridge_depth;
  physics_options->weakzone_2plates_ridge_width =
    weakzone_2pl_ridge_width;
  physics_options->weakzone_2plates_ridge_smoothwidth =
    weakzone_2pl_ridge_smoothwidth;
  physics_options->weakzone_2plates_ridge_weak_factor =
    weakzone_2pl_ridge_weak_factor;

  /* adjust weak zone to temperature for `2plates_poly2` */
  if (physics_options->weakzone_type == SL_WEAKZONE_2PLATES_POLY2) {
    /* compute weak zone paramers from temperature parameters */
    if (physics_options->temperature_type == SL_TEMP_2PLATES_POLY2 ||
        physics_options->temperature_type == SL_TEMP_IMPORT_FILE) {
      slabs_2plates_poly2_set_weakzone_params_from_temp (physics_options);

      YMIR_GLOBAL_INFOF ("%s: New weak zone parameters: longitude %g, "
                         "dip angle %g, zone width %g km\n",
                         this_fn_name,
                         physics_options->weakzone_2plates_subdu_longitude,
                         physics_options->weakzone_2plates_subdu_dip_angle,
                         physics_options->weakzone_2plates_subdu_width / 1.0e3);
    }

    /* overwrite with user input */
    if ( SC_EPS < fabs (weakzone_2pl_subdu_lon -
                        SL_WEAKZONE_2PLATES_SUBDU_LONGITUDE) ) {
      physics_options->weakzone_2plates_subdu_longitude =
        weakzone_2pl_subdu_lon;
      YMIR_GLOBAL_INFOF ("%s: Weak zone parameter overwritten by user: "
                         "longitude %g",
                         this_fn_name,
                         physics_options->weakzone_2plates_subdu_longitude);
    }
    if ( SC_EPS < fabs (weakzone_2pl_subdu_dip_angle -
                        SL_WEAKZONE_2PLATES_SUBDU_DIP_ANGLE) ) {
      physics_options->weakzone_2plates_subdu_dip_angle =
        weakzone_2pl_subdu_dip_angle;
      YMIR_GLOBAL_INFOF ("%s: Weak zone parameter overwritten by user: "
                         "dip angle %g",
                         this_fn_name,
                         physics_options->weakzone_2plates_subdu_dip_angle);
    }
    if ( SC_EPS < fabs (weakzone_2pl_subdu_width -
                        SL_WEAKZONE_2PLATES_SUBDU_WIDTH) ) {
      physics_options->weakzone_2plates_subdu_width = weakzone_2pl_subdu_width;
      YMIR_GLOBAL_INFOF ("%s: Weak zone parameter overwritten by user: "
                         "zone width %g km\n",
                         this_fn_name,
                         physics_options->weakzone_2plates_subdu_width / 1.0e3);
    }
  }

  /* store velocity & pressure options */
  physics_options->velocity_import_filename = slabs_velocity_import_filename;
  physics_options->pressure_import_filename = slabs_pressure_import_filename;

  /* check viscosity bounds */
  YMIR_CHECK_ABORTF (
      visc_min <= 0.0 || visc_max <+ 0.0 || visc_min < visc_max,
      "Invalid viscosity bounds min=%g and max=%g", visc_min, visc_max);

  /* set viscosity model type */
  if (strcmp (visc_model_type_name, "WYUL") == 0) {
    visc_model_type = SL_VISCOSITY_MODEL_WYUL;
  }
  else if (strcmp (visc_model_type_name, "UWYUL") == 0) {
    visc_model_type = SL_VISCOSITY_MODEL_UWYUL;
  }
  else if (strcmp (visc_model_type_name, "UWYL") == 0) {
    visc_model_type = SL_VISCOSITY_MODEL_UWYL;
  }
  else if (strcmp (visc_model_type_name, "UWYL_LREG") == 0) {
    visc_model_type = SL_VISCOSITY_MODEL_UWYL_LREG;
  }
  else if (strcmp (visc_model_type_name, "UWYL_SHIFT_LREG") == 0) {
    visc_model_type = SL_VISCOSITY_MODEL_UWYL_SHIFT_LREG;
  }
  else if (strcmp (visc_model_type_name, "UYWL") == 0) {
    visc_model_type = SL_VISCOSITY_MODEL_UYWL;
  }
  else if (strcmp (visc_model_type_name, "UYWL_SHIFT") == 0) {
    visc_model_type = SL_VISCOSITY_MODEL_UYWL_SHIFT;
  }
  else if (strcmp (visc_model_type_name, "UWL_IIE_REG") == 0) {
    visc_model_type = SL_VISCOSITY_MODEL_UWL_IIE_REG;
  }
  else { /* if unknown viscosity model name */
    YMIR_ABORTF ("Invalid viscosity model `%s`", visc_model_type_name);
  }
  physics_options->viscosity_model_type = visc_model_type;

  /* store viscosity options */
  physics_options->viscosity_type = (slabs_viscosity_t) viscosity_type;
  physics_options->viscosity_type_for_init_nl_stokes =
    (slabs_viscosity_init_nl_stokes_t) viscosity_type_for_init_nl_stokes;
  physics_options->viscosity_IIe_regularization = visc_IIe_reg;
  physics_options->viscosity_min = visc_min;
  physics_options->viscosity_max = visc_max;
  physics_options->viscosity_temp_max = visc_temp_max;
  physics_options->viscosity_scaling = visc_scaling;
  physics_options->viscosity_temp_decay = visc_temp_decay;
  physics_options->viscosity_stress_exponent = visc_stress_exp;
  physics_options->viscosity_stress_yield = visc_stress_yield;
  physics_options->viscosity_yielding_reg = visc_yield_reg;

  if (0.0 < visc_lower_mantle_scaling) {
    physics_options->viscosity_lower_mantle_scaling = visc_lower_mantle_scaling;
  }
  else {
    physics_options->viscosity_lower_mantle_scaling = visc_scaling;
  }
  if (0.0 < visc_lower_mantle_temp_decay) {
    physics_options->viscosity_lower_mantle_temp_decay =
      visc_lower_mantle_temp_decay;
  }
  else {
    physics_options->viscosity_lower_mantle_temp_decay = visc_temp_decay;
  }
  physics_options->viscosity_lower_upper_transition_zone = 0.0;

  physics_options->viscosity_coarsen_eval = visc_coarsen_eval;
  physics_options->viscosity_p_coarsen_eval = visc_p_coarsen_eval;
  physics_options->weakzone_coarsen_type =
    (slabs_weakzone_coarsen_t) weakzone_coarsen_type;
  physics_options->viscosity_coarsen_type =
    (slabs_viscosity_coarsen_t) visc_coarsen_type;
  physics_options->viscosity_lower_upper_transition_zone_incr =
    visc_lower_upper_transition_zone_incr / SL_EARTH_RADIUS;

  /* set additive component of nl. viscous stress coefficient */
  if (physics_options->viscosity_type == SL_VISCOSITY_NONLINEAR) {
    switch (physics_options->viscosity_model_type) {
    case SL_VISCOSITY_MODEL_WYUL:
    case SL_VISCOSITY_MODEL_UWYUL:
    case SL_VISCOSITY_MODEL_UWYL:
    case SL_VISCOSITY_MODEL_UYWL:
    case SL_VISCOSITY_MODEL_UYWL_SHIFT:
    case SL_VISCOSITY_MODEL_UWL_IIE_REG:
      break;
    case SL_VISCOSITY_MODEL_UWYL_LREG:
    case SL_VISCOSITY_MODEL_UWYL_SHIFT_LREG:
      if (0 < physics_options->viscosity_min) {
        ymir_nlstress_op_coeff_tensor_add =
          -2.0 * physics_options->viscosity_min;
        YMIR_GLOBAL_INFOF (
            "%s: Overriding option ymir_nlstress_op_coeff_tensor_add = %g\n",
            this_fn_name, ymir_nlstress_op_coeff_tensor_add);
      }
      break;
    default: /* unknown viscosity model type */
      YMIR_ABORT_NOT_REACHED ();
    }
  }

  /* set right-hand side options */
  physics_options->rhs_scaling = rhs_scaling;
  physics_options->rhs_random = rhs_random;
  physics_options->rhs_multiply_in_weak_zone = rhs_multiply_in_weak_zone;

  /* set plume options */
  physics_options->plume_type = (slabs_plume_t) plume_type;
  switch (domain_shape) {
  case SL_DOMAIN_CUBE:
    physics_options->plume_center_x = SL_CUBE_INNER_PLUME_CENTER_X;
    physics_options->plume_center_y = SL_CUBE_INNER_PLUME_CENTER_Y;
    physics_options->plume_center_z = SL_CUBE_INNER_PLUME_CENTER_Z;
    physics_options->plume_decay = SL_CUBE_INNER_PLUME_DECAY;
    physics_options->plume_scaling = SL_CUBE_INNER_PLUME_SCALING;

  case SL_DOMAIN_BRICK:
    physics_options->plume_center_x = SL_BRICK_INNER_PLUME_CENTER_X
                                     * (double) brick_dx / (double) brick_dz;
    physics_options->plume_center_y = SL_BRICK_INNER_PLUME_CENTER_Y
                                     * (double) brick_dy / (double) brick_dz;
    physics_options->plume_center_z = SL_BRICK_INNER_PLUME_CENTER_Z;
    physics_options->plume_decay = SL_BRICK_INNER_PLUME_DECAY;
    physics_options->plume_scaling = SL_BRICK_INNER_PLUME_SCALING;
    break;

  case SL_DOMAIN_SHELL:
    switch (plume_type) {
    case SL_PLUME_RISING:
      physics_options->plume_center_x = SL_SHELL_RISING_PLUME_CENTER_X;
      physics_options->plume_center_y = SL_SHELL_RISING_PLUME_CENTER_Y;
      physics_options->plume_center_z = SL_SHELL_RISING_PLUME_CENTER_Z;
      physics_options->plume_decay = SL_SHELL_RISING_PLUME_DECAY;
      physics_options->plume_scaling = SL_SHELL_RISING_PLUME_SCALING;
      break;

    default: /* set inner plume */
      physics_options->plume_center_x = SL_SHELL_INNER_PLUME_CENTER_X;
      physics_options->plume_center_y = SL_SHELL_INNER_PLUME_CENTER_Y;
      physics_options->plume_center_z = SL_SHELL_INNER_PLUME_CENTER_Z;
      physics_options->plume_decay = SL_SHELL_INNER_PLUME_DECAY;
      physics_options->plume_scaling = SL_SHELL_INNER_PLUME_SCALING;
    }
    break;

  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    physics_options->plume_center_x = SL_SHELL_SLICE_INNER_PLUME_CENTER_X;
    physics_options->plume_center_y = SL_SHELL_SLICE_INNER_PLUME_CENTER_Y;
    physics_options->plume_center_z = SL_SHELL_SLICE_INNER_PLUME_CENTER_Z;
    physics_options->plume_decay = SL_SHELL_SLICE_INNER_PLUME_DECAY;
    physics_options->plume_scaling = SL_SHELL_SLICE_INNER_PLUME_SCALING;
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
  /* override with (valid) user input for plume */
  if (plume_center_x != 0.0 || plume_center_y != 0.0 || plume_center_z != 0.0) {
    physics_options->plume_center_x = plume_center_x;
    physics_options->plume_center_y = plume_center_y;
    physics_options->plume_center_z = plume_center_z;
  }
  if (0.0 <= plume_decay) {
    physics_options->plume_decay = plume_decay;
  }
  if (plume_scaling != 0.0) {
    physics_options->plume_scaling = plume_scaling;
  }

  /*
   * Discretization Options
   */

  switch (domain_shape) {
  case SL_DOMAIN_CUBE:
  case SL_DOMAIN_BRICK:
    discr_options->X_fn = slabs_discr_identity_X;
    break;

  case SL_DOMAIN_SHELL:
    discr_options->X_fn = slabs_discr_shell_X;
    break;

  case SL_DOMAIN_SHELL_CHUNK:
  case SL_DOMAIN_SHELL_SLICE:
    discr_options->X_fn = slabs_discr_shell_chunk_and_slice_X;
    break;

  default: /* unknown domain type */
    YMIR_ABORT_NOT_REACHED ();
  }
  slabs_discr_set_order (discr_options, slabs_order);
  discr_options->minlevel = (int8_t) slabs_minlevel;
  discr_options->maxlevel = (int8_t) slabs_maxlevel;
  discr_options->refine = slabs_refine;

  /* override polynomial order */
#ifdef YMIR_N
  if (discr_options->order != ymir_n (discr_options->order)) {
    YMIR_GLOBAL_LERRORF (
        "%s: Warning: N = %d overriden with compiled value %d\n",
        this_fn_name, discr_options->order, ymir_n (discr_options->order));
    discr_options->order = ymir_n (discr_options->order);
  }
#endif

  if (refine_radius_str != NULL) {
    int                 i;

    discr_options->refine_radius = slabs_convert_str_to_double (
        refine_radius_str, &(discr_options->refine_n_radii));

    /* check radii to be ascending and in the interval [0,1] */
    YMIR_CHECK_ABORT (0.0 <= discr_options->refine_radius[0],
                      "Invalid mesh refine radii");
    for (i = 1; i < discr_options->refine_n_radii; i++) {
      YMIR_CHECK_ABORT (
          discr_options->refine_radius[i-1] <= discr_options->refine_radius[i]
          && discr_options->refine_radius[i] <= 1.0,
          "Invalid mesh refine radii");
    }
  }
  else {
    discr_options->refine_radius = NULL;
    discr_options->refine_n_radii = 0;
  }

  slabs_discr_compute_domain_size_normalization (discr_options,
                                                 physics_options);

  discr_options->refine_surface_maxdist =
    refine_surface_maxdist / SL_EARTH_RADIUS
    / (SL_SHELL_RADIUS_TOP - SL_SHELL_RADIUS_BOTTOM);
  discr_options->refine_surface_maxlevel =
    (int8_t) round (slabs_discr_resolution_to_level (
        refine_surface_elem_res, physics_options));
  if (0.0 < discr_options->refine_surface_maxdist) {
    YMIR_GLOBAL_INFOF ("%s: Refine at surface down to max level = %i\n",
        this_fn_name, (int) discr_options->refine_surface_maxlevel);
  }

  discr_options->refine_layer_maxdist =
    refine_lm_um_interface_maxdist / SL_EARTH_RADIUS
    / (SL_SHELL_RADIUS_TOP - SL_SHELL_RADIUS_BOTTOM);
  discr_options->refine_layer_maxlevel =
    (int8_t) round (slabs_discr_resolution_to_level (
        refine_lm_um_interface_elem_res, physics_options));
  if (0.0 < discr_options->refine_layer_maxdist) {
    YMIR_GLOBAL_INFOF ("%s: Refine at LM/UM interface down to max level = %i\n",
        this_fn_name, (int) discr_options->refine_layer_maxlevel);
  }
  discr_options->enforce_refinement_at_layer =
    enforce_refinement_lm_um_interface;

  if (0 < slabs_import_order) {
    discr_options->import_mesh_order = slabs_import_order;
  }
  else {
    discr_options->import_mesh_order = discr_options->order;
  }
  if (0 < slabs_import_minlevel) {
    discr_options->import_mesh_minlevel = slabs_import_minlevel;
  }
  else {
    discr_options->import_mesh_minlevel = discr_options->minlevel;
  }

  discr_options->init_amr_max_steps = init_amr_max_steps;
  discr_options->init_amr_rel_threshold = init_amr_rel_threshold;
  discr_options->init_amr_n_elements_max = init_amr_n_elements_max;
  discr_options->init_amr_override_order = init_amr_override_order;
  discr_options->init_amr_lower_mantle = init_amr_lower_mantle;

  discr_options->init_amr_visc_indicator_type =
    (slabs_discr_amr_indicator_type_t) init_amr_visc_indicator_type;
  discr_options->init_amr_visc_tol_min = init_amr_visc_tol_min;
  discr_options->init_amr_visc_tol_max = init_amr_visc_tol_max;
  discr_options->init_amr_visc_inside_plates_tol_min =
    init_amr_visc_in_plates_tol_min;
  discr_options->init_amr_visc_inside_plates_tol_max =
    init_amr_visc_in_plates_tol_max;

  discr_options->init_amr_weak_subdu_indicator_type =
    (slabs_discr_amr_indicator_type_t) init_amr_weak_subdu_indicator_type;
  discr_options->init_amr_weak_subdu_tol_min = init_amr_weak_subdu_tol_min;
  discr_options->init_amr_weak_subdu_tol_max = init_amr_weak_subdu_tol_max;
  discr_options->init_amr_weak_subdu_maxlevel = (int8_t)
    round (slabs_discr_resolution_to_level (init_amr_weak_subdu_elem_res,
                                            physics_options));
  discr_options->init_amr_weak_subdu_maxdist =
    0.5 * weakzone_2pl_subdu_thickness;

  discr_options->init_amr_weak_ridge_indicator_type =
    (slabs_discr_amr_indicator_type_t) init_amr_weak_ridge_indicator_type;
  discr_options->init_amr_weak_ridge_tol_min = init_amr_weak_ridge_tol_min;
  discr_options->init_amr_weak_ridge_tol_max = init_amr_weak_ridge_tol_max;
  discr_options->init_amr_weak_ridge_maxlevel = (int8_t)
    round (slabs_discr_resolution_to_level (init_amr_weak_ridge_elem_res,
                                            physics_options));
  discr_options->init_amr_weak_ridge_maxdist =
    0.5 * weakzone_2pl_ridge_smoothwidth;

  discr_options->init_amr_weak_import_indicator_type =
    (slabs_discr_amr_indicator_type_t) init_amr_weak_import_indicator_type;
  discr_options->init_amr_weak_import_tol_min = init_amr_weak_import_tol_min;
  discr_options->init_amr_weak_import_tol_max = init_amr_weak_import_tol_max;
  discr_options->init_amr_weak_import_maxlevel = (int8_t)
    round (slabs_discr_resolution_to_level (init_amr_weak_import_elem_res,
                                            physics_options));
  discr_options->init_amr_weak_import_maxdist =
    0.5 * weakzone_import_thickness / SL_EARTH_RADIUS;

  discr_options->init_amr_rhs_indicator_type =
    (slabs_discr_amr_indicator_type_t) init_amr_rhs_indicator_type;
  discr_options->init_amr_rhs_tol_min = init_amr_rhs_tol_min;
  discr_options->init_amr_rhs_tol_max = init_amr_rhs_tol_max;
  discr_options->init_amr_rhs_norm_shift = init_amr_rhs_norm_shift;

  discr_options->init_amr_post_uniform_n_steps = init_amr_post_uniform_n_steps;

  discr_options->amr_max_steps = amr_max_steps;
  discr_options->amr_rel_threshold = amr_rel_threshold;
  discr_options->amr_n_elements_max = amr_n_elements_max;
  discr_options->amr_lower_mantle = amr_lower_mantle;

  discr_options->amr_visc_indicator_type =
    (slabs_discr_amr_indicator_type_t) amr_visc_indicator_type;
  discr_options->amr_visc_tol_min = amr_visc_tol_min;
  discr_options->amr_visc_tol_max = amr_visc_tol_max;

  discr_options->amr_visc_dr_indicator_type =
    (slabs_discr_amr_indicator_type_t) amr_visc_dr_indicator_type;
  discr_options->amr_visc_dr_tol_min = amr_visc_dr_tol_min;
  discr_options->amr_visc_dr_tol_max = amr_visc_dr_tol_max;

  discr_options->amr_strain_rate_indicator_type =
    (slabs_discr_amr_indicator_type_t) amr_strain_rate_indicator_type;
  discr_options->amr_strain_rate_tol_min = amr_strain_rate_tol_min;
  discr_options->amr_strain_rate_tol_max = amr_strain_rate_tol_max;

  discr_options->amr_log_maxlevel = amr_log_maxlevel;
  discr_options->mesh_partitioning_type =
    (slabs_mesh_partitioning_type_t) mesh_partitioning_type;
  discr_options->inspect_p4est = 0;

  /*
   * Physics Options That Depend on Discretization
   */

  /* align weak zone in left corner with element boundaries */
  if (   weakzone_2pl_ridge_mesh_align
      && (weakzone_type == SL_WEAKZONE_2PLATES_POLY2 ||
          weakzone_type == SL_WEAKZONE_2PLATES_RHEA1) ) {
    slabs_discr_align_weak_ridge_2plates (physics_options, discr_options);

    YMIR_GLOBAL_INFOF ("%s: Depth and width of weak zone in corner "
                       "(at ridge) aligned with elements, new size: "
                       "depth %g km, width %g km\n",
                       this_fn_name,
                       physics_options->weakzone_2plates_ridge_depth / 1.0e3,
                       physics_options->weakzone_2plates_ridge_width / 1.0e3);
  }

  /* set radius of interface between lower mantle and upper mantle */
  if (visc_lower_mantle) { /* if lower mantle enabled */
    /* find radius for upper mantle s.t. interface is at element boundaries */
    slabs_discr_set_upper_mantle_radius (physics_options, discr_options);

    YMIR_GLOBAL_INFOF ("%s: Depth of upper mantle set to %g km\n", this_fn_name,
        (SL_SHELL_RADIUS_TOP - physics_options->viscosity_upper_mantle_radius)
        * SL_EARTH_RADIUS / 1.0e3);
  }
  else { /* if lower mantle disabled */
    physics_options->viscosity_upper_mantle_radius = -1.0;
    discr_options->refine_layer_radius = -1.0;
  }

  /*
   * Solver Options
   */

  /* set solver type */
  switch (viscosity_type) {
  case SL_VISCOSITY_CONST:
  case SL_VISCOSITY_LINEAR:
    solver_options->nl_solver_type = SL_NL_SOLVER_NONE;
    break;

  case SL_VISCOSITY_NONLINEAR:
    solver_options->nl_solver_type = SL_NL_SOLVER_NEWTON;
    if (nl_solver_type != SL_NL_SOLVER_NONE) {
      solver_options->nl_solver_type = (slabs_nl_solver_type_t) nl_solver_type;
    }
    break;

  default: /* unknown viscosity type */
    YMIR_ABORT_NOT_REACHED ();
  }
  solver_options->nl_solver_primaldual_type =
    (slabs_nl_solver_primaldual_type_t) nl_solver_primaldual_type;
  solver_options->nl_solver_primaldual_scal_type =
    (slabs_nl_solver_primaldual_scal_type_t) nl_solver_primaldual_scal_type;

  solver_options->nl_maxiter = nl_maxiter;
  solver_options->nl_rtol = nl_rtol;
  solver_options->initial_guess_type =
    (slabs_nl_solver_initial_guess_t) nl_initial_guess_type;
  solver_options->nl_resume_at_iter = nl_resume_at_iter;
  solver_options->nl_resume_at_time = nl_resume_at_time;
  solver_options->nl_resume_prev_res = nl_resume_prev_res;
  solver_options->nl_resume_init_res = nl_resume_init_res;
  solver_options->nl_forcing_exponent = nl_forcing_exponent;
  solver_options->nl_forcing_max = nl_forcing_max;
  solver_options->nl_forcing_max_progressive_iter =
    nl_forcing_max_progressive_iter;
  solver_options->nl_forcing_total_min = nl_forcing_total_min;
  solver_options->nl_forcing_saveguard = nl_forcing_saveguard;
  solver_options->nl_forcing_saveguard_threshold = nl_forcing_threshold;
  solver_options->nl_step_length_reduction_type =
    (slabs_nl_solver_step_reduction_type_t) nl_step_length_reduction_type;
  solver_options->nl_step_length_reduction_min = nl_step_length_reduction_min;
  solver_options->nl_step_length_reduction_max = nl_step_length_reduction_max;
  solver_options->nl_step_length_reduction_reg = nl_step_length_reduction_reg;
  solver_options->nl_step_length_descend_cond_relax =
    nl_step_length_descend_relax;
  solver_options->nl_step_length_min = nl_step_length_min;
  solver_options->nl_switch_picard_step_length_min =
    nl_switch_picard_step_length_min;
  solver_options->nl_switch_picard_after_amr = nl_switch_picard_after_amr;
  solver_options->nl_switch_picard_init = nl_switch_picard_init;
  solver_options->nl_switch_picard_maxiter = nl_switch_picard_maxiter;
  solver_options->nl_switch_picard_rtol = nl_switch_picard_rtol;

  solver_options->norm_type = (slabs_norm_type_t) nl_norm_type;
  solver_options->norm_Hminus1_mass_scaling = nl_norm_Hminus1_mass_scaling;

  solver_options->schur_diag_type =
    (slabs_nl_stokes_problem_schur_diag_t) nl_schur_diag_type;
  solver_options->scaling_type =
    (slabs_nl_stokes_problem_scaling_t) nl_scaling_type;
  solver_options->project_out_nullspace =
    (slabs_nl_solver_project_nullspace_t) nl_project_out_nullspace;
  solver_options->enforce_unscaled_reduction = nl_enforce_unscaled_reduction;

  solver_options->grid_continuation_init_amr_threshold =
    nl_grid_cont_init_amr_threshold;
  solver_options->grid_continuation_init_forcing = nl_grid_cont_init_forcing;
  solver_options->grid_continuation_init_forcing_exp =
    nl_grid_cont_init_forcing_exp;
  solver_options->grid_continuation_init_steps = nl_grid_cont_init_steps;
  solver_options->grid_continuation_skipsteps = nl_grid_cont_skipsteps;
  solver_options->grid_continuation_maxsteps = nl_grid_cont_maxsteps;
  solver_options->viscosity_bounds_continuation_min = nl_visc_bounds_cont_min;
  solver_options->viscosity_bounds_continuation_max = nl_visc_bounds_cont_max;
  solver_options->viscosity_bounds_continuation_steps =
    nl_visc_bounds_cont_steps;

  solver_options->nl_check_derivative = nl_check_derivative;
  solver_options->log_physics_stats = log_physics_stats;

  solver_options->krylov_type = (slabs_krylov_type_t) krylov_type;
  solver_options->krylov_gmres_num_vecs = gmres_num_vecs;
  solver_options->krylov_maxiter = krylov_maxiter;
  solver_options->krylov_rtol = krylov_rtol;
  solver_options->krylov_rtol_init = krylov_rtol;
  solver_options->krylov_atol = krylov_atol;

  solver_options->lin_solve_stress_block_only = lin_solve_stress_block_only;
  solver_options->lin_solve_press_bbt_only = lin_solve_press_bbt_only;
  solver_options->lin_solve_cnode_bbt_only = lin_solve_cnode_bbt_only;
}

/**
 * Initializes the weak zone.
 */
static void
slabs_setup_mesh_init_weakzone (slabs_stokes_state_t *state,
                                ymir_mesh_t *mesh,
                                slabs_physics_options_t *physics_options)
{
  MPI_Comm            mpicomm = mesh->ma->mpicomm;

  /* initialize weak zone in Stokes state */
  slabs_stokes_state_init_weakzone (state, mesh,
                                    slabs_physics_compute_weakzone,
                                    physics_options);

  /* initialization of weak zone computation */
  slabs_physics_init_weakzone (mpicomm, physics_options);

  /* fill Stokes state with weak zone factor values */
  slabs_physics_compute_weakzone (state->weak_vec, physics_options);
}

/**
 * Changes the weak zone field to a new mesh and recomputes it.
 */
static void
slabs_setup_mesh_recompute_weakzone (slabs_stokes_state_t *state,
                                     ymir_mesh_t *mesh,
                                     slabs_physics_options_t *physics_options)
{
  /* destroy weak zone in Stokes state */
  slabs_stokes_state_clear_weakzone (state);

  /* initialize weak zone in Stokes state */
  slabs_stokes_state_init_weakzone (state, mesh,
                                    slabs_physics_compute_weakzone,
                                    physics_options);

  /* fill Stokes state with weak zone factor values */
  slabs_physics_compute_weakzone (state->weak_vec, physics_options);
}

/**
 * Sets up the mesh.
 */
void
slabs_setup_mesh (p8est_t **p8est,
                  ymir_mesh_t **mesh,
                  ymir_pressure_elem_t **press_elem,
                  slabs_stokes_state_t **state,
                  slabs_discr_enforce_refinement_data_t
                    **enforce_refinement_data,
                  slabs_physics_coarsen_stokes_coeff_data_t
                    **coarsen_coeff_data,
                  MPI_Comm mpicomm,
                  slabs_physics_options_t *physics_options,
                  slabs_discr_options_t *discr_options,
                  slabs_nl_solver_options_t *solver_options,
                  ymir_perf_counter_t * perf_counter,
                  const char *workload_filepath)
{
  const char         *this_fn_name = "slabs_setup_mesh";

  const int           import_p4est =
    (physics_options->p4est_import_filename != NULL);
  const int           import_temp_data =
    (physics_options->temperature_type == SL_TEMP_IMPORT_FILE);
  const int           import_vel_press_data =
    (physics_options->velocity_import_filename != NULL &&
     physics_options->pressure_import_filename != NULL);
  const double        visc_min = physics_options->viscosity_min;
  const double        visc_max = physics_options->viscosity_max;
  const double        visc_bounds_cont_min =
                        solver_options->viscosity_bounds_continuation_min;
  const double        visc_bounds_cont_max =
                        solver_options->viscosity_bounds_continuation_max;
  const int           visc_bounds_cont_steps =
                        solver_options->viscosity_bounds_continuation_steps;

  const int           order = discr_options->order;
  const int           minlevel = discr_options->minlevel;
  int                 order_import, minlevel_import;
  int                 order_amr;

  const int           init_amr_visc_indicator_type =
                        discr_options->init_amr_visc_indicator_type;
  const double        init_amr_visc_tol_min =
                        discr_options->init_amr_visc_tol_min;
  const double        init_amr_visc_tol_max =
                        discr_options->init_amr_visc_tol_max;
  const int           init_amr_weak_import_indicator_type =
                        discr_options->init_amr_weak_import_indicator_type;
  const double        init_amr_weak_import_tol_min =
                        discr_options->init_amr_weak_import_tol_min;
  const double        init_amr_weak_import_tol_max =
                        discr_options->init_amr_weak_import_tol_max;
  const int           init_amr_rhs_indicator_type =
                        discr_options->init_amr_rhs_indicator_type;
  const double        init_amr_rhs_tol_min =
                        discr_options->init_amr_rhs_tol_min;
  const double        init_amr_rhs_tol_max =
                        discr_options->init_amr_rhs_tol_max;
  const int           init_amr_post_uniform_n_steps =
                        discr_options->init_amr_post_uniform_n_steps;

  mangll_t           *mangll_import, *mangll_amr;
  mangll_cnodes_t    *cnodes_import, *cnodes_amr;
  ymir_mesh_t        *mesh_import, *mesh_amr;
  ymir_pressure_elem_t  *press_elem_import, *press_elem_amr;
  char                path[BUFSIZ];

  YMIR_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* check input */
  YMIR_ASSERT (import_temp_data ||
               physics_options->weakzone_type != SL_WEAKZONE_IMPORT_FILE);

  /* set mesh parameters that differ from the computational mesh */
  if (0 < discr_options->import_mesh_order) {
    order_import = discr_options->import_mesh_order;
  }
  else {
    order_import = order;
  }
  if (0 < discr_options->import_mesh_minlevel) {
    minlevel_import = discr_options->import_mesh_minlevel;
  }
  else {
    minlevel_import = minlevel;
  }
  if (0 < discr_options->init_amr_override_order) {
    order_amr = discr_options->init_amr_override_order;
  }
  else {
    order_amr = order;
  }

  /* start performance counters */
  if (perf_counter != NULL) {
    ymir_perf_counter_start_barrier (perf_counter, mpicomm);
  }

  /*
   * Create p4est
   */

  YMIR_GLOBAL_INFOF ("%s: Generate p4est mesh\n", this_fn_name);

  if (!import_p4est) {
    discr_options->minlevel = minlevel_import;

    /* create p4est mesh and create Stokes state */
    *p8est = slabs_discr_p8est_new (mpicomm, physics_options, discr_options);
    *state = slabs_stokes_state_new (*p8est);

    discr_options->minlevel = minlevel;
  }
  else {
    /* load p4est mesh from file and create Stokes state */
    *state = slabs_stokes_state_load_init (
        p8est, physics_options->p4est_import_filename, mpicomm,
        discr_options->inspect_p4est);
  }

  /* set boundary information in `discr_options` */
  slabs_discr_options_set_boundary (discr_options, *p8est, physics_options);

  /* print memory usage */
  ymir_monitor_print_global_mem_usage (mpicomm);

  /*
   * Write Coordinates If Requested
   */

  if (temp_import_write_coord_path != NULL) {
    slabs_io_coordinate_type_t  coord_type;
    mangll_t           *mangll;
    mangll_cnodes_t    *cnodes;

    /* change mesh parameters */
    if (order_import != order) {
      discr_options->order = order_import;
      YMIR_GLOBAL_INFOF ("%s: Override polynomial discr order for "
                         "coordinates output (%i instead of %i)\n",
                         this_fn_name, order_import, order);
    }
    if (minlevel_import != minlevel) {
      discr_options->minlevel = minlevel_import;
      YMIR_GLOBAL_INFOF ("%s: Override min level "
                         "coordinates output (%i instead of %i)\n",
                         this_fn_name, minlevel_import, minlevel);
    }

    /* refine mesh radially or in z-direction */
    slabs_discr_refine_wrt_depth (*p8est, physics_options, discr_options);

    /* create mangll and cnodes */
    slabs_discr_mangll_and_cnodes_new (&mangll, &cnodes, *p8est, discr_options);

    /* create ymir mesh */
    slabs_discr_ymir_new (mesh, NULL, mangll, cnodes, discr_options);

    /* set domain dependent type of coordinates */
    switch (physics_options->domain_shape) {
    case SL_DOMAIN_CUBE:
    case SL_DOMAIN_BRICK:
      coord_type = SL_CARTESIAN_COORDINATE;
      break;

    case SL_DOMAIN_SHELL_CHUNK:
    case SL_DOMAIN_SHELL_SLICE:
    case SL_DOMAIN_SHELL:
      coord_type = SL_SPHERICAL_COORDINATE_GEO_CONV;
      break;

    default: /* unknown domain type */
      YMIR_ABORT_NOT_REACHED ();
    }

    /* write node coordinates */
    slabs_io_write_node_coordinates_to_textfile (
        temp_import_write_coord_path, *mesh, SL_GLL_CONTINUOUS_NODE,
        coord_type);

    /* write mesh */
#ifdef YMIR_DEBUG
    {
      ymir_vec_t         *dummyvec = ymir_cvec_new_zero (*mesh, 1);

      snprintf (path, BUFSIZ, "%s_mesh", temp_import_write_coord_path);
      ymir_vtk_write (*mesh, path, dummyvec, "dummy", NULL);

      ymir_vec_destroy (dummyvec);
    }
#endif

    /* destroy */
    slabs_discr_ymir_mangll_destroy (*mesh, NULL);

    /* restore mesh parameters for computation */
    discr_options->order = order;
    discr_options->minlevel = minlevel;
  }

  /*
   * Prepare Mesh Generation
   */

  /* set viscosity bounds for bounds continuation */
  if (0 < visc_bounds_cont_steps) {
    if (0.0 < visc_bounds_cont_min && visc_min < visc_bounds_cont_min) {
      physics_options->viscosity_min = visc_bounds_cont_min;
    }
    if (0.0 < visc_bounds_cont_max && visc_bounds_cont_max < visc_max) {
      physics_options->viscosity_max = visc_bounds_cont_max;
    }
  }

  /*
   * Mesh Generation in Case of No Data Input
   */

  if (!import_temp_data) {
    mangll_t           *mangll;
    mangll_cnodes_t    *cnodes;

    /* if no temperature data import, i.e., input comes from functions */
    YMIR_GLOBAL_INFOF ("%s: Generate adaptively refined mesh "
                       "(viscosity from function)\n", this_fn_name);

    /* refine p8est mesh adaptively without interpolating input */
    slabs_discr_initial_amr_no_interp (*p8est, physics_options,
                                       discr_options);

    /* create mangll and cnodes */
    slabs_discr_mangll_and_cnodes_new (&mangll, &cnodes, *p8est, discr_options);

    /* create ymir mesh and pressure element */
    slabs_discr_ymir_new (mesh, press_elem, mangll, cnodes, discr_options);

    YMIR_GLOBAL_INFOF ("%s: Initialize temperature\n", this_fn_name);

    /* initialize temperature in Stokes state */
    slabs_stokes_state_init_temp (*state, cnodes);
    slabs_stokes_state_init_temp_vec (*state, *mesh);

    /* fill Stokes state with temperature values from a function */
    ymir_cvec_set_function ((*state)->temp_vec,
                            slabs_physics_temperature_set_fn, physics_options);
  }

  /*
   * Temperature Data Import
   */

  if (import_temp_data) { /* if temperature data is imported */
    discr_options->order = order_import;
    discr_options->minlevel = minlevel_import;
    YMIR_GLOBAL_INFOF ("%s: Initialize temperature "
                       "(polynomial order %i, min level %i)\n",
                       this_fn_name, order_import, minlevel_import);

    /* refine mesh radially or in z-direction; restore min level */
    if (!import_p4est) {
      slabs_discr_refine_wrt_depth (*p8est, physics_options, discr_options);
    }
    discr_options->minlevel = minlevel;

    /* create mangll, cnodes, ymir mesh associated with import discretization */
    slabs_discr_mangll_and_cnodes_new (&mangll_import, &cnodes_import, *p8est,
                                       discr_options);
    slabs_discr_ymir_new (&mesh_import, &press_elem_import, mangll_import,
                          cnodes_import, discr_options);

    /* load temperature from file */
    slabs_stokes_state_load_temp (*state, mesh_import,
                                  physics_options->temp_import_filename_txt,
                                  physics_options->temp_import_filename_bin);

    /* do not read txt file any more, since content was copied to bin file */
    physics_options->temp_import_filename_txt = NULL;

    /* load velocity & pressure from file */
    if (import_vel_press_data) {
      slabs_stokes_state_load_vel_press (
          *state, mesh_import, press_elem_import,
          physics_options->velocity_import_filename,
          physics_options->pressure_import_filename);
    }

    /* create ymir mesh and interpolate temperature if required */
    //TODO create seperate function for this
    if (order_import != order_amr) { /* if interp. to AMR order required */
      sc_dmatrix_t       *temp_import;
      ymir_vec_t         *temp_import_vec;

      discr_options->order = order_amr;
      YMIR_GLOBAL_INFOF ("%s: Interpolate temperature for AMR to "
                         "polynomial order %i\n", this_fn_name, order_amr);


      /* create mangll, cnodes, ymir mesh associated with AMR discretization */
      slabs_discr_mangll_and_cnodes_new (&mangll_amr, &cnodes_amr, *p8est,
                                         discr_options);
      slabs_discr_ymir_new (&mesh_amr, &press_elem_amr, mangll_amr, cnodes_amr,
                            discr_options);

      /* copy temperature field from Stokes state */
      temp_import = sc_dmatrix_clone ((*state)->temperature);
      temp_import_vec = ymir_cvec_new_data (mesh_import, 1, temp_import);

      /* verify temperature data */
      slabs_physics_verify_temperature (temp_import_vec, physics_options);

      /* initialize new temperature in Stokes state */
      slabs_stokes_state_clear_temp (*state);
      slabs_stokes_state_init_temp (*state, cnodes_amr);
      slabs_stokes_state_init_temp_vec (*state, mesh_amr);

      /* interpolate temperature field */
      if (order_import < order_amr) {
        ymir_gmg_intergrid_p_interpolate_cnode_simple (temp_import_vec,
                                                       (*state)->temp_vec);
      }
      else {
        ymir_gmg_intergrid_p_restrict_cnode_simple (temp_import_vec,
                                                    (*state)->temp_vec);
      }

      /* interpolate velocity & pressure fields */
      if (import_vel_press_data) {
        //TODO interpolate vel & pressure
        YMIR_ABORT_NOT_REACHED ();
      }

      /* destroy import variables */
      slabs_discr_ymir_mangll_destroy (mesh_import, press_elem_import);
      sc_dmatrix_destroy (temp_import);
      ymir_vec_destroy (temp_import_vec);
    }
    else { /* if no change in polynomial order for AMR */
      /* keep mangll and cnodes structures */
      mangll_amr = mangll_import;
      cnodes_amr = cnodes_import;
      mesh_amr = mesh_import;
      press_elem_amr = press_elem_import;

      /* verify temperature data */
      slabs_physics_verify_temperature ((*state)->temp_vec, physics_options);
    }

    /* set mesh for weak zone initialization before AMR */
    *mesh = mesh_amr;

    /* print memory usage */
    ymir_monitor_print_global_mem_usage (mpicomm);

    /* ###DEV### write temp */
#ifdef YMIR_DEBUG
    /*
    slabs_io_write_cvec_to_textfile ("shell_temp_cont", (*state)->temp_vec,
                                     SL_SPHERICAL_COORDINATE_GEO_CONV);
    if (vtk_filepath != NULL && vtk_input) {
      snprintf (path, BUFSIZ, "%s_temp", vtk_filepath);
      ymir_vtk_write (*mesh, path, (*state)->temp_vec, "temp_cont", NULL);
    }
    */
#endif
  }

  /*
   * Weak Zone Initialization
   */

  YMIR_GLOBAL_INFOF ("%s: Initialize weak zone\n", this_fn_name);

  /* initialize weak zone data (required before initial AMR) */
  slabs_setup_mesh_init_weakzone (*state, *mesh, physics_options);

  /* print memory usage */
  ymir_monitor_print_global_mem_usage (mpicomm);

  /*
   * Mesh Generation in Case of Data Input
   */

  if (import_temp_data) { /* if temperature data is imported */
    YMIR_GLOBAL_INFOF ("%s: Generate adaptively refined mesh "
                       "(viscosity from data)\n", this_fn_name);

    /* write global workload to a file */
    if (workload_filepath != NULL) {
      snprintf (path, BUFSIZ, "%s_setup_mesh_begin", workload_filepath);
      ymir_monitor_dump_global_mesh_mem_stats (mesh_amr, press_elem_amr, path);
      ymir_monitor_dump_global_mem_usage_petsc (mpicomm, path);
    }

    /* initial AMR with interpolation of the temperature field */
    if ( (init_amr_visc_indicator_type != SL_AMR_INDICATOR_NONE &&
          init_amr_visc_indicator_type <= SL_AMR_INDICATOR_VISC_PECLET &&
          0.0 < init_amr_visc_tol_max)
         ||
         (init_amr_weak_import_indicator_type != SL_AMR_INDICATOR_NONE)
         ||
         (init_amr_rhs_indicator_type != SL_AMR_INDICATOR_NONE &&
          0.0 < init_amr_rhs_tol_max)
         ||
         (0.0 < discr_options->refine_layer_radius &&
          0.0 < discr_options->refine_layer_maxdist)
         ||
         (0.0 < discr_options->refine_surface_maxdist)
       ) {
      const double        init_amr_n_quadrants_max =
                            discr_options->init_amr_n_elements_max;
      slabs_discr_amr_indicator_params_t  *indicator_params;

      /* create AMR parameters */
      indicator_params = slabs_discr_amr_indicator_params_new (5);

      /* set AMR parameters */
      indicator_params->type[0] =
        (slabs_discr_amr_indicator_type_t) init_amr_visc_indicator_type;
      indicator_params->tol_min[0] = init_amr_visc_tol_min;
      indicator_params->tol_max[0] = init_amr_visc_tol_max;
      indicator_params->level_min[0] = 0;
      indicator_params->level_max[0] = 0;

      indicator_params->type[1] =
        (slabs_discr_amr_indicator_type_t) init_amr_weak_import_indicator_type;
      indicator_params->tol_min[1] = init_amr_weak_import_tol_min;
      indicator_params->tol_max[1] = init_amr_weak_import_tol_max;
      indicator_params->level_min[1] = 0;
      indicator_params->level_max[1] =
        discr_options->init_amr_weak_import_maxlevel;

      indicator_params->type[2] =
        (slabs_discr_amr_indicator_type_t) init_amr_rhs_indicator_type;
      indicator_params->tol_min[2] = init_amr_rhs_tol_min;
      indicator_params->tol_max[2] = init_amr_rhs_tol_max;
      indicator_params->level_min[2] = 0;
      indicator_params->level_max[2] = 0;

      if (0.0 < discr_options->refine_layer_radius &&
          0.0 < discr_options->refine_layer_maxdist) {
        indicator_params->type[3] = SL_AMR_INDICATOR_REFINE_LAYER;
      }
      else {
        indicator_params->type[3] = SL_AMR_INDICATOR_NONE;
      }
      indicator_params->tol_min[3] = 0.0;
      indicator_params->tol_max[3] = 0.5;
      indicator_params->level_min[3] = 0;
      indicator_params->level_max[3] = discr_options->refine_layer_maxlevel;

      if (0.0 < discr_options->refine_surface_maxdist) {
        indicator_params->type[4] = SL_AMR_INDICATOR_REFINE_SURFACE;
      }
      else {
        indicator_params->type[4] = SL_AMR_INDICATOR_NONE;
      }
      indicator_params->tol_min[4] = 0.0;
      indicator_params->tol_max[4] = 0.5;
      indicator_params->level_min[4] = 0;
      indicator_params->level_max[4] = discr_options->refine_surface_maxlevel;

      /* refine mesh */
      if (0 < init_amr_post_uniform_n_steps) {
        discr_options->init_amr_n_elements_max =
          init_amr_n_quadrants_max /
          pow (8.0, (double) init_amr_post_uniform_n_steps);
      }
      slabs_discr_amr (*state, &mesh_amr, &press_elem_amr, *p8est,
                       indicator_params, physics_options, discr_options,
                       !import_vel_press_data /* initial AMR */);
      if (0 < init_amr_post_uniform_n_steps) {
        discr_options->init_amr_n_elements_max = init_amr_n_quadrants_max;
      }

      /* destroy AMR parameters */
      slabs_discr_amr_indicator_params_destroy (indicator_params);
    }
    else { /* if parameters not set for AMR */
      /* postprocess temperature in Stokes state */
      slabs_physics_postprocess_temperature (*state, physics_options);

      /* print mesh statistics */
      ymir_monitor_print_global_mesh_stats (mesh_amr, press_elem_amr);
    }

    /* uniform refinement after AMR */
    if (0 < init_amr_post_uniform_n_steps) {
      const int           init_amr_max_steps =
                            discr_options->init_amr_max_steps;
      const double        init_amr_n_quadrants_max =
                            discr_options->init_amr_n_elements_max;
      const int           amr_max_steps = discr_options->amr_max_steps;
      const double        amr_n_quadrants_max =
                            discr_options->amr_n_elements_max;
      slabs_discr_amr_indicator_params_t  *indicator_params;

      /* create AMR parameters */
      indicator_params = slabs_discr_amr_indicator_params_new (1);

      /* set AMR parameters */
      indicator_params->type[0] =
        (slabs_discr_amr_indicator_type_t) SL_AMR_INDICATOR_REFINE_ALL;
      indicator_params->tol_min[0] = 0.0;
      indicator_params->tol_max[0] = 0.5;
      indicator_params->level_min[0] = 0;
      indicator_params->level_max[0] = 0;
      discr_options->init_amr_max_steps = init_amr_post_uniform_n_steps;
      discr_options->init_amr_n_elements_max = -1.0;
      discr_options->amr_max_steps = init_amr_post_uniform_n_steps;
      discr_options->amr_n_elements_max = -1.0;

      /* refine mesh */
      slabs_discr_amr (*state, &mesh_amr, &press_elem_amr, *p8est,
                       indicator_params, physics_options, discr_options,
                       !import_vel_press_data /* initial AMR */);

      /* destroy AMR parameters; restore old parameters */
      slabs_discr_amr_indicator_params_destroy (indicator_params);
      discr_options->init_amr_max_steps = init_amr_max_steps;
      discr_options->init_amr_n_elements_max = init_amr_n_quadrants_max;
      discr_options->amr_max_steps = amr_max_steps;
      discr_options->amr_n_elements_max = amr_n_quadrants_max;
    }

    /* interpolate temperature if required */
    if (order_amr != order) { /* if interp. to "final" order required */
      mangll_t           *mangll;
      mangll_cnodes_t    *cnodes;
      sc_dmatrix_t       *temp_amr;
      ymir_vec_t         *temp_amr_vec;

      discr_options->order = order;
      YMIR_GLOBAL_INFOF ("%s: Interpolate temperature to "
                         "polynomial order %i\n", this_fn_name, order);

      /* create mangll, cnodes, ymir mesh with "final" order */
      slabs_discr_mangll_and_cnodes_new (&mangll, &cnodes, *p8est,
                                         discr_options);
      slabs_discr_ymir_new (mesh, press_elem, mangll, cnodes, discr_options);

      /* copy temperature field from Stokes state */
      temp_amr = sc_dmatrix_clone ((*state)->temperature);
      temp_amr_vec = ymir_cvec_new_data (mesh_amr, 1, temp_amr);

      /* initialize new temperature in Stokes state */
      slabs_stokes_state_clear_temp (*state);
      slabs_stokes_state_init_temp (*state, cnodes);
      slabs_stokes_state_init_temp_vec (*state, *mesh);

      /* interpolate temperature field */
      if (order_amr < order) {
        ymir_gmg_intergrid_p_interpolate_cnode_simple (temp_amr_vec,
                                                       (*state)->temp_vec);
      }
      else {
        ymir_gmg_intergrid_p_restrict_cnode_simple (temp_amr_vec,
                                                    (*state)->temp_vec);
      }

      /* interpolate velocity & pressure fields */
      if (import_vel_press_data) {
        //TODO interpolate vel & pressure
        YMIR_ABORT_NOT_REACHED ();
      }

      /* destroy AMR variables */
      slabs_discr_ymir_mangll_destroy (mesh_amr, press_elem_amr);
      sc_dmatrix_destroy (temp_amr);
      ymir_vec_destroy (temp_amr_vec);

      /* postprocess temperature in Stokes state */
      slabs_physics_postprocess_temperature (*state, physics_options);

      /* recompute weak zone field */
      slabs_setup_mesh_recompute_weakzone (*state, *mesh, physics_options);

      /* print mesh statistics */
      ymir_monitor_print_global_mesh_stats (*mesh, *press_elem);
    }
    else { /* if no change in polynomial order */
      *mesh = mesh_amr;
      *press_elem = press_elem_amr;
    }
  }
  else { /* if temperature is not imported */
    /* print mesh statistics */
    ymir_monitor_print_global_mesh_stats (*mesh, *press_elem);
  }

  /*
   * Post-Process Mesh Generation
   */

  /* restore original viscosity bounds */
  physics_options->viscosity_min = visc_min;
  physics_options->viscosity_max = visc_max;

  /*
   * Initialize Geometric Multigrid
   */

  YMIR_GLOBAL_INFOF ("%s: Initialize geometric multigrid\n", this_fn_name);

  {
    ymir_gmg_init_mesh (*p8est, discr_options->e_to_fm_fn,
                        discr_options->tree_to_bf);

    *enforce_refinement_data = slabs_discr_enforce_refinement_data_new (
        physics_options, discr_options);
    ymir_gmg_init_enforce_refinement (
        slabs_discr_enforce_refinement_at_layer, *enforce_refinement_data);

    ymir_gmg_init_bc (slabs_set_dirichlet_bc, physics_options, NULL, NULL);

    if (physics_options->viscosity_coarsen_eval ||
        physics_options->viscosity_p_coarsen_eval) {
      *coarsen_coeff_data = slabs_physics_coarsen_stokes_coeff_data_new (
          *state, physics_options);
    }
    else {
      *coarsen_coeff_data = NULL;
    }
    if (physics_options->viscosity_coarsen_eval) {
      ymir_gmg_init_coarsen_coeff (
          slabs_physics_coarsen_stokes_coeff,
          slabs_physics_coarsen_stokes_coeff_data_reset, *coarsen_coeff_data);
    }
    if (physics_options->viscosity_p_coarsen_eval) {
      ymir_gmg_init_p_coarsen_coeff (
          slabs_physics_coarsen_stokes_coeff,
          slabs_physics_coarsen_stokes_coeff_data_reset, *coarsen_coeff_data);
    }
  }

  /* stop performance counters */
  if (perf_counter != NULL) {
    ymir_perf_counter_stop_add (perf_counter);
  }

  /* print memory usage */
  ymir_monitor_print_global_mem_usage (mpicomm);

  /* write global workload to a file */
  if (workload_filepath != NULL) {
    snprintf (path, BUFSIZ, "%s_setup_mesh_end", workload_filepath);
    ymir_monitor_dump_global_mesh_mem_stats (*mesh, *press_elem, path);
    ymir_monitor_dump_global_mem_usage_petsc (mpicomm, path);
  }

  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Sets up a linear or nonlinear Stokes problem.
 */
void
slabs_setup_stokes (slabs_lin_stokes_problem_t **lin_stokes,
                    slabs_nl_stokes_problem_t **nl_stokes,
                    p8est_t *p8est,
                    ymir_mesh_t *mesh,
                    ymir_pressure_elem_t *press_elem,
                    slabs_stokes_state_t *state,
                    slabs_physics_options_t *physics_options,
                    slabs_nl_solver_options_t *solver_options,
                    ymir_perf_counter_t * perf_counter,
                    const char *workload_filepath)
{
  const char         *this_fn_name = "slabs_setup_stokes";
  const int           import_vel_press_data =
    (physics_options->velocity_import_filename != NULL &&
     physics_options->pressure_import_filename != NULL);
  const slabs_nl_solver_type_t             nl_solver_type =
    solver_options->nl_solver_type;
  const slabs_nl_solver_primaldual_type_t  nl_solver_primaldual_type =
    solver_options->nl_solver_primaldual_type;
  MPI_Comm            mpicomm = p8est->mpicomm;
  char                path[BUFSIZ];

  YMIR_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

#if defined(__bgq__)

#if defined(__HAVE_MONEQ)
// MonEQ initialization
  int moneq_status = MonEQ_Initialize();
  if (moneq_status) printf("Error initializing MonEQ\n");

// MonEQ start
MonEQ_StartPowerTag("YMIR_SETUP_STOKES_Block");
  if(p8est->mpirank == 0)
    printf("MONEQ SETUP STOKES START\n");
#endif

#if defined(__HAVE_HPM)
  HPM_Start("YMIR_SETUP_STOKES_Block");
  if(p8est->mpirank == 0)
    printf("HPM SETUP STOKES START\n");
#endif

#endif

  /* start performance counters */
  if (perf_counter != NULL) {
    ymir_perf_counter_start_barrier (perf_counter, mpicomm);
  }

  /* initialize velocity and pressure in Stokes state */
  if (!import_vel_press_data) {
    slabs_stokes_state_init_vel_press (state, mesh, press_elem);
  }

  /* initialize dual tensor in Stokes state for primal-dual nonlinear solver */
  if (   nl_solver_type != SL_NL_SOLVER_NONE
      && nl_solver_primaldual_type != SL_NL_SOLVER_PRIMALDUAL_NONE ) {
    slabs_stokes_state_init_dual_tensor (state, mesh);
  }

  /* create Stokes problem */
  if (nl_solver_type == SL_NL_SOLVER_NONE) { /* if linear solve */
    /* clear allocation for weak zone computation */
    slabs_physics_clear_weakzone (physics_options);

    /* create linear Stokes problem */
    *lin_stokes = slabs_linear_stokes_problem_new (state, mesh, press_elem,
                                                   physics_options);
    *nl_stokes = NULL;
  }
  else {
    /* create nonlinear Stokes problem */
    *nl_stokes = slabs_nonlinear_stokes_problem_new (state, mesh, press_elem,
                                                     physics_options);
    *lin_stokes = NULL;
  }

  /* stop performance counters */
  if (perf_counter != NULL) {
    ymir_perf_counter_stop_add (perf_counter);
  }

#if defined(__bgq__)

#if defined(__HAVE_HPM)
  HPM_Stop("YMIR_SETUP_STOKES_Block");
  if(p8est->mpirank == 0)
    printf("HPM SETUP STOKES STOP\n");
#endif

#if defined(__HAVE_MONEQ)
  // Stop recording power
  MonEQ_EndPowerTag("YMIR_SETUP_STOKES_Block");
  if(p8est->mpirank == 0)
    printf("MONEQ SETUP STOKES STOP\n");
#endif

#endif

  /* print memory usage */
  ymir_monitor_print_global_mem_usage (mpicomm);

  /* write global workload to a file */
  if (workload_filepath != NULL) {
    snprintf (path, BUFSIZ, "%s_setup_stokes", workload_filepath);
    ymir_monitor_dump_global_mesh_mem_stats (mesh, press_elem, path);
    ymir_monitor_dump_global_mem_usage_petsc (mpicomm, path);
  }

  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

/**
 * Cleans Up.
 */
void
slabs_clear (slabs_lin_stokes_problem_t *lin_stokes,
             slabs_nl_stokes_problem_t *nl_stokes,
             p8est_t *p8est,
             ymir_mesh_t *mesh,
             ymir_pressure_elem_t *press_elem,
             slabs_stokes_state_t *state,
             slabs_discr_enforce_refinement_data_t *enforce_refinement_data,
             slabs_physics_coarsen_stokes_coeff_data_t *coarsen_coeff_data,
             slabs_physics_options_t *physics_options,
             slabs_discr_options_t *discr_options)
{
  const char         *this_fn_name = "slabs_clear";

  YMIR_GLOBAL_PRODUCTIONF ("Into %s\n", this_fn_name);

  /* destroy variables for geometric multigrid */
  slabs_discr_enforce_refinement_data_destroy (enforce_refinement_data);
  slabs_physics_coarsen_stokes_coeff_data_destroy (coarsen_coeff_data);

  /* destroy Stokes problem */
  if (lin_stokes != NULL) {
    slabs_linear_stokes_problem_destroy (lin_stokes);
    YMIR_ASSERT (nl_stokes == NULL);
  }
  if (nl_stokes != NULL) {
    slabs_nonlinear_stokes_problem_destroy (nl_stokes);
    YMIR_ASSERT (lin_stokes == NULL);

    /* clear allocation for weak zone computation */
    slabs_physics_clear_weakzone (physics_options);
  }

  /* destroy Stokes state */
  slabs_stokes_state_destroy (state);

  /* destroy mesh */
  slabs_discr_ymir_mangll_destroy (mesh, press_elem);
  slabs_discr_p8est_destroy (p8est);

  /* destroy options */
  slabs_discr_options_clear_boundary (discr_options);
  YMIR_FREE (discr_options->refine_radius);

  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);
}

