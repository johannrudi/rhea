#include <rhea_inversion.h>
#include <rhea_base.h>
#include <rhea_inversion_obs_stress.h>
#include <rhea_inversion_obs_velocity.h>
#include <rhea_inversion_obs_viscosity.h>
#include <rhea_inversion_param.h>
#include <rhea_inversion_qoi.h>
#include <rhea_newton.h>
#include <rhea_newton_check.h>
#include <rhea_velocity.h>
#include <rhea_velocity_pressure.h>
#include <rhea_strainrate.h>
#include <rhea_stress.h>
#include <rhea_stokes_problem_amr.h>
#include <rhea_io_std.h>
#include <rhea_io_mpi.h>
#include <rhea_vtk.h>
#include <ymir_comm.h>
#include <ymir_perf_counter.h>

/* Inverse solver */
#define RHEA_INVERSION_OBJ_N_COMPONENTS 4
#define RHEA_INVERSION_GRAD_NORM_N_COMPONENTS 2

/* I/O labels that are attached at the end of file names */
#define RHEA_INVERSION_IO_LABEL_NL_ITER "_itn"

/* types of Hessians */
typedef enum
{
  RHEA_INVERSION_HESSIAN_NONE             = -1,
  RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND = 0,
  RHEA_INVERSION_HESSIAN_BFGS,
  RHEA_INVERSION_HESSIAN_FIRST_ORDER_APPROX,
  RHEA_INVERSION_HESSIAN_FULL
}
rhea_inversion_hessian_t;

/* types of null space projections */
typedef enum
{
  RHEA_INVERSION_PROJECT_OUT_NULL_NONE     = 0,
  RHEA_INVERSION_PROJECT_OUT_NULL_PRIMAL   = 1,
  RHEA_INVERSION_PROJECT_OUT_NULL_RESIDUAL = 2,
  RHEA_INVERSION_PROJECT_OUT_NULL_BOTH     = 3
}
rhea_inversion_project_out_null_t;

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_INVERSION_DEFAULT_VEL_OBS_TYPE (RHEA_INVERSION_OBS_VELOCITY_NONE)
#define RHEA_INVERSION_DEFAULT_VEL_OBS_WEIGHT_TYPE \
  (RHEA_INVERSION_OBS_VELOCITY_WEIGHT_VALUES)
#define RHEA_INVERSION_DEFAULT_VEL_OBS_STDDEV_MM_YR NULL
#define RHEA_INVERSION_DEFAULT_VISC_OBS_TYPE (RHEA_INVERSION_OBS_VISCOSITY_NONE)
#define RHEA_INVERSION_DEFAULT_VISC_OBS_PAS NULL
#define RHEA_INVERSION_DEFAULT_VISC_OBS_STDDEV_REL NULL
#define RHEA_INVERSION_DEFAULT_STRESS_OBS_TYPE (RHEA_INVERSION_OBS_STRESS_NONE)
#define RHEA_INVERSION_DEFAULT_PARAMETER_PRIOR_STDDEV (1.0)
#define RHEA_INVERSION_DEFAULT_HESSIAN_TYPE (RHEA_INVERSION_HESSIAN_BFGS)
#define RHEA_INVERSION_DEFAULT_HESSIAN_TYPE_WRITE_AFTER_SOLVE \
  (RHEA_INVERSION_HESSIAN_FIRST_ORDER_APPROX)
#define RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_MATRIX (1)
#define RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_ENFORCE_SYMM (0)
#define RHEA_INVERSION_DEFAULT_GRADIENT_DESCENT_DAMPING (NAN)
#define RHEA_INVERSION_DEFAULT_RESTRICT_INIT_STEP_TO_PRIOR_STDDEV (-1.0)
#define RHEA_INVERSION_DEFAULT_RESTRICT_STEP_TO_PRIOR_STDDEV (2.0)
#define RHEA_INVERSION_DEFAULT_HESSIAN_PREV_FILE_PATH_TXT NULL
#define RHEA_INVERSION_DEFAULT_GRADIENT_PREV_FILE_PATH_TXT NULL
#define RHEA_INVERSION_DEFAULT_STEP_PREV_FILE_PATH_TXT NULL
#define RHEA_INVERSION_DEFAULT_INNER_SOLVER_RTOL_ADAPTIVE (0)
#define RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_ALLOW_AMR_FOR_OUTER_N_ITER (1)
#define RHEA_INVERSION_DEFAULT_PROJECT_OUT_NULL \
  (RHEA_INVERSION_PROJECT_OUT_NULL_NONE)
#define RHEA_INVERSION_DEFAULT_CHECK_GRADIENT (0)
#define RHEA_INVERSION_DEFAULT_CHECK_GRADIENT_OF_QOI (0)
#define RHEA_INVERSION_DEFAULT_CHECK_GRADIENT_ELEMENTWISE (0)
#define RHEA_INVERSION_DEFAULT_CHECK_HESSIAN (0)
#define RHEA_INVERSION_DEFAULT_DEBUG (0)
#define RHEA_INVERSION_DEFAULT_MONITOR_PERFORMANCE (0)

#define RHEA_INVERSION_DEFAULT_PRIOR_PRECONDITIONER (0)

/* initialize options */
int                 rhea_inversion_vel_obs_type =
                      RHEA_INVERSION_DEFAULT_VEL_OBS_TYPE;
int                 rhea_inversion_vel_obs_weight_type =
                      RHEA_INVERSION_DEFAULT_VEL_OBS_WEIGHT_TYPE;
char               *rhea_inversion_vel_obs_stddev_mm_yr =
                      RHEA_INVERSION_DEFAULT_VEL_OBS_STDDEV_MM_YR;
int                 rhea_inversion_visc_obs_type =
                      RHEA_INVERSION_DEFAULT_VISC_OBS_TYPE;
char               *rhea_inversion_visc_obs_Pas =
                      RHEA_INVERSION_DEFAULT_VISC_OBS_PAS;
char               *rhea_inversion_visc_obs_stddev_rel =
                      RHEA_INVERSION_DEFAULT_VISC_OBS_STDDEV_REL;
int                 rhea_inversion_stress_obs_type =
                      RHEA_INVERSION_DEFAULT_STRESS_OBS_TYPE;
double              rhea_inversion_parameter_prior_stddev =
                      RHEA_INVERSION_DEFAULT_PARAMETER_PRIOR_STDDEV;
int                 rhea_inversion_hessian_type =
                      RHEA_INVERSION_DEFAULT_HESSIAN_TYPE;
int                 rhea_inversion_hessian_type_write_after_solve =
                      RHEA_INVERSION_DEFAULT_HESSIAN_TYPE_WRITE_AFTER_SOLVE;
int                 rhea_inversion_assemble_hessian_matrix =
                      RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_MATRIX;
int                 rhea_inversion_assemble_hessian_enforce_symm =
                      RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_ENFORCE_SYMM;
double              rhea_inversion_gradient_descent_damping =
                      RHEA_INVERSION_DEFAULT_GRADIENT_DESCENT_DAMPING;
double              rhea_inversion_restrict_init_step_to_prior_stddev =
                      RHEA_INVERSION_DEFAULT_RESTRICT_INIT_STEP_TO_PRIOR_STDDEV;
double              rhea_inversion_restrict_step_to_prior_stddev =
                      RHEA_INVERSION_DEFAULT_RESTRICT_STEP_TO_PRIOR_STDDEV;
char               *rhea_inversion_hessian_prev_file_path_txt =
                      RHEA_INVERSION_DEFAULT_HESSIAN_PREV_FILE_PATH_TXT;
char               *rhea_inversion_gradient_prev_file_path_txt =
                      RHEA_INVERSION_DEFAULT_GRADIENT_PREV_FILE_PATH_TXT;
char               *rhea_inversion_step_prev_file_path_txt =
                      RHEA_INVERSION_DEFAULT_STEP_PREV_FILE_PATH_TXT;
int                 rhea_inversion_forward_solver_iter_max =
                      RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_ITER_MAX;
double              rhea_inversion_forward_solver_rtol =
                      RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_RTOL;
int                 rhea_inversion_adjoint_solver_iter_max =
                      RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_ITER_MAX;
double              rhea_inversion_adjoint_solver_rtol =
                      RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_RTOL;
int                 rhea_inversion_incremental_forward_solver_iter_max =
                      RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_ITER_MAX;
double              rhea_inversion_incremental_forward_solver_rtol =
                      RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_RTOL;
int                 rhea_inversion_incremental_adjoint_solver_iter_max =
                      RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_ITER_MAX;
double              rhea_inversion_incremental_adjoint_solver_rtol =
                      RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_RTOL;
int                 rhea_inversion_inner_solver_rtol_adaptive =
                      RHEA_INVERSION_DEFAULT_INNER_SOLVER_RTOL_ADAPTIVE;
int                 rhea_inversion_allow_amr_for_outer_n_iter =
                      RHEA_INVERSION_DEFAULT_ALLOW_AMR_FOR_OUTER_N_ITER;
int                 rhea_inversion_project_out_null =
                      RHEA_INVERSION_DEFAULT_PROJECT_OUT_NULL;
int                 rhea_inversion_check_gradient =
                      RHEA_INVERSION_DEFAULT_CHECK_GRADIENT;
int                 rhea_inversion_check_gradient_of_qoi =
                      RHEA_INVERSION_DEFAULT_CHECK_GRADIENT_OF_QOI;
int                 rhea_inversion_check_gradient_elementwise =
                      RHEA_INVERSION_DEFAULT_CHECK_GRADIENT_ELEMENTWISE;
int                 rhea_inversion_check_hessian =
                      RHEA_INVERSION_DEFAULT_CHECK_HESSIAN;
int                 rhea_inversion_debug =
                      RHEA_INVERSION_DEFAULT_DEBUG;
int                 rhea_inversion_monitor_performance =
                      RHEA_INVERSION_DEFAULT_MONITOR_PERFORMANCE;

rhea_inversion_param_options_t  rhea_inversion_param_options;
rhea_inversion_qoi_options_t    rhea_inversion_qoi_options;
rhea_newton_options_t           rhea_inversion_newton_options;

void
rhea_inversion_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Inversion";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /* data options */
  YMIR_OPTIONS_I, "velocity-observations-type", '\0',
    &(rhea_inversion_vel_obs_type), RHEA_INVERSION_DEFAULT_VEL_OBS_TYPE,
    "Type of velocity observations at the surface",
  YMIR_OPTIONS_I, "velocity-observations-weight-type", '\0',
    &(rhea_inversion_vel_obs_weight_type),
    RHEA_INVERSION_DEFAULT_VEL_OBS_WEIGHT_TYPE,
    "Velocity obs.: Type of weights",
  YMIR_OPTIONS_S, "velocity-observations-stddev-mm-yr", '\0',
    &(rhea_inversion_vel_obs_stddev_mm_yr),
    RHEA_INVERSION_DEFAULT_VEL_OBS_STDDEV_MM_YR,
    "Velocity obs.: List of standard deviations, one per plate in [mm/yr]",

  YMIR_OPTIONS_I, "viscosity-observations-type", '\0',
    &(rhea_inversion_visc_obs_type), RHEA_INVERSION_DEFAULT_VISC_OBS_TYPE,
    "Type of viscosity observations",
  YMIR_OPTIONS_S, "viscosity-observations-Pas", '\0',
    &(rhea_inversion_visc_obs_Pas), RHEA_INVERSION_DEFAULT_VISC_OBS_PAS,
    "Viscosity obs.: List of (average) viscosities [Pa*s]",
  YMIR_OPTIONS_S, "viscosity-observations-stddev-rel", '\0',
    &(rhea_inversion_visc_obs_stddev_rel),
    RHEA_INVERSION_DEFAULT_VISC_OBS_STDDEV_REL,
    "Viscosity obs.: List of standard deviations for log(avg. viscosity)",
//TODO
//YMIR_OPTIONS_S, "viscosity-observations-polygon-vertices-file-path-txt", '\0',
//  &(rhea_inversion_visc_obs_polygon_vertices_file_path_txt),
//  RHEA_INVERSION_DEFAULT_VISC_OBS_POLYGON_VERTICES_FILE_PATH_TXT,
//  "",
//YMIR_OPTIONS_I, "viscosity-observations-polygon-vertices-num-total", '\0',
//  &(rhea_inversion_visc_obs_polygon_vertices_n_total),
//  RHEA_INVERSION_DEFAULT_VISC_OBS_POLYGON_VERTICES_N_TOTAL,
//  "",

  YMIR_OPTIONS_I, "stress-observations-type", '\0',
    &(rhea_inversion_stress_obs_type), RHEA_INVERSION_DEFAULT_STRESS_OBS_TYPE,
    "Type of stress observations",

  /* parameter options */
  YMIR_OPTIONS_D, "parameter-prior-stddev", '\0',
    &(rhea_inversion_parameter_prior_stddev),
    RHEA_INVERSION_DEFAULT_PARAMETER_PRIOR_STDDEV,
    "Standard deviation of the prior term in the objective functional",

  /* outer solver options (inversion) */
  YMIR_OPTIONS_I, "hessian-type", '\0',
    &(rhea_inversion_hessian_type), RHEA_INVERSION_DEFAULT_HESSIAN_TYPE,
    "Type of Hessian operator: 0: full, 1: Gauss-Newton, 2: BFGS",
  YMIR_OPTIONS_I, "hessian-type-write-after-solve", '\0',
    &(rhea_inversion_hessian_type_write_after_solve),
    RHEA_INVERSION_DEFAULT_HESSIAN_TYPE_WRITE_AFTER_SOLVE,
    "Type of Hessian operator that is written to file after last iteration",
  YMIR_OPTIONS_B, "assemble-hessian", '\0',
    &(rhea_inversion_assemble_hessian_matrix),
    RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_MATRIX,
    "Toggle assembly of the Hessian matrix",
  YMIR_OPTIONS_B, "assemble-hessian-enforce-symmetry", '\0',
    &(rhea_inversion_assemble_hessian_enforce_symm),
    RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_ENFORCE_SYMM,
    "Enforce symmetry of the assembled Hessian matrix",
  YMIR_OPTIONS_D, "gradient-descent-damping", '\0',
    &(rhea_inversion_gradient_descent_damping),
    RHEA_INVERSION_DEFAULT_GRADIENT_DESCENT_DAMPING,
    "Factor by which the steps from gradient descent are damped",
  YMIR_OPTIONS_D, "restrict-init-step-to-prior-stddev", '\0',
    &(rhea_inversion_restrict_init_step_to_prior_stddev),
    RHEA_INVERSION_DEFAULT_RESTRICT_INIT_STEP_TO_PRIOR_STDDEV,
    "Restrict initial Newton step to lie within the range of prior",
  YMIR_OPTIONS_D, "restrict-step-to-prior-stddev", '\0',
    &(rhea_inversion_restrict_step_to_prior_stddev),
    RHEA_INVERSION_DEFAULT_RESTRICT_STEP_TO_PRIOR_STDDEV,
    "Restrict each Newton step to lie within the range of prior",
  YMIR_OPTIONS_S, "hessian-prev-file-path-txt", '\0',
    &(rhea_inversion_hessian_prev_file_path_txt),
    RHEA_INVERSION_DEFAULT_HESSIAN_PREV_FILE_PATH_TXT,
    "Resume inversion: Path to a text file of previous Hessian",
  YMIR_OPTIONS_S, "gradient-prev-file-path-txt", '\0',
    &(rhea_inversion_gradient_prev_file_path_txt),
    RHEA_INVERSION_DEFAULT_GRADIENT_PREV_FILE_PATH_TXT,
    "Resume inversion: Path to a text file of previous gradient",
  YMIR_OPTIONS_S, "step-prev-file-path-txt", '\0',
    &(rhea_inversion_step_prev_file_path_txt),
    RHEA_INVERSION_DEFAULT_STEP_PREV_FILE_PATH_TXT,
    "Resume inversion: Path to a text file of previous step",

  /* inner solver options (i.e., forward/adjoint) */
  YMIR_OPTIONS_B, "inner-solver-rtol-adaptive", '\0',
    &(rhea_inversion_inner_solver_rtol_adaptive),
    RHEA_INVERSION_DEFAULT_INNER_SOLVER_RTOL_ADAPTIVE,
    "Inner solves: Use adaptive relative tolerance provided by Newton's method",
  YMIR_OPTIONS_I, "forward-solver-iter-max", '\0',
    &(rhea_inversion_forward_solver_iter_max),
    RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_ITER_MAX,
    "Forward solver: Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "forward-solver-rtol", '\0',
    &(rhea_inversion_forward_solver_rtol),
    RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_RTOL,
    "Forward solver: Relative tolerance for Stokes solver",
  YMIR_OPTIONS_I, "adjoint-solver-iter-max", '\0',
    &(rhea_inversion_adjoint_solver_iter_max),
    RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_ITER_MAX,
    "Adjoint solver: Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "adjoint-solver-rtol", '\0',
    &(rhea_inversion_adjoint_solver_rtol),
    RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_RTOL,
    "Adjoint solver: Relative tolerance for Stokes solver",
  YMIR_OPTIONS_I, "incremental-forward-solver-iter-max", '\0',
    &(rhea_inversion_incremental_forward_solver_iter_max),
    RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_ITER_MAX,
    "Incr. forward solver: Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "incremental-forward-solver-rtol", '\0',
    &(rhea_inversion_incremental_forward_solver_rtol),
    RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_RTOL,
    "Incr. forward solver: Relative tolerance for Stokes solver",
  YMIR_OPTIONS_I, "incremental-adjoint-solver-iter-max", '\0',
    &(rhea_inversion_incremental_adjoint_solver_iter_max),
    RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_ITER_MAX,
    "Incr. adjoint solver: Maximum number of iterations for Stokes solver",
  YMIR_OPTIONS_D, "incremental-adjoint-solver-rtol", '\0',
    &(rhea_inversion_incremental_adjoint_solver_rtol),
    RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_RTOL,
    "Incr. adjoint solver: Relative tolerance for Stokes solver",
  YMIR_OPTIONS_I, "allow-amr-for-outer-n-iter", '\0',
    &(rhea_inversion_allow_amr_for_outer_n_iter),
    RHEA_INVERSION_DEFAULT_ALLOW_AMR_FOR_OUTER_N_ITER,
    "Allow AMR within inner solves until this number of outer iterations",
  YMIR_OPTIONS_I, "project-out-null", '\0',
    &(rhea_inversion_project_out_null),
    RHEA_INVERSION_DEFAULT_PROJECT_OUT_NULL,
    "Inner solves: Project out null spaces",

  /* derivative checks */
  YMIR_OPTIONS_I, "check-gradient", '\0',
    &(rhea_inversion_check_gradient), RHEA_INVERSION_DEFAULT_CHECK_GRADIENT,
    "Check the gradient during Newton iterations",
  YMIR_OPTIONS_I, "check-gradient-of-qoi", '\0',
    &(rhea_inversion_check_gradient_of_qoi),
    RHEA_INVERSION_DEFAULT_CHECK_GRADIENT_OF_QOI,
    "Check the gradient during QOI computations",
  YMIR_OPTIONS_B, "check-gradient-elementwise", '\0',
    &(rhea_inversion_check_gradient_elementwise),
    RHEA_INVERSION_DEFAULT_CHECK_GRADIENT_ELEMENTWISE,
    "Check the gradient for each entry separately",
  YMIR_OPTIONS_I, "check-hessian", '\0',
    &(rhea_inversion_check_hessian), RHEA_INVERSION_DEFAULT_CHECK_HESSIAN,
    "Check the Hessian during Newton iterations",

  /* debugging */
  YMIR_OPTIONS_B, "debug", '\0',
    &(rhea_inversion_debug), RHEA_INVERSION_DEFAULT_DEBUG,
    "Activate output (text and vtk) for debugging",

  /* monitors */
  YMIR_OPTIONS_B, "monitor-performance", '\0',
    &(rhea_inversion_monitor_performance),
    RHEA_INVERSION_DEFAULT_MONITOR_PERFORMANCE,
    "Measure and print performance statistics (e.g., runtime or flops)",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_inversion_param_add_options (&rhea_inversion_param_options, opt);
  rhea_inversion_qoi_add_options (&rhea_inversion_qoi_options, opt);
  rhea_newton_add_options (&rhea_inversion_newton_options, opt);

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);

  /* initialize (deactivated) performance counters */
  rhea_inversion_perfmon_init (0, 0);
}

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/* perfomance monitor tags and names */
typedef enum
{
  RHEA_INVERSION_PERFMON_SOLVE_FWD,
  RHEA_INVERSION_PERFMON_SOLVE_ADJ,
  RHEA_INVERSION_PERFMON_SOLVE_INCR_FWD,
  RHEA_INVERSION_PERFMON_SOLVE_INCR_ADJ,
  RHEA_INVERSION_PERFMON_N
}
rhea_inversion_perfmon_idx_t;

static const char  *rhea_inversion_perfmon_name[RHEA_INVERSION_PERFMON_N] =
{
  "Solve forward problem",
  "Solve adjoint problem",
  "Solve incremental forward problem",
  "Solve incremental adjoint problem"
};

ymir_perf_counter_t rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_N];
ymir_perf_counter_t rhea_inversion_newton_perfmon[RHEA_NEWTON_PERFMON_N];

void
rhea_inversion_perfmon_init (const int activate, const int skip_if_active)
{
  int                 active = activate && rhea_inversion_monitor_performance;

  ymir_perf_counter_init_all_ext (rhea_inversion_perfmon,
                                  rhea_inversion_perfmon_name,
                                  RHEA_INVERSION_PERFMON_N,
                                  active, skip_if_active);

  /* init sub-monitors */
  active = active && rhea_inversion_newton_options.monitor_performance;
  ymir_perf_counter_init_all_ext (rhea_inversion_newton_perfmon,
                                  rhea_newton_perfmon_name,
                                  RHEA_NEWTON_PERFMON_N,
                                  active, skip_if_active);
}

void
rhea_inversion_perfmon_print (sc_MPI_Comm mpicomm,
                              const int print_wtime,
                              const int print_n_calls,
                              const int print_flops)
{
  const int           active = rhea_inversion_monitor_performance;
  const int           print = (print_wtime || print_n_calls || print_flops);
  int                 n;

  /* exit if nothing to do */
  if (!active || !print) {
    return;
  }

  /* main monitors */
  {
    const int           n_stats = RHEA_INVERSION_PERFMON_N *
                                  YMIR_PERF_COUNTER_N_STATS;
    sc_statinfo_t       stats[n_stats];
    char                stats_name[n_stats][YMIR_PERF_COUNTER_NAME_SIZE];

    /* gather performance statistics */
    n = ymir_perf_counter_gather_stats (
        rhea_inversion_perfmon, RHEA_INVERSION_PERFMON_N, stats, stats_name,
        mpicomm, print_wtime, print_n_calls, print_flops);

    /* print performance statistics */
    ymir_perf_counter_print_stats (stats, n, "Inversion");
  }

  /* sub-monitors */
  {
    const int           n_stats = RHEA_NEWTON_PERFMON_N *
                                  YMIR_PERF_COUNTER_N_STATS;
    sc_statinfo_t       stats[n_stats];
    char                stats_name[n_stats][YMIR_PERF_COUNTER_NAME_SIZE];

    /* gather performance statistics */
    n = ymir_perf_counter_gather_stats (
        rhea_inversion_newton_perfmon, RHEA_NEWTON_PERFMON_N, stats, stats_name,
        mpicomm, print_wtime, print_n_calls, print_flops);

    /* print performance statistics */
    ymir_perf_counter_print_stats (stats, n, "Inversion: Newton");
  }
}

/******************************************************************************
 * Inverse Problem Structure
 *****************************************************************************/

/* inverse problem */
struct rhea_inversion_problem
{
  /* Stokes problem (not owned) */
  rhea_stokes_problem_t  *stokes_problem;

  /* inversion parameters */
  rhea_inversion_param_options_t *inv_param_options;
  rhea_inversion_param_t         *inv_param;

  /* observations: velocity at the surface */
  rhea_inversion_obs_velocity_t         vel_obs_type;
  rhea_inversion_obs_velocity_weight_t  vel_obs_weight_type;
  double                               *vel_obs_weight_values;
  ymir_vec_t                           *vel_obs_surf;
  ymir_vec_t                           *vel_obs_weight_surf;

  /* observations: viscosity inside mantle */
  rhea_inversion_obs_viscosity_t  visc_obs_type;
  int                             visc_obs_n;
  rhea_domain_subset_column_t   **visc_obs_column;
  double                         *visc_obs_value;
  double                         *visc_obs_weight;

  /* observations: stress inside mantle */
  rhea_inversion_obs_stress_t  stress_obs_type;
  ymir_vec_t                  *stress_obs;
  ymir_vec_t                  *stress_obs_weight;

  /* weights for observational data and prior terms */
  double              data_abs_weight;
  double              prior_abs_weight;

  /* quantities of interest for post-processing */
  rhea_inversion_qoi_options_t *inv_qoi_options;
  rhea_inversion_qoi_t         *inv_qoi;

  /* forward and adjoint states */
  ymir_vec_t         *forward_vel_press;
  ymir_vec_t         *adjoint_vel_press;
  int                 forward_is_outdated;
  int                 adjoint_is_outdated;
  int                 forward_nonzero_init;
  int                 adjoint_nonzero_init;

  /* incremental forward and adjoint states */
  ymir_vec_t         *incremental_forward_vel_press;
  ymir_vec_t         *incremental_adjoint_vel_press;
  int                 incremental_forward_is_outdated;
  int                 incremental_adjoint_is_outdated;

  /* variables for inner solves */
  ymir_vec_t                       *rhs_vel_buffer;
  ymir_vec_t                       *rhs_vel_press_buffer;
  double                            inner_solver_adaptive_rtol_comp;
  rhea_inversion_project_out_null_t project_out_null;

  /* assembled Hessian matrix */
  sc_dmatrix_t       *hessian_matrix;

  /* Newton problem */
  rhea_newton_options_t  *newton_options;
  rhea_newton_problem_t  *newton_problem;
  ymir_vec_t             *newton_neg_gradient_vec;
  ymir_vec_t             *newton_step_vec;
  ymir_vec_t             *newton_gradient_adjoint_comp_vec;
  ymir_vec_t             *newton_gradient_prior_comp_vec;
  ymir_vec_t             *newton_gradient_observation_comp_vec;
  ymir_vec_t             *newton_gradient_diff_vec;
  int                     newton_gradient_diff_vec_update_count;
  ymir_vec_t            **check_gradient_perturb_vec;
  int                     check_gradient_perturb_n_vecs;

  /* statistics */
  char              **fwd_adj_solve_stats;
  int                 fwd_adj_solve_stats_length;
  int                 fwd_adj_solve_stats_current_idx;
  int                 fwd_adj_solve_stats_total_num_solves[2];
  int                 fwd_adj_solve_stats_total_iterations[2];
  int                 fwd_adj_solve_stats_total_stop_reason[2];
  char              **ifwd_iadj_solve_stats;
  int                 ifwd_iadj_solve_stats_length;
  int                 ifwd_iadj_solve_stats_current_idx;
  int                 ifwd_iadj_solve_stats_total_num_solves[2];
  int                 ifwd_iadj_solve_stats_total_iterations[2];
  int                 ifwd_iadj_solve_stats_total_stop_reason[2];

  /* statistics for numerical errors */
  double              error_stats_objective;
  double              error_stats_gradient;
  double              error_stats_hessian;

  /* output paths */
  char               *txt_path;
  char               *vtk_path_vol;
  char               *vtk_path_surf;
  char               *vtk_path_debug_gradient;
};

static int
rhea_inversion_mesh_data_exists (rhea_inversion_problem_t *inv_problem)
{
  int                 fwd_adj_states_exist;
  int                 vel_obs_exist;

  fwd_adj_states_exist = (
      NULL != inv_problem->forward_vel_press &&
      NULL != inv_problem->adjoint_vel_press &&
      NULL != inv_problem->incremental_forward_vel_press &&
      NULL != inv_problem->incremental_adjoint_vel_press &&
      NULL != inv_problem->rhs_vel_buffer &&
      NULL != inv_problem->rhs_vel_press_buffer
  );
  vel_obs_exist = (
      inv_problem->vel_obs_type == RHEA_INVERSION_OBS_VELOCITY_NONE ||
      inv_problem->vel_obs_surf != NULL
  );

  return (fwd_adj_states_exist && vel_obs_exist);
}

static int
rhea_inversion_solver_data_exists (rhea_inversion_problem_t *inv_problem)
{
  return rhea_inversion_mesh_data_exists (inv_problem);
}

static void
rhea_inversion_set_weight_data (rhea_inversion_problem_t *inv_problem)
{
  if (RHEA_INVERSION_OBS_VELOCITY_NONE  != inv_problem->vel_obs_type ||
      RHEA_INVERSION_OBS_VISCOSITY_NONE != inv_problem->visc_obs_type ||
      RHEA_INVERSION_OBS_STRESS_NONE    != inv_problem->stress_obs_type) {
    inv_problem->data_abs_weight = 1.0;
  }
  else {
    inv_problem->data_abs_weight = 0.0;
  }
}

static void
rhea_inversion_set_weight_prior (rhea_inversion_problem_t *inv_problem)
{
  const double        prior_stddev = rhea_inversion_parameter_prior_stddev;

  if (isfinite (prior_stddev) && 0.0 < prior_stddev) {
    inv_problem->prior_abs_weight = 1.0/(prior_stddev*prior_stddev);
  }
  else {
    inv_problem->prior_abs_weight = 0.0;
  }
}

/******************************************************************************
 * State Solvers
 *****************************************************************************/

/**
 * Creates mesh dependent data.
 */
static void
rhea_inversion_newton_create_mesh_data (rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t            *ymir_mesh;
  ymir_pressure_elem_t   *press_elem;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* get mesh data */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* create (forward/adjoint) state variables */
  if (inv_problem->forward_vel_press == NULL) {
    inv_problem->forward_vel_press =
      rhea_velocity_pressure_new (ymir_mesh, press_elem);
    inv_problem->forward_is_outdated = 1;
    inv_problem->forward_nonzero_init = 0;
  }
  RHEA_ASSERT (inv_problem->adjoint_vel_press == NULL);
  RHEA_ASSERT (inv_problem->incremental_forward_vel_press == NULL);
  RHEA_ASSERT (inv_problem->incremental_adjoint_vel_press == NULL);
  inv_problem->adjoint_vel_press =
    rhea_velocity_pressure_new (ymir_mesh, press_elem);
  inv_problem->incremental_forward_vel_press =
    rhea_velocity_pressure_new (ymir_mesh, press_elem);
  inv_problem->incremental_adjoint_vel_press =
    rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* reset flags for new state vectors */
  inv_problem->adjoint_is_outdated = 1;
  inv_problem->adjoint_nonzero_init = 0;
  inv_problem->incremental_forward_is_outdated = 1;
  inv_problem->incremental_adjoint_is_outdated = 1;

  /* create variables for inner solves */
  RHEA_ASSERT (inv_problem->rhs_vel_buffer == NULL);
  RHEA_ASSERT (inv_problem->rhs_vel_press_buffer == NULL);
  inv_problem->rhs_vel_buffer       = rhea_velocity_new (ymir_mesh);
  inv_problem->rhs_vel_press_buffer = rhea_velocity_pressure_new (ymir_mesh,
                                                                  press_elem);

  /* create velocity observational data */
  RHEA_ASSERT (inv_problem->vel_obs_surf == NULL);
  RHEA_ASSERT (inv_problem->vel_obs_weight_surf == NULL);
  if (RHEA_INVERSION_OBS_VELOCITY_NONE != inv_problem->vel_obs_type) {
    rhea_plate_options_t *plate_options =
      rhea_stokes_problem_get_plate_options (stokes_problem);
    const int           n_plates = rhea_plate_get_n_plates (plate_options);
    const double        vel_dim_mm_yr = rhea_velocity_get_dim_mm_yr (
        rhea_stokes_problem_get_domain_options (stokes_problem),
        rhea_stokes_problem_get_temperature_options (stokes_problem));
    double             *calculated_weight_values;

    /* create velocities and weights */
    inv_problem->vel_obs_surf = rhea_velocity_surface_new (ymir_mesh);
    inv_problem->vel_obs_weight_surf = ymir_face_cvec_new (
        ymir_mesh, RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
    calculated_weight_values = RHEA_ALLOC (double, SC_MAX (1, n_plates));
    rhea_inversion_obs_velocity_generate (
        inv_problem->vel_obs_surf, inv_problem->vel_obs_weight_surf,
        inv_problem->vel_obs_type, inv_problem->vel_obs_weight_type,
        inv_problem->vel_obs_weight_values, plate_options,
        calculated_weight_values);

    /* print std. dev. of velocity data for each plate */
    if (0 < n_plates) { /* if plates exist */
      int                 pid;

      RHEA_GLOBAL_INFO ("========================================\n");
      RHEA_GLOBAL_INFO ("rhea_inversion: "
                        "obs_velocity_weight, stddev [mm/yr]\n");
      RHEA_GLOBAL_INFO ("----------------------------------------\n");
      for (pid = 0; pid < n_plates; pid++) {
        RHEA_GLOBAL_INFOF ("plate_idx %3i: %.6e, %g\n",
                           pid, calculated_weight_values[pid],
                           vel_dim_mm_yr/calculated_weight_values[pid]);
      }
      RHEA_GLOBAL_INFO ("========================================\n");
    }
    else {
      RHEA_GLOBAL_INFOF_FN_TAG ("rhea_inversion",
                                "obs_velocity_weight=%g, stddev_mm_yr=%g",
                                calculated_weight_values[0],
                                vel_dim_mm_yr/calculated_weight_values[0]);
    }
    RHEA_FREE (calculated_weight_values);
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/**
 * Destroys mesh dependent data.
 */
static void
rhea_inversion_newton_clear_mesh_data (rhea_inversion_problem_t *inv_problem,
                                       const int keep_forward_vel_press)
{
  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* destroy (forward/adjoint) state variables */
  if (!keep_forward_vel_press) {
    rhea_velocity_pressure_destroy (inv_problem->forward_vel_press);
    inv_problem->forward_vel_press = NULL;
  }
  rhea_velocity_pressure_destroy (inv_problem->adjoint_vel_press);
  rhea_velocity_pressure_destroy (inv_problem->incremental_forward_vel_press);
  rhea_velocity_pressure_destroy (inv_problem->incremental_adjoint_vel_press);
  inv_problem->adjoint_vel_press             = NULL;
  inv_problem->incremental_forward_vel_press = NULL;
  inv_problem->incremental_adjoint_vel_press = NULL;

  /* destroy variables for inner solves */
  rhea_velocity_destroy (inv_problem->rhs_vel_buffer);
  rhea_velocity_pressure_destroy (inv_problem->rhs_vel_press_buffer);
  inv_problem->rhs_vel_buffer       = NULL;
  inv_problem->rhs_vel_press_buffer = NULL;

  /* destroy observational data */
  if (RHEA_INVERSION_OBS_VELOCITY_NONE != inv_problem->vel_obs_type) {
    rhea_velocity_surface_destroy (inv_problem->vel_obs_surf);
    ymir_vec_destroy (inv_problem->vel_obs_weight_surf);
    inv_problem->vel_obs_surf        = NULL;
    inv_problem->vel_obs_weight_surf = NULL;
  }
  else {
    RHEA_ASSERT (NULL == inv_problem->vel_obs_surf);
    RHEA_ASSERT (NULL == inv_problem->vel_obs_weight_surf);
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static char        *rhea_inversion_inner_solve_stats_cat (const char *a,
                                                          const char *b);

/**
 * Writes statistics of forward and adjoint solves.
 */
static void
rhea_inversion_fwd_adj_solve_stats_write (rhea_inversion_problem_t *inv_problem,
                                          const int is_adjoint_solve,
                                          const int num_iterations,
                                          const double residual_reduction,
                                          const int stop_reason)
{
  const int           stats_idx = inv_problem->fwd_adj_solve_stats_current_idx;
  char               *stats;
  char                append[64];

  /* check input */
  RHEA_ASSERT (0 <= stats_idx);

  /* create text */
  if (!is_adjoint_solve) {
    snprintf (append, 64, "fwd=(%3d, %.2e, %d) ",
              num_iterations, residual_reduction, stop_reason);
  }
  else {
    snprintf (append, 64, "adj=(%3d, %.2e, %d); ",
              num_iterations, residual_reduction, stop_reason);
  }

  /* append text */
  stats = inv_problem->fwd_adj_solve_stats[stats_idx];
  inv_problem->fwd_adj_solve_stats[stats_idx] =
    rhea_inversion_inner_solve_stats_cat (stats, append);
  RHEA_FREE (stats);

  /* update total values */
  inv_problem->fwd_adj_solve_stats_total_num_solves[is_adjoint_solve]++;
  inv_problem->fwd_adj_solve_stats_total_iterations[is_adjoint_solve] +=
    num_iterations;
  inv_problem->fwd_adj_solve_stats_total_stop_reason[is_adjoint_solve] =
    inv_problem->fwd_adj_solve_stats_total_stop_reason[is_adjoint_solve] &&
    rhea_stokes_problem_has_converged_ext (
      stop_reason, inv_problem->stokes_problem,
      is_adjoint_solve /* force_linear_solve */);
}

/**
 * Runs solver for the forward state.
 */
static void
rhea_inversion_inner_solve_forward (rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  const int           nonzero_init = inv_problem->forward_nonzero_init;
  const int           iter_max = rhea_inversion_forward_solver_iter_max;
  const double        rtol = rhea_inversion_forward_solver_rtol *
                             inv_problem->inner_solver_adaptive_rtol_comp;
  const int           resume = inv_problem->forward_nonzero_init;
  int                 stop_reason, n_iter;
  double              res_reduc;
  int                 mesh_modified_by_solver;

  /* exit if nothing to do */
  if (!inv_problem->forward_is_outdated) {
    return;
  }

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);
  ymir_perf_counter_start (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_SOLVE_FWD]);

  /* update Stokes right-hand-side for forward problem */
  if (!inv_problem->adjoint_is_outdated ||
      !inv_problem->incremental_forward_is_outdated ||
      !inv_problem->incremental_adjoint_is_outdated) {
    rhea_stokes_problem_update_solver (
        stokes_problem, 0 /* !update coeff */, NULL /* since no coeff update */,
        0 /* !override_rhs */, NULL, NULL, NULL, NULL);
  }

  /* run Stokes solver */
  stop_reason = rhea_stokes_problem_solve_ext (
      &(inv_problem->forward_vel_press), nonzero_init, iter_max, rtol,
      stokes_problem, resume, 0 /* !force linear solve */,
      RHEA_INVERSION_KRYLOV_SOLVER_IDX_FORWARD,
      &n_iter, &res_reduc, &mesh_modified_by_solver);
  rhea_stokes_problem_enforce_boundary_conditions (
      inv_problem->forward_vel_press, stokes_problem);
  inv_problem->forward_is_outdated = 0;
  inv_problem->forward_nonzero_init = 1;

  /* write solver statistics */
  rhea_inversion_fwd_adj_solve_stats_write (inv_problem, 0 /* forward */,
                                            n_iter, res_reduc, stop_reason);
  inv_problem->error_stats_objective = res_reduc;

  /* recreate mesh dependent data */
  if (mesh_modified_by_solver) {
    rhea_inversion_newton_clear_mesh_data (inv_problem, 1 /* keep forward */);
    rhea_inversion_newton_create_mesh_data (inv_problem);

    /* deactivate AMR during nonlinear solves (e.g., during backtracking) */
    rhea_stokes_problem_amr_set_nonlinear_type_name ("NONE");
  }

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_SOLVE_FWD]);
  RHEA_GLOBAL_INFO_FN_END (__func__);
}

/**
 * Runs solver for the adjoint state.
 */
static void
rhea_inversion_inner_solve_adjoint (rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  const int           nonzero_init = inv_problem->adjoint_nonzero_init;
  const int           iter_max = rhea_inversion_adjoint_solver_iter_max;
  const double        rtol = rhea_inversion_adjoint_solver_rtol *
                             inv_problem->inner_solver_adaptive_rtol_comp;
  const int           resume = 0;
  int                 stop_reason, n_iter;
  double              res_reduc;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);
  ymir_perf_counter_start (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_SOLVE_ADJ]);

  /* run Stokes solver */
  stop_reason = rhea_stokes_problem_solve_ext (
      &(inv_problem->adjoint_vel_press), nonzero_init, iter_max, rtol,
      stokes_problem, resume, 1 /* force linear solve */,
      RHEA_INVERSION_KRYLOV_SOLVER_IDX_ADJOINT,
      &n_iter, &res_reduc, NULL /* mesh_modified_by_solver */);
  inv_problem->adjoint_is_outdated = 0;
  inv_problem->adjoint_nonzero_init = 1;

  /* write solver statistics */
  rhea_inversion_fwd_adj_solve_stats_write (inv_problem, 1 /* adjoint */,
                                            n_iter, res_reduc, stop_reason);

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_SOLVE_ADJ]);
  RHEA_GLOBAL_INFO_FN_END (__func__);
}

/**
 * Writes statistics of incremental forward and adjoint solves.
 */
static void
rhea_inversion_ifwd_iadj_solve_stats_write (
                                        rhea_inversion_problem_t *inv_problem,
                                        const int is_adjoint_solve,
                                        const int num_iterations,
                                        const double residual_reduction,
                                        const int stop_reason)
{
  const int           stats_idx =
                        inv_problem->ifwd_iadj_solve_stats_current_idx;
  char               *stats;
  char                append[64];

  /* check input */
  RHEA_ASSERT (0 <= stats_idx);

  /* create text */
  if (!is_adjoint_solve) {
    snprintf (append, 64, "(%3d, %.2e, %d)  ",
              num_iterations, residual_reduction, stop_reason);
  }
  else {
    snprintf (append, 64, "(%3d, %.2e, %d); ",
              num_iterations, residual_reduction, stop_reason);
  }

  /* append text */
  stats = inv_problem->ifwd_iadj_solve_stats[stats_idx];
  inv_problem->ifwd_iadj_solve_stats[stats_idx] =
    rhea_inversion_inner_solve_stats_cat (stats, append);
  RHEA_FREE (stats);

  /* update total values */
  inv_problem->ifwd_iadj_solve_stats_total_num_solves[is_adjoint_solve]++;
  inv_problem->ifwd_iadj_solve_stats_total_iterations[is_adjoint_solve] +=
    num_iterations;
  inv_problem->ifwd_iadj_solve_stats_total_stop_reason[is_adjoint_solve] =
    inv_problem->ifwd_iadj_solve_stats_total_stop_reason[is_adjoint_solve] &&
    rhea_stokes_problem_has_converged_ext (
      stop_reason, inv_problem->stokes_problem, 1 /* force_linear_solve */);
}

/**
 * Runs solver for the incremental forward state.
 */
static void
rhea_inversion_inner_solve_incremental_forward (
                                        rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  const int           nonzero_init = 0; //TODO 1?
  const int           iter_max =
                        rhea_inversion_incremental_forward_solver_iter_max;
  const double        rtol =
                        rhea_inversion_incremental_forward_solver_rtol *
                        inv_problem->inner_solver_adaptive_rtol_comp;
  const int           resume = 0;
  int                 stop_reason, n_iter;
  double              res_reduc;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);
  ymir_perf_counter_start (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_SOLVE_INCR_FWD]);

  /* run Stokes solver */
  stop_reason = rhea_stokes_problem_solve_ext (
      &(inv_problem->incremental_forward_vel_press), nonzero_init,
      iter_max, rtol, stokes_problem, resume, 1 /* force linear solve */,
      RHEA_INVERSION_KRYLOV_SOLVER_IDX_INCR_FORWARD,
      &n_iter, &res_reduc, NULL /* mesh_modified_by_solver */);
  inv_problem->incremental_forward_is_outdated = 0;

  /* write solver statistics */
  rhea_inversion_ifwd_iadj_solve_stats_write (inv_problem, 0 /* forward */,
                                              n_iter, res_reduc, stop_reason);

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_SOLVE_INCR_FWD]);
  RHEA_GLOBAL_INFO_FN_END (__func__);
}

/**
 * Runs solver for the incremental adjoint state.
 */
static void
rhea_inversion_inner_solve_incremental_adjoint (
                                        rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  const int           nonzero_init = 0; //TODO 1?
  const int           iter_max =
                        rhea_inversion_incremental_adjoint_solver_iter_max;
  const double        rtol =
                        rhea_inversion_incremental_adjoint_solver_rtol *
                        inv_problem->inner_solver_adaptive_rtol_comp;
  const int           resume = 0;
  int                 stop_reason, n_iter;
  double              res_reduc;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);
  ymir_perf_counter_start (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_SOLVE_INCR_ADJ]);

  /* run Stokes solver */
  stop_reason = rhea_stokes_problem_solve_ext (
      &(inv_problem->incremental_adjoint_vel_press), nonzero_init,
      iter_max, rtol, stokes_problem, resume, 1 /* force linear solve */,
      RHEA_INVERSION_KRYLOV_SOLVER_IDX_INCR_ADJOINT,
      &n_iter, &res_reduc, NULL /* mesh_modified_by_solver */);
  inv_problem->incremental_adjoint_is_outdated = 0;

  /* write solver statistics */
  rhea_inversion_ifwd_iadj_solve_stats_write (inv_problem, 1 /* adjoint */,
                                              n_iter, res_reduc, stop_reason);

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_SOLVE_INCR_ADJ]);
  RHEA_GLOBAL_INFO_FN_END (__func__);
}

/******************************************************************************
 * Helper Functions, Pre- and Post-Processing
 *****************************************************************************/

static void         rhea_inversion_assemble_hessian (
                                        sc_dmatrix_t *hessian_matrix,
                                        rhea_inversion_problem_t *inv_problem,
                                        const rhea_inversion_hessian_t type,
                                        const int enforce_symmetry);

static void         rhea_inversion_apply_hessian (
                                        ymir_vec_t *param_vec_out,
                                        ymir_vec_t *param_vec_in,
                                        rhea_inversion_problem_t *inv_problem,
                                        const rhea_inversion_hessian_t type);

static double       rhea_inversion_newton_evaluate_objective_fn (
                                            ymir_vec_t *solution, void *data,
                                            double *obj_comp);

static void         rhea_inversion_newton_compute_negative_gradient_fn (
                                            ymir_vec_t *neg_gradient,
                                            ymir_vec_t *solution, void *data);

static int
rhea_inversion_read_prev_newton (sc_dmatrix_t *hessian_mat,
                                 ymir_vec_t *neg_gradient_vec,
                                 ymir_vec_t *step_vec,
                                 const char *hessian_prev_file_path_txt,
                                 const char *gradient_prev_file_path_txt,
                                 const char *step_prev_file_path_txt,
                                 rhea_inversion_problem_t *inv_problem)
{
  rhea_inversion_param_t *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t        *ymir_mesh =
                        rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  const int           n_active = rhea_inversion_param_get_n_active (inv_param);
  int                 success = 0;

  /* check input */
  RHEA_ASSERT (NULL != hessian_prev_file_path_txt);
  RHEA_ASSERT (NULL != gradient_prev_file_path_txt);
  RHEA_ASSERT (NULL != step_prev_file_path_txt);

  /* read hessian */
  if (NULL != hessian_mat && NULL != hessian_prev_file_path_txt) {
    int                 success_n;

    success_n = rhea_io_mpi_read_broadcast_double (
        hessian_mat->e[0], hessian_mat->m * hessian_mat->n,
        NULL /* path bin */, hessian_prev_file_path_txt, mpicomm);
    RHEA_ASSERT (success_n == n_active*n_active);
    RHEA_ASSERT (sc_dmatrix_is_valid (hessian_mat));
    success += (success_n == n_active*n_active);
  }

  /* read gradient */
  if (NULL != neg_gradient_vec && gradient_prev_file_path_txt) {
    sc_dmatrix_t       *grad_read;
    sc_dmatrix_t       *grad_reduced;
    int                 success_n, row;

    grad_read = sc_dmatrix_new (n_active, 3);
    grad_reduced = rhea_inversion_param_vec_reduced_new (
        neg_gradient_vec, inv_param);

    success_n = rhea_io_mpi_read_broadcast_double (
        grad_read->e[0], grad_read->m * grad_read->n,
        NULL /* path bin */, gradient_prev_file_path_txt, mpicomm);
    RHEA_ASSERT (success_n == 3*n_active);
    RHEA_ASSERT (sc_dmatrix_is_valid (grad_read));
    success += (success_n == 3*n_active);

    for (row = 0; row < n_active; row++) {
      grad_reduced->e[row][0] = grad_read->e[row][0];
    }
    rhea_inversion_param_vec_reduced_copy (
        neg_gradient_vec, grad_reduced, inv_param);
    ymir_vec_scale (-1.0, neg_gradient_vec);
    RHEA_ASSERT (rhea_inversion_param_vec_is_valid (neg_gradient_vec,
                                                    inv_param));

    sc_dmatrix_destroy (grad_read);
    rhea_inversion_param_vec_reduced_destroy (grad_reduced);
  }

  /* read step */
  if (NULL != step_vec && step_prev_file_path_txt) {
    sc_dmatrix_t       *step_reduced;
    int                 success_n;

    step_reduced = rhea_inversion_param_vec_reduced_new (
        step_vec, inv_param);

    success_n = rhea_io_mpi_read_broadcast_double (
        step_reduced->e[0], step_reduced->m * step_reduced->n,
        NULL /* path bin */, step_prev_file_path_txt, mpicomm);
    RHEA_ASSERT (success_n == n_active);
    RHEA_ASSERT (sc_dmatrix_is_valid (step_reduced));
    success += (success_n == n_active);

    rhea_inversion_param_vec_reduced_copy (
        step_vec, step_reduced, inv_param);
    rhea_inversion_param_vec_reduced_destroy (step_reduced);
  }

  return success;
}

static void
rhea_inversion_write_posterior_hessian (
                            const char *hessian_path,
                            const rhea_inversion_hessian_t write_hessian_type,
                            rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t            *ymir_mesh;
  sc_MPI_Comm         mpicomm;
  int                 mpirank, mpiret;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "write_hessian_type=%i",
                              write_hessian_type);

  /* exit if nothing to do */
  if (RHEA_INVERSION_HESSIAN_NONE == write_hessian_type) {
    RHEA_GLOBAL_INFO_FN_END (__func__);
    return;
  }

  /* check input */
  RHEA_ASSERT (NULL != hessian_path);
  RHEA_ASSERT (NULL != inv_problem->hessian_matrix);

  switch (write_hessian_type) {
  case RHEA_INVERSION_HESSIAN_FIRST_ORDER_APPROX:
  case RHEA_INVERSION_HESSIAN_FULL:
    /* assemble Hessian matrix */
    rhea_inversion_assemble_hessian (
        inv_problem->hessian_matrix, inv_problem, write_hessian_type,
        rhea_inversion_assemble_hessian_enforce_symm);

    /* get parallel environment */
    ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
    mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
    mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

    /* write to file */
    if (mpirank == 0) {
      rhea_io_std_write_double_to_txt (
          hessian_path, inv_problem->hessian_matrix->e[0],
          inv_problem->hessian_matrix->m *
          inv_problem->hessian_matrix->n,
          inv_problem->hessian_matrix->n);
    }
    break;
  default: /* unsupported Hessian type */
    RHEA_ABORT_NOT_REACHED ();
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

static void
rhea_inversion_propagate_solution_to_qoi (const char *qoi_path,
                                          ymir_vec_t *parameter_vec,
                                          rhea_inversion_problem_t *inv_problem)
{
  /* backup of settings */
  const rhea_inversion_obs_velocity_t  vel_obs_type    = inv_problem->vel_obs_type;
  const rhea_inversion_obs_viscosity_t visc_obs_type   = inv_problem->visc_obs_type;
  const rhea_inversion_obs_stress_t    stress_obs_type = inv_problem->stress_obs_type;
  ymir_vec_t const   *stress_obs_weight = inv_problem->stress_obs_weight;
  const double        prior_abs_weight  = inv_problem->prior_abs_weight;
  /* QOI variables */
  rhea_inversion_qoi_t *inv_qoi = inv_problem->inv_qoi;
  const int           n_parameters = rhea_inversion_param_get_n_active (
                          inv_problem->inv_param);
  sc_dmatrix_t       *stress_qoi_vec, *stress_qoi_jac;
  size_t              idx_stress_qoi, idx_param;
  char                path[BUFSIZ];

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "#qoi=%i, #parameters=%i, path=%s",
                              rhea_inversion_qoi_get_count_all (inv_qoi),
                              n_parameters, qoi_path);

  /* nothing to do if
   * - no QOI,
   * - no parameter vector,
   * - Hessian matrix does not exist */
  if (!rhea_inversion_qoi_exists (inv_qoi) ||
      NULL == parameter_vec ||
      NULL == inv_problem->hessian_matrix) {
    return;
  }

  /* create stress QOI */
  stress_qoi_vec = rhea_inversion_qoi_stress_vec_reduced_new (inv_qoi);
  stress_qoi_jac = rhea_inversion_qoi_stress_to_param_matrix_new (
      n_parameters, inv_qoi);
  RHEA_ASSERT (NULL != stress_qoi_vec);
  RHEA_ASSERT (NULL != stress_qoi_jac);
  RHEA_ASSERT (stress_qoi_jac->m == stress_qoi_vec->m);
  RHEA_ASSERT (stress_qoi_jac->n == n_parameters);

  /* deactivate all non-stress observations and prior */
  inv_problem->vel_obs_type     = RHEA_INVERSION_OBS_VELOCITY_NONE;
  inv_problem->visc_obs_type    = RHEA_INVERSION_OBS_VISCOSITY_NONE;
  inv_problem->prior_abs_weight = 0.0;

  /* compute stress QOI and their gradients */
  {
    ymir_vec_t         *neg_gradient_vec = rhea_inversion_param_vec_new (
                            inv_problem->inv_param);

    for (idx_stress_qoi = 0; idx_stress_qoi < stress_qoi_vec->m; idx_stress_qoi++) {
      sc_dmatrix_t       *neg_gradient_reduced;

      RHEA_GLOBAL_INFOF_FN_TAG (__func__, "stress_qoi_id=%u, #stress_qoi=%u",
                                (unsigned int) idx_stress_qoi,
                                (unsigned int) stress_qoi_vec->m);

      /* set observation types for QOI */
      rhea_inversion_qoi_stress_create_obs (&inv_problem->stress_obs_type,
                                            &inv_problem->stress_obs_weight,
                                            inv_qoi, idx_stress_qoi);

      /* compute QOI */
      stress_qoi_vec->e[idx_stress_qoi][0] =
        rhea_inversion_newton_evaluate_objective_fn (parameter_vec,
                                                     inv_problem, NULL);

      /* compute gradient */
      if (rhea_inversion_debug) { /* if set vtk debug path */
        snprintf (path, BUFSIZ, "%s_debug_stress_qoi%02u",
                  inv_problem->vtk_path_vol, (unsigned int) idx_stress_qoi);
        inv_problem->vtk_path_debug_gradient = path;
      }
      rhea_inversion_newton_compute_negative_gradient_fn (
          neg_gradient_vec, parameter_vec, inv_problem);
      if (rhea_inversion_debug) { /* if remove vtk debug path */
        inv_problem->vtk_path_debug_gradient = NULL;
      }
      rhea_newton_check_gradient (parameter_vec, neg_gradient_vec,
                                  rhea_inversion_check_gradient_of_qoi,
                                  inv_problem->check_gradient_perturb_vec,
                                  inv_problem->check_gradient_perturb_n_vecs,
                                  1000 + (int) idx_stress_qoi,
                                  inv_problem->newton_problem);

      /* fill gradient into QOI Jacobian */
      neg_gradient_reduced = rhea_inversion_param_vec_reduced_new (
          neg_gradient_vec, inv_problem->inv_param);
      RHEA_ASSERT (neg_gradient_reduced->m == stress_qoi_jac->n);
      RHEA_ASSERT (neg_gradient_reduced->m == n_parameters);
      for (idx_param = 0; idx_param < n_parameters; idx_param++) {
        stress_qoi_jac->e[idx_stress_qoi][idx_param] =
          - neg_gradient_reduced->e[idx_param][0];
      }

      /* destroy */
      rhea_inversion_param_vec_reduced_destroy (neg_gradient_reduced);
      rhea_inversion_qoi_stress_clear_obs (inv_problem->stress_obs_weight,
                                           inv_qoi);
    }

    rhea_inversion_param_vec_destroy (neg_gradient_vec);
  }

  /* restore settings from backup */
  inv_problem->vel_obs_type      = vel_obs_type;
  inv_problem->visc_obs_type     = visc_obs_type;
  inv_problem->stress_obs_type   = stress_obs_type;
  inv_problem->stress_obs_weight = (ymir_vec_t *) stress_obs_weight;
  inv_problem->prior_abs_weight  = prior_abs_weight;

  /* write to file */
  {
    ymir_mesh_t        *ymir_mesh;
    sc_MPI_Comm         mpicomm;
    int                 mpirank, mpiret;
    char                path[BUFSIZ];
    /* matrices */
    const int           n_qoi       = stress_qoi_jac->m;
    sc_dmatrix_t       *Cov_qoi_mat = sc_dmatrix_new (n_qoi, n_qoi);
    sc_dmatrix_t       *Id_mat  = sc_dmatrix_new (n_parameters, n_parameters);
    sc_dmatrix_t       *Cov_mat = sc_dmatrix_new (n_parameters, n_parameters);
    sc_dmatrix_t       *CJ_mat  = sc_dmatrix_new (n_parameters, n_qoi);

    /* get parallel environment */
    ymir_mesh = rhea_stokes_problem_get_ymir_mesh (
        inv_problem->stokes_problem);
    mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
    mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

    /* set identity matrix */
    sc_dmatrix_set_value (Id_mat, 0.0);
    for (idx_param = 0; idx_param < n_parameters; idx_param++) {
      Id_mat->e[idx_param][idx_param] = 1.0;
    }

    /* invert Hessian to compute covariance matrix */
    RHEA_ASSERT (inv_problem->hessian_matrix->m == n_parameters);
    RHEA_ASSERT (inv_problem->hessian_matrix->n == n_parameters);
    sc_dmatrix_rdivide (SC_NO_TRANS, Id_mat, inv_problem->hessian_matrix,
                        Cov_mat);
    sc_dmatrix_destroy (Id_mat);

    /* multipy covariance by Jacobian^T (from right) and
     * Jacobian (from left) */
    sc_dmatrix_multiply (SC_NO_TRANS /* Cov */, SC_TRANS /* Jac^T */, 1.0,
                         Cov_mat, stress_qoi_jac, 0.0, CJ_mat);
    sc_dmatrix_multiply (SC_NO_TRANS /* Jac */, SC_NO_TRANS /* CJ */, 1.0,
                         stress_qoi_jac, CJ_mat, 0.0, Cov_qoi_mat);
    sc_dmatrix_destroy (Cov_mat);
    sc_dmatrix_destroy (CJ_mat);

    /* write QOI covariance matrix to file */
    snprintf (path, BUFSIZ, "%s_cov.txt", qoi_path);
    if (mpirank == 0) {
      rhea_io_std_write_double_to_txt (path, Cov_qoi_mat->e[0],
                                       Cov_qoi_mat->m * Cov_qoi_mat->n,
                                       Cov_qoi_mat->n);
    }
    sc_dmatrix_destroy (Cov_qoi_mat);

    /* write QOI Jacobian matrix to file */
    snprintf (path, BUFSIZ, "%s_jac.txt", qoi_path);
    if (mpirank == 0) {
      rhea_io_std_write_double_to_txt (path, stress_qoi_jac->e[0],
                                       stress_qoi_jac->m * stress_qoi_jac->n,
                                       stress_qoi_jac->n);
    }

    /* write QOI mean vector to file */
    snprintf (path, BUFSIZ, "%s_mean.txt", qoi_path);
    if (mpirank == 0) {
      rhea_io_std_write_double_to_txt (path, stress_qoi_vec->e[0],
                                       stress_qoi_vec->m * stress_qoi_vec->n,
                                       stress_qoi_vec->n);
    }
  }

  /* destroy */
  rhea_inversion_qoi_stress_to_param_matrix_destroy (stress_qoi_jac);
  rhea_inversion_qoi_stress_vec_reduced_destroy (stress_qoi_vec);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

/**
 * Updates inversion parameters in the Stokes model.
 */
static void
rhea_inversion_update_param (ymir_vec_t *solution,
                             rhea_inversion_param_t *inv_param)
{
  const int           solution_exists = (solution != NULL);
  ymir_vec_t         *parameter_vec;
#ifdef RHEA_ENABLE_DEBUG
  ymir_vec_t         *parameter_vec_prev =
                        rhea_inversion_param_vec_new (inv_param);

  /* get parameters before update */
  rhea_inversion_param_pull_from_model (parameter_vec_prev, inv_param);
#endif

  /* get parameter vector */
  if (solution_exists) { /* if parameters are given */
    parameter_vec = solution;
  }
  else { /* otherwise use zeros as parameters */
    parameter_vec = rhea_inversion_param_vec_new (inv_param);
    ymir_vec_set_zero (parameter_vec);
  }

  /* update parameters in Stokes model */
  rhea_inversion_param_push_to_model (parameter_vec, inv_param);

  /* print previous and new parameters */
#ifdef RHEA_ENABLE_DEBUG
  {
    const int          *active = rhea_inversion_param_get_active (inv_param);
    const double       *param = parameter_vec->meshfree->e[0];
    const double       *param_prev = parameter_vec_prev->meshfree->e[0];
    int                 i;

    RHEA_GLOBAL_VERBOSE ("========================================\n");
    RHEA_GLOBAL_VERBOSEF ("%s\n", __func__);
    RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
    for (i = 0; i < rhea_inversion_param_get_n_parameters (inv_param); i++) {
      if (active[i]) {
        RHEA_GLOBAL_VERBOSEF ("param# %3i: %.15e -> %.15e\n", i,
                              param_prev[i], param[i]);
      }
    }
    RHEA_GLOBAL_VERBOSE ("========================================\n");
  }
  rhea_inversion_param_vec_destroy (parameter_vec_prev);
#endif

  /* destroy */
  if (!solution_exists) {
    rhea_inversion_param_vec_destroy (parameter_vec);
  }
}

/******************************************************************************
 * Callback Functions for Newton Solver
 *****************************************************************************/

static void
rhea_inversion_newton_create_solver_data_fn (ymir_vec_t *solution, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;

  /* exit if nothing to do */
  if (rhea_inversion_solver_data_exists (inv_problem)) {
    return;
  }

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (solution == NULL ||
               rhea_inversion_param_vec_check_type (solution, inv_param));
  RHEA_ASSERT (solution == NULL ||
               rhea_inversion_param_vec_is_valid (solution, inv_param));

  /* create mesh dependent data */
  rhea_inversion_newton_create_mesh_data (inv_problem);

  /* create Hessian matrix */
  if (rhea_inversion_assemble_hessian_matrix) {
    const int           n = rhea_inversion_param_get_n_active (inv_param);

    RHEA_ASSERT (inv_problem->hessian_matrix == NULL);
    inv_problem->hessian_matrix = sc_dmatrix_new (n, n);
  }

  /* read previous Hessian and gradient from file */
  if (rhea_inversion_assemble_hessian_matrix &&
      NULL != rhea_inversion_hessian_prev_file_path_txt &&
      NULL != rhea_inversion_gradient_prev_file_path_txt &&
      NULL != rhea_inversion_step_prev_file_path_txt) {
    int                 success;

    success = rhea_inversion_read_prev_newton (
        inv_problem->hessian_matrix,
        inv_problem->newton_neg_gradient_vec,
        inv_problem->newton_step_vec,
        rhea_inversion_hessian_prev_file_path_txt,
        rhea_inversion_gradient_prev_file_path_txt,
        rhea_inversion_step_prev_file_path_txt,
        inv_problem);
    RHEA_ASSERT (3 == success);
    inv_problem->newton_gradient_diff_vec_update_count = 2;
  }

  /* update inversion parameters */
  rhea_inversion_update_param (solution, inv_param);

  /* set up Stokes solver for forward problem */
  rhea_stokes_problem_setup_solver (stokes_problem);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

static void
rhea_inversion_newton_clear_solver_data_fn (ymir_vec_t *solution, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  const rhea_inversion_hessian_t  write_hessian_type =
    (rhea_inversion_hessian_t) rhea_inversion_hessian_type_write_after_solve;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "write_hessian_type=%i",
                              write_hessian_type);

  /* check input */
  RHEA_ASSERT (rhea_inversion_solver_data_exists (inv_problem));

  /* write Hessian matrix to file */
  if (NULL != inv_problem->txt_path && NULL != inv_problem->hessian_matrix) {
    char                path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s_posterior_hessian.txt", inv_problem->txt_path);
    rhea_inversion_write_posterior_hessian (path, write_hessian_type,
                                            inv_problem);
  }

  /* propagate results to QOI */
  if (NULL != inv_problem->txt_path) {
    char                path[BUFSIZ];

    snprintf (path, BUFSIZ, "%s_posterior_qoi", inv_problem->txt_path);
    rhea_inversion_propagate_solution_to_qoi (path, solution, inv_problem);
  }

  /* destroy Hessian matrix */
  if (NULL != inv_problem->hessian_matrix) {
    sc_dmatrix_destroy (inv_problem->hessian_matrix);
    inv_problem->hessian_matrix = NULL;
  }

  /* destroy mesh dependent data */
  rhea_inversion_newton_clear_mesh_data (inv_problem, 0 /* destroy forward */);

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

static void
rhea_inversion_newton_update_operator_fn (ymir_vec_t *solution, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (solution == NULL ||
               rhea_inversion_param_vec_check_type (solution, inv_param));
  RHEA_ASSERT (solution == NULL ||
               rhea_inversion_param_vec_is_valid (solution, inv_param));

  /* ensure the forward state is up-to-date */
  rhea_inversion_inner_solve_forward (inv_problem);

  /* update inversion parameters */
  rhea_inversion_update_param (solution, inv_param);

  /* print updated parameters */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFO ("Inversion candidate parameters\n");
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  rhea_inversion_param_print (
      solution, RHEA_INVERSION_PARAM_VERBOSE_REAL_NONDIM_DIM, inv_param);
  RHEA_GLOBAL_INFO ("========================================\n");

  /* update Stokes solver for forward problem
   * Note: The Stokes coefficient needs only to be updated for linear Stokes,
   *       since the nonlinear solver will re-compute the coefficient later. */
  RHEA_ASSERT (!inv_problem->forward_is_outdated);
  rhea_stokes_problem_update_solver (
      stokes_problem,
      !rhea_stokes_problem_is_nonlinear (stokes_problem) /* update if lin. */,
      inv_problem->forward_vel_press /* velocity-pressure for coeff update */,
      0 /* !override_rhs */, NULL, NULL, NULL, NULL);

  /* reset flags for all states */
  inv_problem->forward_is_outdated = 1;
  inv_problem->adjoint_is_outdated = 1;
  inv_problem->incremental_forward_is_outdated = 1;
  inv_problem->incremental_adjoint_is_outdated = 1;

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void
rhea_inversion_assemble_hessian_bfgs_init (
                                        sc_dmatrix_t *inv_hessian_matrix,
                                        ymir_vec_t *gradient_vec,
                                        const double damping,
                                        rhea_inversion_problem_t *inv_problem)
{
  rhea_inversion_param_t *inv_param = inv_problem->inv_param;
  ymir_vec_t         *param_vec_in, *param_vec_out;
  sc_dmatrix_t       *param_vec_reduced;
  double              inv_hessian_scaling = NAN;
  const int           hessian_dim = inv_hessian_matrix->m;
  int                 k;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "damping=%g", damping);

  /* check input */
  RHEA_ASSERT (inv_hessian_matrix->m == inv_hessian_matrix->n);

  /* create work vectors */
  param_vec_in  = rhea_inversion_param_vec_new (inv_param);
  param_vec_out = rhea_inversion_param_vec_new (inv_param);

  /* compute lumped Hessian to extract diagonal entries of its prior term */
  ymir_vec_set_value (param_vec_in, 1.0);
  rhea_inversion_apply_hessian (param_vec_out, param_vec_in, inv_problem,
                                RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND);
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (param_vec_out, inv_param));
  param_vec_reduced = rhea_inversion_param_vec_reduced_new (param_vec_out,
                                                            inv_param);
  RHEA_ASSERT (param_vec_reduced->m == hessian_dim);
  RHEA_ASSERT (param_vec_reduced->n == 1);

  /* fill diagonal of inverse Hessian matrix */
  sc_dmatrix_set_zero (inv_hessian_matrix);
  for (k = 0; k < hessian_dim; k++) {
    RHEA_ASSERT (0.0 < param_vec_reduced->e[k][0]);
    inv_hessian_matrix->e[k][k] = 1.0/param_vec_reduced->e[k][0];
  }

  /* set scaling factor dependent on gradient norm */
  if (NULL != gradient_vec) {
    double              hess_norm, grad_norm;

    /* compute the norm of the Hessian and gradient */
    hess_norm = sqrt (rhea_inversion_param_vec_reduced_ip (
        param_vec_reduced, param_vec_reduced, NULL /* weight */,
        0 /* !normalization wrt. size */));
    grad_norm = sqrt (rhea_inversion_param_vec_reduced_inner_product (
        gradient_vec, gradient_vec, NULL /* weight */, inv_param,
        0 /* !normalization wrt. size */));
    RHEA_GLOBAL_VERBOSEF_FN_TAG (
        __func__, "hessian_norm=%g, gradient_norm=%g, normalization=%g",
        hess_norm, grad_norm, hess_norm/grad_norm);

    /* set scaling for inverse of Hessian to be bounded by ||gradient|| */
    if (isfinite (inv_hessian_scaling)) {
      inv_hessian_scaling *= hess_norm/grad_norm;
    }
    else {
      inv_hessian_scaling = hess_norm/grad_norm;
    }
  }

  /* apply damping to scaling factor */
  if (isfinite (damping) && 0.0 < damping) {
    if (isfinite (inv_hessian_scaling)) {
      inv_hessian_scaling *= damping;
    }
    else {
      inv_hessian_scaling = damping;
    }
  }

  /* scale inverse Hessian matrix */
  if (isfinite (inv_hessian_scaling) && 0.0 < inv_hessian_scaling) {
    RHEA_GLOBAL_INFOF_FN_TAG (__func__, "inv_hessian_scaling=%g",
                              inv_hessian_scaling);
    for (k = 0; k < hessian_dim; k++) {
      inv_hessian_matrix->e[k][k] *= inv_hessian_scaling;
    }
  }

  /* destroy */
  rhea_inversion_param_vec_reduced_destroy (param_vec_reduced);
  rhea_inversion_param_vec_destroy (param_vec_in);
  rhea_inversion_param_vec_destroy (param_vec_out);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static double
rhea_inversion_assemble_hessian_bfgs_update (
                                        sc_dmatrix_t *inv_hessian_matrix,
                                        ymir_vec_t *gradient_diff_vec,
                                        ymir_vec_t *step_vec,
                                        const double step_length,
                                        rhea_inversion_problem_t *inv_problem,
                                        const int reconstruct_init_matrix)
{
  rhea_inversion_param_t *inv_param = inv_problem->inv_param;
  ymir_vec_t         *hessian_step_prod_vec, *gradient_diff_mod_vec;
  sc_dmatrix_t       *vec_reduced;
  sc_dmatrix_t       *grad_diff_vec_reduced, *step_vec_reduced;
  const double        curvature_threshold = 0.2;
  double              curvature_grad, curvature_hess, curvature_effective;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "curvature_threshold=%g, reconstruct_init_matrix=%i",
      curvature_threshold, reconstruct_init_matrix);

  /* compute curvature based on gradient difference */
  curvature_grad =
    step_length * rhea_inversion_param_vec_reduced_inner_product (
      step_vec, gradient_diff_vec, NULL /* weight */, inv_param,
      0 /* !normalization wrt. size */);

  /* compute matvec: H*s */
  hessian_step_prod_vec = rhea_inversion_param_vec_new (inv_param);
  rhea_inversion_param_vec_reduced_solve_matrix (
      hessian_step_prod_vec, inv_hessian_matrix, step_vec, inv_param, NULL);
  ymir_vec_scale (step_length, hessian_step_prod_vec);

  /* compute curvature based on previous Hessian: s^T*H*s */
  curvature_hess =
    step_length * rhea_inversion_param_vec_reduced_inner_product (
      step_vec, hessian_step_prod_vec, NULL /* weight */, inv_param,
      0 /* !normalization wrt. size */);

  /* set modified gradient difference */
  gradient_diff_mod_vec = rhea_inversion_param_vec_new (inv_param);
  if (curvature_threshold*curvature_hess <= curvature_grad) {
    ymir_vec_copy (gradient_diff_vec, gradient_diff_mod_vec);
    curvature_effective = curvature_grad;
  }
  else { /* otherwise enforce sufficiently positive curvature */
    double              theta;

    /* set correction factor for convex combination */
    theta = (1.0 - curvature_threshold) * curvature_hess /
            (curvature_hess - curvature_grad);
    RHEA_ASSERT (0.0 < theta && theta <= 1.0);
    /* set modified gradient difference: theta*gd + (1-theta)*H*s */
    ymir_vec_set_zero (gradient_diff_mod_vec);
    ymir_vec_add (theta, gradient_diff_vec, gradient_diff_mod_vec);
    ymir_vec_add (1.0 - theta, hessian_step_prod_vec, gradient_diff_mod_vec);
    /* set modified curvature: theta*s^T*gd + (1-theta)*s^T*H*s */
    curvature_effective = theta*curvature_grad + (1.0 - theta)*curvature_hess;

    RHEA_GLOBAL_VERBOSEF_FN_TAG (__func__, "curvature_correction=%g", theta);
  }

  RHEA_GLOBAL_INFOF_FN_TAG (
      __func__, "curvature_grad=%g, curvature_hess=%g, curvature_effective=%g",
      curvature_grad, curvature_hess, curvature_effective);

  /* re-construct the initial BFGS approximation of inverse Hessian */
  if (reconstruct_init_matrix) {
    double              grad_diff_norm;

    grad_diff_norm = sqrt (rhea_inversion_param_vec_reduced_inner_product (
        gradient_diff_mod_vec, gradient_diff_mod_vec, NULL /* weight */,
        inv_param, 0 /* !normalization wrt. size */));
    rhea_inversion_assemble_hessian_bfgs_init (
        inv_hessian_matrix, NULL /* !gradient_vec */,
        curvature_effective/grad_diff_norm /* damping */, inv_problem);
  }

  /* create work vectors */
  vec_reduced = rhea_inversion_param_vec_reduced_new (NULL, inv_param);
  RHEA_ASSERT (vec_reduced->m == inv_hessian_matrix->m);
  RHEA_ASSERT (vec_reduced->m == inv_hessian_matrix->n);
  RHEA_ASSERT (vec_reduced->n == 1);
  grad_diff_vec_reduced = rhea_inversion_param_vec_reduced_new (
      gradient_diff_mod_vec, inv_param);
  step_vec_reduced = rhea_inversion_param_vec_reduced_new (
      step_vec, inv_param);

  /* update the BFGS approximation of the inverse Hessian:
   *   H^-1 <- (I - scal*s*gd^T) * H^-1 * (I - scal*gd*s^T) + scal*s*s^T
   * where
   *   s = step_length * step
   *   gd = grad - grad_prev
   *   scal = 1 / s^T*gd
   */
  {
    const double        scal = 1.0/curvature_effective;
    const double       *s = step_vec_reduced->e[0];
    const double       *v = vec_reduced->e[0];
    double            **H = inv_hessian_matrix->e;
    int                 row, col;

    /* compute: M <- (I - scal*s*gd^T) * M */
    sc_dmatrix_vector (
        SC_TRANS, SC_NO_TRANS, SC_NO_TRANS,
        1.0 /* factor for in vec */, inv_hessian_matrix, grad_diff_vec_reduced,
        0.0 /* factor for out vec */, vec_reduced);
    for (row = 0; row < inv_hessian_matrix->m; row++) {
      for (col = 0; col < inv_hessian_matrix->n; col++) {
        H[row][col] -= scal * step_length*s[row] * v[col];
      }
    }

    /* compute: M <- M * (I - scal*gd*s^T) */
    sc_dmatrix_vector (
        SC_NO_TRANS, SC_NO_TRANS, SC_NO_TRANS,
        1.0 /* factor for in vec */, inv_hessian_matrix, grad_diff_vec_reduced,
        0.0 /* factor for out vec */, vec_reduced);
    for (row = 0; row < inv_hessian_matrix->m; row++) {
      for (col = 0; col < inv_hessian_matrix->n; col++) {
        H[row][col] -= scal * step_length*s[col] * v[row];
      }
    }

    /* compute: M <- M + scal*s*s^T */
    for (row = 0; row < inv_hessian_matrix->m; row++) {
      for (col = 0; col < inv_hessian_matrix->n; col++) {
        H[row][col] += scal * step_length*s[col] * step_length*s[row];
      }
    }
  }

  /* destroy */
  rhea_inversion_param_vec_reduced_destroy (vec_reduced);
  rhea_inversion_param_vec_reduced_destroy (grad_diff_vec_reduced);
  rhea_inversion_param_vec_reduced_destroy (step_vec_reduced);
  rhea_inversion_param_vec_destroy (hessian_step_prod_vec);
  rhea_inversion_param_vec_destroy (gradient_diff_mod_vec);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);

  /* return curvature from BFGS method */
  return curvature_effective;
}

static void
rhea_inversion_assemble_hessian_gradient_descend (
                                        sc_dmatrix_t *inv_hessian_matrix,
                                        ymir_vec_t *gradient_vec,
                                        const double damping,
                                        rhea_inversion_problem_t *inv_problem)
{
  const int           hessian_dim = inv_hessian_matrix->m;
  int                 k;
  double              max_val;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "damping=%g", damping);

  /* use the initial BFGS approximation of the inverse Hessian matrix as an
   * estimate for the Hessian's scale */
  rhea_inversion_assemble_hessian_bfgs_init (
      inv_hessian_matrix, gradient_vec, damping, inv_problem);

  /* compute harmonic mean of diagonal entries of inverse Hessian matrix */
  max_val = 0.0;
  for (k = 0; k < hessian_dim; k++) {
    RHEA_ASSERT (0.0 < inv_hessian_matrix->e[k][k]);
    max_val = SC_MAX (max_val, inv_hessian_matrix->e[k][k]);
  }

  /* fill diagonal of inverse Hessian matrix: `H = mean * Id` */
  sc_dmatrix_set_zero (inv_hessian_matrix);
  for (k = 0; k < hessian_dim; k++) {
    inv_hessian_matrix->e[k][k] = max_val;
  }

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/**
 * Assembles the Hessian matrix.
 */
static void
rhea_inversion_assemble_hessian (sc_dmatrix_t *hessian_matrix,
                                 rhea_inversion_problem_t *inv_problem,
                                 const rhea_inversion_hessian_t type,
                                 const int enforce_symmetry)
{
  rhea_inversion_param_t *inv_param = inv_problem->inv_param;
  const int           n_parameters =
                        rhea_inversion_param_get_n_parameters (inv_param);
  const int          *active = rhea_inversion_param_get_active (inv_param);
  ymir_vec_t         *param_vec_in, *param_vec_out;
  sc_dmatrix_t       *param_vec_reduced;
  int                 i, row, col;

  /* create work vectors */
  param_vec_in = rhea_inversion_param_vec_new (inv_param);
  param_vec_out = rhea_inversion_param_vec_new (inv_param);

  /* assemble Hessian matrix */
  row = 0;
  col = 0;
  for (i = 0; i < n_parameters; i++) { /* loop over all parameters */
    if (active[i]) {
      /* apply unit vector to Hessian */
      ymir_vec_set_zero (param_vec_in);
      param_vec_in->meshfree->e[0][i] = 1.0;
      rhea_inversion_apply_hessian (param_vec_out, param_vec_in, inv_problem,
                                    type);
      RHEA_ASSERT (rhea_inversion_param_vec_is_valid (param_vec_out,
                                                      inv_param));
      param_vec_reduced = rhea_inversion_param_vec_reduced_new (param_vec_out,
                                                                inv_param);
      RHEA_ASSERT (param_vec_reduced->m == hessian_matrix->m);
      RHEA_ASSERT (param_vec_reduced->n == 1);
      RHEA_ASSERT (col < hessian_matrix->n);

      /* fill one column of the Hessian matrix */
      for (row = 0; row < param_vec_reduced->m; row++) {
        hessian_matrix->e[row][col] = param_vec_reduced->e[row][0];
      }
      rhea_inversion_param_vec_reduced_destroy (param_vec_reduced);

      /* increment column index of matrix */
      col++;
    }
  }

  /* enforce symmetry of Hessian matrix */
  if (rhea_inversion_assemble_hessian_enforce_symm) {
    sc_dmatrix_t       *hessian_transp;

    hessian_transp = sc_dmatrix_new (hessian_matrix->n, hessian_matrix->m);
    sc_dmatrix_transpose (hessian_matrix, hessian_transp);
    RHEA_ASSERT (hessian_matrix->m == hessian_matrix->n);
    sc_dmatrix_add (1.0, hessian_transp, hessian_matrix);
    sc_dmatrix_scale (0.5, hessian_matrix);

    sc_dmatrix_destroy (hessian_transp);
  }

  /* destroy */
  rhea_inversion_param_vec_destroy (param_vec_in);
  rhea_inversion_param_vec_destroy (param_vec_out);
}

static void
rhea_inversion_newton_update_hessian_fn (ymir_vec_t *solution,
                                         ymir_vec_t *step_vec,
                                         const double step_length,
                                         void *data)
{
  const rhea_inversion_hessian_t  type = (rhea_inversion_hessian_t)
                                         rhea_inversion_hessian_type;
  const double              gradient_descent_damping =
                              rhea_inversion_gradient_descent_damping;
  rhea_inversion_problem_t *inv_problem = data;

  RHEA_GLOBAL_INFOF_FN_BEGIN (
      __func__, "type=%i, assemble_matrix=%i", type,
      rhea_inversion_assemble_hessian_matrix);

  /* check input */
  RHEA_ASSERT (!rhea_inversion_assemble_hessian_matrix ||
               inv_problem->hessian_matrix != NULL);

  /* update Hessian */
  switch (type) {
  case RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND:
    RHEA_ASSERT (rhea_inversion_assemble_hessian_matrix);
    RHEA_ASSERT (NULL != inv_problem->newton_gradient_diff_vec);
    if (1 == inv_problem->newton_gradient_diff_vec_update_count) {
      /* set inverse Hessian matrix */
      rhea_inversion_assemble_hessian_gradient_descend (
          inv_problem->hessian_matrix, inv_problem->newton_gradient_diff_vec,
          gradient_descent_damping, inv_problem);
    }
    else if (step_length < 0.5) { /* if adjust scaling of Hessian */
      sc_dmatrix_scale (step_length, inv_problem->hessian_matrix);
    }
    break;
  case RHEA_INVERSION_HESSIAN_BFGS:
    RHEA_ASSERT (rhea_inversion_assemble_hessian_matrix);
    RHEA_ASSERT (NULL != inv_problem->newton_gradient_diff_vec);
    if (1 == inv_problem->newton_gradient_diff_vec_update_count) {
      /* compute initial BFGS approximation of the inverse Hessian matrix */
      rhea_inversion_assemble_hessian_bfgs_init (
          inv_problem->hessian_matrix, inv_problem->newton_gradient_diff_vec,
          gradient_descent_damping, inv_problem);
    }
    else { /* otherwise at 2nd nonlinear iteration and above */
      const int           reconstruct_init =
        (2 == inv_problem->newton_gradient_diff_vec_update_count);

      /* update BFGS approximation of the inverse Hessian matrix */
      if (NULL != step_vec) { /* if step given as argument */
        rhea_inversion_assemble_hessian_bfgs_update (
            inv_problem->hessian_matrix, inv_problem->newton_gradient_diff_vec,
            step_vec, step_length, inv_problem, reconstruct_init);
      }
      else if (NULL != rhea_inversion_step_prev_file_path_txt) { /* if read */
        rhea_inversion_assemble_hessian_bfgs_update (
            inv_problem->hessian_matrix, inv_problem->newton_gradient_diff_vec,
            inv_problem->newton_step_vec, 1.0 /* step_length */, inv_problem,
            reconstruct_init);
      }
      else { /* otherwise unknown state */
        RHEA_ABORT_NOT_REACHED ();
      }
    }
    break;
  case RHEA_INVERSION_HESSIAN_FIRST_ORDER_APPROX:
  case RHEA_INVERSION_HESSIAN_FULL:
    if (rhea_inversion_assemble_hessian_matrix) { /* if Hessian is assembled */
      rhea_inversion_assemble_hessian (
          inv_problem->hessian_matrix, inv_problem, type,
          rhea_inversion_assemble_hessian_enforce_symm);
    }
    break;
  default: /* unknown Hessian type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* print & check output */
#ifdef RHEA_ENABLE_DEBUG
  if (inv_problem->hessian_matrix != NULL) {
    sc_dmatrix_t       *mat = inv_problem->hessian_matrix;
    int                 row, col;

    RHEA_GLOBAL_VERBOSE ("========================================\n");
    RHEA_GLOBAL_VERBOSEF ("Inversion Hessian (type=%i)\n", type);
    RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
    for (row = 0; row < mat->m; row++) {
      char                line[BUFSIZ] = "";
      char               *pos = line;
      const char         *end = line + BUFSIZ;

      for (col = 0; col < mat->n; col++) {
        pos += snprintf (pos, end - pos, "  %+.6e", mat->e[row][col]);
      }
      RHEA_GLOBAL_VERBOSEF ("%s\n", line);
    }
    RHEA_GLOBAL_VERBOSE ("========================================\n");
    RHEA_ASSERT (sc_dmatrix_is_valid (mat));
  }
#endif

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

static double
rhea_inversion_newton_evaluate_objective_fn (ymir_vec_t *solution, void *data,
                                             double *obj_comp)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;
  double              data_abs_weight;
  double              prior_abs_weight;
  ymir_vec_t         *vel;
  double              obj_data_misfit[3], obj_prior_misfit, obj_val;

  /* set weights if they do not exist */
  if (!isfinite (inv_problem->data_abs_weight)) {
    rhea_inversion_set_weight_data (inv_problem);
  }
  if (!isfinite (inv_problem->prior_abs_weight)) {
    rhea_inversion_set_weight_prior (inv_problem);
  }
  data_abs_weight = inv_problem->data_abs_weight;
  prior_abs_weight = inv_problem->prior_abs_weight;

  RHEA_GLOBAL_INFOF_FN_BEGIN (__func__, "weights=[%.6e,%.6e]",
                              data_abs_weight, prior_abs_weight);

  /* solve for forward state */
  rhea_inversion_inner_solve_forward (inv_problem);

  /* retrieve velocity of the forward state */
  {
    ymir_mesh_t              *ymir_mesh =
      rhea_stokes_problem_get_ymir_mesh (stokes_problem);
    ymir_pressure_elem_t     *press_elem =
      rhea_stokes_problem_get_press_elem (stokes_problem);

    vel = rhea_velocity_new (ymir_mesh);
    rhea_velocity_pressure_copy_components (
        vel, NULL, inv_problem->forward_vel_press, press_elem);
    rhea_stokes_problem_velocity_enforce_boundary_conditions (
        vel, stokes_problem);
    RHEA_ASSERT (rhea_velocity_is_valid (vel));
  }

  /* compute the observation data misfit term */
  if (0.0 < data_abs_weight) {
    obj_data_misfit[0] = rhea_inversion_obs_velocity_misfit (
        vel, inv_problem->vel_obs_surf, inv_problem->vel_obs_weight_surf,
        inv_problem->vel_obs_type,
        rhea_stokes_problem_get_domain_options (stokes_problem));
    obj_data_misfit[1] = rhea_inversion_obs_viscosity_misfit (
        vel, inv_problem->visc_obs_type, inv_problem->visc_obs_n,
        inv_problem->visc_obs_column, inv_problem->visc_obs_value,
        inv_problem->visc_obs_weight, stokes_problem);
    obj_data_misfit[2] = rhea_inversion_obs_stress_misfit (
        inv_problem->forward_vel_press, inv_problem->stress_obs,
        inv_problem->stress_obs_weight, inv_problem->stress_obs_type,
        stokes_problem);
    obj_data_misfit[0] *= data_abs_weight;
    obj_data_misfit[1] *= data_abs_weight;
    obj_data_misfit[2] *= data_abs_weight;
  }
  else {
    obj_data_misfit[0] = 0.0;
    obj_data_misfit[1] = 0.0;
    obj_data_misfit[2] = 0.0;
  }
  rhea_velocity_destroy (vel);

  /* compute prior term */
  if (0.0 < prior_abs_weight) {
    obj_prior_misfit =
      prior_abs_weight * rhea_inversion_param_prior (solution, inv_param);
  }
  else {
    obj_prior_misfit = 0.0;
  }

  /* add terms to get full objective value */
  obj_val = obj_data_misfit[0] + obj_data_misfit[1] + obj_data_misfit[2] +
            obj_prior_misfit;

  RHEA_GLOBAL_INFOF_FN_END (
      __func__, "value=%.15e, components[%.15e, %.15e, %.15e; %.15e]",
      obj_val, obj_data_misfit[0], obj_data_misfit[1], obj_data_misfit[2],
      obj_prior_misfit);

  /* return value of objective functional */
  if (NULL != obj_comp) {
    obj_comp[0] = obj_data_misfit[0];
    obj_comp[1] = obj_data_misfit[1];
    obj_comp[2] = obj_data_misfit[2];
    obj_comp[3] = obj_prior_misfit;
  }
  return obj_val;
}

static int
rhea_inversion_newton_evaluate_objective_err_fn (double *error_mean,
                                                 double *error_stddev,
                                                 void *data)
{
  rhea_inversion_problem_t *inv_problem = data;

  *error_mean = inv_problem->error_stats_objective;
  *error_stddev = NAN;
  return 1;
}

static double
_param_derivative_from_observations (
                                  ymir_vec_t *forward_vel_press,
                                  ymir_vec_t *adjoint_vel_press,
                                  const int parameter_idx,
                                  const int derivative_type,
                                  const rhea_weakzone_label_t weak_label,
                                  ymir_vec_t *coeff_param_derivative,
                                  ymir_stress_op_t *stress_op_param_derivative,
                                  void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  double              val = 0.0;

  /* gather derivatives from observation operators */
  val += rhea_inversion_obs_velocity_misfit_param_derivative (
      inv_problem->vel_obs_type);
  val += rhea_inversion_obs_viscosity_misfit_param_derivative (
      inv_problem->visc_obs_type);
  val += rhea_inversion_obs_stress_misfit_param_derivative (
      forward_vel_press, inv_problem->stress_obs,
      inv_problem->stress_obs_weight, inv_problem->stress_obs_type,
      inv_problem->stokes_problem, stress_op_param_derivative);

  return val;
}

static void
rhea_inversion_newton_compute_negative_gradient_fn (
                                            ymir_vec_t *neg_gradient,
                                            ymir_vec_t *solution, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t              *ymir_mesh;
  ymir_pressure_elem_t     *press_elem;
  const rhea_inversion_project_out_null_t project_out_null =
                                            inv_problem->project_out_null;
  const double        data_abs_weight = inv_problem->data_abs_weight;
  const double        prior_abs_weight = inv_problem->prior_abs_weight;
  ymir_vec_t         *gradient_adj_comp =
                        inv_problem->newton_gradient_adjoint_comp_vec;
  ymir_vec_t         *gradient_prior_comp =
                        inv_problem->newton_gradient_prior_comp_vec;
  ymir_vec_t         *gradient_obs_comp =
                        inv_problem->newton_gradient_observation_comp_vec;
  ymir_vec_t         *rhs_vel_mass  = inv_problem->rhs_vel_buffer;
  ymir_vec_t         *rhs_vel_press = inv_problem->rhs_vel_press_buffer;
  ymir_vec_t         *vel;

  RHEA_GLOBAL_INFO_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (neg_gradient, inv_param));

  /* begin the computation of gradient difference: g - g_prev */
  switch (rhea_inversion_hessian_type) {
  case RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND:
  case RHEA_INVERSION_HESSIAN_BFGS:
    RHEA_ASSERT (NULL != inv_problem->newton_gradient_diff_vec);
    if (0 < inv_problem->newton_gradient_diff_vec_update_count) {
      ymir_vec_copy (neg_gradient, inv_problem->newton_gradient_diff_vec);
    }
    else {
      ymir_vec_set_zero (inv_problem->newton_gradient_diff_vec);
    }
    break;
  case RHEA_INVERSION_HESSIAN_FIRST_ORDER_APPROX:
  case RHEA_INVERSION_HESSIAN_FULL:
    break;
  default: /* unknown Hessian type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* ensure the forward state is up-to-date */
  rhea_inversion_inner_solve_forward (inv_problem);

  /* get mesh data */
  ymir_mesh  = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* retrieve velocity of the forward state */
  vel = rhea_velocity_new (ymir_mesh);
  rhea_velocity_pressure_copy_components (
      vel, NULL, inv_problem->forward_vel_press, press_elem);
  rhea_stokes_problem_velocity_enforce_boundary_conditions (
      vel, stokes_problem);
  RHEA_ASSERT (rhea_velocity_is_valid (vel));

  /* set the right-hand side of the momentum eq. for the adjoint problem */
  ymir_vec_set_zero (rhs_vel_mass);
  rhea_inversion_obs_velocity_add_adjoint_rhs (
      rhs_vel_mass, vel, inv_problem->vel_obs_surf,
      inv_problem->vel_obs_weight_surf, inv_problem->vel_obs_type,
      rhea_stokes_problem_get_domain_options (stokes_problem));
  rhea_inversion_obs_viscosity_add_adjoint_rhs (
      rhs_vel_mass, vel, inv_problem->visc_obs_type, inv_problem->visc_obs_n,
      inv_problem->visc_obs_column, inv_problem->visc_obs_value,
      inv_problem->visc_obs_weight, stokes_problem);
  rhea_inversion_obs_stress_add_adjoint_rhs (
      rhs_vel_mass, inv_problem->forward_vel_press,
      inv_problem->stress_obs, inv_problem->stress_obs_weight,
      inv_problem->stress_obs_type, stokes_problem);
  ymir_vec_scale (data_abs_weight, rhs_vel_mass);
  rhea_velocity_destroy (vel);

  /* project out null spaces and enforce Dirichlet BC's on right-hand side */
  if (project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_RESIDUAL ||
      project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_BOTH) {
    rhea_stokes_problem_velocity_project_out_mean_rotation (
        rhs_vel_mass, 1 /* residual_space */, stokes_problem);
  }
  rhea_stokes_problem_velocity_set_boundary_zero (
      rhs_vel_mass, stokes_problem);

  /* set Stokes right-hand side for the adjoint problem */
  ymir_vec_set_zero (rhs_vel_press);
  rhea_velocity_pressure_set_components (rhs_vel_press, rhs_vel_mass, NULL,
                                         press_elem);

  /* update Stokes solver for adjoint problem
   * Note: The Stokes coefficient is updated only for nonlinear Stokes, since
   *       the coefficient depends on the forward velocity. */
  rhea_stokes_problem_update_solver (
      stokes_problem,
      rhea_stokes_problem_is_nonlinear (stokes_problem) /* update if nl. */,
      inv_problem->forward_vel_press /* velocity-pressure for coeff update */,
      1 /* override_rhs */, rhs_vel_press, NULL, NULL, NULL);

  /* solve for adjoint state */
  rhea_inversion_inner_solve_adjoint (inv_problem);
  if (project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_PRIMAL ||
      project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_BOTH) {
    rhea_stokes_problem_project_out_nullspace (
        inv_problem->adjoint_vel_press, stokes_problem);
  }

  /* write vtk output pertaining to adjoint */
  if (inv_problem->vtk_path_debug_gradient != NULL) {
    ymir_vec_t         *adj_vel, *adj_press;

    adj_vel   = rhea_velocity_new (ymir_mesh);
    adj_press = rhea_pressure_new (ymir_mesh, press_elem);
    rhea_velocity_pressure_copy_components (
        adj_vel, adj_press, inv_problem->adjoint_vel_press, press_elem);

    rhea_velocity_remove_mass (rhs_vel_mass);
    ymir_vec_share_owned (rhs_vel_mass);

    rhea_vtk_write_inversion_adjoint (
        inv_problem->vtk_path_debug_gradient,
        adj_vel, adj_press, rhs_vel_mass, NULL /* rhs_pressure_adj */);

    rhea_velocity_destroy (adj_vel);
    rhea_pressure_destroy (adj_press);
  }

  /* compute the (positive) gradient */
  rhea_inversion_param_compute_gradient (
      neg_gradient, solution, inv_problem->forward_vel_press,
      inv_problem->adjoint_vel_press, prior_abs_weight, inv_param,
      gradient_adj_comp, gradient_prior_comp, gradient_obs_comp,
      _param_derivative_from_observations, inv_problem);

  /* print gradient */
  {
    const int          *active = rhea_inversion_param_get_active (inv_param);
    const double       *grad       = neg_gradient->meshfree->e[0];
    const double       *grad_adj   = gradient_adj_comp->meshfree->e[0];
    const double       *grad_prior = gradient_prior_comp->meshfree->e[0];
    const double       *grad_obs   = gradient_obs_comp->meshfree->e[0];
    int                 i;

    RHEA_GLOBAL_INFO ("========================================\n");
    RHEA_GLOBAL_INFO ("Inversion gradient\n");
    RHEA_GLOBAL_INFO ("----------------------------------------\n");
    for (i = 0; i < rhea_inversion_param_get_n_parameters (inv_param); i++) {
      if (active[i]) {
        RHEA_GLOBAL_INFOF ("param# %3i: %+.6e [%+.6e, %+.6e, %+.6e]\n", i,
                           grad[i], grad_adj[i], grad_prior[i], grad_obs[i]);
      }
    }
    RHEA_GLOBAL_INFO ("========================================\n");
  }

  /* finalize the computation of gradient difference: g - g_prev */
  switch (rhea_inversion_hessian_type) {
  case RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND:
  case RHEA_INVERSION_HESSIAN_BFGS:
    RHEA_ASSERT (NULL != inv_problem->newton_gradient_diff_vec);
    ymir_vec_add (1.0, neg_gradient, inv_problem->newton_gradient_diff_vec);
    inv_problem->newton_gradient_diff_vec_update_count++;
    break;
  case RHEA_INVERSION_HESSIAN_FIRST_ORDER_APPROX:
  case RHEA_INVERSION_HESSIAN_FULL:
    break;
  default: /* unknown Hessian type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* flip the sign of the (positive) gradient */
  ymir_vec_scale (-1.0, neg_gradient);

  /* check output */
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (neg_gradient, inv_param));

  /* set up perturbations for the (elementwise) gradient check */
  if ((rhea_inversion_check_gradient || rhea_inversion_check_gradient_of_qoi) &&
      rhea_inversion_check_gradient_elementwise) {
    const double        restrict_to_prior_stddev =
                          rhea_inversion_restrict_step_to_prior_stddev;
    const int          *active = rhea_inversion_param_get_active (inv_param);
    ymir_vec_t        **perturb_vec = inv_problem->check_gradient_perturb_vec;
    int                 vidx, i;

    RHEA_ASSERT (NULL != inv_problem->check_gradient_perturb_vec);
    RHEA_ASSERT (NULL != solution);

    for (vidx = 0; vidx < inv_problem->check_gradient_perturb_n_vecs; vidx++) {
      ymir_vec_set_zero (perturb_vec[vidx]);
    }
    vidx = 0;
    for (i = 0; i < rhea_inversion_param_get_n_parameters (inv_param); i++) {
      if (active[i]) {
        RHEA_ASSERT (vidx < inv_problem->check_gradient_perturb_n_vecs);
        /* set perturbation for an entry in the solution vector;
         * assume: `solution->meshfree->e[0][i]` is of order 1 for all i */
        perturb_vec[vidx]->meshfree->e[0][i] = 0.1;
        vidx++;
      }
    }
    RHEA_ASSERT (vidx == inv_problem->check_gradient_perturb_n_vecs);
    for (vidx = 0; vidx < inv_problem->check_gradient_perturb_n_vecs; vidx++) {
      rhea_inversion_param_restrict_step_length_to_feasible (
          perturb_vec[vidx], solution, restrict_to_prior_stddev, inv_param,
          NULL /* step_length_new */);
    }
  }

  RHEA_GLOBAL_INFO_FN_END (__func__);
}

static double
rhea_inversion_newton_compute_gradient_norm_fn (ymir_vec_t *neg_gradient,
                                                void *data, double *norm_comp)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  double              grad_norm, grad_adj_norm, grad_prior_norm, grad_obs_norm;

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (neg_gradient, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (neg_gradient, inv_param));

  /* compute norms */
  grad_norm = rhea_inversion_param_compute_gradient_norm (
      neg_gradient, inv_param);
  grad_adj_norm = rhea_inversion_param_compute_gradient_norm (
      inv_problem->newton_gradient_adjoint_comp_vec, inv_param);
  grad_prior_norm = rhea_inversion_param_compute_gradient_norm (
      inv_problem->newton_gradient_prior_comp_vec, inv_param);
  grad_obs_norm = rhea_inversion_param_compute_gradient_norm (
      inv_problem->newton_gradient_observation_comp_vec, inv_param);
  RHEA_GLOBAL_VERBOSEF_FN_TAG (
      __func__, "gradient_norm=%.6e, gradient_adjoint_comp=%.6e, "
      "gradient_prior_comp=%.6e, gradient_observation_comp=%.6e",
      grad_norm, grad_adj_norm, grad_prior_norm, grad_obs_norm);

  /* return gradient norms */
  if (norm_comp != NULL) {
    norm_comp[0] = grad_adj_norm + grad_obs_norm;
    norm_comp[1] = grad_prior_norm;
  }
  return grad_norm;
}

static void
rhea_inversion_apply_hessian (ymir_vec_t *param_vec_out,
                              ymir_vec_t *param_vec_in,
                              rhea_inversion_problem_t *inv_problem,
                              const rhea_inversion_hessian_t type)
{
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t              *ymir_mesh;
  ymir_pressure_elem_t     *press_elem;
  const rhea_inversion_project_out_null_t project_out_null =
                                            inv_problem->project_out_null;
  const double        data_abs_weight  = inv_problem->data_abs_weight;
  const double        prior_abs_weight = inv_problem->prior_abs_weight;
  int                 compute_incr_fwd, compute_incr_adj;
  ymir_vec_t         *rhs_vel_mass  = inv_problem->rhs_vel_buffer;
  ymir_vec_t         *rhs_vel_press = inv_problem->rhs_vel_press_buffer;

  /* set which computations to perform */
  switch (type) {
  case RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND:
  case RHEA_INVERSION_HESSIAN_BFGS:
    compute_incr_fwd = 0;
    compute_incr_adj = 0;
    break;
  case RHEA_INVERSION_HESSIAN_FIRST_ORDER_APPROX:
  case RHEA_INVERSION_HESSIAN_FULL:
    compute_incr_fwd = 1;
    compute_incr_adj = 1;
    break;
  default: /* unknown Hessian type */
    compute_incr_fwd = 0;
    compute_incr_adj = 0;
    RHEA_ABORT_NOT_REACHED ();
  }

  /* set up incremental solves */
  if (compute_incr_fwd || compute_incr_adj) {
    /* ensure the forward state is up-to-date */
    rhea_inversion_inner_solve_forward (inv_problem);

    /* get mesh data */
    ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
    press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);
  }

  if (compute_incr_fwd) { /* BEGIN: incremental forward */
    /* set right-hand side of the momentum eq. for the incr. forward problem */
    RHEA_ASSERT (!inv_problem->forward_is_outdated);
    rhea_inversion_param_incremental_forward_rhs (
        rhs_vel_mass, param_vec_in, inv_problem->forward_vel_press, inv_param);

    /* project out null spaces and enforce Dirichlet BC's on right-hand side */
    if (project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_RESIDUAL ||
        project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_BOTH) {
      rhea_stokes_problem_velocity_project_out_mean_rotation (
          rhs_vel_mass, 1 /* residual_space */, stokes_problem);
    }
    rhea_stokes_problem_velocity_set_boundary_zero (
        rhs_vel_mass, stokes_problem);

    /* set Stokes right-hand side for the incremental forward problem */
    ymir_vec_set_zero (rhs_vel_press);
    rhea_velocity_pressure_set_components (rhs_vel_press, rhs_vel_mass, NULL,
                                           press_elem);

    /* update Stokes right-hand side for incremental forward problem */
    rhea_stokes_problem_update_solver (
        stokes_problem, 0 /* !update coeff */, NULL /* since no coeff update */,
        1 /* override_rhs */, rhs_vel_press, NULL, NULL, NULL);

    /* solve for incremental forward state */
    rhea_inversion_inner_solve_incremental_forward (inv_problem);
    if (project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_PRIMAL ||
        project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_BOTH) {
      rhea_stokes_problem_project_out_nullspace (
          inv_problem->incremental_forward_vel_press, stokes_problem);
    }
  } /* END: incremental forward */

  if (compute_incr_adj) { /* BEGIN: incremental adjoint */
    ymir_vec_t         *vel;

    /* retrieve velocity of the incremental forward state */
    vel = rhea_velocity_new (ymir_mesh);
    rhea_velocity_pressure_copy_components (
        vel, NULL, inv_problem->incremental_forward_vel_press, press_elem);
    rhea_stokes_problem_velocity_set_boundary_zero (vel, stokes_problem);
    RHEA_ASSERT (rhea_velocity_is_valid (vel));

    /* set right-hand side of the momentum eq. for the incr. adjoint problem */
    rhea_inversion_obs_velocity_incremental_adjoint_rhs (
        rhs_vel_mass, vel, inv_problem->vel_obs_weight_surf,
        inv_problem->vel_obs_type,
        rhea_stokes_problem_get_domain_options (stokes_problem));
    //TODO add rhs from viscosity obs.
    ymir_vec_scale (data_abs_weight, rhs_vel_mass);

    /* add 2nd-order terms to right-hand side */
    switch (type) {
    case RHEA_INVERSION_HESSIAN_FIRST_ORDER_APPROX:
      break;
    case RHEA_INVERSION_HESSIAN_FULL:
      /* retrieve velocity of the adjont state */
      RHEA_ASSERT (!inv_problem->adjoint_is_outdated);
      rhea_velocity_pressure_copy_components (
          vel, NULL, inv_problem->adjoint_vel_press, press_elem);
      rhea_stokes_problem_velocity_set_boundary_zero (vel, stokes_problem);
      RHEA_ASSERT (rhea_velocity_is_valid (vel));

      RHEA_ABORT_NOT_REACHED (); //TODO missing implementation
      break;
    default: /* unknown Hessian type */
      RHEA_ABORT_NOT_REACHED ();
    }
    rhea_velocity_destroy (vel);

    /* project out null spaces and enforce Dirichlet BC's on right-hand side */
    if (project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_RESIDUAL ||
        project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_BOTH) {
      rhea_stokes_problem_velocity_project_out_mean_rotation (
          rhs_vel_mass, 1 /* residual_space */, stokes_problem);
    }
    rhea_stokes_problem_velocity_set_boundary_zero (
        rhs_vel_mass, stokes_problem);

    /* set Stokes right-hand side for incremental adjoint problem */
    ymir_vec_set_zero (rhs_vel_press);
    rhea_velocity_pressure_set_components (rhs_vel_press, rhs_vel_mass, NULL,
                                           press_elem);

    /* update Stokes right-hand side for incremental adjoint problem */
    rhea_stokes_problem_update_solver (
        stokes_problem, 0 /* !update coeff */, NULL /* since no coeff update */,
        1 /* override_rhs */, rhs_vel_press, NULL, NULL, NULL);

    /* solve for incremental adjoint state */
    rhea_inversion_inner_solve_incremental_adjoint (inv_problem);
    if (project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_PRIMAL ||
        project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_BOTH) {
      rhea_stokes_problem_project_out_nullspace (
          inv_problem->incremental_adjoint_vel_press, stokes_problem);
    }
  } /* END: incremental adjoint */

  /* compute the Hessian application */
  switch (type) {
  case RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND:
  case RHEA_INVERSION_HESSIAN_BFGS:
    //TODO this is not the proper implementation of matrix-free BFGS
    RHEA_ASSERT (rhea_inversion_assemble_hessian_matrix);
    rhea_inversion_param_apply_hessian (
        param_vec_out, param_vec_in, NULL, NULL, NULL, NULL,
        prior_abs_weight, inv_param);
    break;
  case RHEA_INVERSION_HESSIAN_FIRST_ORDER_APPROX:
    rhea_inversion_param_apply_hessian (
        param_vec_out, param_vec_in,
        inv_problem->forward_vel_press, NULL /* adjoint */,
        NULL /* incr. forward */, inv_problem->incremental_adjoint_vel_press,
        prior_abs_weight, inv_param);
    break;
  case RHEA_INVERSION_HESSIAN_FULL:
    rhea_inversion_param_apply_hessian (
        param_vec_out, param_vec_in,
        inv_problem->forward_vel_press,
        inv_problem->adjoint_vel_press,
        inv_problem->incremental_forward_vel_press,
        inv_problem->incremental_adjoint_vel_press,
        prior_abs_weight, inv_param);
    break;
  default: /* unknown Hessian type */
    RHEA_ABORT_NOT_REACHED ();
  }
}

static void
rhea_inversion_newton_apply_hessian_fn (ymir_vec_t *out, ymir_vec_t *in,
                                        void *data)
{
  const rhea_inversion_hessian_t  type = (rhea_inversion_hessian_t)
                                         rhea_inversion_hessian_type;
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "type=%i, assemble_matrix=%i", type,
      rhea_inversion_assemble_hessian_matrix);

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (out, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (in, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (in, inv_param));
  RHEA_ASSERT (!rhea_inversion_assemble_hessian_matrix ||
               inv_problem->hessian_matrix != NULL);

  /* apply the Hessian */
  switch (type) {
  case RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND:
    RHEA_ASSERT (rhea_inversion_assemble_hessian_matrix);
    ymir_vec_copy (in, out);
    ymir_vec_scale (1.0/inv_problem->hessian_matrix->e[0][0], out);
    break;
  case RHEA_INVERSION_HESSIAN_BFGS:
    //TODO this factorization should not be performed for large matrices
    RHEA_ASSERT (rhea_inversion_assemble_hessian_matrix);
    rhea_inversion_param_vec_reduced_solve_matrix (
        out, inv_problem->hessian_matrix, in, inv_param, NULL);
    break;
  case RHEA_INVERSION_HESSIAN_FIRST_ORDER_APPROX:
  case RHEA_INVERSION_HESSIAN_FULL:
    if (!rhea_inversion_assemble_hessian_matrix) { /* if call apply function */
      rhea_inversion_apply_hessian (out, in, inv_problem, type);
    }
    else { /* otherwise use assembled Hessian matrix */
      rhea_inversion_param_vec_reduced_apply_matrix (
          out, inv_problem->hessian_matrix, in, inv_param);
    }
    break;
  default: /* unknown Hessian type */
    RHEA_ABORT_NOT_REACHED ();
  }
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (out, inv_param));

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

#if 0 //TODO unused
/**
 * Applies prior preconditioner to Hessian--gradient system.
 */
static void
rhea_inversion_apply_prior_preconditioner (sc_dmatrix_t *hessian_matrix,
                                           sc_dmatrix_t *neg_gradient_rhs,
                                           const sc_dmatrix_t *prior_icov)
{
  const size_t        n_rows = hessian_matrix->m;
  const size_t        n_cols = hessian_matrix->n;
  const double       *icov = prior_icov->e[0];
  double            **hessian = hessian_matrix->e;
  int                 row, col;

  /* check input */
  RHEA_ASSERT (hessian_matrix->m == prior_icov->m * prior_icov->n);

  /* apply covariance of prior to the right-hand side vector */
  sc_dmatrix_dotdivide (prior_icov, neg_gradient_rhs);

  /* apply covariance of prior to Hessian matrix */
  for (row = 0; row < n_rows; row++) {
    for (col = 0; col < n_cols; col++) {
      hessian[row][col] *= 1.0/icov[row];
    }
  }
}
#endif

static int
rhea_inversion_newton_solve_hessian_system_fn (
                                            ymir_vec_t *step,
                                            ymir_vec_t *neg_gradient,
                                            const int lin_iter_max,
                                            const double lin_res_norm_rtol,
                                            const int nonzero_initial_guess,
                                            void *data, int *lin_iter_count)
{
  const rhea_inversion_hessian_t  type = (rhea_inversion_hessian_t)
                                         rhea_inversion_hessian_type;
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  int                 stop_reason;

  RHEA_GLOBAL_INFOF_FN_BEGIN (
      __func__, "lin_iter_max=%i, lin_rtol=%.1e, nonzero_init_guess=%i, "
      "type=%i, assemble_matrix=%i", lin_iter_max, lin_res_norm_rtol,
      nonzero_initial_guess, type, rhea_inversion_assemble_hessian_matrix);

  /* check input */
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (step, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_check_type (neg_gradient, inv_param));
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (neg_gradient, inv_param));
  RHEA_ASSERT (!rhea_inversion_assemble_hessian_matrix ||
               inv_problem->hessian_matrix != NULL);

  /* update relative tolerance for inner solves */
  if (rhea_inversion_inner_solver_rtol_adaptive) {
    inv_problem->inner_solver_adaptive_rtol_comp = lin_res_norm_rtol;
  }
  else {
    inv_problem->inner_solver_adaptive_rtol_comp = 1.0;
  }

  /* (approximately) invert the Hessian */
  switch (type) {
  case RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND:
    RHEA_ASSERT (rhea_inversion_assemble_hessian_matrix);
    ymir_vec_copy (neg_gradient, step);
    ymir_vec_scale (inv_problem->hessian_matrix->e[0][0], step);
    if (lin_iter_count != NULL) {
      *lin_iter_count = 1;
    }
    stop_reason = 1;
    break;
  case RHEA_INVERSION_HESSIAN_BFGS:
    RHEA_ASSERT (rhea_inversion_assemble_hessian_matrix);
    rhea_inversion_param_vec_reduced_apply_matrix (
        step, inv_problem->hessian_matrix, neg_gradient, inv_param);
    if (lin_iter_count != NULL) {
      *lin_iter_count = 1;
    }
    stop_reason = 1;
    break;
  case RHEA_INVERSION_HESSIAN_FULL:
  case RHEA_INVERSION_HESSIAN_FIRST_ORDER_APPROX:
    if (!rhea_inversion_assemble_hessian_matrix) { /* if matrix-free */
      RHEA_ABORT_NOT_REACHED (); //TODO implement non-assembled case
      stop_reason = 1;
    }
    else { /* if Hessian matrix is assembled */
      stop_reason = rhea_inversion_param_vec_reduced_solve_matrix (
          step, inv_problem->hessian_matrix, neg_gradient, inv_param,
          lin_iter_count);
    }
    break;
  default: /* unknown Hessian type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* print step */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFO ("Inversion step\n");
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  rhea_inversion_param_vec_print (step, inv_param);
  RHEA_GLOBAL_INFO ("========================================\n");

  /* check output */
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (step, inv_param));

  /* increment index for statistics of inner solves */
  inv_problem->fwd_adj_solve_stats_current_idx++;
  inv_problem->ifwd_iadj_solve_stats_current_idx++;

  RHEA_GLOBAL_INFO_FN_END (__func__);

  /* return iteraton count and "stopping" reason */
  return stop_reason;
}

static double
rhea_inversion_newton_modify_step_fn (ymir_vec_t *step, ymir_vec_t *solution,
                                      const int iter, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  double              restrict_to_prior_stddev;
  int                 step_modified;
  double              step_length_new;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* get value */
  if (0 == iter &&
      0.0 < rhea_inversion_restrict_init_step_to_prior_stddev) {
    restrict_to_prior_stddev =
      rhea_inversion_restrict_init_step_to_prior_stddev;
  }
  else {
    restrict_to_prior_stddev = rhea_inversion_restrict_step_to_prior_stddev;
  }

  /* modify step */
  step_modified = rhea_inversion_param_restrict_step_length_to_feasible (
      step, solution, restrict_to_prior_stddev, inv_param, &step_length_new);
  if (step_modified) {
    RHEA_GLOBAL_INFOF_FN_TAG (__func__, "step_modified=%i, step_length_new=%g",
                              step_modified, step_length_new);

    /* print step */
    RHEA_GLOBAL_INFO ("========================================\n");
    RHEA_GLOBAL_INFOF ("Inversion feasible step, length reduced by %.15e\n",
                       step_length_new);
    RHEA_GLOBAL_INFO ("----------------------------------------\n");
    rhea_inversion_param_vec_print (step, inv_param);
    RHEA_GLOBAL_INFO ("========================================\n");
  }

  /* print relative step */
#ifdef RHEA_ENABLE_DEBUG
  {
    ymir_vec_t         *rel_step = rhea_inversion_param_vec_new (inv_param);

    ymir_vec_fabs (solution, rel_step);
    ymir_vec_reciprocal (rel_step);
    ymir_vec_multiply_in (step, rel_step);

    RHEA_GLOBAL_VERBOSE ("========================================\n");
    RHEA_GLOBAL_VERBOSE ("Inversion relative step\n");
    RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
    rhea_inversion_param_vec_print (rel_step, inv_param);
    RHEA_GLOBAL_VERBOSE ("========================================\n");

    rhea_inversion_param_vec_destroy (rel_step);
  }
#endif

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
  return step_length_new;
}

static void
rhea_inversion_write_txt (ymir_vec_t *solution,
                          const int iter_start, const int iter,
                          rhea_inversion_problem_t *inv_problem)
{
  rhea_inversion_param_t *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t            *ymir_mesh;
  const char         *txt_path = inv_problem->txt_path;
  sc_MPI_Comm         mpicomm;
  int                 mpirank, mpiret;
  sc_dmatrix_t       *sol, *step;
  sc_dmatrix_t       *grad_combined, *neg_grad;
  sc_dmatrix_t       *grad_adj, *grad_prior, *grad_obs;
  int                 i;
  char                path[BUFSIZ];

  /* get parallel environment */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

  /* exit if nothing to do */
  if (txt_path == NULL || mpirank != 0) {
    return;
  }

  /* write parameters */
  sol = rhea_inversion_param_vec_reduced_new (solution, inv_param);
  snprintf (path, BUFSIZ, "%s%s%02i_parameters.txt", txt_path,
            RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
  rhea_io_std_write_double_to_txt (
      path, sol->e[0], sol->m * sol->n, sol->n);
  rhea_inversion_param_vec_reduced_destroy (sol);

  /* write gradient */
  neg_grad = rhea_inversion_param_vec_reduced_new (
      inv_problem->newton_neg_gradient_vec, inv_param);
  grad_adj = rhea_inversion_param_vec_reduced_new (
      inv_problem->newton_gradient_adjoint_comp_vec, inv_param);
  grad_prior = rhea_inversion_param_vec_reduced_new (
      inv_problem->newton_gradient_prior_comp_vec, inv_param);
  grad_obs = rhea_inversion_param_vec_reduced_new (
      inv_problem->newton_gradient_observation_comp_vec, inv_param);
  grad_combined = sc_dmatrix_new (neg_grad->m, 4);
  RHEA_ASSERT (neg_grad->m == grad_adj->m &&
               neg_grad->m == grad_prior->m &&
               neg_grad->m == grad_obs->m);
  for (i = 0; i < grad_combined->m; i++) {
    grad_combined->e[i][0] = -neg_grad->e[i][0];
    grad_combined->e[i][1] = grad_adj->e[i][0];
    grad_combined->e[i][2] = grad_prior->e[i][0];
    grad_combined->e[i][3] = grad_obs->e[i][0];
  }
  snprintf (path, BUFSIZ, "%s%s%02i_gradient.txt", txt_path,
            RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
  rhea_io_std_write_double_to_txt (
      path, grad_combined->e[0], grad_combined->m * grad_combined->n,
      grad_combined->n);
  rhea_inversion_param_vec_reduced_destroy (neg_grad);
  rhea_inversion_param_vec_reduced_destroy (grad_adj);
  rhea_inversion_param_vec_reduced_destroy (grad_prior);
  rhea_inversion_param_vec_reduced_destroy (grad_obs);
  sc_dmatrix_destroy (grad_combined);

  /* write step */
  if (iter_start < iter) {
    step = rhea_inversion_param_vec_reduced_new (
        inv_problem->newton_step_vec, inv_param);
    snprintf (path, BUFSIZ, "%s%s%02i_step.txt", txt_path,
              RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
    rhea_io_std_write_double_to_txt (
        path, step->e[0], step->m * step->n, step->n);
    rhea_inversion_param_vec_reduced_destroy (step);
  }

  /* write Hessian/covariance */
  if (rhea_inversion_assemble_hessian_matrix) {
    if (iter_start == iter) { /* at first iteration */
      ymir_vec_t         *param_vec_in, *param_vec_out;
      sc_dmatrix_t       *param_vec_reduced;
      double             *C;
      const size_t        n = inv_problem->hessian_matrix->n;
      size_t              k;

      /* extract prior covariance */
      param_vec_in = rhea_inversion_param_vec_new (inv_param);
      param_vec_out = rhea_inversion_param_vec_new (inv_param);
      ymir_vec_set_value (param_vec_in, 1.0);
      rhea_inversion_apply_hessian (param_vec_out, param_vec_in, inv_problem,
                                    RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND);
      param_vec_reduced = rhea_inversion_param_vec_reduced_new (param_vec_out,
                                                                inv_param);
      RHEA_ASSERT (n == param_vec_reduced->m &&
                   n == inv_problem->hessian_matrix->m &&
                   n == inv_problem->hessian_matrix->n);

      /* create diagonal covariance matrix */
      C = RHEA_ALLOC_ZERO (double, n*n);
      for (k = 0; k < n; k++) {
        C[n*k + k] = 1.0/param_vec_reduced->e[k][0];
      }

      /* write covariance */
      snprintf (path, BUFSIZ, "%s%s%02i_prior_covariance.txt", txt_path,
                RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
      rhea_io_std_write_double_to_txt (path, C, n*n, n);

      /* destroy */
      RHEA_FREE (C);
      rhea_inversion_param_vec_reduced_destroy (param_vec_reduced);
      rhea_inversion_param_vec_destroy (param_vec_in);
      rhea_inversion_param_vec_destroy (param_vec_out);
    }
    else { /* otherwise at iterations after first */
      snprintf (path, BUFSIZ, "%s%s%02i_hessian.txt", txt_path,
                RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
      rhea_io_std_write_double_to_txt (
          path, inv_problem->hessian_matrix->e[0],
          inv_problem->hessian_matrix->m *
          inv_problem->hessian_matrix->n,
          inv_problem->hessian_matrix->n);
    }
  }
}

static void
rhea_inversion_write_vis (const int iter,
                          rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t      *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t                *ymir_mesh;
  ymir_pressure_elem_t       *press_elem;
  rhea_domain_options_t      *domain_options;
  rhea_temperature_options_t *temp_options;
  rhea_viscosity_options_t   *visc_options;
  rhea_weakzone_options_t    *weak_options;
  const char         *vtk_path_vol = inv_problem->vtk_path_vol;
  const char         *vtk_path_surf = inv_problem->vtk_path_surf;
  ymir_vec_t         *vel_fwd_vol, *press_fwd_vol;
  ymir_vec_t         *vel_adj_vol, *press_adj_vol;
  ymir_vec_t         *viscosity, *marker;
  double              strainrate_dim_1_s;
  char                path[BUFSIZ];

  /* exit if nothing to do */
  if (vtk_path_vol == NULL && vtk_path_surf == NULL) {
    return;
  }

  /* get mesh data */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

   /* get options */
  domain_options = rhea_stokes_problem_get_domain_options (stokes_problem);
  temp_options = rhea_stokes_problem_get_temperature_options (stokes_problem);
  visc_options = rhea_stokes_problem_get_viscosity_options (stokes_problem);
  weak_options = rhea_stokes_problem_get_weakzone_options (stokes_problem);

  /* get volume fields */
  vel_fwd_vol   = rhea_velocity_new (ymir_mesh);
  press_fwd_vol = rhea_pressure_new (ymir_mesh, press_elem);
  vel_adj_vol   = rhea_velocity_new (ymir_mesh);
  press_adj_vol = rhea_pressure_new (ymir_mesh, press_elem);
  viscosity     = rhea_viscosity_new (ymir_mesh);
  marker        = rhea_viscosity_new (ymir_mesh);
  rhea_velocity_pressure_copy_components (
      vel_fwd_vol, press_fwd_vol, inv_problem->forward_vel_press, press_elem);
  rhea_velocity_pressure_copy_components (
      vel_adj_vol, press_adj_vol, inv_problem->adjoint_vel_press, press_elem);
  rhea_stokes_problem_copy_viscosity (viscosity, stokes_problem);
  rhea_stokes_problem_copy_marker (marker, stokes_problem);

  /* process surface fields */
  if (vtk_path_surf != NULL) {
    ymir_vec_t         *vel_fwd_surf, *vel_adj_surf, *vel_obs_surf;
    ymir_vec_t         *misfit_surf;

    /* get surface fields */
    vel_fwd_surf = rhea_velocity_surface_new (ymir_mesh);
    vel_adj_surf = rhea_velocity_surface_new (ymir_mesh);
    rhea_velocity_surface_interpolate (vel_fwd_surf, vel_fwd_vol);
    rhea_velocity_surface_interpolate (vel_adj_surf, vel_adj_vol);

    /* get observations */
    vel_obs_surf = rhea_velocity_surface_new (ymir_mesh);
    ymir_vec_copy (inv_problem->vel_obs_surf, vel_obs_surf);

    /* compute the data misfit term */
    misfit_surf = rhea_velocity_surface_new (ymir_mesh);
    rhea_inversion_obs_velocity_diff (
        misfit_surf, vel_fwd_vol, inv_problem->vel_obs_surf,
        inv_problem->vel_obs_weight_surf, inv_problem->vel_obs_type,
        rhea_stokes_problem_get_domain_options (stokes_problem));

    /* convert to physical dimensions */
    rhea_velocity_convert_to_dimensional_mm_yr (
        vel_fwd_surf, domain_options, temp_options);
    rhea_velocity_convert_to_dimensional_mm_yr (
        vel_adj_surf, domain_options, temp_options);
    rhea_velocity_convert_to_dimensional_mm_yr (
        vel_obs_surf, domain_options, temp_options);
    rhea_velocity_convert_to_dimensional_mm_yr (
        misfit_surf, domain_options, temp_options);

    /* write VTK of surface fields*/
    snprintf (path, BUFSIZ, "%s%s%02i", vtk_path_surf,
              RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
    rhea_vtk_write_inversion_iteration_surf (
        path, vel_fwd_surf, vel_adj_surf, vel_obs_surf,
        inv_problem->vel_obs_weight_surf, misfit_surf);

    /* destroy */
    rhea_velocity_surface_destroy (vel_fwd_surf);
    rhea_velocity_surface_destroy (vel_adj_surf);
    rhea_velocity_surface_destroy (vel_obs_surf);
    rhea_velocity_surface_destroy (misfit_surf);
  }

  /* process volume fields */
  if (vtk_path_vol != NULL) {
    ymir_vec_t         *weak_normal, *stress;

    /* compute normal direction at plate boundaries */
    weak_normal = rhea_weakzone_normal_new (ymir_mesh);
    rhea_weakzone_compute_normal (weak_normal, weak_options);

    /* compute stress tensor */
    stress = rhea_stress_new (ymir_mesh);
    rhea_stokes_problem_stress_compute (
        stress, inv_problem->forward_vel_press, stokes_problem,
        NULL /* !override_stress_op */, 0 /* !linearized */,
        0 /* !skip_pressure */);

    /* convert to physical dimensions */
    rhea_velocity_convert_to_dimensional_mm_yr (
        vel_fwd_vol, domain_options, temp_options);
    rhea_velocity_convert_to_dimensional_mm_yr (
        vel_adj_vol, domain_options, temp_options);
    rhea_pressure_convert_to_dimensional_Pa (
        press_fwd_vol, domain_options, temp_options, visc_options);
    rhea_pressure_convert_to_dimensional_Pa (
        press_adj_vol, domain_options, temp_options, visc_options);
    rhea_viscosity_convert_to_dimensional_Pas (viscosity, visc_options);
    rhea_stress_convert_to_dimensional_Pa (
        stress, domain_options, temp_options, visc_options);
    strainrate_dim_1_s =
      rhea_strainrate_get_dim_1_s (domain_options, temp_options) /
      rhea_velocity_get_dim_mm_yr (domain_options, temp_options);

    /* write VTK of volume fields*/
    snprintf (path, BUFSIZ, "%s%s%02i", vtk_path_vol,
              RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
    rhea_vtk_write_inversion_iteration (
        path, vel_fwd_vol, press_fwd_vol, vel_adj_vol, press_adj_vol,
        viscosity, marker, stress, weak_normal, strainrate_dim_1_s);

    /* destroy */
    rhea_stress_destroy (stress);
    rhea_weakzone_normal_destroy (weak_normal);
  }

  /* destroy */
  rhea_velocity_destroy (vel_fwd_vol);
  rhea_velocity_destroy (vel_adj_vol);
  rhea_pressure_destroy (press_fwd_vol);
  rhea_pressure_destroy (press_adj_vol);
  rhea_viscosity_destroy (viscosity);
  rhea_viscosity_destroy (marker);
}

static void
rhea_inversion_newton_output_prestep_fn (ymir_vec_t *solution,
                                         ymir_vec_t *neg_residual,
                                         const int iter, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  const int           iter_start = inv_problem->newton_options->iter_start;
  const int           amr_n_iter = rhea_inversion_allow_amr_for_outer_n_iter;
  const int           amr_active = (iter + 1) < amr_n_iter;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "newton_iter=%i", iter);

  /* print active inversion parameters */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFOF ("Inversion parameters, newton_iter=%i\n", iter);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  rhea_inversion_param_print (
      solution, RHEA_INVERSION_PARAM_VERBOSE_REAL_NONDIM_DIM,
      inv_problem->inv_param);
  RHEA_GLOBAL_INFO ("========================================\n");

  /* write text files */
  rhea_inversion_write_txt (solution, iter_start, iter, inv_problem);

  /* write visualization output */
  rhea_inversion_write_vis (iter, inv_problem);

  /* activate/deactivate AMR for next nonlinear iteration */
  RHEA_GLOBAL_INFOF_FN_TAG (__func__, "nonlinear_amr=%i", amr_active);
  if (amr_active) {
    rhea_stokes_problem_amr_process_options ();
  }
  else {
    rhea_stokes_problem_amr_set_nonlinear_type_name ("NONE");
  }

  RHEA_GLOBAL_VERBOSEF_FN_END (__func__, "newton_iter=%i", iter);
}

//static int
//rhea_inversion_newton_setup_poststep_fn (ymir_vec_t **solution, const int iter,
//                                         void *data)
//{
//}

static void
rhea_inversion_newton_problem_create (rhea_inversion_problem_t *inv_problem)
{
  rhea_inversion_param_t *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t            *ymir_mesh =
                            rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  sc_MPI_Comm             mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  rhea_newton_problem_t  *newton_problem;

  /* create vectors */
  inv_problem->newton_neg_gradient_vec =
    rhea_inversion_param_vec_new (inv_param);
  inv_problem->newton_gradient_adjoint_comp_vec =
    rhea_inversion_param_vec_new (inv_param);
  inv_problem->newton_gradient_prior_comp_vec =
    rhea_inversion_param_vec_new (inv_param);
  inv_problem->newton_gradient_observation_comp_vec =
    rhea_inversion_param_vec_new (inv_param);
  inv_problem->newton_step_vec =
    rhea_inversion_param_vec_new (inv_param);
  switch (rhea_inversion_hessian_type) {
  case RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND:
  case RHEA_INVERSION_HESSIAN_BFGS:
    inv_problem->newton_gradient_diff_vec =
      rhea_inversion_param_vec_new (inv_param);
    inv_problem->newton_gradient_diff_vec_update_count = 0;
    break;
  case RHEA_INVERSION_HESSIAN_FIRST_ORDER_APPROX:
  case RHEA_INVERSION_HESSIAN_FULL:
    inv_problem->newton_gradient_diff_vec = NULL;
    inv_problem->newton_gradient_diff_vec_update_count = 0;
    break;
  default: /* unknown Hessian type */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* create Newton problem */
  newton_problem = rhea_newton_problem_new (
      rhea_inversion_newton_compute_negative_gradient_fn,
      rhea_inversion_newton_compute_gradient_norm_fn,
      RHEA_INVERSION_GRAD_NORM_N_COMPONENTS,
      rhea_inversion_newton_solve_hessian_system_fn);
  inv_problem->newton_problem = newton_problem;
  inv_problem->newton_options = &rhea_inversion_newton_options;

  /* set vectors */
  rhea_newton_problem_set_vectors (
      inv_problem->newton_neg_gradient_vec,
      inv_problem->newton_step_vec, newton_problem);

  /* set functions */
  rhea_newton_problem_set_data_fn (
      inv_problem,
      rhea_inversion_newton_create_solver_data_fn,
      rhea_inversion_newton_clear_solver_data_fn, newton_problem);
  rhea_newton_problem_set_evaluate_objective_fn (
      rhea_inversion_newton_evaluate_objective_fn,
      RHEA_INVERSION_OBJ_N_COMPONENTS,
      rhea_inversion_newton_evaluate_objective_err_fn, newton_problem);
  rhea_newton_problem_set_apply_hessian_fn (
      rhea_inversion_newton_apply_hessian_fn, newton_problem);
  rhea_newton_problem_set_modify_step_fn (
      rhea_inversion_newton_modify_step_fn, newton_problem);
  rhea_newton_problem_set_update_fn (
      rhea_inversion_newton_update_operator_fn,
      rhea_inversion_newton_update_hessian_fn,
      NULL /*rhea_inversion_newton_modify_hessian_system_fn*/, newton_problem);
//rhea_newton_problem_set_setup_poststep_fn (
//    rhea_inversion_newton_setup_poststep_fn, newton_problem);
  rhea_newton_problem_set_output_fn (
      rhea_inversion_newton_output_prestep_fn, newton_problem);
  rhea_newton_problem_set_perfmon (
      rhea_inversion_newton_perfmon, newton_problem);
  rhea_newton_problem_set_mpicomm (
      mpicomm, newton_problem);

  /* set up derivative checks */
  rhea_newton_problem_set_checks (
      rhea_inversion_check_gradient, rhea_inversion_check_hessian,
      newton_problem);
  if ((rhea_inversion_check_gradient || rhea_inversion_check_gradient_of_qoi) &&
      rhea_inversion_check_gradient_elementwise) {
    const int          *active = rhea_inversion_param_get_active (inv_param);
    const int           n_active =
                          rhea_inversion_param_get_n_active (inv_param);
    ymir_vec_t        **perturb_vec;
    int                 vidx, i;

    inv_problem->check_gradient_perturb_vec = perturb_vec =
      RHEA_ALLOC (ymir_vec_t *, n_active);
    inv_problem->check_gradient_perturb_n_vecs = n_active;
    for (vidx = 0; vidx < inv_problem->check_gradient_perturb_n_vecs; vidx++) {
      perturb_vec[vidx] = rhea_inversion_param_vec_new (inv_param);
      ymir_vec_set_zero (perturb_vec[vidx]);
    }

    vidx = 0;
    for (i = 0; i < rhea_inversion_param_get_n_parameters (inv_param); i++) {
      if (active[i]) {
        RHEA_ASSERT (vidx < inv_problem->check_gradient_perturb_n_vecs);
        perturb_vec[vidx]->meshfree->e[0][i] = 1.0;
        vidx++;
      }
    }
    RHEA_ASSERT (vidx == inv_problem->check_gradient_perturb_n_vecs);
  }
  else {
    inv_problem->check_gradient_perturb_vec = NULL;
    inv_problem->check_gradient_perturb_n_vecs = 0;
  }
  rhea_newton_problem_check_gradient_set_perturbations (
      inv_problem->check_gradient_perturb_vec,
      inv_problem->check_gradient_perturb_n_vecs, newton_problem);
}

static void
rhea_inversion_newton_problem_clear (rhea_inversion_problem_t *inv_problem)
{
  /* destroy Newton problem */
  rhea_newton_problem_destroy (inv_problem->newton_problem);

  /* destroy vectors */
  rhea_inversion_param_vec_destroy (inv_problem->newton_neg_gradient_vec);
  rhea_inversion_param_vec_destroy (inv_problem->newton_step_vec);
  rhea_inversion_param_vec_destroy (
      inv_problem->newton_gradient_adjoint_comp_vec);
  rhea_inversion_param_vec_destroy (
      inv_problem->newton_gradient_prior_comp_vec);
  rhea_inversion_param_vec_destroy (
      inv_problem->newton_gradient_observation_comp_vec);
  if (NULL != inv_problem->newton_gradient_diff_vec) {
    rhea_inversion_param_vec_destroy (inv_problem->newton_gradient_diff_vec);
  }
  if (NULL != inv_problem->check_gradient_perturb_vec) {
    int                 vidx;

    RHEA_ASSERT (0 < inv_problem->check_gradient_perturb_n_vecs);
    for (vidx = 0; vidx < inv_problem->check_gradient_perturb_n_vecs; vidx++) {
      rhea_inversion_param_vec_destroy (
          inv_problem->check_gradient_perturb_vec[vidx]);
    }
    RHEA_FREE (inv_problem->check_gradient_perturb_vec);
  }
}

/******************************************************************************
 * Inversion Statistics
 *****************************************************************************/

/**
 * Creates a list to store statistics of inner solves.
 */
static void
rhea_inversion_inner_solve_stats_create (rhea_inversion_problem_t *inv_problem,
                                         const int max_num_iterations)
{
  int                 k;

  /* create text storage for forward and adjoint solves */
  inv_problem->fwd_adj_solve_stats = RHEA_ALLOC (char *, max_num_iterations);
  for (k = 0; k < max_num_iterations; k++) {
    inv_problem->fwd_adj_solve_stats[k] = NULL;
  }
  inv_problem->fwd_adj_solve_stats_length = max_num_iterations;
  inv_problem->fwd_adj_solve_stats_current_idx = 0;

  /* create text storage for incremental forward and adjoint solves */
  inv_problem->ifwd_iadj_solve_stats = RHEA_ALLOC (char *, max_num_iterations);
  for (k = 0; k < max_num_iterations; k++) {
    inv_problem->ifwd_iadj_solve_stats[k] = NULL;
  }
  inv_problem->ifwd_iadj_solve_stats_length = max_num_iterations;
  inv_problem->ifwd_iadj_solve_stats_current_idx = 0;

  /* init total values */
  inv_problem->fwd_adj_solve_stats_total_num_solves[0] = 0;
  inv_problem->fwd_adj_solve_stats_total_num_solves[1] = 0;
  inv_problem->fwd_adj_solve_stats_total_iterations[0] = 0;
  inv_problem->fwd_adj_solve_stats_total_iterations[1] = 0;
  inv_problem->fwd_adj_solve_stats_total_stop_reason[0] = 1;
  inv_problem->fwd_adj_solve_stats_total_stop_reason[1] = 1;
  inv_problem->ifwd_iadj_solve_stats_total_num_solves[0] = 0;
  inv_problem->ifwd_iadj_solve_stats_total_num_solves[1] = 0;
  inv_problem->ifwd_iadj_solve_stats_total_iterations[0] = 0;
  inv_problem->ifwd_iadj_solve_stats_total_iterations[1] = 0;
  inv_problem->ifwd_iadj_solve_stats_total_stop_reason[0] = 1;
  inv_problem->ifwd_iadj_solve_stats_total_stop_reason[1] = 1;
}

/**
 * Destroys a list with statistics of inner solves.
 */
static void
rhea_inversion_inner_solve_stats_clear (rhea_inversion_problem_t *inv_problem)
{
  int                 k;

  /* destroy text storage of forward and adjoint solves */
  for (k = 0; k < inv_problem->fwd_adj_solve_stats_length; k++) {
    if (NULL != inv_problem->fwd_adj_solve_stats[k]) {
      RHEA_FREE (inv_problem->fwd_adj_solve_stats[k]);
    }
  }
  RHEA_FREE (inv_problem->fwd_adj_solve_stats);

  /* destroy text storage of incremental forward and adjoint solves */
  for (k = 0; k < inv_problem->ifwd_iadj_solve_stats_length; k++) {
    if (NULL != inv_problem->ifwd_iadj_solve_stats[k]) {
      RHEA_FREE (inv_problem->ifwd_iadj_solve_stats[k]);
    }
  }
  RHEA_FREE (inv_problem->ifwd_iadj_solve_stats);
}

/**
 * Concatenates two strings.
 */
static char *
rhea_inversion_inner_solve_stats_cat (const char *a, const char *b)
{
  size_t              length;
  char               *cat;

  /* create string to store contatenation */
  length = 1 + (NULL != a ? strlen(a) : 0) + (NULL != b ? strlen(b) : 0);
  cat = RHEA_ALLOC (char, length);
  cat[0] = '\0';

  /* append strings */
  if (NULL != a) {
    strcat(cat, a);
  }
  if (NULL != b) {
    strcat(cat, b);
  }

  /* return concatenation */
  return cat;
}

/**
 * Print statistics of inner solves.
 */
static void
rhea_inversion_inner_solve_stats_print (rhea_inversion_problem_t *inv_problem)
{
  int                 stats_count, k;

  /* print statistics of forward and adjoint solves */
  stats_count = inv_problem->fwd_adj_solve_stats_current_idx + 1;
  RHEA_ASSERT (stats_count <= inv_problem->fwd_adj_solve_stats_length);
  if (0 < stats_count) {
    RHEA_GLOBAL_INFO ("========================================\n");
    RHEA_GLOBAL_INFO ("Inversion summary: forward & adjoint solves\n");
    RHEA_GLOBAL_INFO ("----------------------------------------\n");
    for (k = 0; k < stats_count; k++) {
      if (NULL == inv_problem->fwd_adj_solve_stats[k]) {
        break;
      }
      RHEA_GLOBAL_INFOF ("%3d; %s\n", k, inv_problem->fwd_adj_solve_stats[k]);
    }
    RHEA_GLOBAL_INFO ("----------------------------------------\n");
    RHEA_GLOBAL_INFOF ("fwd=(%d, %d, %d) adj=(%d, %d, %d)\n",
        inv_problem->fwd_adj_solve_stats_total_num_solves[0],
        inv_problem->fwd_adj_solve_stats_total_iterations[0],
        inv_problem->fwd_adj_solve_stats_total_stop_reason[0],
        inv_problem->fwd_adj_solve_stats_total_num_solves[1],
        inv_problem->fwd_adj_solve_stats_total_iterations[1],
        inv_problem->fwd_adj_solve_stats_total_stop_reason[1]);
    RHEA_GLOBAL_INFO ("========================================\n");
  }

  /* print statistics of incremental forward and adjoint solves */
  stats_count = inv_problem->ifwd_iadj_solve_stats_current_idx;
  RHEA_ASSERT (stats_count <= inv_problem->ifwd_iadj_solve_stats_length);
  if (0 < stats_count) {
    RHEA_GLOBAL_INFO ("========================================\n");
    RHEA_GLOBAL_INFO ("Inversion summary: "
                      "incremental forward & adjoint solves\n");
    RHEA_GLOBAL_INFO ("----------------------------------------\n");
    for (k = 0; k < stats_count; k++) {
      if (NULL == inv_problem->ifwd_iadj_solve_stats[k]) {
        break;
      }
      RHEA_GLOBAL_INFOF ("%3d; %s\n", k, inv_problem->ifwd_iadj_solve_stats[k]);
    }
    RHEA_GLOBAL_INFO ("----------------------------------------\n");
    RHEA_GLOBAL_INFOF ("incr_fwd=(%d, %d, %d) incr_adj=(%d, %d, %d)\n",
        inv_problem->ifwd_iadj_solve_stats_total_num_solves[0],
        inv_problem->ifwd_iadj_solve_stats_total_iterations[0],
        inv_problem->ifwd_iadj_solve_stats_total_stop_reason[0],
        inv_problem->ifwd_iadj_solve_stats_total_num_solves[1],
        inv_problem->ifwd_iadj_solve_stats_total_iterations[1],
        inv_problem->ifwd_iadj_solve_stats_total_stop_reason[1]);
    RHEA_GLOBAL_INFO ("========================================\n");
  }
}

/******************************************************************************
 * Inverse Problem
 *****************************************************************************/

rhea_inversion_problem_t *
rhea_inversion_new (rhea_stokes_problem_t *stokes_problem)
{
  rhea_inversion_problem_t *inv_problem;

  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (RHEA_INVERSION_KRYLOV_SOLVER_N <=
               rhea_stokes_problem_get_num_krylov_solvers (stokes_problem));

  /* initialize inverse problem */
  inv_problem = RHEA_ALLOC (rhea_inversion_problem_t, 1);
  inv_problem->stokes_problem = stokes_problem;
  inv_problem->forward_vel_press               = NULL;
  inv_problem->adjoint_vel_press               = NULL;
  inv_problem->incremental_forward_vel_press   = NULL;
  inv_problem->incremental_adjoint_vel_press   = NULL;
  inv_problem->forward_is_outdated             = 1;
  inv_problem->adjoint_is_outdated             = 1;
  inv_problem->forward_nonzero_init            = 0;
  inv_problem->adjoint_nonzero_init            = 0;
  inv_problem->incremental_forward_is_outdated = 1;
  inv_problem->incremental_adjoint_is_outdated = 1;
  inv_problem->rhs_vel_buffer                  = NULL;
  inv_problem->rhs_vel_press_buffer            = NULL;
  inv_problem->project_out_null =
    (rhea_inversion_project_out_null_t) rhea_inversion_project_out_null;
  inv_problem->hessian_matrix = NULL;

  /* initialize velocity data */
  inv_problem->vel_obs_type =
    (rhea_inversion_obs_velocity_t) rhea_inversion_vel_obs_type;
  inv_problem->vel_obs_weight_type =
    (rhea_inversion_obs_velocity_weight_t) rhea_inversion_vel_obs_weight_type;
  inv_problem->vel_obs_surf        = NULL;
  inv_problem->vel_obs_weight_surf = NULL;
  { /* get weight values from options */
    rhea_plate_options_t *plate_options =
      rhea_stokes_problem_get_plate_options (stokes_problem);
    const int           n_plates = rhea_plate_get_n_plates (plate_options);
    const double        vel_dim_mm_yr = rhea_velocity_get_dim_mm_yr (
        rhea_stokes_problem_get_domain_options (stokes_problem),
        rhea_stokes_problem_get_temperature_options (stokes_problem));
    int                 pid;

    inv_problem->vel_obs_weight_values =
      RHEA_ALLOC (double, SC_MAX (1, n_plates));
    if (NULL != rhea_inversion_vel_obs_stddev_mm_yr) {
      double              stddev, *stddev_vals = NULL;
      int                 n_entries = 0;

      n_entries = ymir_options_convert_string_to_double (
          rhea_inversion_vel_obs_stddev_mm_yr, &stddev_vals);
      if (0 < n_plates && n_plates == n_entries) {
        for (pid = 0; pid < n_plates; pid++) {
          stddev = stddev_vals[pid]/vel_dim_mm_yr;
          inv_problem->vel_obs_weight_values[pid] = 1.0/stddev;
        }
      }
      else if (0 < n_plates && 1 == n_entries) {
        stddev = stddev_vals[0]/vel_dim_mm_yr;
        for (pid = 0; pid < n_plates; pid++) {
          inv_problem->vel_obs_weight_values[pid] = 1.0/stddev;
        }
      }
      else if (1 == n_entries) {
        stddev = stddev_vals[0]/vel_dim_mm_yr;
        inv_problem->vel_obs_weight_values[0] = 1.0/stddev;
      }
      else {
        RHEA_ABORT_NOT_REACHED ();
      }
      YMIR_FREE (stddev_vals); /* was allocated in ymir */
    }
    else { /* if no input given */
      if (0 < n_plates) {
        for (pid = 0; pid < n_plates; pid++) {
          inv_problem->vel_obs_weight_values[pid] = 1.0;
        }
      }
      else {
        inv_problem->vel_obs_weight_values[0] = 1.0;
      }
    }
  }

  /* initialize viscosity data */
  inv_problem->visc_obs_type =
    (rhea_inversion_obs_viscosity_t) rhea_inversion_visc_obs_type;
  inv_problem->visc_obs_column = rhea_inversion_obs_viscosity_new (
      &inv_problem->visc_obs_n, &inv_problem->visc_obs_value,
      &inv_problem->visc_obs_weight, inv_problem->visc_obs_type,
      rhea_inversion_visc_obs_Pas, rhea_inversion_visc_obs_stddev_rel,
      rhea_stokes_problem_get_plate_options (stokes_problem),
      rhea_stokes_problem_get_viscosity_options (stokes_problem));

  /* initialize stress data */
  inv_problem->stress_obs_type =
    (rhea_inversion_obs_stress_t) rhea_inversion_stress_obs_type;
  inv_problem->stress_obs = NULL;               //TODO
  inv_problem->stress_obs_weight = NULL;        //TODO

  /* create parameters */
  inv_problem->inv_param_options = &rhea_inversion_param_options;
  inv_problem->inv_param = rhea_inversion_param_new (
      stokes_problem, inv_problem->inv_param_options);

  /* initialize weights */
  inv_problem->data_abs_weight  = NAN;
  inv_problem->prior_abs_weight = NAN;

  /* create quantities of interest */
  inv_problem->inv_qoi_options = &rhea_inversion_qoi_options;
  inv_problem->inv_qoi = rhea_inversion_qoi_new (
      stokes_problem, inv_problem->inv_qoi_options);

  /* create Newton problem */
  rhea_inversion_newton_problem_create (inv_problem);

  /* initialize adaptive relative tolerance for inner Stokes solves */
  if (rhea_inversion_inner_solver_rtol_adaptive) {
    inv_problem->inner_solver_adaptive_rtol_comp =
      inv_problem->newton_options->lin_rtol_init;
  }
  else {
    inv_problem->inner_solver_adaptive_rtol_comp = 1.0;
  }

  /* create statistics */
  {
    const int           max_n_iter = inv_problem->newton_options->iter_max -
                                     inv_problem->newton_options->iter_start+1;

    rhea_inversion_inner_solve_stats_create (inv_problem, max_n_iter);
  }

  /* initialize error statistics */
  inv_problem->error_stats_objective = NAN;
  inv_problem->error_stats_gradient  = NAN;
  inv_problem->error_stats_hessian   = NAN;

  /* initialize output paths */
  inv_problem->txt_path                = NULL;
  inv_problem->vtk_path_vol            = NULL;
  inv_problem->vtk_path_surf           = NULL;
  inv_problem->vtk_path_debug_gradient = NULL;

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);

  return inv_problem;
}

void
rhea_inversion_destroy (rhea_inversion_problem_t *inv_problem)
{
  RHEA_GLOBAL_PRODUCTION_FN_BEGIN (__func__);

  /* clear statistics */
  rhea_inversion_inner_solve_stats_clear (inv_problem);

  /* destroy solver data */
  if (rhea_inversion_solver_data_exists (inv_problem)) {
    rhea_inversion_newton_clear_solver_data_fn (NULL, inv_problem);
  }

  /* destroy Newton problem */
  rhea_inversion_newton_problem_clear (inv_problem);

  /* destroy quantities of interest */
  rhea_inversion_qoi_destroy (inv_problem->inv_qoi);

  /* destroy parameters */
  rhea_inversion_param_destroy (inv_problem->inv_param);

  /* destroy viscosity data */
  rhea_inversion_obs_viscosity_destroy (
      inv_problem->visc_obs_column, inv_problem->visc_obs_n,
      inv_problem->visc_obs_value, inv_problem->visc_obs_weight);

  /* destroy velocity data */
  if (NULL != inv_problem->vel_obs_weight_values) {
    RHEA_FREE (inv_problem->vel_obs_weight_values);
  }

  /* destroy structure */
  RHEA_FREE (inv_problem);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}

void
rhea_inversion_set_txt_output (rhea_inversion_problem_t *inv_problem,
                               char *txt_path)
{
  inv_problem->txt_path = txt_path;
}

void
rhea_inversion_set_vtk_output (rhea_inversion_problem_t *inv_problem,
                               char *vtk_path_vol, char *vtk_path_surf)
{
  inv_problem->vtk_path_vol = vtk_path_vol;
  inv_problem->vtk_path_surf = vtk_path_surf;
}

/******************************************************************************
 * Inverse Problem Solver
 *****************************************************************************/

void
rhea_inversion_solve (rhea_inversion_problem_t *inv_problem,
                      const int use_initial_guess,
                      ymir_vec_t *inv_parameter_vec)
{
  const int           init_guess = use_initial_guess &&
                                   (inv_parameter_vec != NULL);
  rhea_inversion_param_t *inv_param = inv_problem->inv_param;
  rhea_newton_options_t  *newton_options = inv_problem->newton_options;
  rhea_newton_problem_t  *newton_problem = inv_problem->newton_problem;
  ymir_vec_t         *parameter_vec;

  RHEA_GLOBAL_PRODUCTIONF_FN_BEGIN (__func__, "use_given_initial_guess=%i",
                                    init_guess);

  /* create vector for inversion parameters */
  parameter_vec = rhea_inversion_param_vec_new (inv_param);
  if (init_guess) {
    ymir_vec_copy (inv_parameter_vec, parameter_vec);
  }
  else {
    rhea_inversion_param_set_initial_guess (
        parameter_vec, inv_param, inv_problem->inv_param_options);
  }

  /* run Newton solver */
  newton_options->nonzero_initial_guess = 1;
  rhea_newton_solve (&parameter_vec, newton_problem, newton_options);

  /* print inner solver statistics */
  rhea_inversion_inner_solve_stats_print (inv_problem);

  /* destroy */
  if (inv_parameter_vec != NULL) {
    ymir_vec_copy (parameter_vec, inv_parameter_vec);
  }
  rhea_inversion_param_vec_destroy (parameter_vec);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}

void
rhea_inversion_solve_with_vel_obs (rhea_inversion_problem_t *inv_problem,
                                   const int use_initial_guess,
                                   ymir_vec_t *inv_parameter_vec,
                                   ymir_vec_t *vel_obs_surf,
                                   ymir_vec_t *vel_obs_weight_surf,
                                   const double vel_obs_add_noise_stddev)
{
  const int           init_guess = use_initial_guess &&
                                   (inv_parameter_vec != NULL);
  rhea_inversion_param_t *inv_param = inv_problem->inv_param;
  ymir_vec_t         *parameter_vec;

  RHEA_GLOBAL_PRODUCTIONF_FN_BEGIN (
      __func__, "use_given_initial_guess=%i, vel_obs_add_noise_stddev=%g",
      init_guess, vel_obs_add_noise_stddev);

  /* create vector for inversion parameters */
  parameter_vec = rhea_inversion_param_vec_new (inv_param);
  if (init_guess) {
    ymir_vec_copy (inv_parameter_vec, parameter_vec);
  }
  else {
    rhea_inversion_param_set_initial_guess (
        parameter_vec, inv_param, inv_problem->inv_param_options);
  }

  /* create solver data */
  rhea_inversion_newton_create_solver_data_fn (parameter_vec, inv_problem);

  /* check given velocity observations */
  RHEA_ASSERT (vel_obs_surf != NULL);
  RHEA_ASSERT (rhea_velocity_surface_check_vec_type (vel_obs_surf));
  RHEA_ASSERT (rhea_velocity_surface_is_valid (vel_obs_surf));
  RHEA_ASSERT (ymir_vec_get_mesh (vel_obs_surf) ==
               ymir_vec_get_mesh (inv_problem->vel_obs_surf));
  RHEA_ASSERT (vel_obs_weight_surf == NULL ||
               ymir_vec_get_mesh (vel_obs_weight_surf) ==
               ymir_vec_get_mesh (inv_problem->vel_obs_surf));

  /* override observational data */
  ymir_vec_copy (vel_obs_surf, inv_problem->vel_obs_surf);
  if (isfinite (vel_obs_add_noise_stddev) &&
      0.0 < fabs (vel_obs_add_noise_stddev)) {
    double              mean, variance, size, stddev;
    ymir_vec_t         *noise;

    /* calculate standard deviation for noise */
    ymir_vec_compute_mean_and_variance (inv_problem->vel_obs_surf, &mean,
                                        &variance, &size);
    stddev = sqrt (variance) * fabs (vel_obs_add_noise_stddev);
    RHEA_GLOBAL_VERBOSEF_FN_TAG (
        __func__, "vel_mean=%g, vel_variance=%g, noise_stddev=%g",
        mean, variance, stddev);

    /* add noise to observations; weighted by the l2-norm of the observations */
    //TODO make data noise independent of mesh size
    noise = ymir_vec_template (inv_problem->vel_obs_surf);
    ymir_vec_set_random_normal (noise, stddev, 0.0 /* mean */);
    ymir_vec_add (1.0, noise, inv_problem->vel_obs_surf);
    ymir_vec_destroy (noise);

    /* possibly remove normal/tangential components */
    switch (inv_problem->vel_obs_type) {
    case RHEA_INVERSION_OBS_VELOCITY_NORMAL:
      rhea_inversion_obs_velocity_remove_tangential (inv_problem->vel_obs_surf);
      break;
    case RHEA_INVERSION_OBS_VELOCITY_TANGENTIAL:
    case RHEA_INVERSION_OBS_VELOCITY_TANGENTIAL_ROTFREE:
      rhea_inversion_obs_velocity_remove_normal (inv_problem->vel_obs_surf);
      break;
    case RHEA_INVERSION_OBS_VELOCITY_ALL:
    case RHEA_INVERSION_OBS_VELOCITY_ALL_ROTFREE:
      break;
    default: /* unknown observation type */
      RHEA_ABORT_NOT_REACHED ();
    }
  }
  if (vel_obs_weight_surf != NULL) {
    RHEA_ASSERT (inv_problem->vel_obs_weight_surf != NULL);
    ymir_vec_copy (vel_obs_weight_surf, inv_problem->vel_obs_weight_surf);
  }
  RHEA_ASSERT (rhea_velocity_surface_is_valid (inv_problem->vel_obs_surf));

  /* run solver */
  rhea_inversion_solve (inv_problem, 1 /* use_initial_guess */, parameter_vec);

  /* destroy */
  if (inv_parameter_vec != NULL) {
    ymir_vec_copy (parameter_vec, inv_parameter_vec);
  }
  rhea_inversion_param_vec_destroy (parameter_vec);

  RHEA_GLOBAL_PRODUCTION_FN_END (__func__);
}
