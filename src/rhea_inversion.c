#include <rhea_inversion.h>
#include <rhea_base.h>
#include <rhea_inversion_param.h>
#include <rhea_inversion_obs_velocity.h>
#include <rhea_newton.h>
#include <rhea_velocity.h>
#include <rhea_velocity_pressure.h>
#include <rhea_stokes_problem_amr.h>
#include <rhea_io_std.h>
#include <rhea_vtk.h>
#include <ymir_perf_counter.h>

/* Inverse solver */
#define RHEA_INVERSION_GRAD_NORM_N_COMPONENTS 2

/* I/O labels that are attached at the end of file names */
#define RHEA_INVERSION_IO_LABEL_NL_ITER "_itn"

/* types of Hessians */
typedef enum
{
  RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND,
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
#define RHEA_INVERSION_DEFAULT_PARAMETER_PRIOR_REL_WEIGHT (NAN)
#define RHEA_INVERSION_DEFAULT_VEL_OBS_TYPE (RHEA_INVERSION_OBS_VELOCITY_NONE)
#define RHEA_INVERSION_DEFAULT_HESSIAN_TYPE (RHEA_INVERSION_HESSIAN_BFGS)
#define RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_MATRIX (1)
#define RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_ENFORCE_SYMM (0)
#define RHEA_INVERSION_DEFAULT_RESTRICT_STEP_TO_PRIOR_STDDEV (2.0)
#define RHEA_INVERSION_DEFAULT_INNER_SOLVER_RTOL_ADAPTIVE (0)
#define RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_FORWARD_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_ADJOINT_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_INCREMENTAL_FORWARD_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_ITER_MAX (100)
#define RHEA_INVERSION_DEFAULT_INCREMENTAL_ADJOINT_SOLVER_RTOL (1.0e-6)
#define RHEA_INVERSION_DEFAULT_PROJECT_OUT_NULL \
  (RHEA_INVERSION_PROJECT_OUT_NULL_NONE)
#define RHEA_INVERSION_DEFAULT_CHECK_GRADIENT (0)
#define RHEA_INVERSION_DEFAULT_CHECK_GRADIENT_ELEMENTWISE (0)
#define RHEA_INVERSION_DEFAULT_CHECK_HESSIAN (0)
#define RHEA_INVERSION_DEFAULT_MONITOR_PERFORMANCE (0)

#define RHEA_INVERSION_DEFAULT_PRIOR_PRECONDITIONER (0)

/* initialize options */
double              rhea_inversion_parameter_prior_rel_weight =
                      RHEA_INVERSION_DEFAULT_PARAMETER_PRIOR_REL_WEIGHT;
int                 rhea_inversion_vel_obs_type =
                      RHEA_INVERSION_DEFAULT_VEL_OBS_TYPE;
int                 rhea_inversion_hessian_type =
                      RHEA_INVERSION_DEFAULT_HESSIAN_TYPE;
int                 rhea_inversion_assemble_hessian_matrix =
                      RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_MATRIX;
int                 rhea_inversion_assemble_hessian_enforce_symm =
                      RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_ENFORCE_SYMM;
double              rhea_inversion_restrict_step_to_prior_stddev =
                      RHEA_INVERSION_DEFAULT_RESTRICT_STEP_TO_PRIOR_STDDEV;
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
int                 rhea_inversion_project_out_null =
                      RHEA_INVERSION_DEFAULT_PROJECT_OUT_NULL;
int                 rhea_inversion_check_gradient =
                      RHEA_INVERSION_DEFAULT_CHECK_GRADIENT;
int                 rhea_inversion_check_gradient_elementwise =
                      RHEA_INVERSION_DEFAULT_CHECK_GRADIENT_ELEMENTWISE;
int                 rhea_inversion_check_hessian =
                      RHEA_INVERSION_DEFAULT_CHECK_HESSIAN;
int                 rhea_inversion_monitor_performance =
                      RHEA_INVERSION_DEFAULT_MONITOR_PERFORMANCE;

rhea_inversion_param_options_t  rhea_inversion_param_options;
rhea_newton_options_t           rhea_inversion_newton_options;

void
rhea_inversion_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Inversion";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /* parameter options */
  YMIR_OPTIONS_D, "parameter-prior-weight", '\0',
    &(rhea_inversion_parameter_prior_rel_weight),
    RHEA_INVERSION_DEFAULT_PARAMETER_PRIOR_REL_WEIGHT,
    "Relative weight of the prior term in the objective functional",

  /* observational data options */
  YMIR_OPTIONS_I, "velocity-observations-type", '\0',
    &(rhea_inversion_vel_obs_type), RHEA_INVERSION_DEFAULT_VEL_OBS_TYPE,
    "Type of velocity observations at the surface",

  /* outer solver options (inversion) */
  YMIR_OPTIONS_I, "hessian-type", '\0',
    &(rhea_inversion_hessian_type), RHEA_INVERSION_DEFAULT_HESSIAN_TYPE,
    "Type of Hessian operator: 0: full, 1: Gauss-Newton, 2: BFGS",
  YMIR_OPTIONS_B, "assemble-hessian", '\0',
    &(rhea_inversion_assemble_hessian_matrix),
    RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_MATRIX,
    "Toggle assembly of the Hessian matrix",
  YMIR_OPTIONS_B, "assemble-hessian-enforce-symmetry", '\0',
    &(rhea_inversion_assemble_hessian_enforce_symm),
    RHEA_INVERSION_DEFAULT_ASSEMBLE_HESSIAN_ENFORCE_SYMM,
    "Enforce symmetry of the assembled Hessian matrix",
  YMIR_OPTIONS_D, "restrict-step-to-prior-stddev", '\0',
    &(rhea_inversion_restrict_step_to_prior_stddev),
    RHEA_INVERSION_DEFAULT_RESTRICT_STEP_TO_PRIOR_STDDEV,
    "Restrict each Newton step to lie within the range of prior",

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
  YMIR_OPTIONS_I, "project-out-null", '\0',
    &(rhea_inversion_project_out_null),
    RHEA_INVERSION_DEFAULT_PROJECT_OUT_NULL,
    "Inner solves: Project out null spaces",

  /* derivative checks */
  YMIR_OPTIONS_I, "check-gradient", '\0',
    &(rhea_inversion_check_gradient), RHEA_INVERSION_DEFAULT_CHECK_GRADIENT,
    "Check the gradient during Newton iterations",
  YMIR_OPTIONS_B, "check-gradient-elementwise", '\0',
    &(rhea_inversion_check_gradient_elementwise),
    RHEA_INVERSION_DEFAULT_CHECK_GRADIENT_ELEMENTWISE,
    "Check the gradient for each entry separately",
  YMIR_OPTIONS_I, "check-hessian", '\0',
    &(rhea_inversion_check_hessian), RHEA_INVERSION_DEFAULT_CHECK_HESSIAN,
    "Check the Hessian during Newton iterations",

  /* monitors */
  YMIR_OPTIONS_B, "monitor-performance", '\0',
    &(rhea_inversion_monitor_performance),
    RHEA_INVERSION_DEFAULT_MONITOR_PERFORMANCE,
    "Measure and print performance statistics (e.g., runtime or flops)",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  rhea_inversion_param_add_options (&rhea_inversion_param_options, opt);
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
  RHEA_INVERSION_PERFMON_NEWTON_SETUP_INV_SOLVER,
  RHEA_INVERSION_PERFMON_NEWTON_UPDATE_OPERATOR,
  RHEA_INVERSION_PERFMON_NEWTON_UPDATE_HESSIAN,
  RHEA_INVERSION_PERFMON_NEWTON_OBJ,
  RHEA_INVERSION_PERFMON_NEWTON_GRAD,
  RHEA_INVERSION_PERFMON_NEWTON_HESSIAN_APPLY,
  RHEA_INVERSION_PERFMON_NEWTON_HESSIAN_SOLVE,
  RHEA_INVERSION_PERFMON_SOLVE_FWD,
  RHEA_INVERSION_PERFMON_SOLVE_ADJ,
  RHEA_INVERSION_PERFMON_SOLVE_INCR_FWD,
  RHEA_INVERSION_PERFMON_SOLVE_INCR_ADJ,
  RHEA_INVERSION_PERFMON_N
}
rhea_inversion_perfmon_idx_t;

static const char  *rhea_inversion_perfmon_name[RHEA_INVERSION_PERFMON_N] =
{
  "Newton: setup solver",
  "Newton: update operator",
  "Newton: update Hessian",
  "Newton: evaluate objective",
  "Newton: compute gradient",
  "Newton: apply Hessian matrix",
  "Newton: solve Hessian system",
  "Solve forward problem",
  "Solve adjoint problem",
  "Solve incremental forward problem",
  "Solve incremental adjoint problem"
};
ymir_perf_counter_t rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_N];

void
rhea_inversion_perfmon_init (const int activate, const int skip_if_active)
{
  const int           active = activate && rhea_inversion_monitor_performance;

  ymir_perf_counter_init_all_ext (rhea_inversion_perfmon,
                                  rhea_inversion_perfmon_name,
                                  RHEA_INVERSION_PERFMON_N,
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
  int                 n_stats = RHEA_INVERSION_PERFMON_N *
                                YMIR_PERF_COUNTER_N_STATS;
  sc_statinfo_t       stats[n_stats];
  char                stats_name[n_stats][YMIR_PERF_COUNTER_NAME_SIZE];

  /* exit if nothing to do */
  if (!active || !print) {
    return;
  }

  /* gather performance statistics */
  n_stats = ymir_perf_counter_gather_stats (
      rhea_inversion_perfmon, RHEA_INVERSION_PERFMON_N, stats, stats_name,
      mpicomm, print_wtime, print_n_calls, print_flops);

  /* print performance statistics */
  ymir_perf_counter_print_stats (stats, n_stats, "Inversion");
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
  rhea_inversion_obs_velocity_t vel_obs_type;
  ymir_vec_t         *vel_obs_surf;
  ymir_vec_t         *vel_obs_weight_surf;

  /* weights for observational data and prior terms */
  double              data_abs_weight;
  double              prior_abs_weight;

  /* forward and adjoint states */
  ymir_vec_t         *forward_vel_press;
  ymir_vec_t         *adjoint_vel_press;
  int                 forward_is_outdated;
  int                 adjoint_is_outdated;
  int                 forward_nonzero_init;

  /* incremental forward and adjoint states */
  ymir_vec_t         *incremental_forward_vel_press;
  ymir_vec_t         *incremental_adjoint_vel_press;
  int                 incremental_forward_is_outdated;
  int                 incremental_adjoint_is_outdated;

  /* options for inner solves */
  double              inner_solver_adaptive_rtol_comp;
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
  ymir_vec_t             *newton_gradient_diff_vec;
  int                     newton_gradient_diff_vec_update_count;
  ymir_vec_t        **check_gradient_perturb_vec;
  int                 check_gradient_perturb_n_vecs;

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
};

static int
rhea_inversion_mesh_data_exists (rhea_inversion_problem_t *inv_problem)
{
  int                 fwd_adj_states_exist;
  int                 vel_obs_exist;

  fwd_adj_states_exist = (inv_problem->forward_vel_press != NULL &&
                          inv_problem->adjoint_vel_press != NULL &&
                          inv_problem->incremental_forward_vel_press != NULL &&
                          inv_problem->incremental_adjoint_vel_press != NULL);
  vel_obs_exist = (inv_problem->vel_obs_surf != NULL);

  return (fwd_adj_states_exist && vel_obs_exist);
}

static int
rhea_inversion_solver_data_exists (rhea_inversion_problem_t *inv_problem)
{
  return rhea_inversion_mesh_data_exists (inv_problem);
}

static void
rhea_inversion_set_weight_obs_misfit (rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;

  inv_problem->data_abs_weight = 1.0 / rhea_inversion_obs_velocity_misfit (
        NULL /* vel_fwd_vol */, inv_problem->vel_obs_surf,
        inv_problem->vel_obs_weight_surf, inv_problem->vel_obs_type,
        rhea_stokes_problem_get_domain_options (stokes_problem));
}

static void
rhea_inversion_set_weight_prior (rhea_inversion_problem_t *inv_problem)
{
  const double        prior_rel_weight =
                        rhea_inversion_parameter_prior_rel_weight;
  rhea_inversion_param_t *inv_param = inv_problem->inv_param;

  if (isfinite (prior_rel_weight) && 0.0 < prior_rel_weight) {
    inv_problem->prior_abs_weight =
      prior_rel_weight / rhea_inversion_param_prior (NULL /* !param */,
                                                     inv_param);
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
  inv_problem->incremental_forward_is_outdated = 1;
  inv_problem->incremental_adjoint_is_outdated = 1;

  /* create and generate observational data */
  RHEA_ASSERT (inv_problem->vel_obs_surf == NULL);
  RHEA_ASSERT (inv_problem->vel_obs_weight_surf == NULL);
  inv_problem->vel_obs_surf = rhea_velocity_surface_new (ymir_mesh);
  inv_problem->vel_obs_weight_surf = ymir_face_cvec_new (
      ymir_mesh, RHEA_DOMAIN_BOUNDARY_FACE_TOP, 1);
  rhea_inversion_obs_velocity_generate (
      inv_problem->vel_obs_surf, inv_problem->vel_obs_weight_surf,
      inv_problem->vel_obs_type,
      rhea_stokes_problem_get_plate_options (stokes_problem));

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
  inv_problem->adjoint_vel_press = NULL;
  inv_problem->incremental_forward_vel_press = NULL;
  inv_problem->incremental_adjoint_vel_press = NULL;

  /* destroy observational data */
  rhea_velocity_surface_destroy (inv_problem->vel_obs_surf);
  ymir_vec_destroy (inv_problem->vel_obs_weight_surf);
  inv_problem->vel_obs_surf = NULL;
  inv_problem->vel_obs_weight_surf = NULL;

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

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);
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
  inv_problem->forward_is_outdated = 0;
  //inv_problem->forward_nonzero_init = 1; //TODO

  /* write solver statistics */
  rhea_inversion_fwd_adj_solve_stats_write (inv_problem, 0 /* forward */,
                                            n_iter, res_reduc, stop_reason);
  inv_problem->error_stats_objective = res_reduc;

  /* recreate mesh dependent data */
  if (mesh_modified_by_solver) {
    rhea_inversion_newton_clear_mesh_data (inv_problem, 1 /* keep forward */);
    rhea_inversion_newton_create_mesh_data (inv_problem);

    /* deactivate AMR during future nonlinear solves */
    rhea_stokes_problem_amr_set_nonlinear_type_name ("NONE");
  }

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_SOLVE_FWD]);
  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/**
 * Runs solver for the adjoint state.
 */
static void
rhea_inversion_inner_solve_adjoint (rhea_inversion_problem_t *inv_problem)
{
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  const int           nonzero_init = 0; //TODO 1?
  const int           iter_max = rhea_inversion_adjoint_solver_iter_max;
  const double        rtol = rhea_inversion_adjoint_solver_rtol *
                             inv_problem->inner_solver_adaptive_rtol_comp;
  const int           resume = 0;
  int                 stop_reason, n_iter;
  double              res_reduc;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);
  ymir_perf_counter_start (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_SOLVE_ADJ]);

  //###DEV###
#if 0
  if (nonzero_init) {
    ymir_mesh_t          *ymir_mesh =
      rhea_stokes_problem_get_ymir_mesh (stokes_problem);
    ymir_pressure_elem_t *press_elem =
      rhea_stokes_problem_get_press_elem (stokes_problem);
    ymir_vec_t         *vel = rhea_velocity_new (ymir_mesh);

    ymir_vec_set_random (vel);
    ymir_vec_set_zero (inv_problem->adjoint_vel_press);
    rhea_stokes_problem_velocity_boundary_set_zero (vel, stokes_problem);
    rhea_velocity_pressure_set_components (inv_problem->adjoint_vel_press, vel,
                                           NULL, press_elem);
    rhea_velocity_destroy (vel);
  }
#endif

  /* run Stokes solver */
  stop_reason = rhea_stokes_problem_solve_ext (
      &(inv_problem->adjoint_vel_press), nonzero_init, iter_max, rtol,
      stokes_problem, resume, 1 /* force linear solve */,
      RHEA_INVERSION_KRYLOV_SOLVER_IDX_ADJOINT,
      &n_iter, &res_reduc, NULL /* mesh_modified_by_solver */);
  inv_problem->adjoint_is_outdated = 0;

  /* write solver statistics */
  rhea_inversion_fwd_adj_solve_stats_write (inv_problem, 1 /* adjoint */,
                                            n_iter, res_reduc, stop_reason);

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_SOLVE_ADJ]);
  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
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

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);
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
  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
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

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);
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
  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

/******************************************************************************
 * Callback Functions for Newton Solver
 *****************************************************************************/

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

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);
  ymir_perf_counter_start (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_SETUP_INV_SOLVER]);

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

  /* update inversion parameters */
  rhea_inversion_update_param (solution, inv_param);

  /* set up Stokes solver for forward problem */
  rhea_stokes_problem_setup_solver (stokes_problem);

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_SETUP_INV_SOLVER]);
  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void
rhea_inversion_newton_clear_solver_data_fn (void *data)
{
  rhea_inversion_problem_t *inv_problem = data;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  /* check input */
  RHEA_ASSERT (rhea_inversion_solver_data_exists (inv_problem));

  /* destroy Hessian matrix */
  if (inv_problem->hessian_matrix != NULL) {
    sc_dmatrix_destroy (inv_problem->hessian_matrix);
    inv_problem->hessian_matrix = NULL;
  }

  /* destroy mesh dependent data */
  rhea_inversion_newton_clear_mesh_data (inv_problem, 0 /* destroy forward */);

  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void
rhea_inversion_newton_update_operator_fn (ymir_vec_t *solution, void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);
  ymir_perf_counter_start (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_UPDATE_OPERATOR]);

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
#ifdef RHEA_ENABLE_DEBUG
  RHEA_GLOBAL_VERBOSE ("========================================\n");
  RHEA_GLOBAL_VERBOSE ("Inversion candidate parameters\n");
  RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
  rhea_inversion_param_print (
      solution, RHEA_INVERSION_PARAM_VERBOSE_REAL_NONDIM, inv_param);
  RHEA_GLOBAL_VERBOSE ("========================================\n");
#endif

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

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_UPDATE_OPERATOR]);
  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static void         rhea_inversion_apply_hessian (
                                        ymir_vec_t *param_vec_out,
                                        ymir_vec_t *param_vec_in,
                                        rhea_inversion_problem_t *inv_problem,
                                        const rhea_inversion_hessian_t type);

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
  double              scaling = NAN;
  const int           hessian_dim = inv_hessian_matrix->m;
  int                 k;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "damping=%g", damping);

  /* check input */
  RHEA_ASSERT (inv_hessian_matrix->m == inv_hessian_matrix->n);

  /* create work vectors */
  param_vec_in = rhea_inversion_param_vec_new (inv_param);
  param_vec_out = rhea_inversion_param_vec_new (inv_param);

  /* compute lumped Hessian to extract diagonal entries of its prior term */
  ymir_vec_set_value (param_vec_in, 1.0);
  rhea_inversion_apply_hessian (param_vec_out, param_vec_in, inv_problem,
                                RHEA_INVERSION_HESSIAN_BFGS);
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
    if (isfinite (scaling)) scaling *= hess_norm/grad_norm;
    else                    scaling = hess_norm/grad_norm;
  }

  /* apply damping to scaling factor */
  if (isfinite (damping) && 0.0 < damping) {
    if (isfinite (scaling)) scaling *= damping;
    else                    scaling = damping;
  }

  /* scale inverse Hessian matrix */
  if (isfinite (scaling) && 0.0 < scaling && scaling < 1.0) {
    RHEA_GLOBAL_VERBOSEF_FN_TAG (__func__, "scaling=%g", scaling);
    for (k = 0; k < hessian_dim; k++) {
      inv_hessian_matrix->e[k][k] *= scaling;
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

  RHEA_GLOBAL_VERBOSEF_FN_TAG (
      __func__, "curvature_grad=%g, curvature_hess=%g, curvature_effective=%g",
      curvature_grad, curvature_hess, curvature_effective);

  /* re-construct the initial BFGS approximation of inverse Hessian */
  if (reconstruct_init_matrix) {
    double              grad_diff_norm;

    grad_diff_norm = sqrt (rhea_inversion_param_vec_reduced_inner_product (
        gradient_diff_mod_vec, gradient_diff_mod_vec, NULL /* weight */,
        inv_param, 0 /* !normalization wrt. size */));
    rhea_inversion_assemble_hessian_bfgs_init (
        inv_hessian_matrix, gradient_diff_mod_vec,
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
  const rhea_inversion_hessian_t  type = rhea_inversion_hessian_type;
  const double        prior_rel_weight =
                        rhea_inversion_parameter_prior_rel_weight;
  rhea_inversion_problem_t *inv_problem = data;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "type=%i, assemble_matrix=%i", type,
      rhea_inversion_assemble_hessian_matrix);
  ymir_perf_counter_start (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_UPDATE_HESSIAN]);

  /* check input */
  RHEA_ASSERT (!rhea_inversion_assemble_hessian_matrix ||
               inv_problem->hessian_matrix != NULL);

  /* update Hessian */
  switch (type) {
  case RHEA_INVERSION_HESSIAN_GRADIENT_DESCEND:
    RHEA_ASSERT (rhea_inversion_assemble_hessian_matrix);
    RHEA_ASSERT (NULL != inv_problem->newton_gradient_diff_vec);
    if (1 == inv_problem->newton_gradient_diff_vec_update_count) {
      double              damping;

      /* set damping factor in case of rel. small prior term */
      if (isfinite (prior_rel_weight) &&
          0.0 < prior_rel_weight && prior_rel_weight < 1.0) {
        damping = prior_rel_weight;
      }
      else {
        damping = NAN;
      }

      /* set inverse Hessian matrix */
      rhea_inversion_assemble_hessian_gradient_descend (
          inv_problem->hessian_matrix, inv_problem->newton_gradient_diff_vec,
          damping, inv_problem);
    }
    else if (step_length < 0.5) { /* if adjust scaling of Hessian */
      sc_dmatrix_scale (step_length, inv_problem->hessian_matrix);
    }
    break;
  case RHEA_INVERSION_HESSIAN_BFGS:
    RHEA_ASSERT (rhea_inversion_assemble_hessian_matrix);
    RHEA_ASSERT (NULL != inv_problem->newton_gradient_diff_vec);
    if (1 == inv_problem->newton_gradient_diff_vec_update_count) {
      double              damping;

      /* set damping factor for Hessian in case of rel. small prior term */
      if (isfinite (prior_rel_weight) &&
          0.0 < prior_rel_weight && prior_rel_weight < 1.0) {
        damping = prior_rel_weight;
      }
      else {
        damping = NAN;
      }

      /* compute initial BFGS approximation of the inverse Hessian matrix */
      rhea_inversion_assemble_hessian_bfgs_init (
          inv_problem->hessian_matrix, inv_problem->newton_gradient_diff_vec,
          damping, inv_problem);
    }
    else { /* otherwise at 2nd nonlinear iteration and above */
      const int           reconstruct_init =
        (2 == inv_problem->newton_gradient_diff_vec_update_count);

      /* update BFGS approximation of the inverse Hessian matrix */
      RHEA_ASSERT (NULL != step_vec);
      rhea_inversion_assemble_hessian_bfgs_update (
          inv_problem->hessian_matrix, inv_problem->newton_gradient_diff_vec,
          step_vec, step_length, inv_problem, reconstruct_init);
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

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_UPDATE_HESSIAN]);
  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
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
  double              obj_data_misfit, obj_prior_misfit, obj_val;

  /* set weights if they do not exist */
  if (!isfinite (inv_problem->data_abs_weight)) {
    rhea_inversion_set_weight_obs_misfit (inv_problem);
  }
  if (!isfinite (inv_problem->prior_abs_weight)) {
    rhea_inversion_set_weight_prior (inv_problem);
  }
  data_abs_weight = inv_problem->data_abs_weight;
  prior_abs_weight = inv_problem->prior_abs_weight;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "weights=[%.6e,%.6e]",
                                 data_abs_weight, prior_abs_weight);
  ymir_perf_counter_start (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_OBJ]);

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
    rhea_stokes_problem_velocity_boundary_set_zero (vel, stokes_problem);
    RHEA_ASSERT (rhea_velocity_is_valid (vel));
  }

  /* compute the observation data misfit term */
  if (0.0 < data_abs_weight) {
    obj_data_misfit =
      0.5 * data_abs_weight * rhea_inversion_obs_velocity_misfit (
        vel, inv_problem->vel_obs_surf, inv_problem->vel_obs_weight_surf,
        inv_problem->vel_obs_type,
        rhea_stokes_problem_get_domain_options (stokes_problem));
  }
  else {
    obj_data_misfit = 0.0;
  }
  rhea_velocity_destroy (vel);

  /* compute prior term */
  if (0.0 < prior_abs_weight) {
    obj_prior_misfit =
      0.5 * prior_abs_weight * rhea_inversion_param_prior (solution, inv_param);
  }
  else {
    obj_prior_misfit = 0.0;
  }

  /* add terms to get full objective value */
  obj_val = obj_data_misfit + obj_prior_misfit;

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_OBJ]);
  RHEA_GLOBAL_VERBOSEF_FN_END (
      __func__, "value=%.15e, components[%.15e, %.15e]",
      obj_val, obj_data_misfit, obj_prior_misfit);

  /* return value of objective functional */
  if (NULL != obj_comp) {
    obj_comp[0] = obj_data_misfit;
    obj_comp[1] = obj_prior_misfit;
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
  ymir_vec_t         *gradient_adjoint_comp =
                        inv_problem->newton_gradient_adjoint_comp_vec;
  ymir_vec_t         *gradient_prior_comp =
                        inv_problem->newton_gradient_prior_comp_vec;
  ymir_vec_t         *vel, *rhs_vel_mass, *rhs_vel_press;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);
  ymir_perf_counter_start (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_GRAD]);

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
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* create work variables */
  vel = rhea_velocity_new (ymir_mesh);
  rhs_vel_mass = rhea_velocity_new (ymir_mesh);
  rhs_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);

  /* retrieve velocity of the forward state */
  rhea_velocity_pressure_copy_components (
      vel, NULL, inv_problem->forward_vel_press, press_elem);
  rhea_stokes_problem_velocity_boundary_set_zero (vel, stokes_problem);
  RHEA_ASSERT (rhea_velocity_is_valid (vel));

  /* set the right-hand side of the momentum eq. for the adjoint problem */
  rhea_inversion_obs_velocity_adjoint_rhs (
      rhs_vel_mass, vel, inv_problem->vel_obs_surf,
      inv_problem->vel_obs_weight_surf, inv_problem->vel_obs_type,
      rhea_stokes_problem_get_domain_options (stokes_problem));
  ymir_vec_scale (data_abs_weight, rhs_vel_mass);

  /* project out null spaces and enforce Dirichlet BC's on right-hand side */
  if (project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_RESIDUAL ||
      project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_BOTH) {
    rhea_stokes_problem_velocity_project_out_mean_rotation (
        rhs_vel_mass, 1 /* residual_space */, stokes_problem);
  }
  rhea_stokes_problem_velocity_boundary_set_zero (
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

  /* destroy */
  rhea_velocity_destroy (vel);
  rhea_velocity_destroy (rhs_vel_mass);
  rhea_velocity_pressure_destroy (rhs_vel_press);

  /* solve for adjoint state */
  rhea_inversion_inner_solve_adjoint (inv_problem);
  if (project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_PRIMAL ||
      project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_BOTH) {
    rhea_stokes_problem_project_out_nullspace (
        inv_problem->adjoint_vel_press, stokes_problem);
  }

  /* compute the (positive) gradient */
  rhea_inversion_param_compute_gradient (
      neg_gradient, solution, inv_problem->forward_vel_press,
      inv_problem->adjoint_vel_press, prior_abs_weight, inv_param,
      gradient_adjoint_comp, gradient_prior_comp);

  /* print gradient */
#ifdef RHEA_ENABLE_DEBUG
  {
    const int          *active = rhea_inversion_param_get_active (inv_param);
    const double       *grad = neg_gradient->meshfree->e[0];
    const double       *grad_adj = gradient_adjoint_comp->meshfree->e[0];
    const double       *grad_prior = gradient_prior_comp->meshfree->e[0];
    int                 i;

    RHEA_GLOBAL_VERBOSE ("========================================\n");
    RHEA_GLOBAL_VERBOSE ("Inversion gradient\n");
    RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
    for (i = 0; i < rhea_inversion_param_get_n_parameters (inv_param); i++) {
      if (active[i]) {
        RHEA_GLOBAL_VERBOSEF ("param# %3i: %+.6e [%+.6e, %+.6e]\n", i,
                              grad[i], grad_adj[i], grad_prior[i]);
      }
    }
    RHEA_GLOBAL_VERBOSE ("========================================\n");
  }
#endif

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
  if (rhea_inversion_check_gradient &&
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
        perturb_vec[vidx]->meshfree->e[0][i] =
            0.5 * fabs (solution->meshfree->e[0][i]);
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

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_GRAD]);
  RHEA_GLOBAL_VERBOSE_FN_END (__func__);
}

static double
rhea_inversion_newton_compute_gradient_norm_fn (ymir_vec_t *neg_gradient,
                                                void *data, double *norm_comp)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  double              grad_norm, grad_adj_norm, grad_prior_norm;

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
  RHEA_GLOBAL_VERBOSEF_FN_TAG (
      __func__, "gradient_norm=%.6e, gradient_adjoint_comp=%.6e, "
      "gradient_prior_comp=%.6e", grad_norm, grad_adj_norm, grad_prior_norm);

  /* return gradient norms */
  if (norm_comp != NULL) {
    norm_comp[0] = grad_adj_norm;
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
  const double        data_abs_weight = inv_problem->data_abs_weight;
  const double        prior_abs_weight = inv_problem->prior_abs_weight;
  int                 compute_incr_fwd, compute_incr_adj;
  ymir_vec_t         *vel, *rhs_vel_mass, *rhs_vel_press;

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

    /* create work variables */
    vel = rhea_velocity_new (ymir_mesh);
    rhs_vel_mass = rhea_velocity_new (ymir_mesh);
    rhs_vel_press = rhea_velocity_pressure_new (ymir_mesh, press_elem);
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
    rhea_stokes_problem_velocity_boundary_set_zero (
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
    /* retrieve velocity of the incremental forward state */
    rhea_velocity_pressure_copy_components (
        vel, NULL, inv_problem->incremental_forward_vel_press, press_elem);
    rhea_stokes_problem_velocity_boundary_set_zero (vel, stokes_problem);
    RHEA_ASSERT (rhea_velocity_is_valid (vel));

    /* set right-hand side of the momentum eq. for the incr. adjoint problem */
    rhea_inversion_obs_velocity_incremental_adjoint_rhs (
        rhs_vel_mass, vel, inv_problem->vel_obs_type,
        rhea_stokes_problem_get_domain_options (stokes_problem));
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
      rhea_stokes_problem_velocity_boundary_set_zero (vel, stokes_problem);
      RHEA_ASSERT (rhea_velocity_is_valid (vel));

      RHEA_ABORT_NOT_REACHED (); //TODO
      break;
    default: /* unknown Hessian type */
      RHEA_ABORT_NOT_REACHED ();
    }

    /* project out null spaces and enforce Dirichlet BC's on right-hand side */
    if (project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_RESIDUAL ||
        project_out_null == RHEA_INVERSION_PROJECT_OUT_NULL_BOTH) {
      rhea_stokes_problem_velocity_project_out_mean_rotation (
          rhs_vel_mass, 1 /* residual_space */, stokes_problem);
    }
    rhea_stokes_problem_velocity_boundary_set_zero (
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

  /* destroy after incremental solves */
  if (compute_incr_fwd || compute_incr_adj) {
    rhea_velocity_destroy (vel);
    rhea_velocity_destroy (rhs_vel_mass);
    rhea_velocity_pressure_destroy (rhs_vel_press);
  }

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
  const rhea_inversion_hessian_t  type = rhea_inversion_hessian_type;
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "type=%i, assemble_matrix=%i", type,
      rhea_inversion_assemble_hessian_matrix);
  ymir_perf_counter_start (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_HESSIAN_APPLY]);

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

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_HESSIAN_APPLY]);
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
  const rhea_inversion_hessian_t  type = rhea_inversion_hessian_type;
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  int                 stop_reason;

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (
      __func__, "lin_iter_max=%i, lin_rtol=%.1e, nonzero_init_guess=%i, "
      "type=%i, assemble_matrix=%i", lin_iter_max, lin_res_norm_rtol,
      nonzero_initial_guess, type, rhea_inversion_assemble_hessian_matrix);
  ymir_perf_counter_start (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_HESSIAN_SOLVE]);

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
  switch (rhea_inversion_hessian_type) {
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
#ifdef RHEA_ENABLE_DEBUG
  RHEA_GLOBAL_VERBOSE ("========================================\n");
  RHEA_GLOBAL_VERBOSE ("Inversion step\n");
  RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
  rhea_inversion_param_vec_print (step, inv_param);
  RHEA_GLOBAL_VERBOSE ("========================================\n");
#endif

  /* check output */
  RHEA_ASSERT (rhea_inversion_param_vec_is_valid (step, inv_param));

  /* increment index for statistics of inner solves */
  inv_problem->fwd_adj_solve_stats_current_idx++;
  inv_problem->ifwd_iadj_solve_stats_current_idx++;

  ymir_perf_counter_stop_add (
      &rhea_inversion_perfmon[RHEA_INVERSION_PERFMON_NEWTON_HESSIAN_SOLVE]);
  RHEA_GLOBAL_VERBOSE_FN_END (__func__);

  /* return iteraton count and "stopping" reason */
  return stop_reason;
}

static void
rhea_inversion_newton_modify_step_fn (ymir_vec_t *step, ymir_vec_t *solution,
                                      void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_inversion_param_t   *inv_param = inv_problem->inv_param;
  const double        restrict_to_prior_stddev =
                        rhea_inversion_restrict_step_to_prior_stddev;
  int                 step_modified;
  double              step_length_new;

  RHEA_GLOBAL_VERBOSE_FN_BEGIN (__func__);

  step_modified = rhea_inversion_param_restrict_step_length_to_feasible (
      step, solution, restrict_to_prior_stddev, inv_param, &step_length_new);
  if (step_modified) {
    RHEA_GLOBAL_INFOF_FN_TAG (__func__, "step_modified=%i, step_length_new=%g",
                              step_modified, step_length_new);

    /* print step */
#ifdef RHEA_ENABLE_DEBUG
    RHEA_GLOBAL_VERBOSE ("========================================\n");
    RHEA_GLOBAL_VERBOSEF ("Inversion feasible step, length reduced by %.15e\n",
                          step_length_new);
    RHEA_GLOBAL_VERBOSE ("----------------------------------------\n");
    rhea_inversion_param_vec_print (step, inv_param);
    RHEA_GLOBAL_VERBOSE ("========================================\n");
#endif
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
}

static void
rhea_inversion_newton_output_prestep_fn (ymir_vec_t *solution, const int iter,
                                         void *data)
{
  rhea_inversion_problem_t *inv_problem = data;
  rhea_stokes_problem_t    *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t              *ymir_mesh;
  ymir_pressure_elem_t     *press_elem;
  sc_MPI_Comm         mpicomm;
  int                 mpirank, mpiret;
  const int           iter_start = inv_problem->newton_options->iter_start;
  const char         *txt_path = inv_problem->txt_path;
  const char         *vtk_path_vol = inv_problem->vtk_path_vol;
  const char         *vtk_path_surf = inv_problem->vtk_path_surf;
  char                path[BUFSIZ];

  RHEA_GLOBAL_VERBOSEF_FN_BEGIN (__func__, "newton_iter=%i", iter);

  /* print active inversion parameters */
  RHEA_GLOBAL_INFO ("========================================\n");
  RHEA_GLOBAL_INFOF ("Inversion parameters, newton_iter=%i\n", iter);
  RHEA_GLOBAL_INFO ("----------------------------------------\n");
  rhea_inversion_param_print (
      solution, RHEA_INVERSION_PARAM_VERBOSE_REAL_NONDIM,
      inv_problem->inv_param);
  RHEA_GLOBAL_INFO ("========================================\n");

  /* get mesh data */
  ymir_mesh = rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  press_elem = rhea_stokes_problem_get_press_elem (stokes_problem);

  /* get parallel environment */
  mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

  /* write text files */
  if (txt_path != NULL && mpirank == 0) {
    rhea_inversion_param_t *inv_param = inv_problem->inv_param;
    sc_dmatrix_t       *sol, *step;
    sc_dmatrix_t       *grad_combined, *neg_grad, *grad_adj, *grad_prior;
    int                 i;

    /* write parameters */
    sol = rhea_inversion_param_vec_reduced_new (solution, inv_param);
    snprintf (path, BUFSIZ, "%s%s%02i_parameters.txt", txt_path,
              RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
    rhea_io_std_write_double_to_txt (
        path, sol->e[0], (size_t) sol->m * sol->n, sol->n);
    rhea_inversion_param_vec_reduced_destroy (sol);

    /* write gradient */
    neg_grad = rhea_inversion_param_vec_reduced_new (
        inv_problem->newton_neg_gradient_vec, inv_param);
    grad_adj = rhea_inversion_param_vec_reduced_new (
        inv_problem->newton_gradient_adjoint_comp_vec, inv_param);
    grad_prior = rhea_inversion_param_vec_reduced_new (
        inv_problem->newton_gradient_prior_comp_vec, inv_param);
    grad_combined = sc_dmatrix_new (neg_grad->m, 3);
    RHEA_ASSERT (neg_grad->m == grad_adj->m && neg_grad->m == grad_prior->m);
    for (i = 0; i < grad_combined->m; i++) {
      grad_combined->e[i][0] = -neg_grad->e[i][0];
      grad_combined->e[i][1] = grad_adj->e[i][0];
      grad_combined->e[i][2] = grad_prior->e[i][0];
    }
    snprintf (path, BUFSIZ, "%s%s%02i_gradient.txt", txt_path,
              RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
    rhea_io_std_write_double_to_txt (
        path, grad_combined->e[0], (size_t) grad_combined->m * grad_combined->n,
        grad_combined->n);
    rhea_inversion_param_vec_reduced_destroy (neg_grad);
    rhea_inversion_param_vec_reduced_destroy (grad_adj);
    rhea_inversion_param_vec_reduced_destroy (grad_prior);
    sc_dmatrix_destroy (grad_combined);

    /* write step */
    if (iter_start < iter) {
      step = rhea_inversion_param_vec_reduced_new (
          inv_problem->newton_step_vec, inv_param);
      snprintf (path, BUFSIZ, "%s%s%02i_step.txt", txt_path,
                RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
      rhea_io_std_write_double_to_txt (
          path, step->e[0], (size_t) step->m * step->n, step->n);
      rhea_inversion_param_vec_reduced_destroy (step);
    }

    /* write Hessian */
    if (iter_start < iter && rhea_inversion_assemble_hessian_matrix) {
      snprintf (path, BUFSIZ, "%s%s%02i_hessian.txt", txt_path,
                RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
      rhea_io_std_write_double_to_txt (
          path, inv_problem->hessian_matrix->e[0],
          (size_t) inv_problem->hessian_matrix->m *
                   inv_problem->hessian_matrix->n,
          inv_problem->hessian_matrix->n);
    }
  }

  /* create visualization */
  if (vtk_path_vol != NULL || vtk_path_surf != NULL) {
    ymir_vec_t         *vel_fwd_vol, *press_fwd_vol;
    ymir_vec_t         *vel_adj_vol, *press_adj_vol;

    /* get volume fields */
    rhea_velocity_pressure_create_components (
        &vel_fwd_vol, &press_fwd_vol, inv_problem->forward_vel_press,
        press_elem);
    rhea_velocity_pressure_create_components (
        &vel_adj_vol, &press_adj_vol, inv_problem->adjoint_vel_press,
        press_elem);

    /* write VTK of volume fields */
    if (vtk_path_vol != NULL) {
      snprintf (path, BUFSIZ, "%s%s%02i", vtk_path_vol,
                RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
      rhea_vtk_write_inversion_iteration (
          path, vel_fwd_vol, press_fwd_vol, vel_adj_vol, press_adj_vol);
    }

    /* destroy */
    rhea_velocity_destroy (press_fwd_vol);
    rhea_velocity_destroy (press_adj_vol);

    if (vtk_path_surf != NULL) {
      ymir_vec_t         *vel_fwd_surf, *vel_adj_surf, *misfit_surf;

      /* get surface fields */
      vel_fwd_surf = rhea_velocity_surface_new (ymir_mesh);
      vel_adj_surf = rhea_velocity_surface_new (ymir_mesh);
      rhea_velocity_surface_interpolate (vel_fwd_surf, vel_fwd_vol);
      rhea_velocity_surface_interpolate (vel_adj_surf, vel_adj_vol);

      /* compute the observation data misfit term */
      misfit_surf = rhea_velocity_surface_new (ymir_mesh);
      rhea_inversion_obs_velocity_misfit_vec (
          misfit_surf, vel_fwd_vol, inv_problem->vel_obs_surf,
          inv_problem->vel_obs_weight_surf, inv_problem->vel_obs_type,
          rhea_stokes_problem_get_domain_options (stokes_problem));

      /* write VTK of surface fields*/
      snprintf (path, BUFSIZ, "%s_surf%s%02i", vtk_path_surf,
                RHEA_INVERSION_IO_LABEL_NL_ITER, iter);
      rhea_vtk_write_inversion_iteration_surf (
          path, vel_fwd_surf, vel_adj_surf, inv_problem->vel_obs_surf,
          misfit_surf);

      /* destroy */
      rhea_velocity_surface_destroy (vel_fwd_surf);
      rhea_velocity_surface_destroy (vel_adj_surf);
      rhea_velocity_surface_destroy (misfit_surf);
    }

    /* destroy */
    rhea_velocity_destroy (vel_fwd_vol);
    rhea_velocity_destroy (vel_adj_vol);
  }

  RHEA_GLOBAL_VERBOSEF_FN_END (__func__, "newton_iter=%i", iter);
}

static void
rhea_inversion_newton_problem_create (rhea_inversion_problem_t *inv_problem)
{
  rhea_inversion_param_t *inv_param = inv_problem->inv_param;
  rhea_stokes_problem_t  *stokes_problem = inv_problem->stokes_problem;
  ymir_mesh_t        *ymir_mesh =
                        rhea_stokes_problem_get_ymir_mesh (stokes_problem);
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  const int           obj_n_components = 2;
  rhea_newton_problem_t *newton_problem;


  /* create vectors */
  inv_problem->newton_neg_gradient_vec =
    rhea_inversion_param_vec_new (inv_param);
  inv_problem->newton_gradient_adjoint_comp_vec =
    rhea_inversion_param_vec_new (inv_param);
  inv_problem->newton_gradient_prior_comp_vec =
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
      rhea_inversion_newton_evaluate_objective_fn, obj_n_components,
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
//    rhea_inversion_newton_amr_fn, newton_problem);
  rhea_newton_problem_set_output_fn (
      rhea_inversion_newton_output_prestep_fn, newton_problem);
  rhea_newton_problem_set_mpicomm (
      mpicomm, newton_problem);

  /* set up derivative checks */
  rhea_newton_problem_set_checks (
      rhea_inversion_check_gradient, rhea_inversion_check_hessian,
      newton_problem);
  if (rhea_inversion_check_gradient &&
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
      RHEA_ASSERT (NULL != inv_problem->fwd_adj_solve_stats[k]);
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
  inv_problem->vel_obs_type =
    (rhea_inversion_obs_velocity_t) rhea_inversion_vel_obs_type;
  inv_problem->vel_obs_surf = NULL;
  inv_problem->vel_obs_weight_surf = NULL;
  inv_problem->forward_vel_press = NULL;
  inv_problem->adjoint_vel_press = NULL;
  inv_problem->incremental_forward_vel_press = NULL;
  inv_problem->incremental_adjoint_vel_press = NULL;
  inv_problem->forward_is_outdated = 1;
  inv_problem->adjoint_is_outdated = 1;
  inv_problem->forward_nonzero_init = 0;
  inv_problem->incremental_forward_is_outdated = 1;
  inv_problem->incremental_adjoint_is_outdated = 1;
  inv_problem->project_out_null =
      (rhea_inversion_project_out_null_t) rhea_inversion_project_out_null;
  inv_problem->hessian_matrix = NULL;

  /* create parameters */
  inv_problem->inv_param_options = &rhea_inversion_param_options;
  inv_problem->inv_param = rhea_inversion_param_new (
      stokes_problem, inv_problem->inv_param_options);

  /* initialize weights */
  inv_problem->data_abs_weight = NAN;
  inv_problem->prior_abs_weight = NAN;

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
  inv_problem->error_stats_gradient = NAN;
  inv_problem->error_stats_hessian = NAN;

  /* initialize output paths */
  inv_problem->txt_path = NULL;
  inv_problem->vtk_path_vol = NULL;
  inv_problem->vtk_path_surf = NULL;

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
    rhea_inversion_newton_clear_solver_data_fn (inv_problem);
  }

  /* destroy Newton problem */
  rhea_inversion_newton_problem_clear (inv_problem);

  /* destroy parameters */
  rhea_inversion_param_destroy (inv_problem->inv_param);

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
