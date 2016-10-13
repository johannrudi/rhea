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

  This example runs perfomace tests of matvecs.

*/

#include <slabs_setup.h>
#include <ymir.h>
#include <ymir_comm.h>
#include <ymir_stress_op_optimized.h>
#include <ymir_stokes_op_optimized.h>
#include <ymir_bbt_optimized.h>
#include <ymir_stiff_op_optimized.h>
#include <ymir_stokes_vec.h>
#include <ymir_perf_counter.h>

#define SLABS_PERF_TENSOR_REF (0)

#include <mangll_tensor_optimized.h>
#include <mangll_tensor_test.h>
#if (1 == SLABS_PERF_TENSOR_REF)
# include <mangll_tensor.h>
#endif

#if defined(__bgq__)
#include <ymir_bgq.h>
#endif

/* for performance statistics */
//#include <ymir_gmg_hierarchy_mesh.h>
//#include <ymir_gmg_hierarchy_stress.h>
//#include <ymir_gmg_hierarchy_stiff.h>
//#include <ymir_gmg_hierarchy_bbt.h>

/* perfomance counters */
typedef enum
{
  SLABS_PERF_COUNTER_SETUP_MESH,
  SLABS_PERF_COUNTER_SETUP_STOKES,

  SLABS_PERF_COUNTER_TENSOR_OPT,
  SLABS_PERF_COUNTER_TENSOR_REF,
  SLABS_PERF_COUNTER_TENSOR_GRAD_VEL_OPT,

  SLABS_PERF_COUNTER_STRESS_MATVEC,
  SLABS_PERF_COUNTER_STOKES_MATVEC,
  SLABS_PERF_COUNTER_BBT_MATVEC,
  SLABS_PERF_COUNTER_STIFF_MATVEC,
  SLABS_PERF_COUNTER_TOTAL,
  SLABS_PERF_COUNTER_N
}
slabs_perf_counter_idx_t;
ymir_perf_counter_t slabs_perf_counter[SLABS_PERF_COUNTER_N];
const char         *slabs_perf_counter_name[SLABS_PERF_COUNTER_N] =
{
  "Setup mesh",
  "Setup Stokes",

  "Tensor all optimized",
  "Tensor all reference",
  "Tensor grad, velocity, optimized",

  "Stress matvec",
  "Stokes matvec",
  "BB^T matvec",
  "Stiffness matvec",
  "Total"
};
sc_statinfo_t       slabs_perf_stats[
                      SLABS_PERF_COUNTER_N * YMIR_PERF_COUNTER_N_STATS];
char                slabs_perf_stats_name[
                      SLABS_PERF_COUNTER_N * YMIR_PERF_COUNTER_N_STATS][
                      YMIR_PERF_COUNTER_NAME_SIZE];

/* set default options */
#define SLABS_PERF_MATVEC_DEFAULT_TENSOR_APPLY_NUM (0)
#define SLABS_PERF_MATVEC_DEFAULT_STRESS_APPLY_NUM (0)
#define SLABS_PERF_MATVEC_DEFAULT_STRESS_APPLY_NL (0)
#define SLABS_PERF_MATVEC_DEFAULT_STRESS_APPLY_DIRTY (0)
#define SLABS_PERF_MATVEC_DEFAULT_STRESS_APPLY_OPTIMIZED (0)
#define SLABS_PERF_MATVEC_DEFAULT_STOKES_APPLY_NUM (0)
#define SLABS_PERF_MATVEC_DEFAULT_STOKES_APPLY_NL (0)
#define SLABS_PERF_MATVEC_DEFAULT_STOKES_APPLY_DIRTY (0)
#define SLABS_PERF_MATVEC_DEFAULT_STOKES_APPLY_OPTIMIZED (0)
#define SLABS_PERF_MATVEC_DEFAULT_BBT_APPLY_NUM (0)
#define SLABS_PERF_MATVEC_DEFAULT_BBT_APPLY_OPTIMIZED (0)
#define SLABS_PERF_MATVEC_DEFAULT_STIFF_APPLY_NUM (0)
#define SLABS_PERF_MATVEC_DEFAULT_STIFF_APPLY_DIRTY (0)
#define SLABS_PERF_MATVEC_DEFAULT_STIFF_APPLY_OPTIMIZED (0)

/**
 * Sets random input vector.
 */
static void
slabs_perf_matvec_set_random_vec (ymir_vec_t *vec)
{
#ifdef YMIR_PETSC
  slabs_petsc_vec_set_random (vec, YMIR_MESH_PETSCLAYOUT_NONE);
#else
  YMIR_ABORT_NOT_REACHED ();
#endif
}

#if (1 == SLABS_PERF_TENSOR_REF)
/**
 * Tests apply functions for stress operator.
 */
static double
slabs_perf_tensor_test_get_error (const double * chk, const double * ref,
                                  const int n_entries)
{
  double              rel_error_max = 0.0;
  int                 k;

  for (k = 0; k < n_entries; k++) {
    rel_error_max = SC_MAX ( rel_error_max, fabs ((chk[k] - ref[k]) / ref[k]) );
  }

  return rel_error_max;
}
#endif

/**
 * Tests apply functions for stress operator.
 */
static void
slabs_perf_matvec_test_stress (ymir_linearop_fn fn_chk,
                               ymir_linearop_fn fn_ref,
                               ymir_stress_op_t *stress_op)
{
  const char         *this_fn_name = "slabs_perf_matvec_test_stress";
  ymir_mesh_t        *mesh = stress_op->viscosity->mesh;
  ymir_vec_t         *in = ymir_cvec_new (mesh, 3);
  ymir_vec_t         *out_ref = ymir_cvec_new (mesh, 3);
  ymir_vec_t         *out_chk = ymir_cvec_new (mesh, 3);
  double              abs_error, rel_error;

  /* set random input vector */
  slabs_perf_matvec_set_random_vec (in);
  ymir_vel_dir_separate (in, in, NULL, NULL, stress_op->vel_dir);

  /* apply function that's being tested */
  fn_chk (in, out_chk, stress_op);

  /* apply reference function */
  fn_ref (in, out_ref, stress_op);

  /* calculate error in output vectors */
  ymir_vec_add (-1.0, out_ref, out_chk);
  abs_error = ymir_vec_norm (out_chk);
  rel_error = abs_error / ymir_vec_norm (out_ref);

  YMIR_GLOBAL_INFOF ("%s: abs error %1.3e rel error %1.3e\n",
                     this_fn_name, abs_error, rel_error);

  /* destroy */
  ymir_vec_destroy (in);
  ymir_vec_destroy (out_ref);
  ymir_vec_destroy (out_chk);
}

/**
 * Tests apply functions for Stokes operator.
 */
static void
slabs_perf_matvec_test_stokes (ymir_linearop_fn fn_chk,
                               ymir_linearop_fn fn_ref,
                               ymir_stokes_op_t *stokes_op)
{
  const char         *this_fn_name = "slabs_perf_matvec_test_stokes";
  ymir_mesh_t        *mesh = stokes_op->stress_op->viscosity->mesh;
  ymir_pressure_elem_t  *press_elem = stokes_op->press_elem;
  ymir_vec_t         *in = ymir_stokes_vec_new (mesh, press_elem);
  ymir_vec_t         *out_ref = ymir_stokes_vec_new (mesh, press_elem);
  ymir_vec_t         *out_chk = ymir_stokes_vec_new (mesh, press_elem);
  double              abs_error, rel_error;

  /* set random input vector */
  slabs_perf_matvec_set_random_vec (in);
  ymir_vel_dir_separate (in, in, NULL, NULL, stokes_op->stress_op->vel_dir);

  /* apply function that's being tested */
  fn_chk (in, out_chk, stokes_op);

  /* apply reference function */
  fn_ref (in, out_ref, stokes_op);

  /* calculate error in output vectors */
  ymir_vec_add (-1.0, out_ref, out_chk);
  abs_error = ymir_vec_norm (out_chk);
  rel_error = abs_error / ymir_vec_norm (out_ref);

  YMIR_GLOBAL_INFOF ("%s: abs error %1.3e rel error %1.3e\n",
                     this_fn_name, abs_error, rel_error);

  /* destroy */
  ymir_vec_destroy (in);
  ymir_vec_destroy (out_ref);
  ymir_vec_destroy (out_chk);
}

/**
 * Tests apply functions for BB^T operator.
 */
static void
slabs_perf_matvec_test_bbt (ymir_linearop_fn fn_chk,
                            ymir_linearop_fn fn_ref,
                            ymir_bbt_t *bbt)
{
  const char         *this_fn_name = "slabs_perf_matvec_test_bbt";
  ymir_mesh_t        *mesh = bbt->stress_op->viscosity->mesh;
  ymir_pressure_elem_t  *press_elem = bbt->press_elem;
  ymir_vec_t         *in = ymir_pressure_vec_new (mesh, press_elem);
  ymir_vec_t         *out_ref = ymir_pressure_vec_new (mesh, press_elem);
  ymir_vec_t         *out_chk = ymir_pressure_vec_new (mesh, press_elem);
  double              abs_error, rel_error;

  /* set random input vector */
  slabs_perf_matvec_set_random_vec (in);

  /* apply function that's being tested */
  fn_chk (in, out_chk, bbt);

  /* apply reference function */
  fn_ref (in, out_ref, bbt);

  /* calculate error in output vectors */
  ymir_vec_add (-1.0, out_ref, out_chk);
  abs_error = ymir_vec_norm (out_chk);
  rel_error = abs_error / ymir_vec_norm (out_ref);

  YMIR_GLOBAL_INFOF ("%s: abs error %1.3e rel error %1.3e\n",
                     this_fn_name, abs_error, rel_error);

  /* destroy */
  ymir_vec_destroy (in);
  ymir_vec_destroy (out_ref);
  ymir_vec_destroy (out_chk);
}

/**
 * Tests apply functions for stiffness operator.
 */
static void
slabs_perf_matvec_test_stiff (ymir_linearop_fn fn_chk,
                              ymir_linearop_fn fn_ref,
                              ymir_stiff_op_t *stiff_op)
{
  const char         *this_fn_name = "slabs_perf_matvec_test_stiff";
  ymir_mesh_t        *mesh = stiff_op->invec->mesh;
  ymir_vec_t         *in = ymir_cvec_new (mesh, 1);
  ymir_vec_t         *out_ref = ymir_cvec_new (mesh, 1);
  ymir_vec_t         *out_chk = ymir_cvec_new (mesh, 1);
  double              abs_error, rel_error;

  /* set random input vector */
  slabs_perf_matvec_set_random_vec (in);
  if (stiff_op->dirfaces != NULL) {
    ymir_stiff_op_zero_dirichlet_simple (in, in, NULL, stiff_op->dirscale);
  }

  /* apply function that's being tested */
  fn_chk (in, out_chk, stiff_op);

  /* apply reference function */
  fn_ref (in, out_ref, stiff_op);

  /* calculate error in output vectors */
  ymir_vec_add (-1.0, out_ref, out_chk);
  abs_error = ymir_vec_norm (out_chk);
  rel_error = abs_error / ymir_vec_norm (out_ref);

  YMIR_GLOBAL_INFOF ("%s: abs error %1.3e rel error %1.3e\n",
                     this_fn_name, abs_error, rel_error);

  /* destroy */
  ymir_vec_destroy (in);
  ymir_vec_destroy (out_ref);
  ymir_vec_destroy (out_chk);
}

/******************************************************************************
 * Main program
 *****************************************************************************/

int
main (int argc, char **argv)
{
  const char         *this_fn_name = "slabs_perf_matvec:main";
  /* MPI */
  MPI_Comm            mpicomm = MPI_COMM_WORLD;
  int                 mpisize, mpirank, ompsize;
  int                 mpiret;
  /* options */
  ymir_options_t     *opt;
  slabs_physics_options_t   physics_options;
  slabs_discr_options_t     discr_options;
  slabs_nl_solver_options_t solver_options;
  /* options local to this function */
  int                 tensor_apply_num;
  int                 stress_apply_num;
  int                 stress_apply_nl;
  int                 stress_apply_dirty;
  int                 stress_apply_optimized;
  int                 stokes_apply_num;
  int                 stokes_apply_nl;
  int                 stokes_apply_dirty;
  int                 stokes_apply_optimized;
  int                 bbt_apply_num;
  int                 bbt_apply_optimized;
  int                 stiff_apply_num;
  int                 stiff_apply_dirty;
  int                 stiff_apply_optimized;
  int                 production_run;
  int                 monitor_performance;
  char               *workload_out_path;
  /* mesh */
  p8est_t            *p8est;
  ymir_mesh_t        *mesh;
  ymir_pressure_elem_t  *press_elem;
  /* Stokes problem */
  slabs_stokes_state_t  *state;
  slabs_discr_enforce_refinement_data_t     *enforce_refinement_data;
  slabs_physics_coarsen_stokes_coeff_data_t *coarsen_coeff_data;
  slabs_lin_stokes_problem_t  *lin_stokes;
  slabs_nl_stokes_problem_t   *nl_stokes;
  /* other */
  ymir_perf_counter_t perf_counter_single_opt;
  ymir_perf_counter_t perf_counter_single_ref;

  /*
   * Initialize Libraries
   */

  /* initialize ymir and dependent libraries */
  ymir_initialize (argc, argv, mpicomm, NULL, SC_LP_DEFAULT, SC_LP_INFO);

  /* get parallel environment */
  mpiret = MPI_Comm_size (mpicomm, &mpisize); YMIR_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank); YMIR_CHECK_MPI (mpiret);

#ifdef YMIR_ENABLE_OPENMP
  #pragma omp parallel
  {
    #pragma omp single
    {
      ompsize = omp_get_num_threads ();
    }
  }
#else
  ompsize = 1;
#endif

  /*
   * Define & Parse Options
   */

  opt = ymir_options_global_new (argv[0] /* program path */);

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  /* basic options */
  YMIR_OPTIONS_CALLBACK, "help", 'h', 0 /* no callback fn args */,
    ymir_options_print_usage_and_exit_fn, NULL /* no arg usage */,
    "Print usage and exit",
  YMIR_OPTIONS_INIFILE, "options-file", 'f',
    ".ini file with option values",

  /* test of matvec */
  YMIR_OPTIONS_I, "tensor-test-num", '\0', &tensor_apply_num,
    SLABS_PERF_MATVEC_DEFAULT_TENSOR_APPLY_NUM,
    "Number of applications of the element tensor",

  /* test of matvec */
  YMIR_OPTIONS_I, "matvec-test-stress-num", '\0', &stress_apply_num,
    SLABS_PERF_MATVEC_DEFAULT_STRESS_APPLY_NUM,
    "Number of applications of the viscous stress operator",
  YMIR_OPTIONS_I, "matvec-test-stress-nl", '\0', &stress_apply_nl,
    SLABS_PERF_MATVEC_DEFAULT_STRESS_APPLY_NL,
    "Version of viscous stress operator (linear or nonlinear)",
  YMIR_OPTIONS_I, "matvec-test-stress-dirty", '\0', &stress_apply_dirty,
    SLABS_PERF_MATVEC_DEFAULT_STRESS_APPLY_DIRTY,
    "Version of viscous stress operator (dirty or not)",
  YMIR_OPTIONS_I, "matvec-test-stress-optimized", '\0', &stress_apply_optimized,
    SLABS_PERF_MATVEC_DEFAULT_STRESS_APPLY_OPTIMIZED,
    "Version of viscous stress operator (optimized or not)",

  YMIR_OPTIONS_I, "matvec-test-stokes-num", '\0', &stokes_apply_num,
    SLABS_PERF_MATVEC_DEFAULT_STOKES_APPLY_NUM,
    "Number of applications of the Stokes operator",
  YMIR_OPTIONS_I, "matvec-test-stokes-nl", '\0', &stokes_apply_nl,
    SLABS_PERF_MATVEC_DEFAULT_STOKES_APPLY_NL,
    "Version of Stokes operator (linear or nonlinear)",
  YMIR_OPTIONS_I, "matvec-test-stokes-dirty", '\0', &stokes_apply_dirty,
    SLABS_PERF_MATVEC_DEFAULT_STOKES_APPLY_DIRTY,
    "Version of Stokes operator (dirty or not)",
  YMIR_OPTIONS_I, "matvec-test-stokes-optimized", '\0', &stokes_apply_optimized,
    SLABS_PERF_MATVEC_DEFAULT_STOKES_APPLY_OPTIMIZED,
    "Version of Stokes operator (optimized or not)",

  YMIR_OPTIONS_I, "matvec-test-bbt-num", '\0', &bbt_apply_num,
    SLABS_PERF_MATVEC_DEFAULT_BBT_APPLY_NUM,
    "Number of applications of the BB^T operator",
  YMIR_OPTIONS_I, "matvec-test-bbt-optimized", '\0', &bbt_apply_optimized,
    SLABS_PERF_MATVEC_DEFAULT_BBT_APPLY_OPTIMIZED,
    "Version of BB^T operator (optimized or not)",

  YMIR_OPTIONS_I, "matvec-test-stiff-num", '\0', &stiff_apply_num,
    SLABS_PERF_MATVEC_DEFAULT_STIFF_APPLY_NUM,
    "Number of applications of the stiffness operator",
  YMIR_OPTIONS_I, "matvec-test-stiff-dirty", '\0', &stiff_apply_dirty,
    SLABS_PERF_MATVEC_DEFAULT_STIFF_APPLY_DIRTY,
    "Version of stiffness operator (dirty or not)",
  YMIR_OPTIONS_I, "matvec-test-stiff-optimized", '\0', &stiff_apply_optimized,
    SLABS_PERF_MATVEC_DEFAULT_STIFF_APPLY_OPTIMIZED,
    "Version of stiffness operator (optimized or not)",

  /* performance & monitoring options */
  YMIR_OPTIONS_B, "production-run", '\0',
                  &production_run, 0,
    "Execute as a production run (to reduce some overhead and checks)",

  YMIR_OPTIONS_B, "monitor-performance", '\0',
                  &monitor_performance, 0,
    "Measure and print performance statistics, e.g., runtime or flops",

  YMIR_OPTIONS_S, "workload-out-path", '\0',
                  &workload_out_path, NULL,
    "File path for output of global workload for each processor",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add sub-options */
  slabs_setup_add_suboptions (opt);
  ymir_options_add_suboptions_solver_stokes (opt);

  /* parse options */
  {
    int                 optret;

    optret = ymir_options_parse (SC_LP_INFO, opt, argc, argv);
    if (optret < 0) { /* if parsing was not successful */
      ymir_options_print_usage (SC_LP_INFO, opt, NULL /* args usage */);
      YMIR_GLOBAL_INFO ("Option parsing failed\n");
      exit (0);
    }
  }

  /*
   * Initialize Main Program
   */

  YMIR_GLOBAL_PRODUCTIONF ("Into %s (production %i)\n",
      this_fn_name, production_run);
  YMIR_GLOBAL_PRODUCTIONF (
      "Parallel environment: MPI size %i, OpenMP size %i\n",
      mpisize, ompsize);
  ymir_set_up (argc, argv, mpicomm, production_run);

  /* initialize performance counters */
  ymir_perf_counter_init_all (slabs_perf_counter, slabs_perf_counter_name,
                              SLABS_PERF_COUNTER_N, monitor_performance);

  /* start performance counters */
  ymir_perf_counter_start_barrier (
      &slabs_perf_counter[SLABS_PERF_COUNTER_TOTAL], mpicomm);

  /* print & process options */
  ymir_options_print_summary (SC_LP_INFO, opt);
  slabs_setup_process_options (&physics_options, &discr_options,
                               &solver_options);

  /*
   * Setup Mesh
   */

  /* create mesh and Stokes state */
  slabs_setup_mesh (&p8est, &mesh, &press_elem, &state,
                    &enforce_refinement_data, &coarsen_coeff_data,
                    mpicomm, &physics_options, &discr_options, &solver_options,
                    &slabs_perf_counter[SLABS_PERF_COUNTER_SETUP_MESH],
                    workload_out_path);

  /*
   * Setup Stokes Problem
   */

  /* setup Stokes problem */
  slabs_setup_stokes (&lin_stokes, &nl_stokes, p8est, mesh, press_elem, state,
                      &physics_options, &solver_options,
                      &slabs_perf_counter[SLABS_PERF_COUNTER_SETUP_STOKES],
                      workload_out_path);

#if defined(__bgq__)

#if defined(__HAVE_HPM)
  HPM_Start("YMIR_LOOP_Block");
  if(mpirank == 0)
    printf("HPM MPI START\n");
#endif

  trace_start();
  summary_start();
#endif

  /*
   * ###DEV### Performance Test for Element Tensor Products ###DEV###
   */

  if (0 < tensor_apply_num) {
    const ymir_locidx_t n_elements = 10000;
  //const int           N = mesh->ma->N;
  //const int           Nrp = mesh->ma->Nrp;
    const int           Np = mesh->ma->Np;
    const double       *_sc_restrict D = mesh->drst->e[0];
  //const double       *_sc_restrict Dt = mesh->drst_trans->e[0];
  //const double       *_sc_restrict B = mesh->brst->e[0];
  //const double       *_sc_restrict Bt = mesh->brst_trans->e[0];
    double             *_sc_restrict scal_in =
                          YMIR_ALLOC (double, n_elements * Np);
    double             *_sc_restrict scal_out =
                          YMIR_ALLOC (double, n_elements * Np);
    double             *_sc_restrict vel_in =
                          YMIR_ALLOC (double, n_elements * Np * 3);
    double             *_sc_restrict vel_out =
                          YMIR_ALLOC (double, n_elements * Np * 3);
    double             *_sc_restrict vec4_in =
                          YMIR_ALLOC (double, n_elements * Np * 4);
    double             *_sc_restrict vec4_out =
                          YMIR_ALLOC (double, n_elements * Np * 4);
    ymir_locidx_t       elid;
    int                 k;

    ymir_perf_counter_init (&perf_counter_single_opt, "perf_counter_single_opt",
                            monitor_performance);

    for (k = 0; k < n_elements * Np; k++) {
      scal_in[k] = (double) k;
      vel_in[3*k    ] = 0.2 + (double) k;
      vel_in[3*k + 1] = 0.4 + (double) k;
      vel_in[3*k + 2] = 0.6 + (double) k;
      vec4_in[4*k    ] = 0.1 + (double) k;
      vec4_in[4*k + 1] = 0.3 + (double) k;
      vec4_in[4*k + 2] = 0.5 + (double) k;
      vec4_in[4*k + 3] = 0.7 + (double) k;
    }

    YMIR_GLOBAL_INFOF (
        "%s: ###DEV### Run tensor apply %i times on %i elements\n",
        this_fn_name, tensor_apply_num, n_elements);

    /* velocity IIA */

    ymir_perf_counter_start (&perf_counter_single_opt);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_optimized_IIA_vel_N2_opt1 (
            &vel_out[3*Np*elid], D, &vel_in[3*Np*elid]);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_opt);
    YMIR_GLOBAL_INFOF (
        "%s: mangll_tensor_optimized_IIA_vel_N2_opt1 time: %g\n",
        this_fn_name, perf_counter_single_opt.wtime_intv);

    /* vector[4] IIA */

    ymir_perf_counter_start (&perf_counter_single_opt);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_optimized_IIA_vec4_N2_opt1 (
            &vec4_out[4*Np*elid], D, &vec4_in[4*Np*elid]);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_opt);
    YMIR_GLOBAL_INFOF (
        "%s: mangll_tensor_optimized_IIA_vec4_N2_opt1 time: %g\n",
        this_fn_name, perf_counter_single_opt.wtime_intv);

    ymir_perf_counter_start (&perf_counter_single_opt);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_optimized_IIA_vec4_N2_opt2 (
            &vec4_out[4*Np*elid], D, &vec4_in[4*Np*elid]);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_opt);
    YMIR_GLOBAL_INFOF (
        "%s: mangll_tensor_optimized_IIA_vec4_N2_opt2 time: %g\n",
        this_fn_name, perf_counter_single_opt.wtime_intv);

    /* destroy */
    YMIR_FREE (scal_in);
    YMIR_FREE (scal_out);
    YMIR_FREE (vel_in);
    YMIR_FREE (vel_out);
    YMIR_FREE (vec4_in);
    YMIR_FREE (vec4_out);
  }

  /*
   * Performance Test for Element Tensor Products
   */

  if (0 < tensor_apply_num) {
    const ymir_locidx_t n_elements = 10000;
    const int           N = mesh->ma->N;
    const int           Nrp = mesh->ma->Nrp;
    const int           Np = mesh->ma->Np;
    double             *_sc_restrict D = mesh->drst->e[0];
    double             *_sc_restrict Dt = mesh->drst_trans->e[0];
    double             *_sc_restrict B = mesh->brst->e[0];
    double             *_sc_restrict Bt = mesh->brst_trans->e[0];
    double             *_sc_restrict scal_in =
                          YMIR_ALLOC (double, n_elements * Np);
    double             *scal_out_opt, *scal_out_ref;
    double             *_sc_restrict vel_in =
                          YMIR_ALLOC (double, n_elements * Np * 3);
    double             *vel_out_opt, *vel_out_ref;
    double             *vel_dr, *vel_ds, *vel_dt;
    ymir_locidx_t       elid;
    int                 k;

    ymir_perf_counter_init (&perf_counter_single_opt, "perf_counter_single_opt",
                            monitor_performance);
    ymir_perf_counter_init (&perf_counter_single_ref, "perf_counter_single_ref",
                            monitor_performance);

    YMIR_GLOBAL_INFOF ("%s: Run tensor apply %i times on %i elements\n",
                       this_fn_name, tensor_apply_num, n_elements);

    for (k = 0; k < n_elements * Np; k++) {
      scal_in[k] = (double) k;
      vel_in[3*k    ] = 0.2 + (double) k;
      vel_in[3*k + 1] = 0.4 + (double) k;
      vel_in[3*k + 2] = 0.6 + (double) k;
    }

    /*** scalar ***/
    scal_out_opt = YMIR_ALLOC (double, n_elements * Np);
    scal_out_ref = YMIR_ALLOC (double, n_elements * Np);

    /* apply tensor IIA */
    ymir_perf_counter_start (&perf_counter_single_opt);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_optimized_IIA_scal (&scal_out_opt[Np*elid], D,
                                          &scal_in[Np*elid], N);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_opt);

    YMIR_GLOBAL_INFOF (
        "%s: IIA scalar optimized time: %g\n",
        this_fn_name, perf_counter_single_opt.wtime_intv);

#if (1 == SLABS_PERF_TENSOR_REF)
    ymir_perf_counter_start (&perf_counter_single_ref);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_IIAx_apply_elem (Nrp, 1, D,
                                       &scal_in[Np*elid],
                                       &scal_out_ref[Np*elid]);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_ref);

    YMIR_GLOBAL_INFOF (
        "%s: IIA scalar reference time: %g\n",
        this_fn_name, perf_counter_single_ref.wtime_intv);
    YMIR_GLOBAL_INFOF (
        "%s: IIA scalar rel max-error %.3e\n", this_fn_name,
        slabs_perf_tensor_test_get_error (scal_out_opt, scal_out_ref,
                                          n_elements * Np) );
#endif

    /* apply tensor IAI */
    ymir_perf_counter_start (&perf_counter_single_opt);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_optimized_IAI_scal (&scal_out_opt[Np*elid], Dt,
                                          &scal_in[Np*elid], N);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_opt);

    YMIR_GLOBAL_INFOF (
        "%s: IAI scalar optimized time: %g\n",
        this_fn_name, perf_counter_single_opt.wtime_intv);

#if (1 == SLABS_PERF_TENSOR_REF)
    ymir_perf_counter_start (&perf_counter_single_ref);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_IAIx_apply_elem (Nrp, 1, D,
                                       &scal_in[Np*elid],
                                       &scal_out_ref[Np*elid]);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_ref);

    YMIR_GLOBAL_INFOF (
        "%s: IAI scalar reference time: %g\n",
        this_fn_name, perf_counter_single_ref.wtime_intv);
    YMIR_GLOBAL_INFOF (
        "%s: IAI scalar rel max-error %.3e\n", this_fn_name,
        slabs_perf_tensor_test_get_error (scal_out_opt, scal_out_ref,
                                          n_elements * Np) );
#endif

    /* apply tensor AII */
    ymir_perf_counter_start (&perf_counter_single_opt);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_optimized_AII_scal (&scal_out_opt[Np*elid], Dt,
                                          &scal_in[Np*elid], N);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_opt);

    YMIR_GLOBAL_INFOF (
        "%s: AII scalar optimized time: %g\n",
        this_fn_name, perf_counter_single_opt.wtime_intv);

#if (1 == SLABS_PERF_TENSOR_REF)
    ymir_perf_counter_start (&perf_counter_single_ref);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_AIIx_apply_elem (Nrp, 1, D,
                                       &scal_in[Np*elid],
                                       &scal_out_ref[Np*elid]);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_ref);

    YMIR_GLOBAL_INFOF (
        "%s: AII scalar reference time: %g\n",
        this_fn_name, perf_counter_single_ref.wtime_intv);
    YMIR_GLOBAL_INFOF (
        "%s: AII scalar rel max-error %.3e\n", this_fn_name,
        slabs_perf_tensor_test_get_error (scal_out_opt, scal_out_ref,
                                          n_elements * Np) );
#endif

    /* destroy */
    YMIR_FREE (scal_in);
    YMIR_FREE (scal_out_opt);
    YMIR_FREE (scal_out_ref);

    /*** velocity ***/
    vel_out_opt = YMIR_ALLOC (double, n_elements * Np * 3);
    vel_out_ref = YMIR_ALLOC (double, n_elements * Np * 3);
    vel_dr = YMIR_ALLOC (double, n_elements * Np * 3);
    vel_ds = YMIR_ALLOC (double, n_elements * Np * 3);
    vel_dt = YMIR_ALLOC (double, n_elements * Np * 3);

    /* apply tensor IIA */
    ymir_perf_counter_start (&perf_counter_single_opt);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_optimized_IIA_vel (&vel_out_opt[3*Np*elid], D,
                                         &vel_in[3*Np*elid], N);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_opt);

    YMIR_GLOBAL_INFOF (
        "%s: IIA velocity optimized time: %g\n",
        this_fn_name, perf_counter_single_opt.wtime_intv);

#if (1 == SLABS_PERF_TENSOR_REF)
    ymir_perf_counter_start (&perf_counter_single_ref);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_IIAx_apply_elem (Nrp, 3, D,
                                       &vel_in[3*Np*elid],
                                       &vel_out_ref[3*Np*elid]);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_ref);

    YMIR_GLOBAL_INFOF (
        "%s: IIA velocity reference time: %g\n",
        this_fn_name, perf_counter_single_ref.wtime_intv);
    YMIR_GLOBAL_INFOF (
        "%s: IIA velocity rel max-error %.3e\n", this_fn_name,
        slabs_perf_tensor_test_get_error (vel_out_opt, vel_out_ref,
                                          n_elements * Np * 3) );
#endif

    /* apply tensor IAI */
    ymir_perf_counter_start (&perf_counter_single_opt);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_optimized_IAI_vel (&vel_out_opt[3*Np*elid], Dt,
                                         &vel_in[3*Np*elid], N);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_opt);

    YMIR_GLOBAL_INFOF (
        "%s: IAI velocity optimized time: %g\n",
        this_fn_name, perf_counter_single_opt.wtime_intv);

#if (1 == SLABS_PERF_TENSOR_REF)
    ymir_perf_counter_start (&perf_counter_single_ref);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_IAIx_apply_elem (Nrp, 3, D,
                                       &vel_in[3*Np*elid],
                                       &vel_out_ref[3*Np*elid]);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_ref);

    YMIR_GLOBAL_INFOF (
        "%s: IAI velocity reference time: %g\n",
        this_fn_name, perf_counter_single_ref.wtime_intv);
    YMIR_GLOBAL_INFOF (
        "%s: IAI velocity rel max-error %.3e\n", this_fn_name,
        slabs_perf_tensor_test_get_error (vel_out_opt, vel_out_ref,
                                          n_elements * Np * 3) );
#endif

    /* apply tensor AII */
    ymir_perf_counter_start (&perf_counter_single_opt);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_optimized_AII_vel (&vel_out_opt[3*Np*elid], Dt,
                                         &vel_in[3*Np*elid], N);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_opt);

    YMIR_GLOBAL_INFOF (
        "%s: AII velocity optimized time: %g\n",
        this_fn_name, perf_counter_single_opt.wtime_intv);

#if (1 == SLABS_PERF_TENSOR_REF)
    ymir_perf_counter_start (&perf_counter_single_ref);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_AIIx_apply_elem (Nrp, 3, D,
                                       &vel_in[3*Np*elid],
                                       &vel_out_ref[3*Np*elid]);
      }
    }
    ymir_perf_counter_stop_add (&perf_counter_single_ref);

    YMIR_GLOBAL_INFOF (
        "%s: AII velocity reference time: %g\n",
        this_fn_name, perf_counter_single_ref.wtime_intv);
    YMIR_GLOBAL_INFOF (
        "%s: AII velocity rel max-error %.3e\n", this_fn_name,
        slabs_perf_tensor_test_get_error (vel_out_opt, vel_out_ref,
                                          n_elements * Np * 3) );
#endif

    ymir_perf_counter_copy (&slabs_perf_counter[SLABS_PERF_COUNTER_TENSOR_OPT],
                            &perf_counter_single_opt);
    ymir_perf_counter_copy (&slabs_perf_counter[SLABS_PERF_COUNTER_TENSOR_REF],
                            &perf_counter_single_ref);

    /* apply tensor gradient */
    for (elid = 0; elid < n_elements; elid++) {
      mangll_tensor_optimized_grad_vel (
          vel_dr, vel_ds, vel_dt, D, Dt, B, Bt, &vel_in[3*Np*elid], N);
    }
    for (k = 0; k < n_elements * Np; k++) {
      vel_in[3*k    ] = 0.2 + (double) k;
      vel_in[3*k + 1] = 0.4 + (double) k;
      vel_in[3*k + 2] = 0.6 + (double) k;
    }
    ymir_perf_counter_start (
        &slabs_perf_counter[SLABS_PERF_COUNTER_TENSOR_GRAD_VEL_OPT]);
    for (k = 0; k < tensor_apply_num; k++ ){
      for (elid = 0; elid < n_elements; elid++) {
        mangll_tensor_optimized_grad_vel (
            vel_dr, vel_ds, vel_dt, D, Dt, B, Bt, &vel_in[3*Np*elid], N);
      }
    }
    ymir_perf_counter_stop_add (
        &slabs_perf_counter[SLABS_PERF_COUNTER_TENSOR_GRAD_VEL_OPT]);

    /* destroy */
    YMIR_FREE (vel_in);
    YMIR_FREE (vel_out_opt);
    YMIR_FREE (vel_out_ref);
    YMIR_FREE (vel_dr);
    YMIR_FREE (vel_ds);
    YMIR_FREE (vel_dt);
  }

  /*
   * Performance Test for Stress Matvec
   */

  if (0 < stress_apply_num) {
    const int           nl = stress_apply_nl;
    const int           dirty = stress_apply_dirty;
    const int           optimized_matvec = stress_apply_optimized;
    ymir_stress_op_t   *stress_op;
    ymir_linearop_fn    fn;
    ymir_linearop_fn    fn_test_against = NULL;
    ymir_vec_t         *in = ymir_cvec_new (mesh, 3);
    ymir_vec_t         *out = ymir_cvec_new (mesh, 3);
    int                 k;

    YMIR_GLOBAL_INFOF (
        "%s: Run Stress apply (nl %i, dirty %i, opt %i) matvec %i times\n",
        this_fn_name, nl, dirty, optimized_matvec, stress_apply_num);

    /* get stress operator */
    if (stress_apply_nl) {
      YMIR_ASSERT (nl_stokes != NULL);
      if (nl_stokes->stokes_op == NULL) {
        slabs_nonlinear_stokes_op_init (nl_stokes, state, &physics_options, 1,
                                        SL_NL_SOLVER_NEWTON,
                                        SL_NL_SOLVER_PRIMALDUAL_NONE,
                                        SL_NL_SOLVER_PRIMALDUAL_SCAL_NONE, 0);
      }
      stress_op = nl_stokes->stokes_op->stress_op;
    }
    else {
      YMIR_ASSERT (lin_stokes != NULL);
      YMIR_ASSERT (lin_stokes->stokes_op != NULL);
      stress_op = lin_stokes->stokes_op->stress_op;
    }

    /* set random input vector */
    slabs_perf_matvec_set_random_vec (in);

    /* get stress operator function */
    if (!nl && !dirty) {
      if (optimized_matvec) {
        /* linear, not dirty, optimized */
        fn = (ymir_linearop_fn) ymir_stress_op_optimized_apply;
        fn_test_against = (ymir_linearop_fn) ymir_stress_op_apply;
      }
      else {
        /* linear, not dirty, not optimized */
        fn = (ymir_linearop_fn) ymir_stress_op_apply;
      }
    }
    else if (nl && !dirty) {
      if (optimized_matvec) {
        /* nonlinear, not dirty, optimized */
        fn = (ymir_linearop_fn) ymir_nlstress_op_optimized_apply;
        fn_test_against = (ymir_linearop_fn) ymir_nlstress_op_apply;
      }
      else {
        /* nonlinear, not dirty, not optimized */
        fn = (ymir_linearop_fn) ymir_nlstress_op_apply;
      }
    }
    else if (!nl && dirty) {
      if (optimized_matvec) {
        /* linear, dirty, optimized */
        fn = (ymir_linearop_fn) ymir_stress_op_optimized_apply_dirty;
        fn_test_against = (ymir_linearop_fn) ymir_stress_op_apply_dirty;
      }
      else {
        /* linear, dirty, not optimized */
        fn = (ymir_linearop_fn) ymir_stress_op_apply_dirty;
      }
    }
    else if (nl && dirty) {
      if (optimized_matvec) {
        /* nonlinear, dirty, optimized */
        fn = (ymir_linearop_fn) ymir_nlstress_op_optimized_apply_dirty;
        fn_test_against = (ymir_linearop_fn) ymir_nlstress_op_apply_dirty;
      }
      else {
        /* nonlinear, dirty, not optimized */
        fn = (ymir_linearop_fn) ymir_nlstress_op_apply_dirty;
      }
    }
    else {
      YMIR_ABORT_NOT_REACHED ();
    }

    /* test function */
    if (fn_test_against != NULL) {
      slabs_perf_matvec_test_stress (fn, fn_test_against, stress_op);
    }

    /* start performance counters */
    ymir_perf_counter_start_barrier (
        &slabs_perf_counter[SLABS_PERF_COUNTER_STRESS_MATVEC], mpicomm);

    /* apply stress operator */
    for (k = 0; k < stress_apply_num; k++) {
      fn (in, out, stress_op);
    }

    /* stop performance counters */
    ymir_perf_counter_stop_add_barrier (
        &slabs_perf_counter[SLABS_PERF_COUNTER_STRESS_MATVEC], mpicomm);

    /* destroy */
    ymir_vec_destroy (in);
    ymir_vec_destroy (out);
  }

  /*
   * Performance Test for Stokes Matvec
   */

  if (0 < stokes_apply_num) {
    const int           nl = stokes_apply_nl;
    const int           dirty = stokes_apply_dirty;
    const int           optimized_matvec = stokes_apply_optimized;
    ymir_stokes_op_t   *stokes_op;
    ymir_linearop_fn    fn;
    ymir_linearop_fn    fn_test_against = NULL;
    ymir_vec_t         *in = ymir_stokes_vec_new (mesh, press_elem);
    ymir_vec_t         *out = ymir_stokes_vec_new (mesh, press_elem);
    int                 k;

    YMIR_GLOBAL_INFOF (
        "%s: Run Stokes apply (nl %i, dirty %i, opt %i) matvec %i times\n",
        this_fn_name, nl, dirty, optimized_matvec, stokes_apply_num);

    /* get Stokes operator */
    if (stokes_apply_nl) {
      YMIR_ASSERT (nl_stokes != NULL);
      if (nl_stokes->stokes_op == NULL) {
        slabs_nonlinear_stokes_op_init (nl_stokes, state, &physics_options, 1,
                                        SL_NL_SOLVER_NEWTON,
                                        SL_NL_SOLVER_PRIMALDUAL_NONE,
                                        SL_NL_SOLVER_PRIMALDUAL_SCAL_NONE, 0);
      }
      stokes_op = nl_stokes->stokes_op;
    }
    else {
      YMIR_ASSERT (lin_stokes != NULL);
      YMIR_ASSERT (lin_stokes->stokes_op != NULL);
      stokes_op = lin_stokes->stokes_op;
    }

    /* set random input vector */
    slabs_perf_matvec_set_random_vec (in);

    /* get Stokes operator function */
    if (!nl && !dirty) {
      if (optimized_matvec) {
        /* linear, not dirty, optimized */
        fn = (ymir_linearop_fn) ymir_stokes_op_optimized_apply;
        fn_test_against = (ymir_linearop_fn) ymir_stokes_op_apply;
      }
      else {
        /* linear, not dirty, not optimized */
        fn = (ymir_linearop_fn) ymir_stokes_op_apply;
      }
    }
    else if (nl && !dirty) {
      if (optimized_matvec) {
        /* nonlinear, not dirty, optimized */
        fn = (ymir_linearop_fn) ymir_nlstokes_op_optimized_apply;
        fn_test_against = (ymir_linearop_fn) ymir_nlstokes_op_apply;
      }
      else {
        /* nonlinear, not dirty, not optimized */
        fn = (ymir_linearop_fn) ymir_nlstokes_op_apply;
      }
    }
    else if (!nl && dirty) {
      if (optimized_matvec) {
        /* linear, dirty, optimized */
        fn = (ymir_linearop_fn) ymir_stokes_op_optimized_apply_dirty;
        fn_test_against = (ymir_linearop_fn) ymir_stokes_op_apply_dirty;
      }
      else {
        /* linear, dirty, not optimized */
        fn = (ymir_linearop_fn) ymir_stokes_op_apply_dirty;
      }
    }
    else if (nl && dirty) {
      if (optimized_matvec) {
        /* nonlinear, dirty, optimized */
        fn = (ymir_linearop_fn) ymir_nlstokes_op_optimized_apply_dirty;
        fn_test_against = (ymir_linearop_fn) ymir_nlstokes_op_apply_dirty;
      }
      else {
        /* nonlinear, dirty, not optimized */
        fn = (ymir_linearop_fn) ymir_nlstokes_op_apply_dirty;
      }
    }
    else {
      YMIR_ABORT_NOT_REACHED ();
    }

    /* test function */
    if (fn_test_against != NULL) {
      slabs_perf_matvec_test_stokes (fn, fn_test_against, stokes_op);
    }

    /* start performance counters */
    ymir_perf_counter_start_barrier (
        &slabs_perf_counter[SLABS_PERF_COUNTER_STOKES_MATVEC], mpicomm);

    /* apply Stokes operator */
    for (k = 0; k < stokes_apply_num; k++) {
      fn (in, out, stokes_op);
    }

    /* stop performance counters */
    ymir_perf_counter_stop_add_barrier (
        &slabs_perf_counter[SLABS_PERF_COUNTER_STOKES_MATVEC], mpicomm);

    /* destroy */
    ymir_vec_destroy (in);
    ymir_vec_destroy (out);
  }

  /*
   * Performance Test for BB^T Matvec
   */

  if (0 < bbt_apply_num) {
    const int           optimized_matvec = bbt_apply_optimized;
    ymir_bbt_t         *bbt;
    ymir_linearop_fn    fn;
    ymir_linearop_fn    fn_test_against = NULL;
    ymir_vec_t         *in = ymir_pressure_vec_new (mesh, press_elem);
    ymir_vec_t         *out = ymir_pressure_vec_new (mesh, press_elem);
    int                 k;

    YMIR_GLOBAL_INFOF (
        "%s: Run BB^T apply (opt %i) matvec %i times\n",
        this_fn_name, optimized_matvec, bbt_apply_num);

    /* get BB^T operator */
    if (nl_stokes != NULL) {
      if (nl_stokes->stokes_pc == NULL) {
        slabs_nonlinear_stokes_pc_init (
            nl_stokes, state, SL_NL_STOKES_PROB_SCHUR_DIAG_INV_VISC_PMASS,
            SL_NL_STOKES_PROB_SCALING_NONE, 0);
      }
      YMIR_ASSERT (nl_stokes->stokes_pc->bfbt != NULL);
      YMIR_ASSERT (nl_stokes->stokes_pc->bfbt->right_bbt != NULL);
      bbt = nl_stokes->stokes_pc->bfbt->right_bbt;
    }
    else {
      YMIR_ASSERT (lin_stokes != NULL);
      YMIR_ASSERT (lin_stokes->stokes_pc->bfbt != NULL);
      YMIR_ASSERT (lin_stokes->stokes_pc->bfbt->right_bbt != NULL);
      bbt = lin_stokes->stokes_pc->bfbt->right_bbt;
    }

    /* set random input vector */
    slabs_perf_matvec_set_random_vec (in);

    /* get BB^T operator function */
    if (optimized_matvec) {
      /* not dirty, optimized */
      fn = (ymir_linearop_fn) ymir_bbt_optimized_apply;
      fn_test_against = (ymir_linearop_fn) ymir_bbt_apply;
    }
    else {
      /* not dirty, not optimized */
      fn = (ymir_linearop_fn) ymir_bbt_apply;
    }

    /* test function */
    if (fn_test_against != NULL) {
      slabs_perf_matvec_test_bbt (fn, fn_test_against, bbt);
    }

    /* start performance counters */
    ymir_perf_counter_start_barrier (
        &slabs_perf_counter[SLABS_PERF_COUNTER_BBT_MATVEC], mpicomm);

    /* apply BB^T operator */
    for (k = 0; k < bbt_apply_num; k++) {
      fn (in, out, bbt);
    }

    /* stop performance counters */
    ymir_perf_counter_stop_add_barrier (
        &slabs_perf_counter[SLABS_PERF_COUNTER_BBT_MATVEC], mpicomm);

    /* destroy */
    ymir_vec_destroy (in);
    ymir_vec_destroy (out);
  }

  /*
   * Performance Test for Stiffness Matvec
   */

  if (0 < stiff_apply_num) {
    const int           dirty = stiff_apply_dirty;
    const int           optimized_matvec = stiff_apply_optimized;
    ymir_stiff_op_t    *stiff_op;
    ymir_linearop_fn    fn;
    ymir_linearop_fn    fn_test_against = NULL;
    ymir_vec_t         *in = ymir_cvec_new (mesh, 1);
    ymir_vec_t         *out = ymir_cvec_new (mesh, 1);
    int                 k;

    YMIR_GLOBAL_INFOF (
        "%s: Run stiffness apply (dirty %i, opt %i) matvec %i times\n",
        this_fn_name, dirty, optimized_matvec, stiff_apply_num);

    /* get stiffness operator */
    if (nl_stokes != NULL) {
      YMIR_ASSERT (nl_stokes->stokes_pc->bfbt != NULL);
      YMIR_ASSERT (nl_stokes->stokes_pc->bfbt->right_bbt != NULL);
      YMIR_ASSERT (nl_stokes->stokes_pc->bfbt->right_bbt->stiff_op != NULL);
      stiff_op = nl_stokes->stokes_pc->bfbt->right_bbt->stiff_op;
    }
    else {
      YMIR_ASSERT (lin_stokes != NULL);
      YMIR_ASSERT (lin_stokes->stokes_pc->bfbt != NULL);
      YMIR_ASSERT (lin_stokes->stokes_pc->bfbt->right_bbt != NULL);
      YMIR_ASSERT (lin_stokes->stokes_pc->bfbt->right_bbt->stiff_op != NULL);
      stiff_op = lin_stokes->stokes_pc->bfbt->right_bbt->stiff_op;
    }

    /* set random input vector */
    slabs_perf_matvec_set_random_vec (in);

    /* get stiffness operator function */
    if (!dirty) {
      if (optimized_matvec) {
        /* not dirty, optimized */
        fn = (ymir_linearop_fn) ymir_stiff_op_optimized_apply;
        fn_test_against = (ymir_linearop_fn) ymir_stiff_op_apply;
      }
      else {
        /* not dirty, not optimized */
        fn = (ymir_linearop_fn) ymir_stiff_op_apply;
      }
    }
    else {
      if (optimized_matvec) {
        /* dirty, optimized */
        fn = (ymir_linearop_fn) ymir_stiff_op_optimized_apply_dirty;
        fn_test_against = (ymir_linearop_fn) ymir_stiff_op_apply_dirty;
      }
      else {
        /* dirty, not optimized */
        fn = (ymir_linearop_fn) ymir_stiff_op_apply_dirty;
      }
    }

    /* test function */
    if (fn_test_against != NULL) {
      slabs_perf_matvec_test_stiff (fn, fn_test_against, stiff_op);
    }

    /* start performance counters */
    ymir_perf_counter_start_barrier (
        &slabs_perf_counter[SLABS_PERF_COUNTER_STIFF_MATVEC], mpicomm);

    /* apply stiffness operator */
    for (k = 0; k < stiff_apply_num; k++) {
      fn (in, out, stiff_op);
    }

    /* stop performance counters */
    ymir_perf_counter_stop_add_barrier (
        &slabs_perf_counter[SLABS_PERF_COUNTER_STIFF_MATVEC], mpicomm);

    /* destroy */
    ymir_vec_destroy (in);
    ymir_vec_destroy (out);
  }

#if defined(__bgq__)
  summary_stop();
  trace_stop();

#if defined(__HAVE_HPM)
  HPM_Stop("YMIR_LOOP_Block");
  if(mpirank == 0)
    printf("HPM MPI STOP\n");
#endif

#endif

  /*
   * Finalize
   */

  /* stop performance counters */
  ymir_perf_counter_stop_add_barrier (
      &slabs_perf_counter[SLABS_PERF_COUNTER_TOTAL], mpicomm);

  /* destroy Stokes problem, Stokes state, and mesh */
  slabs_clear (lin_stokes, nl_stokes, p8est, mesh, press_elem, state,
               enforce_refinement_data, coarsen_coeff_data,
               &physics_options, &discr_options);

  /* destroy options */
  YMIR_FREE (discr_options.refine_radius);

  /* print performance statistics */
  if (0 < stress_apply_num) {
    ymir_stress_op_perf_counter_print ();            /* Stress Op */
  }
  if (0 < stokes_apply_num) {
    ymir_stokes_op_perf_counter_print ();            /* Stokes Op */
  }
  if (0 < bbt_apply_num) {
    ymir_bbt_perf_counter_print ();                  /* BBT */
  }
  if (0 < stiff_apply_num) {
    ymir_stiff_op_perf_counter_print ();             /* Stiffness Op */
  }
  //ymir_gmg_hierarchy_mesh_perf_counter_print ();   /* GMG mesh */
  //ymir_stress_pc_perf_counter_print ();            /* Stress PC */
  //ymir_gmg_hierarchy_stress_perf_counter_print (); /* GMG stress */
  //ymir_stiff_pc_perf_counter_print ();             /* Stiffness PC */
  //ymir_gmg_hierarchy_stiff_perf_counter_print ();  /* GMG stiffness */
  //ymir_gmg_hierarchy_bbt_perf_counter_print ();    /* GMG BB^T */
  //ymir_stokes_pc_perf_counter_print ();            /* Stokes PC */

  /* print slabs performance statistics */
  if (monitor_performance) {
    const int           print_wtime = 1;
    const int           print_n_calls = 0;
    const int           print_flops = 1;
    int                 n_stats;

    /* gather and print slabs statistics */
    n_stats = ymir_perf_counter_gather_stats (
        slabs_perf_counter, SLABS_PERF_COUNTER_N,
        slabs_perf_stats, slabs_perf_stats_name, mpicomm,
        print_wtime, print_n_calls, print_flops);
    ymir_perf_counter_print_stats (
        slabs_perf_stats, n_stats, "Slabs performance test, matvec");
  }

  /* destroy options */
  ymir_options_global_destroy ();

  /* print that this function is ending */
  YMIR_GLOBAL_PRODUCTIONF ("Done %s\n", this_fn_name);

  /* finalize ymir */
  ymir_finalize ();

  return 0;
}

