/*
 */

#include <rhea_io_mpi.h>
#include <rhea_io_std.h>
#include <rhea_base.h>
#include <ymir_perf_counter.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_IO_MPI_DEFAULT_MONITOR_PERFORMANCE (0)

int                 rhea_io_mpi_monitor_performance =
                      RHEA_IO_MPI_DEFAULT_MONITOR_PERFORMANCE;

void
rhea_io_mpi_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "IO-MPI";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_B, "monitor-performance", '\0',
    &(rhea_io_mpi_monitor_performance), RHEA_IO_MPI_DEFAULT_MONITOR_PERFORMANCE,
    "Measure and print performance statistics (e.g., runtime or flops)",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);

  /* initialize (deactivated) performance counters */
  rhea_io_mpi_perfmon_init (0, 0);
}

/******************************************************************************
 * Monitoring
 *****************************************************************************/

/* perfomance monitor tags and names */
typedef enum
{
  RHEA_IO_MPI_PERFMON_READ_BCAST,
  RHEA_IO_MPI_PERFMON_READ_SCATTER,
  RHEA_IO_MPI_PERFMON_GATHER_WRITE,
  RHEA_IO_MPI_PERFMON_N
}
rhea_io_mpi_perfmon_idx_t;

static const char  *rhea_io_mpi_perfmon_name[RHEA_IO_MPI_PERFMON_N] =
{
  "Read & Broadcast",
  "Read & Scatter",
  "Gather & Write"
};
ymir_perf_counter_t rhea_io_mpi_perfmon[RHEA_IO_MPI_PERFMON_N];

void
rhea_io_mpi_perfmon_init (const int activate, const int skip_if_active)
{
  const int           active = activate && rhea_io_mpi_monitor_performance;

  ymir_perf_counter_init_all_ext (rhea_io_mpi_perfmon,
                                  rhea_io_mpi_perfmon_name,
                                  RHEA_IO_MPI_PERFMON_N,
                                  active, skip_if_active);
}

void
rhea_io_mpi_perfmon_print (sc_MPI_Comm mpicomm,
                           const int print_wtime,
                           const int print_n_calls,
                           const int print_flops)
{
  const int           active = rhea_io_mpi_monitor_performance;
  const int           print = (print_wtime || print_n_calls || print_flops);
  int                 n_stats = RHEA_IO_MPI_PERFMON_N *
                                YMIR_PERF_COUNTER_N_STATS;
  sc_statinfo_t       stats[n_stats];
  char                stats_name[n_stats][YMIR_PERF_COUNTER_NAME_SIZE];

  /* exit if nothing to do */
  if (!active || !print) {
    return;
  }

  /* gather performance statistics */
  n_stats = ymir_perf_counter_gather_stats (
      rhea_io_mpi_perfmon, RHEA_IO_MPI_PERFMON_N, stats, stats_name,
      mpicomm, print_wtime, print_n_calls, print_flops);

  /* print performance statistics */
  ymir_perf_counter_print_stats (stats, n_stats, "I/O MPI");
}

/******************************************************************************
 * Read
 *****************************************************************************/

int
rhea_io_mpi_read_broadcast_double (double *values_all,
                                   int n_entries,
                                   const char *file_path_bin,
                                   const char *file_path_txt,
                                   sc_MPI_Comm mpicomm)
{
  int                 mpirank, mpiret;
  size_t              total_size = (size_t) n_entries;

  /* start performance monitors */
  ymir_perf_counter_start (
      &rhea_io_mpi_perfmon[RHEA_IO_MPI_PERFMON_READ_BCAST]);

  if (file_path_txt == NULL) { /* if read from binary file */
    RHEA_GLOBAL_INFOF ("Into %s (%s, #entries %i)\n", __func__,
                       file_path_bin, n_entries);
  }
  else {
    RHEA_GLOBAL_INFOF ("Into %s (bin: %s, txt: %s)\n", __func__,
                       file_path_bin, file_path_txt);
  }

  /* check input */
  RHEA_ASSERT (file_path_txt != NULL ||
               (file_path_bin != NULL && 0 < n_entries));

  /* read file with one processor */
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);
  if (mpirank == 0) {
    if (file_path_txt == NULL) { /* if read from binary file */
      rhea_io_std_read_double (values_all, total_size, file_path_bin);
    }
    else { /* otherwise read from text file */
      total_size = rhea_io_std_read_double_from_txt (values_all, total_size,
                                                     file_path_txt);
      if (file_path_bin != NULL) {
        rhea_io_std_write_double (file_path_bin, values_all, total_size);
      }
    }
  }

  /* broadcast values to all processors */
  if (n_entries <= 0) { /* if #entries is not known */
    n_entries = (int) total_size;
    mpiret = sc_MPI_Bcast (&n_entries, 1, sc_MPI_INT, 0, mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  mpiret = sc_MPI_Bcast (values_all, n_entries, sc_MPI_DOUBLE, 0, mpicomm);
  SC_CHECK_MPI (mpiret);

  if (file_path_txt == NULL) { /* if read from binary file */
    RHEA_GLOBAL_INFOF ("Done %s (%s)\n", __func__, file_path_bin);
  }
  else {
    RHEA_GLOBAL_INFOF ("Done %s (bin: %s, txt: %s)\n", __func__,
                       file_path_bin, file_path_txt);
  }

  /* stop performance monitors */
  ymir_perf_counter_stop_add (
      &rhea_io_mpi_perfmon[RHEA_IO_MPI_PERFMON_READ_BCAST]);

  /* return number of entries read */
  return n_entries;
}

int
rhea_io_mpi_read_broadcast_int (int *values_all,
                                int n_entries,
                                const char *file_path_bin,
                                const char *file_path_txt,
                                sc_MPI_Comm mpicomm)
{
  int                 mpirank, mpiret;
  size_t              total_size = (size_t) n_entries;

  /* start performance monitors */
  ymir_perf_counter_start (
      &rhea_io_mpi_perfmon[RHEA_IO_MPI_PERFMON_READ_BCAST]);

  if (file_path_txt == NULL) { /* if read from binary file */
    RHEA_GLOBAL_INFOF ("Into %s (%s, #entries %i)\n", __func__,
                       file_path_bin, n_entries);
  }
  else {
    RHEA_GLOBAL_INFOF ("Into %s (bin: %s, txt: %s)\n", __func__,
                       file_path_bin, file_path_txt);
  }

  /* check input */
  RHEA_ASSERT (file_path_txt != NULL ||
               (file_path_bin != NULL && 0 < n_entries));

  /* read file with one processor */
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);
  if (mpirank == 0) {
    if (file_path_txt == NULL) { /* if read from binary file */
      rhea_io_std_read_int (values_all, total_size, file_path_bin);
    }
    else { /* otherwise read from text file */
      total_size = rhea_io_std_read_int_from_txt (values_all, total_size,
                                                  file_path_txt);
      if (file_path_bin != NULL) {
        rhea_io_std_write_int (file_path_bin, values_all, total_size);
      }
    }
  }

  /* broadcast values to all processors */
  if (n_entries <= 0) { /* if #entries is not known */
    n_entries = (int) total_size;
    mpiret = sc_MPI_Bcast (&n_entries, 1, sc_MPI_INT, 0, mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  mpiret = sc_MPI_Bcast (values_all, n_entries, sc_MPI_INT, 0, mpicomm);
  SC_CHECK_MPI (mpiret);

  if (file_path_txt == NULL) { /* if read from binary file */
    RHEA_GLOBAL_INFOF ("Done %s (%s)\n", __func__, file_path_bin);
  }
  else {
    RHEA_GLOBAL_INFOF ("Done %s (bin: %s, txt: %s)\n", __func__,
                       file_path_bin, file_path_txt);
  }

  /* stop performance monitors */
  ymir_perf_counter_stop_add (
      &rhea_io_mpi_perfmon[RHEA_IO_MPI_PERFMON_READ_BCAST]);

  /* return number of entries read */
  return n_entries;
}

void
rhea_io_mpi_read_scatter_double (double *values_segment,
                                 const int *segment_offset,
                                 const char *file_path_bin,
                                 const char *file_path_txt,
                                 sc_MPI_Comm mpicomm)
{
  int                 mpisize, mpirank, mpiret;
  double             *values_all;
  int                *segment_size;
  int                 r;

  /* start performance monitors */
  ymir_perf_counter_start (
      &rhea_io_mpi_perfmon[RHEA_IO_MPI_PERFMON_READ_SCATTER]);

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

  if (file_path_txt == NULL) { /* if read from binary file */
    RHEA_GLOBAL_INFOF ("Into %s (%s, #entries %i)\n", __func__,
                       file_path_bin, segment_offset[mpisize]);
  }
  else {
    RHEA_GLOBAL_INFOF ("Into %s (bin: %s, txt: %s)\n", __func__,
                       file_path_bin, file_path_txt);
  }

  /* check input */
  RHEA_ASSERT (file_path_txt != NULL || file_path_bin != NULL);
  RHEA_ASSERT (0 < segment_offset[mpisize]);

  /* read file with one processor */
  if (mpirank == 0) {
    size_t              total_size = (size_t) segment_offset[mpisize];

    values_all = RHEA_ALLOC (double, total_size);
    if (file_path_txt == NULL) { /* if read from binary file */
      rhea_io_std_read_double (values_all, total_size, file_path_bin);
    }
    else { /* otherwise read from text file */
      total_size = rhea_io_std_read_double_from_txt (values_all, total_size,
                                                     file_path_txt);
      if (file_path_bin != NULL) {
        rhea_io_std_write_double (file_path_bin, values_all, total_size);
      }
    }
  }
  else {
    values_all = NULL;
  }

  /* create segment size array */
  segment_size = RHEA_ALLOC (int, mpisize);
  for (r = 0; r < mpisize; r++) {
    segment_size[r] = segment_offset[r+1] - segment_offset[r];
  }

  /* scatter segments to corresponding processors */
  mpiret = MPI_Scatterv (
      values_all, segment_size, segment_offset, sc_MPI_DOUBLE,
      values_segment, segment_size[mpirank], sc_MPI_DOUBLE, 0, mpicomm);
  SC_CHECK_MPI (mpiret);
  //TODO `sc_MPI_Scatterv` is not defined

  /* destroy */
  if (mpirank == 0) {
    RHEA_FREE (values_all);
  }
  RHEA_FREE (segment_size);

  if (file_path_txt == NULL) { /* if read from binary file */
    RHEA_GLOBAL_INFOF ("Done %s (%s)\n", __func__, file_path_bin);
  }
  else {
    RHEA_GLOBAL_INFOF ("Done %s (bin: %s, txt: %s)\n", __func__,
                       file_path_bin, file_path_txt);
  }

  /* stop performance monitors */
  ymir_perf_counter_stop_add (
      &rhea_io_mpi_perfmon[RHEA_IO_MPI_PERFMON_READ_SCATTER]);
}

/******************************************************************************
 * Write
 *****************************************************************************/

void
rhea_io_mpi_gather_write_double_to_txt (const char *file_path_txt,
                                        const double *values_segment,
                                        const int *segment_offset,
                                        int n_entries_per_line,
                                        sc_MPI_Comm mpicomm)
{
  int                 mpisize, mpirank, mpiret;
  double             *values_all;
  size_t              total_size;
  int                *segment_size;
  int                 r;

  /* start performance monitors */
  ymir_perf_counter_start (
      &rhea_io_mpi_perfmon[RHEA_IO_MPI_PERFMON_GATHER_WRITE]);

  RHEA_GLOBAL_INFOF ("Into %s (%s)\n", __func__, file_path_txt);

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

  /* set total (=global) size */
  RHEA_ASSERT (0 < segment_offset[mpisize]);
  total_size = (size_t) segment_offset[mpisize];

  /* create segment size array */
  segment_size = RHEA_ALLOC (int, mpisize);
  for (r = 0; r < mpisize; r++) {
    segment_size[r] = segment_offset[r+1] - segment_offset[r];
  }

  /* gather segments from all processors */
  if (mpirank == 0) {
    values_all = RHEA_ALLOC (double, total_size);
  }
  else {
    values_all = NULL;
  }
  mpiret = MPI_Gatherv (
      values_segment, segment_size[mpirank], sc_MPI_DOUBLE,
      values_all, segment_size, segment_offset, sc_MPI_DOUBLE, 0, mpicomm);
  SC_CHECK_MPI (mpiret);
  //TODO `sc_MPI_Gatherv` is not defined

  /* write file with one processor */
  if (mpirank == 0) {
    rhea_io_std_write_double_to_txt (file_path_txt, values_all, total_size,
                                     n_entries_per_line);
  }

  /* destroy */
  if (mpirank == 0) {
    RHEA_FREE (values_all);
  }
  RHEA_FREE (segment_size);

  RHEA_GLOBAL_INFOF ("Done %s (%s)\n", __func__, file_path_txt);

  /* stop performance monitors */
  ymir_perf_counter_stop_add (
      &rhea_io_mpi_perfmon[RHEA_IO_MPI_PERFMON_GATHER_WRITE]);
}
