/*
 */

#include <rhea_io_std.h>
#include <rhea_base.h>

/**
 * Opens a file.
 */
static FILE *
rhea_io_std_file_open (const char *file_path, const char *access_type,
                       const char *this_fn_name)
{
  FILE               *file_ptr;

  /* open file */
  file_ptr = fopen (file_path, access_type);

  /* check status */
  if (!file_ptr) {
    if (this_fn_name != NULL) {
      RHEA_LERRORF ("%s: Could not open file \"%s\"\n",
                    this_fn_name, file_path);
    }
    else {
      RHEA_LERRORF ("Could not open file \"%s\"\n", file_path);
    }
  }

  /* return file pointer */
  return file_ptr;
}

/**
 * Closes a file.
 */
static void
rhea_io_std_file_close (FILE *file_ptr, const char *file_path,
                        const char *this_fn_name)
{
  int                 status;

  /* close file */
  status = fclose (file_ptr);

  /* check status */
  if (0 != status) {
    if (this_fn_name != NULL && file_path != NULL) {
      RHEA_LERRORF ("%s: Could not close file \"%s\"\n",
                    this_fn_name, file_path);
    }
    else if (file_path != NULL) {
      RHEA_LERRORF ("Could not close file \"%s\"\n", file_path);
    }
    else {
      RHEA_LERROR ("Could not close file\n");
    }
  }
}

/******************************************************************************
 * Binary Files
 *****************************************************************************/

size_t
rhea_io_std_read_double (double *values, const size_t n_entries,
                         const char *file_path)
{
  const char         *this_fn_name = "rhea_io_std_read_double";
  FILE               *file_ptr;
  size_t              n_read;

  RHEA_INFOF ("Into %s (%s, #entries %i)\n", this_fn_name, file_path,
              (int) n_entries);

  /* exit if nothing to do */
  if (0 == n_entries) {
    return 0;
  }

  /* open binary file */
  file_ptr = rhea_io_std_file_open (file_path, "rb", this_fn_name);
  if (!file_ptr) {
    return 0;
  }

  /* read values */
  n_read = fread (values, sizeof (double), n_entries, file_ptr);

  /* check #entries read vs. requested */
  if (n_read != n_entries) {
    RHEA_LERRORF (
        "%s: Mismatch of #entries read: "
        "Requested #entries %lli, actually read %lli, file path %s\n",
        this_fn_name, (long long int) n_entries, (long long int) n_read,
        file_path);
  }

  /* close file */
  rhea_io_std_file_close (file_ptr, file_path, this_fn_name);

  RHEA_INFOF ("Done %s (%s)\n", this_fn_name, file_path);

  /* return number of read entries */
  return n_read;
}

size_t
rhea_io_std_write_double (const char *file_path, const double *values,
                          const size_t n_entries)
{
  const char         *this_fn_name = "rhea_io_std_write_double";
  FILE               *file_ptr;
  size_t              n_written;

  RHEA_INFOF ("Into %s (%s, #entries %i)\n", this_fn_name, file_path,
              (int) n_entries);

  /* exit if nothing to do */
  if (0 == n_entries) {
    return 0;
  }

  /* open binary file */
  file_ptr = rhea_io_std_file_open (file_path, "wb", this_fn_name);
  if (!file_ptr) {
    return 0;
  }

  /* write values */
  n_written = fwrite (values, sizeof (double), n_entries, file_ptr);

  /* check #entries written vs. requested */
  if (n_written != n_entries) {
    RHEA_LERRORF (
        "%s: Mismatch of #entries written: "
        "Requested #entries %lli, actually written %lli, file path %s\n",
        this_fn_name, (long long int) n_entries, (long long int) n_written,
        file_path);
  }

  /* close file */
  rhea_io_std_file_close (file_ptr, file_path, this_fn_name);

  RHEA_INFOF ("Done %s (%s)\n", this_fn_name, file_path);

  /* return number of written entries */
  return n_written;
}

/******************************************************************************
 * Text Files
 *****************************************************************************/

size_t
rhea_io_std_read_txt (void *data, rhea_io_std_fscanf_fn_t fscanf_fn,
                      void *fscanf_params, const char *file_path)
{
  const char         *this_fn_name = "rhea_io_std_read_txt";
  FILE               *file_ptr;
  size_t              n_read;
  int                 status;

  RHEA_INFOF ("Into %s (%s)\n", this_fn_name, file_path);

  /* open text file */
  file_ptr = rhea_io_std_file_open (file_path, "r", this_fn_name);
  if (!file_ptr) {
    return 0;
  }

  /* read data */
  n_read = 0;
  while (!feof (file_ptr)) { /* while not at the end of file */
    status = fscanf_fn (data, fscanf_params, file_ptr, n_read);
    if (0 < status) {
      n_read += (size_t) status;
    }
  }

  /* close file */
  rhea_io_std_file_close (file_ptr, file_path, this_fn_name);

  RHEA_INFOF ("Done %s (%s)\n", this_fn_name, file_path);

  /* return number of read items */
  return n_read;
}

size_t
rhea_io_std_write_txt (const char *file_path, const void *data,
                       rhea_io_std_fprintf_fn_t fprintf_fn,
                       void *fprintf_params)
{
  const char         *this_fn_name = "rhea_io_std_write_txt";
  FILE               *file_ptr;
  size_t              n_written;
  int                 status;

  RHEA_INFOF ("Into %s (%s)\n", this_fn_name, file_path);

  /* open binary file */
  file_ptr = rhea_io_std_file_open (file_path, "w", this_fn_name);
  if (!file_ptr) {
    return 0;
  }

  /* write data */
  n_written = 0;
  status = 1;
  while (0 < status) { /* while items are written */
    status = fprintf_fn (file_ptr, data, fprintf_params, n_written);
    if (0 < status) {
      n_written += (size_t) status;
    }
  }

  /* close file */
  rhea_io_std_file_close (file_ptr, file_path, this_fn_name);

  RHEA_INFOF ("Done %s (%s)\n", this_fn_name, file_path);

  /* return number of written items */
  return n_written;
}

static int
rhea_io_std_fscanf_double_fn (void *data, void *params, FILE *file_ptr,
                              const size_t n_read)
{
  double             *values = data;
  const size_t        n_entries = *((size_t *) params);
  int                 status;

  if (0 < n_entries) {
    RHEA_ASSERT (n_read <= n_entries);
    status = fscanf (file_ptr, "%lf", &values[SC_MIN (n_read, n_entries-1)]);
  }
  else {
    status = fscanf (file_ptr, "%lf", &values[n_read]);
  }

  return status;
}

size_t
rhea_io_std_read_double_from_txt (double *values,
                                  size_t n_entries,
                                  const char *file_path)
{
  const char         *this_fn_name = "rhea_io_std_read_double_from_txt";
  size_t              n_read;

  /* read double values */
  n_read = rhea_io_std_read_txt (values, rhea_io_std_fscanf_double_fn,
                                 &n_entries, file_path);

  /* check #entries requested vs. read */
  if (0 < n_entries && n_entries != n_read) {
    RHEA_INFOF (
        "%s: Mismatch of #entries read: "
        "Requested #entries %lli, actually read %lli, file path %s\n",
        this_fn_name, (long long int) n_entries, (long long int) n_read,
        file_path);
    RHEA_ASSERT (n_read < n_entries);
  }

  /* return number of read entries */
  return n_read;
}

static int
rhea_io_std_fprintf_double_fn (FILE *file_ptr, const void *data, void *params,
                               const size_t n_written)
{
  const double       *values = data;
  int                *p = params;
  const int           n_entries = p[0];
  const int           n_entries_per_line = p[1];
  int                 status;

  if (n_written >= n_entries) {
    return 0;
  }

  status = fprintf (file_ptr, "%+.16e ", values[n_written]);

  if (0 < status && 0 < n_entries_per_line &&
      0 == (n_written + 1) % n_entries_per_line){
    fprintf (file_ptr, "\n");
  }

  return (0 < status);
}

size_t
rhea_io_std_write_double_to_txt (const char *file_path,
                                 const double *values,
                                 size_t n_entries,
                                 int n_entries_per_line)
{
  const char         *this_fn_name = "rhea_io_std_write_double_to_txt";
  int                 params[2] = {(int) n_entries, n_entries_per_line};
  size_t              n_written;

  /* check input */
  RHEA_ASSERT (0 < n_entries);

  /* write double values */
  n_written = rhea_io_std_write_txt (file_path, values,
                                     rhea_io_std_fprintf_double_fn, params);

  /* check #entries written vs. requested */
  if (n_entries != n_written) {
    RHEA_LERRORF (
        "%s: Mismatch of #entries written: "
        "Requested #entries %lli, actually written %lli, file path %s\n",
        this_fn_name, (long long int) n_entries, (long long int) n_written,
        file_path);
  }

  /* return number of written entries */
  return n_written;
}
