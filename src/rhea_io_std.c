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
  if (!status) {
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

size_t
rhea_io_std_read_double (double *values, const size_t n_entries,
                         const char *file_path)
{
  const char         *this_fn_name = "rhea_io_std_read_double";
  FILE               *file_ptr;
  size_t              n_read;

  RHEA_INFOF ("Into %s (%s)\n", this_fn_name, file_path);

  /* open binary file */
  file_ptr = rhea_io_std_file_open (file_path, "rb", this_fn_name);

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

  RHEA_INFOF ("Into %s (%s)\n", this_fn_name, file_path);

  /* open binary file */
  file_ptr = rhea_io_std_file_open (file_path, "wb", this_fn_name);

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
