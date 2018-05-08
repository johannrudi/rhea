/*
 */

#include <rhea_velocity_pressure.h>
#include <rhea_base.h>
#include <rhea_velocity.h>
#include <rhea_pressure.h>
#include <rhea_io_mpi.h>
#include <ymir_pressure_vec.h>
#include <ymir_stokes_vec.h>
#include <ymir_comm.h>

ymir_vec_t *
rhea_velocity_pressure_new (ymir_mesh_t *ymir_mesh,
                            ymir_pressure_elem_t *press_elem)
{
  return ymir_stokes_vec_new (ymir_mesh, press_elem);
}

void
rhea_velocity_pressure_destroy (ymir_vec_t *velocity_pressure)
{
  ymir_vec_destroy (velocity_pressure);
}

int
rhea_velocity_pressure_check_vec_type (ymir_vec_t *vec)
{
  return ymir_stokes_vec_is_stokes_vec (vec);
}

int
rhea_velocity_pressure_is_valid (ymir_vec_t *vec)
{
  return sc_dmatrix_is_valid (vec->dataown) && sc_dmatrix_is_valid (vec->coff);
}

int
rhea_velocity_pressure_create_components (ymir_vec_t **vel, ymir_vec_t **press,
                                          ymir_vec_t *vel_press,
                                          ymir_pressure_elem_t *press_elem)
{
  ymir_mesh_t        *mesh = ymir_vec_get_mesh (vel_press);
  int                 is_view;

  /* check input */
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (vel_press));

  if (press_elem->space == YMIR_PRESSURE_SPACE_STAB) { /* if copy data */
    if (vel != NULL && press != NULL) {
      *vel = ymir_cvec_new (mesh, 3);
      *press = ymir_pressure_vec_new (mesh, press_elem);
      ymir_stokes_vec_get_components (vel_press, *vel, *press, press_elem);
    }
    else if (vel != NULL) {
      *vel = ymir_cvec_new (mesh, 3);
      ymir_stokes_vec_get_velocity (vel_press, *vel, press_elem);
    }
    else if (press != NULL) {
      *press = ymir_pressure_vec_new (mesh, press_elem);
      ymir_stokes_vec_get_pressure (vel_press, *press, press_elem);
    }
    is_view = 0;
  }
  else { /* otherwise create view onto data */
    RHEA_ASSERT (press_elem->space == YMIR_PRESSURE_SPACE_POLY ||
                 press_elem->space == YMIR_PRESSURE_SPACE_TENS);
    if (vel != NULL) {
      *vel = ymir_cvec_new_data (mesh, vel_press->ncfields, vel_press->cvec);
    }
    if (press != NULL) {
      *press = ymir_evec_new_data (mesh, vel_press->nefields,
                                   vel_press->e_to_d_fn,
                                   vel_press->e_to_d_data,
                                   vel_press->evec);
    }
    is_view = 1;
  }

  return is_view;
}

void
rhea_velocity_pressure_copy_components (ymir_vec_t *vel, ymir_vec_t *press,
                                        ymir_vec_t *vel_press,
                                        ymir_pressure_elem_t *press_elem)
{
  if (vel != NULL && press != NULL) {
    ymir_stokes_vec_get_components (vel_press, vel, press, press_elem);
  }
  else if (vel != NULL) {
    ymir_stokes_vec_get_velocity (vel_press, vel, press_elem);
  }
  else if (press != NULL) {
    ymir_stokes_vec_get_pressure (vel_press, press, press_elem);
  }
}

void
rhea_velocity_pressure_set_components (ymir_vec_t *vel_press,
                                       ymir_vec_t *vel, ymir_vec_t *press,
                                       ymir_pressure_elem_t *press_elem)
{
  if (vel != NULL && press != NULL) {
    ymir_stokes_vec_set_components (vel, press, vel_press, press_elem);
  }
  else if (vel != NULL) {
    ymir_stokes_vec_set_velocity (vel, vel_press, press_elem);
  }
  else if (press != NULL) {
    ymir_stokes_vec_set_pressure (press, vel_press, press_elem);
  }
}

int
rhea_velocity_pressure_read (ymir_vec_t *vel_press,
                             char *vel_file_path_bin,
                             char *press_file_path_bin,
                             ymir_pressure_elem_t *press_elem,
                             sc_MPI_Comm mpicomm)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vel_press);
  ymir_vec_t         *velocity, *pressure;
  MPI_Offset          segment_offset;
  int                 segment_size;
  int                 n_read;
  double             *data;
  int                 success = 0;

  /* check input */
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (vel_press));

  /* create velocity & pressure */
  velocity = rhea_velocity_new (ymir_mesh);
  pressure = rhea_pressure_new (ymir_mesh, press_elem);

  /* read velocity */
  {
    data = ymir_cvec_index (velocity, 0, 0);
    segment_offset = rhea_velocity_segment_offset_get (velocity);
    segment_size = rhea_velocity_segment_size_get (velocity);

    n_read = rhea_io_mpi_read_segment_double (
          data, segment_offset, segment_size, vel_file_path_bin, mpicomm);
    RHEA_ASSERT (n_read == segment_size);
    success += (n_read == segment_size);
  }

  /* read pressure */
  {
    if (press_elem->space == YMIR_PRESSURE_SPACE_STAB) {
      RHEA_ABORT_NOT_REACHED ();
      //TODO
    }
    else {
      RHEA_ASSERT (press_elem->space == YMIR_PRESSURE_SPACE_POLY ||
                   press_elem->space == YMIR_PRESSURE_SPACE_TENS);
      data = ymir_evec_index (pressure, 0, 0);
      segment_offset = rhea_pressure_segment_offset_get (pressure);
      segment_size = rhea_pressure_segment_size_get (pressure);
    }

    n_read = rhea_io_mpi_read_segment_double (
          data, segment_offset, segment_size, vel_file_path_bin, mpicomm);
    RHEA_ASSERT (n_read == segment_size);
    success += (n_read == segment_size);
  }

  /* set velocity & pressure */
  rhea_velocity_pressure_set_components (vel_press, velocity, pressure,
                                         press_elem);

  /* destroy */
  rhea_velocity_destroy (velocity);
  rhea_pressure_destroy (pressure);

  /* communicate shared node values */
  ymir_vec_share_owned (vel_press);
  RHEA_ASSERT (rhea_velocity_pressure_is_valid (vel_press));

  return success;
}

int
rhea_velocity_pressure_write (char *vel_file_path_bin,
                              char *press_file_path_bin,
                              ymir_vec_t *vel_press,
                              ymir_pressure_elem_t *press_elem,
                              sc_MPI_Comm mpicomm)
{
  ymir_vec_t         *velocity, *pressure;
  MPI_Offset          segment_offset;
  int                 segment_size;
  int                 n_written;
  double             *data;
  int                 success = 0;

  /* check input */
  RHEA_ASSERT (rhea_velocity_pressure_check_vec_type (vel_press));
  RHEA_ASSERT (sc_dmatrix_is_valid (vel_press->dataown));

  /* get velocity & pressure */
  rhea_velocity_pressure_create_components (&velocity, &pressure, vel_press,
                                            press_elem);

  /* write velocity */
  {
    data = ymir_cvec_index (velocity, 0, 0);
    segment_offset = rhea_velocity_segment_offset_get (velocity);
    segment_size = rhea_velocity_segment_size_get (velocity);

    n_written = rhea_io_mpi_write_segment_double (
          vel_file_path_bin, data, segment_offset, segment_size, mpicomm);
    RHEA_ASSERT (n_written == segment_size);

    success += (n_written == segment_size);
  }

  /* write pressure */
  {
    if (press_elem->space == YMIR_PRESSURE_SPACE_STAB) {
      RHEA_ABORT_NOT_REACHED ();
      //TODO
    }
    else {
      RHEA_ASSERT (press_elem->space == YMIR_PRESSURE_SPACE_POLY ||
                   press_elem->space == YMIR_PRESSURE_SPACE_TENS);
      data = ymir_evec_index (pressure, 0, 0);
      segment_offset = rhea_pressure_segment_offset_get (pressure);
      segment_size = rhea_pressure_segment_size_get (pressure);
    }

    n_written = rhea_io_mpi_write_segment_double (
          press_file_path_bin, data, segment_offset, segment_size, mpicomm);
    RHEA_ASSERT (n_written == segment_size);

    success += (n_written == segment_size);
  }

  /* destroy */
  rhea_velocity_destroy (velocity);
  rhea_pressure_destroy (pressure);

  return success;
}
