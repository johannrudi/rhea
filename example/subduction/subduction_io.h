#ifndef SUBDUCTION_IO_H
#define SUBDUCTION_IO_H

#include <subduction_options.h>
#include <subduction_physics.h>
#include <subduction_postp.h>
#include <subduction_TI.h>
#include <rhea_domain.h>
#include <rhea_io_mpi.h>
#include <rhea_io_std.h>
#include <rhea_base.h>
#include <ymir_vec_ops.h>
#include <ymir_vec_getset.h>
#include <ymir_stress_op.h>
#include <ymir_velocity_vec.h>
#include <ymir_comm.h>
#include <ymir_stokes_vec.h>
#include <ymir_derivative_elem.h>
#include <ymir_stokes_pc.h>
#include <ymir_interp_vec.h>
#include <ymir_mass_vec.h>

void
subduction_add_txt_options (ymir_options_t * opt);

void
subd_visual_mesh_write (ymir_mesh_t * mesh, const char *file_path_txt);

void
subd_gauss_coord_pressure_write (ymir_vec_t *pressure, const char *file_path_txt);

void
subd_coord_topo_write (ymir_vec_t *topography, const char *file_path_txt);

/*write out volume ymir_vec in txt*/
void
subd_vol_cvec_write (ymir_vec_t *volvec, const char *file_path_txt);

void
subd_vol_dvec_write (ymir_vec_t *volvec, const char *file_path_txt);

/*write out face mesh ymir_vec in txt. This is more general, can replace subd_volvec_write*/
void
subd_facevec_write (ymir_vec_t *facevec, const char *file_path_txt);

void
subd_vol_dvec_read (ymir_vec_t *facevec, const char *file_path_txt);

/* read face mesh ymir_vec data from txt.*/
void
subd_facevec_read (ymir_vec_t *facevec, const char *file_path_txt);

void
subd_txt_io (rhea_stokes_problem_t *stokes_problem, subd_options_t *subd_options);

#endif
