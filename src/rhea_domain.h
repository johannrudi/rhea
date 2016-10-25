/*
 */

#ifndef RHEA_DOMAIN_H
#define RHEA_DOMAIN_H

#include <p4est_to_p8est.h>
#include <p8est.h>
#include <ymir_options.h>
#include <ymir_mesh.h>
#include <ymir_velocity_dirichlet.h>

/* enumerator for domain shapes */
typedef enum
{
  RHEA_DOMAIN_CUBE,
  RHEA_DOMAIN_BOX,
  RHEA_DOMAIN_SHELL,
  RHEA_DOMAIN_CUBE_SPHERICAL,
  RHEA_DOMAIN_BOX_SPHERICAL
}
rhea_domain_shape_t;

/* enumerator for velocity boundray conditions */
typedef enum
{
  /* Dirichlet all */
  RHEA_DOMAIN_VELOCITY_BC_DIRICHLET_ALL = 0,

  /* Dirchlet in normal direction */
  RHEA_DOMAIN_VELOCITY_BC_DIRICHLET_NORM = 1,

  /* Dirchlet norm (if shell: fix two points on inner sphere for rotation
   * invariance) */
  RHEA_DOMAIN_VELOCITY_BC_DIRICHLET_NORM_FIXDOF = 2,

  /* Dirichlet all on inner sphere and Dirichlet norm on outer sphere of a shell
   * (for shell domain only) */
  RHEA_DOMAIN_VELOCITY_BC_DIRICHLET_ALL_INNER = 3,

  /* Dirichlet all on side faces and Neumann on top and bottom faces
   * (for rectangular domains only) */
  RHEA_DOMAIN_VELOCITY_BC_DIR_SIDES_NEU_TB = 4,

  /* Dirichlet norm on side faces and Dirichlet all on top and bottom faces
   * (for rectangular domains only) */
  RHEA_DOMAIN_VELOCITY_BC_DIR_NORM_SIDES_DIR_ALL_TB = 5,

  /* Neumann */
  RHEA_DOMAIN_VELOCITY_BC_NEUMANN_ALL = 6
}
rhea_domain_velocity_bc_t;

/* boundary information */
typedef struct rhea_domain_boundary
{
  ymir_mesh_e_to_fm_t e_to_fm_fn;
  ymir_topidx_t      *tree_to_bf; /* data for `e_to_fm_fn` */
}
rhea_domain_boundary_t;

/* options & properties of a computational domain */
typedef struct rhea_domain_options
{
  /* shape of the domain */
  rhea_domain_shape_t shape;

  /* extension in each Cartesian direction for `box` domain */
  double              box_x_extension;
  double              box_y_extension;
  double              box_z_extension;

  /* the domain knows the location of the lower-upper mantle interface,
   * which causes discontinuous material properties */
  double              lm_um_interface_radius;
  double              lm_um_interface_smooth_transition_width;

  /* velocity boundary conditions */
  rhea_domain_velocity_bc_t  velocity_bc_type;

  /* properties of the domain (reference coordinate system) */
  double              x_min;
  double              x_max;
  double              y_min;
  double              y_max;
  double              z_min;
  double              z_max;
  double              lon_min;
  double              lon_max;
  double              radius_min;
  double              radius_max;
  double              volume;
  double              center[3];
  double              moment_of_inertia[3];
}
rhea_domain_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_domain_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_domain_process_options (rhea_domain_options_t *opt);

/**
 * Computes the radius of a shell domain or the corresponding value for a
 * rectangular domain.
 */
double              rhea_domain_compute_radius (const double x, const double y,
                                                const double z,
                                                rhea_domain_options_t *opt);

/**
 * Returns whether the element's center is located in the upper mantle.
 */
int                 rhea_domain_elem_is_in_upper_mantle (
                                                  const double *x,
                                                  const double *y,
                                                  const double *z,
                                                  const int *Vmask,
                                                  rhea_domain_options_t *opt);

/**
 * Creates a new boundary object.
 */
rhea_domain_boundary_t *rhea_domain_boundary_new (p4est_t *p4est,
                                                  rhea_domain_options_t *opt);

/**
 * Destroys a boundary object.
 */
void                rhea_domain_boundary_destroy (
                                            rhea_domain_boundary_t *boundary);

/**
 * Creates Dirichlet boundary conditions.
 */
ymir_vel_dir_t     *rhea_domain_create_velocity_dirichlet_bc (
                                                        ymir_mesh_t *ymir_mesh,
                                                        ymir_vec_t *dirscal,
                                                        void *data);

#endif /* RHEA_DOMAIN_H */
