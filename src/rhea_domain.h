/* RHEA_DOMAIN  Physical domain properties, independent of discretization. */

#ifndef RHEA_DOMAIN_H
#define RHEA_DOMAIN_H

#include <p4est_to_p8est.h>
#include <p8est.h>
#include <ymir_options.h>
#include <ymir_mesh.h>
#include <ymir_velocity_dirichlet.h>

/******************************************************************************
 * Options
 *****************************************************************************/

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
  /* user provided function */
  RHEA_DOMAIN_VELOCITY_BC_USER = -1,

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

  /* Dirichlet norm on side faces and bottom face; Neumann on top face
   * (for rectangular domains only) */
  RHEA_DOMAIN_VELOCITY_BC_DIR_NORM_SIDES_B_NEU_T = 6,

  /* Neumann */
  RHEA_DOMAIN_VELOCITY_BC_NEUMANN_ALL = 7
}
rhea_domain_velocity_bc_t;

/* enumerator for boundary faces */
typedef enum
{
  RHEA_DOMAIN_BOUNDARY_FACE_NONE = -1,
  RHEA_DOMAIN_BOUNDARY_FACE_BASE = 0,
  RHEA_DOMAIN_BOUNDARY_FACE_TOP,
  RHEA_DOMAIN_BOUNDARY_FACE_SIDE1,
  RHEA_DOMAIN_BOUNDARY_FACE_SIDE2,
  RHEA_DOMAIN_BOUNDARY_FACE_SIDE3,
  RHEA_DOMAIN_BOUNDARY_FACE_SIDE4
}
rhea_domain_boundary_face_t;

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
  /* shape of the domain (input) */
  rhea_domain_shape_t shape;

  /* length in each Cartesian direction for `box` domain (input) */
  double              box_length_x;
  double              box_length_y;
  double              box_length_z;

  /* subdivision of box into cubes along each Cartesian direction (input) */
  int                 box_subdivision_x;
  int                 box_subdivision_y;
  int                 box_subdivision_z;

  /* correct a pillow-type distortion of domain shape `box_spherical` */
  int                 box_spherical_distortion_corr;

  /* the domain knows the location of the lower-upper mantle interface,
   * which causes discontinuous material properties (input) */
  double              lm_um_interface_radius;
  double              lm_um_interface_smoothing_width;

  /* velocity boundary conditions (input) */
  rhea_domain_velocity_bc_t  velocity_bc_type;

  /* nondimensional properties of the domain (computed) */
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
  double              depth;
  double              volume;
  double              center[3];
  double              moment_of_inertia[3];
  double              moment_of_inertia_surface[3];

  /* dimensional properties of the domain (input) */
  double              radius_min_m;
  double              radius_max_m;
  double              density_kg_m3;
  double              gravity_m_s2;
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
 * Prints domain options that are assumed to remain constant.
 */
void                rhea_domain_print_const_options (
                                                  rhea_domain_options_t *opt);

/******************************************************************************
 * Basic Calculations
 *****************************************************************************/

/* types of coordinates */
typedef enum
{
  RHEA_DOMAIN_COORDINATE_CARTESIAN,         /* Cartesian non-dimensional
                                             * (x,y,z) coordinates
                                             */
  RHEA_DOMAIN_COORDINATE_SPHERICAL_MATH,    /* Spherical non-dimensional
                                             * (r,theta,phi) coordinates with
                                             * - radius:  0<=r<=1
                                             * - azimuthal angle: -pi<theta<=pi
                                             * - polar angle:        0<=phi<=pi
                                             */
  RHEA_DOMAIN_COORDINATE_SPHERICAL_GEO,     /* Spherical non-dimensional
                                             * (r,phi,theta) coordinates with
                                             * - radius:  0<=r<=1
                                             * - azimuthal angle:   -pi<phi<=pi
                                             * - polar angle:      0<=theta<=pi
                                             */
  RHEA_DOMAIN_COORDINATE_SPHERICAL_GEO_DIM  /* Spherical dimensional
                                             * (r,phi,theta) coordinates with
                                             * - radius: 0<=r<=6371.0e3
                                             * - east longitude: 0<phi<=360
                                             * - colatitude:  0<=theta<=180
                                             */
}
rhea_domain_coordinate_type_t;

/**
 * Converts nondimensional Cartesian coordinates (x,y,z) into coordinates
 * (coord1,coord2,coord3) of the type `coord_type`.
 */
void                rhea_domain_convert_coordinates (
                                double *coord1, double *coord2, double *coord3,
                                const double x, const double y, const double z,
                                rhea_domain_coordinate_type_t coord_type,
                                rhea_domain_options_t *opt);

/**
 * Converts dimensional depth [m] to nondimensional depth.
 */
double              rhea_domain_depth_m_to_depth (const double depth_m,
                                                  rhea_domain_options_t *opt);

/**
 * Converts dimensional depth [m] to nondimensional radius.
 */
double              rhea_domain_depth_m_to_radius (const double depth_m,
                                                   rhea_domain_options_t *opt);

/**
 * Converts nondimensional radius to dimensional radius [m].
 */
double              rhea_domain_radius_to_radius_m (const double radius,
                                                    rhea_domain_options_t *opt);

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
 * Maps coordinates onto the surface parallel to the radial direction.
 */
void                rhea_domain_project_to_surface (double *x, double *y,
                                                    double *z,
                                                    rhea_domain_options_t *opt);

/**
 * Extracts lateral 2-dimensional coordinates of type `coord_type`.
 */
void                rhea_domain_extract_lateral (
                                double *coord1, double *coord2,
                                const double x, const double y, const double z,
                                rhea_domain_coordinate_type_t coord_type,
                                rhea_domain_options_t *opt);

/**
 * Rotates Cartesian coordinates about x-, y-, or z-axis by a given angle.
 */
void                rhea_domain_rotate_x_axis (double c[3],
                                               const double angle);

void                rhea_domain_rotate_y_axis (double c[3],
                                               const double angle);

void                rhea_domain_rotate_z_axis (double c[3],
                                               const double angle);

/******************************************************************************
 * Domain Boundary
 *****************************************************************************/

/**
 * Returns the number of boundary faces.
 */
int                 rhea_domain_boundary_get_num (rhea_domain_options_t *opt);

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

/******************************************************************************
 * Boundary Conditions
 *****************************************************************************/

/**
 * Sets function that defines velocity Dirichlet BC's.
 */
void                rhea_domain_set_user_velocity_dirichlet_bc (
                                            ymir_vel_dir_fn_t vel_dir_bc_fn,
                                            void *vel_dir_bc_data,
                                            const int vel_dir_bc_nonzero);

/**
 * Creates Dirichlet boundary conditions.
 */
ymir_vel_dir_t     *rhea_domain_create_velocity_dirichlet_bc (
                                            ymir_mesh_t *ymir_mesh,
                                            ymir_vec_t *dirscal,
                                            void *data);

#endif /* RHEA_DOMAIN_H */
