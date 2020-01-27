/*
 */

#include <rhea_composition.h>
#include <rhea_base.h>
#include <rhea_io_mpi.h>
#include <ymir_vec_getset.h>
#include <ymir_comm.h>

/******************************************************************************
 * Options
 *****************************************************************************/

/* default options */
#define RHEA_COMPOSITION_DEFAULT_TYPE_NAME "NONE"
#define RHEA_COMPOSITION_DEFAULT_DATA_FILE_PATH_BIN NULL
#define RHEA_COMPOSITION_DEFAULT_DATA_FILE_PATH_TXT NULL
#define RHEA_COMPOSITION_DEFAULT_WRITE_DATA_FILE_PATH_BIN NULL
#define RHEA_COMPOSITION_DEFAULT_RHS_SCALING (1.0)
#define RHEA_COMPOSITION_DEFAULT_FLAVOR_NUMBER (0)
#define RHEA_COMPOSITION_DEFAULT_FLAVOR_DENSITY_ANOMALY NULL
#define RHEA_COMPOSITION_DEFAULT_FLAVOR_VISC_EFFECT NULL

/* initialize options */
char               *rhea_composition_type_name =
  RHEA_COMPOSITION_DEFAULT_TYPE_NAME;
char               *rhea_composition_data_file_path_bin =
  RHEA_COMPOSITION_DEFAULT_DATA_FILE_PATH_BIN;
char               *rhea_composition_data_file_path_txt =
  RHEA_COMPOSITION_DEFAULT_DATA_FILE_PATH_TXT;
char               *rhea_composition_write_data_file_path_bin =
  RHEA_COMPOSITION_DEFAULT_WRITE_DATA_FILE_PATH_BIN;
double              rhea_composition_rhs_scaling =
  RHEA_COMPOSITION_DEFAULT_RHS_SCALING;
int				   rhea_composition_flavor_number =
  RHEA_COMPOSITION_DEFAULT_FLAVOR_NUMBER;
char			   *rhea_composition_flavor_density_anomaly =
  RHEA_COMPOSITION_DEFAULT_FLAVOR_DENSITY_ANOMALY;
char			   *rhea_composition_flavor_visc_effect =
  RHEA_COMPOSITION_DEFAULT_FLAVOR_VISC_EFFECT;

void
rhea_composition_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Composition";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

    YMIR_OPTIONS_S, "type", '\0',
	&(rhea_composition_type_name), RHEA_COMPOSITION_DEFAULT_TYPE_NAME,
	"Type of composition: NONE, flavor, density",

	YMIR_OPTIONS_I, "number-of-flavor", '\0',
	&(rhea_composition_flavor_number), RHEA_COMPOSITION_DEFAULT_FLAVOR_NUMBER,
	"number of flavors for composition when type is flavor",

	YMIR_OPTIONS_S, "density-anomaly-per-flavor", '\0',
	&(rhea_composition_flavor_density_anomaly),
	RHEA_COMPOSITION_DEFAULT_FLAVOR_DENSITY_ANOMALY,
	"density anomaly for each flavor (comma separated list of doubles)",

	YMIR_OPTIONS_S, "visc-effect-per-flavor", '\0',
	&(rhea_composition_flavor_visc_effect),
	RHEA_COMPOSITION_DEFAULT_FLAVOR_VISC_EFFECT,
	"viscosity effect for each flavor (comma separated list of doubles)",

	YMIR_OPTIONS_D, "composition-right-hand-side-scaling", '\0',
	&(rhea_composition_rhs_scaling), RHEA_COMPOSITION_DEFAULT_RHS_SCALING,
	"Scaling factor for velocity right-hand side from composition",

	YMIR_OPTIONS_S, "data-file-path-bin", '\0',
	&(rhea_composition_data_file_path_bin),
	RHEA_COMPOSITION_DEFAULT_DATA_FILE_PATH_BIN,
	"Path to a binary file that contains a composition field",
	YMIR_OPTIONS_S, "data-file-path-txt", '\0',
	&(rhea_composition_data_file_path_txt),
	RHEA_COMPOSITION_DEFAULT_DATA_FILE_PATH_TXT,
	"Path to a text file that contains a composition field",
	YMIR_OPTIONS_S, "write-data-file-path-bin", '\0',
	&(rhea_composition_write_data_file_path_bin),
	RHEA_COMPOSITION_DEFAULT_WRITE_DATA_FILE_PATH_BIN,
	"Output path for a binary file containing a composition field",

  YMIR_OPTIONS_END_OF_LIST);

  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);

}

void
rhea_composition_process_options (rhea_composition_options_t *opt,
                                  rhea_domain_options_t *domain_options)
{

  /* set composition type */
  if (strcmp (rhea_composition_type_name, "NONE") == 0) {
    opt->type = RHEA_COMPOSITION_NONE;
  }
  else if (strcmp (rhea_composition_type_name, "flavor") == 0) {
    opt->type = RHEA_COMPOSITION_FLAVOR;
    opt->n_flavor = rhea_composition_flavor_number;

    /* set density anomaly and viscosity effect of flavors */
    int                 flavors, k;
	double             *n_density_labels = NULL;
	double			   *n_visc_labels = NULL;

	flavors = ymir_options_convert_string_to_double (
		rhea_composition_flavor_density_anomaly, &n_density_labels);
	RHEA_CHECK_ABORT (rhea_composition_flavor_number == flavors,
					  "Mismatch with provided number of flavors");
	flavors = ymir_options_convert_string_to_double (
		rhea_composition_flavor_visc_effect, &n_visc_labels);
	RHEA_CHECK_ABORT (rhea_composition_flavor_number == flavors,
					  "Mismatch with provided number of flavors");
	for (k = 0; k < flavors; k++) {
	  opt->density_flavor[k] = (double) n_density_labels[k];
	  opt->visc_flavor[k] = (double) n_visc_labels[k];
	}
	YMIR_FREE (n_density_labels); /* was allocated in ymir */
	YMIR_FREE (n_visc_labels); /* was allocated in ymir */
  }
  else if (strcmp (rhea_composition_type_name, "density") == 0) {
    opt->type = RHEA_COMPOSITION_DENSITY;
  }
  else { /* unknown composition type */
    RHEA_ABORT ("Unknown composition type");
  }

  /* set right-hand side options */
  opt->rhs_scaling_comp = rhea_composition_rhs_scaling;

  /* store domain options */
  opt->domain_options = domain_options;
}

ymir_vec_t *
rhea_composition_new (ymir_mesh_t *ymir_mesh)
{
  return ymir_cvec_new (ymir_mesh, 1);
}

void
rhea_composition_destroy (ymir_vec_t *composition)
{
  ymir_vec_destroy (composition);
}

int
rhea_composition_check_vec_type (ymir_vec_t *vec)
{
  return (
      ymir_vec_is_cvec (vec) &&
      vec->ncfields == 1 &&
      vec->node_type == YMIR_GLL_NODE
  );
}

MPI_Offset *
rhea_composition_segment_offset_create (ymir_vec_t *vec)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
  const mangll_locidx_t *n_nodes = ymir_mesh->cnodes->Ngo;
  const int           n_fields = vec->ncfields;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpisize, mpiret;
  MPI_Offset         *segment_offset;
  int                 r;

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);

  /* create segment offsets */
  segment_offset = RHEA_ALLOC (MPI_Offset, mpisize + 1);
  segment_offset[0] = 0;
  for (r = 0; r < mpisize; r++) {
    segment_offset[r+1] = segment_offset[r] +
                          (MPI_Offset) (n_fields * n_nodes[r]);
  }

  return segment_offset;
}

int
rhea_composition_segment_size_get (ymir_vec_t *vec)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (vec);
  const mangll_locidx_t *n_nodes = ymir_mesh->cnodes->Ngo;
  const int           n_fields = vec->ncfields;
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpirank, mpiret;

  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);
  return (int) (n_fields * n_nodes[mpirank]);
}

void rhea_composition_read (ymir_vec_t *composition,
							rhea_composition_options_t *opt)
{
  switch (opt->type) {
  case RHEA_COMPOSITION_FLAVOR:
	  rhea_composition_read_internal (composition);
	  break;
  case RHEA_COMPOSITION_DENSITY:
	  rhea_composition_read_internal (composition);
	  break;
  case RHEA_COMPOSITION_NONE:
	  break;
  default: /* unknown composition type */
	  RHEA_ABORT_NOT_REACHED ();
  }
}

static void
rhea_composition_read_internal (ymir_vec_t *composition)
{
  ymir_mesh_t        *ymir_mesh = ymir_vec_get_mesh (composition);
  sc_MPI_Comm         mpicomm = ymir_mesh_get_MPI_Comm (ymir_mesh);
  int                 mpisize, mpirank, mpiret;
  MPI_Offset         *segment_offset;
  char               *file_path_bin;
  char               *file_path_txt;
  double             *comp_data = ymir_cvec_index (composition, 0, 0);

  /* set file paths */
  if (rhea_composition_data_file_path_bin != NULL) {
    file_path_bin = rhea_composition_data_file_path_bin;
    file_path_txt = NULL;
  }
  else if (rhea_composition_data_file_path_txt != NULL) {
    file_path_bin = rhea_composition_write_data_file_path_bin;
    file_path_txt = rhea_composition_data_file_path_txt;
  }
  else { /* unknown file path */
    RHEA_ABORT_NOT_REACHED ();
  }

  /* get parallel environment */
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize); SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank); SC_CHECK_MPI (mpiret);

  /* create segment offsets */
  segment_offset = rhea_composition_segment_offset_create (composition);

  /* read composition */
  if (file_path_txt != NULL || segment_offset[mpisize] < INT_MAX) {
    int                *segment_offset_int = RHEA_ALLOC (int, mpisize + 1);
    int                 r;

    for (r = 0; r <= mpisize; r++) {
      segment_offset_int[r] = (int) segment_offset[r];
    }
    rhea_io_mpi_read_scatter_double (comp_data, segment_offset_int,
                                     file_path_bin, file_path_txt, mpicomm);
    RHEA_FREE (segment_offset_int);
  }
  else {
    const int           segment_size = (int) ymir_mesh->cnodes->Ngo[mpirank];
    int                 n_read;

    RHEA_ASSERT (file_path_bin != NULL);
    n_read = rhea_io_mpi_read_segment_double (
        comp_data, segment_offset[mpirank], segment_size,
        file_path_bin, mpicomm);
    if (!(n_read == segment_size)) {
      RHEA_LERRORF ("%s: Mismatch of #entries read: "
                    "#requested %lli, #read %lli, file path %s\n",
                    __func__, (long long int) segment_size,
                    (long long int) n_read, file_path_bin);
    }
  }

  /* destroy */
  RHEA_FREE (segment_offset);

  /* communicate shared node values */
  ymir_vec_share_owned (composition);
}

int
rhea_composition_write (char *file_path_bin,
                        ymir_vec_t *composition,
                        sc_MPI_Comm mpicomm)
{
  MPI_Offset          segment_offset;
  int                 segment_size;
  int                 n_written;
  double             *data;

  /* check input */
  RHEA_ASSERT (rhea_composition_check_vec_type (composition));
  RHEA_ASSERT (sc_dmatrix_is_valid (composition->dataown));

  /* set parameters for I/O */
  data = ymir_cvec_index (composition, 0, 0);
  segment_offset = rhea_composition_segment_offset_get (composition);
  segment_size = rhea_composition_segment_size_get (composition);

  /* write */
  n_written = rhea_io_mpi_write_segment_double (
        file_path_bin, data, segment_offset, segment_size, mpicomm);
  RHEA_ASSERT (n_written == segment_size);

  /* return whether writing was success */
  return (n_written == segment_size);
}

/**
 * Computes velocity right-hand side at one node from composition.
 */
static void
rhea_composition_rhs_vel_node (double *rhs, const double x, const double y,
                               const double z, const double comp_density,
                               rhea_composition_options_t *opt)
{
  rhea_domain_options_t  *domain_options = opt->domain_options;
  const double        scaling = opt->rhs_scaling_comp;

  switch (domain_options->shape) {
  case RHEA_DOMAIN_CUBE:
  case RHEA_DOMAIN_BOX:
    /**
     * Computes right-hand side from composition
     * on the cube domain:
     *
     *   f(x) = e_z * C
     *
     * where `e_z` is unit vector in z-direction, `C` is compositional density.
     */
    rhs[0] = 0.0;
    rhs[1] = 0.0;
    rhs[2] = scaling * comp_density;
    break;
  case RHEA_DOMAIN_SHELL:
  case RHEA_DOMAIN_CUBE_SPHERICAL:
  case RHEA_DOMAIN_BOX_SPHERICAL:
    /**
     * Computes right-hand side from composition
     * on the shell domain:
     *
     *   f(x) = e_r * C
     *
     * where `e_r` is normalized spherical position vector, `C` is compositional density.
     */
    {
      const double        radius = rhea_domain_compute_radius (x, y, z,
                                                               domain_options);

      rhs[0] = scaling * x / radius * comp_density;
      rhs[1] = scaling * y / radius * comp_density;
      rhs[2] = scaling * z / radius * comp_density;
    }
    break;
  default: /* unknown domain shape */
    RHEA_ABORT_NOT_REACHED ();
  }
}

/* data for callback function to compute the velocity right-hand side */
typedef struct rhea_composition_compute_rhs_vel_fn_data
{
  ymir_vec_t          *comp_density;
  rhea_composition_options_t  *comp_options;
}
rhea_composition_compute_rhs_vel_fn_data_t;

/**
 * Callback function to compute the velocity right-hand side.
 */
static void
rhea_composition_compute_rhs_vel_fn (double *rhs, double x, double y, double z,
                                     ymir_locidx_t nodeid, void *data)
{
  rhea_composition_compute_rhs_vel_fn_data_t  *d = data;
  rhea_composition_options_t  *opt = d->comp_options;
  const double        comp_density = *ymir_cvec_index (d->comp_density, nodeid, 0);

  /* compute right-hand side */
  rhea_composition_rhs_vel_node (rhs, x, y, z, comp_density, opt);
}

void
rhea_composition_compute_rhs_vel (ymir_vec_t *rhs_vel,
                                  ymir_vec_t *comp_density,
                                  rhea_composition_options_t *opt)
{
  rhea_composition_compute_rhs_vel_fn_data_t  data;

  /* set right-hand side */
  data.comp_density = comp_density;
  data.comp_options = opt;
  ymir_cvec_set_function (rhs_vel, rhea_composition_compute_rhs_vel_fn, &data);
}

void rhea_composition_convert (ymir_vec_t *composition,
							ymir_vec_t *comp_density, ymir_vec_t *comp_visc,
							ymir_mesh_t *ymir_mesh,
							rhea_composition_options_t *opt)
{
  switch (opt->type) {
  case RHEA_COMPOSITION_FLAVOR:
	  comp_density = rhea_composition_to_density(composition, ymir_mesh, opt);
	  comp_visc = rhea_composition_to_viscosity(composition, ymir_mesh, opt);;
	  break;
  case RHEA_COMPOSITION_DENSITY:
	  comp_density = rhea_composition_new(ymir_mesh);
	  ymir_vec_copy (composition, comp_density);
	  comp_visc = NULL;
	  break;
  case RHEA_COMPOSITION_NONE:
	  comp_density = NULL;
	  comp_visc = NULL;
	  break;
  default: /* unknown composition type */
  	  RHEA_ABORT_NOT_REACHED ();
  }
}

static ymir_vec_t *
rhea_composition_to_density(ymir_vec_t *composition,
							ymir_mesh_t *ymir_mesh,
							rhea_composition_options_t *opt)
{
	ymir_vec_t		*comp_density;
	rhea_composition_compute_rhs_vel_fn_data_t  data;

	/* set right-hand side */
	data.comp_density = composition;
	data.comp_options = opt;

	comp_density = rhea_composition_new(ymir_mesh);

	ymir_cvec_set_function (comp_density, rhea_composition_to_density_fn, &data);

	return comp_density;
}

/**
 * Callback function to compute the density based on composition.
 */
static void
rhea_composition_to_density_fn (double *comp_density, double x, double y, double z,
                                     ymir_locidx_t nodeid, void *data)
{
  rhea_composition_compute_rhs_vel_fn_data_t  *d = data;
  rhea_composition_options_t  *opt = d->comp_options;
  const double  composition = *ymir_cvec_index (d->comp_density, nodeid, 0);
  int			index;

  /* compute right-hand side */
  index = composition;
  *comp_density = opt->density_flavor[index];
}

static ymir_vec_t *
rhea_composition_to_viscosity(ymir_vec_t *composition,
							ymir_mesh_t *ymir_mesh,
							rhea_composition_options_t *opt)
{
	ymir_vec_t		*comp_visc;
	rhea_composition_compute_rhs_vel_fn_data_t  data;

	/* set right-hand side */
	data.comp_density = composition;
	data.comp_options = opt;

	comp_visc = rhea_composition_new(ymir_mesh);

	ymir_cvec_set_function (comp_visc, rhea_composition_to_viscosity_fn, &data);

	return comp_visc;
}

/**
 * Callback function to compute the density based on composition.
 */
static void
rhea_composition_to_viscosity_fn (double *comp_visc, double x, double y, double z,
                                     ymir_locidx_t nodeid, void *data)
{
  rhea_composition_compute_rhs_vel_fn_data_t  *d = data;
  rhea_composition_options_t  *opt = d->comp_options;
  const double      composition = *ymir_cvec_index (d->comp_density, nodeid, 0);
  int				index;

  /* compute right-hand side */
  index = composition;
  *comp_visc = opt->visc_flavor[index];
}
