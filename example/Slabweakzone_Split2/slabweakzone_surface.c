/**************************************
 * Geometry Transformations
 *************************************/

static void
slabs_X_fn_identity (mangll_tag_t tag, mangll_locidx_t np,
                       const double *_sc_restrict EX,
                       const double *_sc_restrict EY,
                       const double *_sc_restrict EZ,
                       double *_sc_restrict X,
                       double *_sc_restrict Y,
                       double *_sc_restrict Z, void *data)
{
  mangll_locidx_t     il;

  for (il = 0; il < np; ++il) {
    X[il] = EX[il];
    Y[il] = EY[il];
    Z[il] = EZ[il];
  }
}

static void
slabs_X_fn_sine (mangll_tag_t tag, mangll_locidx_t np,
                       const double *_sc_restrict EX,
                       const double *_sc_restrict EY,
                       const double *_sc_restrict EZ,
                       double *_sc_restrict X,
                       double *_sc_restrict Y,
                       double *_sc_restrict Z, void *data)
{
  mangll_locidx_t     il;
  const double        slope = 0.5;

  for (il = 0; il < np; ++il) {
    X[il] = EX[il];
    Y[il] = EY[il];
    Z[il] = EZ[il] * (1 + slope * sin(M_PI * EY[il]));
  }
}

static void
slabs_X_fn_profile (mangll_tag_t tag, mangll_locidx_t np,
                   const double *_sc_restrict EX,
                   const double *_sc_restrict EY,
                   const double *_sc_restrict EZ,
                   double *_sc_restrict X,
                   double *_sc_restrict Y,
                   double *_sc_restrict Z, void *data)
{
  mangll_locidx_t     il;
  double factor;
  slabs_topo_profile_t *topo = (slabs_topo_profile_t *) data;
  double *tX = topo->tX;
  double *tY = topo->tY;
  double *tZ = topo->tZ;
  int     m, nsurf = topo->nsurf;

  for (il = 0; il < np; ++il) {
    X[il] = EX[il];
    Y[il] = EY[il];
    for (m = 0; m < nsurf; m++)  {
      if (fabs(EX[il] - tX[m]) < SC_1000_EPS &&
          fabs(EY[il] - tY[m]) < SC_1000_EPS)  {
        factor = tZ[m];
      }
    }
    Z[il] = EZ[il] * factor;
  }
}

void
slabs_surface_location (slabs_options_t *slabs_options,
                       rhea_discretization_options_t *discr_options)
{
 slabs_surf_options_t  *slabs_surf_options = slabs_options->slabs_surf_options;
 slabs_topo_profile_t *topo = slabs_surf_options->topo_profile;

    /* set custom X-function */
  switch (slabs_surf_options->x_func) {
    case SLABS_X_FUNCTION_IDENTITY:
      rhea_discretization_set_user_X_fn (discr_options,
                                     slabs_X_fn_identity, NULL);
      break;

    case SLABS_X_FUNCTION_SINE:
      rhea_discretization_set_user_X_fn (discr_options,
                                       slabs_X_fn_sine, NULL);
      break;

    case SLABS_X_FUNCTION_PROFILE:
      rhea_discretization_set_user_X_fn (discr_options,
                                       slabs_X_fn_profile, topo);
      break;

    default:
      RHEA_ABORT_NOT_REACHED ();
    break;
  }
}

