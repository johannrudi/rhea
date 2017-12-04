
/****************************************************
 * Curvature Slabs 2plate_poly2
*****************************************************/

/* Computes the radius of a shell domain or the corresponding value for a
 * rectangular domain. */
double
slabs_compute_radius (const double x, const double y, const double z,
                      slabs_options_t *slabs_options)
{
      slabs_domain_options_t    * dom_opt = slabs_options->slabs_domain_options;
      const double    z_max = dom_opt->z_max;
      const double    radius_min = dom_opt->radius_min;
      const double    radius_max = dom_opt->radius_max;

      return z / z_max * (radius_max - radius_min) + radius_min;
}

/* Computes the longitude of a shell domain or the corresponding value for a
 * rectangular domain. */
double
slabs_compute_longitude (const double x, const double y, const double z,
                         slabs_options_t *slabs_options)
{
      slabs_domain_options_t    * dom_opt = slabs_options->slabs_domain_options;
      const double        lon_min = dom_opt->lon_min;
      const double        lon_max = dom_opt->lon_min;
      const double        y_max = dom_opt->y_max;

      return y / y_max * (lon_max - lon_min) + lon_min;
}


/* Evaluates linear polynomial. */
#define slabs_poly1(x,poly_coeff) ( (poly_coeff)[0] + (poly_coeff)[1]*(x) )

/* Evaluates quadratic polynomial.*/
#define slabs_poly2(x,poly_coeff) \
  ( (poly_coeff)[0] + (poly_coeff)[1]*(x) + (poly_coeff)[2]*(x)*(x) )

/* Evaluates derivative of quadratic polynomial.*/
#define slabs_poly2_deriv(x,poly_coeff) \
  ( (poly_coeff)[1] + 2*(poly_coeff)[2]*(x) )

/* Computes interpolating polynomial in 2D via Hermite interpolation.*/
static double *
slabs_compute_poly2_interpolation (double start_node, double start_val,
                                   double start_deriv,
                                   double end_node, double end_val)
{
  double             *poly_coeff;

  /* allocate coefficients */
  poly_coeff = RHEA_ALLOC (double, 3);

  /* calculate coefficients via Hermite interpolation */
  poly_coeff[2] = ( (end_val - start_val) / (end_node - start_node)
                    - start_deriv ) / (end_node - start_node);
  poly_coeff[1] = start_deriv - 2 * poly_coeff[2] * start_node;
  poly_coeff[0] = start_val - (start_deriv * start_node)
                  + (poly_coeff[2] * start_node * start_node);

  /* return coefficients */
  return poly_coeff;
}

/* Runs Newton's method to find the closest point on the curve of a
 * piecewise linear/quadratic polynomial to a given point.*/
static int
slabs_closest_pt_newton (double *x, double point_x, double point_y,
                         double rtol, double maxiter,
                         double *poly1_left_coeff, double stitch_node_left,
                         double *poly2_coeff, double stitch_node_right,
                         double *poly1_right_coeff) {
  double              x_prev;
  int                 k;

  /* run Newton iterations */
  for (k = 0; k < maxiter; k++) {
    /* store previous step */
    x_prev = *x;

    /* compute one Newton step */
    if (*x <= stitch_node_left) { /* if at left lin. poly. */
      *x = *x - ( (*x - point_x) + poly1_left_coeff[1]
                  * (slabs_poly1 (*x, poly1_left_coeff) - point_y) )
              / ( 1.0 + poly1_left_coeff[1] * poly1_left_coeff[1] );
    }
    else if (stitch_node_right <= *x) { /* if at right lin. poly. */
      *x = *x - ( (*x - point_x) + poly1_right_coeff[1]
                  * (slabs_poly1 (*x, poly1_right_coeff) - point_y) )
              / ( 1.0 + poly1_right_coeff[1] * poly1_right_coeff[1] );
    }
    else { /* if in middle at quadratic poly. */
      *x = *x - ( (*x - point_x) + (slabs_poly2 (*x, poly2_coeff) - point_y)
                                 * slabs_poly2_deriv (*x, poly2_coeff) )
              / ( 1.0 + slabs_poly2_deriv (*x, poly2_coeff)
                      * slabs_poly2_deriv (*x, poly2_coeff)
                      + (slabs_poly2 (*x, poly2_coeff) - point_y)
                      * 2.0 * poly2_coeff[2] );
    }

    /* check for convergence */
    RHEA_ASSERT (*x != 0.0);
    if (fabs ((x_prev - *x) / *x) < rtol) {
      break;
    }
  }

  /* return number of iterations */
  return k;
}

/* Computes distance and orientation of a point to a quadratic polynomial in 2D.*/
static double *
slabs_compute_closest_pt_on_poly2 (double point_x, double point_y,
                                   double *poly2_coeff, double start_node,
                                   double start_val, double start_deriv,
                                   double end_node, double end_val,
                                   int *orientation_wrt_curve)
{
  double              poly1_left_coeff[2], poly1_right_coeff[2];
  double              tangent_start[2], tangent_end[2];
  double              orth_start[2], orth_end[2];
  double              line_orientation_start, line_orientation_end;
  double              proj;
  double              dist_sq_mid, dist_sq_start, dist_sq_end;
  //int                 newton_num_iter;
  double             *closest_pt;

  /* check input parameters */
  RHEA_ASSERT (orientation_wrt_curve != NULL);

  /* compute coefficients of linear polynomial extending the quadratic
   * polynomial to the left and right */
  poly1_left_coeff[1] = slabs_poly2_deriv (start_node, poly2_coeff);
  poly1_left_coeff[0] = slabs_poly2 (start_node, poly2_coeff)
                        - poly1_left_coeff[1] * start_node;
  poly1_right_coeff[1] = slabs_poly2_deriv (end_node, poly2_coeff);
  poly1_right_coeff[0] = slabs_poly2 (end_node, poly2_coeff)
                         - poly1_right_coeff[1] * end_node;

  /* compute normalized tangent vectors of curve at start and end points */
  tangent_start[1] = slabs_poly2_deriv (start_node, poly2_coeff);
  tangent_start[0] = 1.0 / sqrt (1.0 + tangent_start[1] * tangent_start[1]);
  tangent_start[1] *= tangent_start[0];
  tangent_end[1] = slabs_poly2_deriv (end_node, poly2_coeff);
  tangent_end[0] = 1.0 / sqrt (1.0 + tangent_end[1] * tangent_end[1]);
  tangent_end[1] *= tangent_end[0];

  /* compute normal vectors to slopes at start and end points */
  orth_start[0] = slabs_poly2_deriv (start_node, poly2_coeff);
  orth_start[1] = -1.0;
  orth_end[0] = slabs_poly2_deriv (end_node, poly2_coeff);
  orth_end[1] = -1.0;

  /* flip signs if x-value is negative*/
  if (orth_start[0] < 0.0) {
    orth_start[0] = -orth_start[0];
    orth_start[1] = -orth_start[1];
  }
  if (orth_end[0] < 0.0) {
    orth_end[0] = -orth_end[0];
    orth_end[1] = -orth_end[1];
  }

  /* compute orientation of target point w.r.t. lines through start and
   * end points along the resp. normals (use cross product) */
  line_orientation_start = orth_start[0] * (point_y - start_val) -
                           orth_start[1] * (point_x - start_node);
  line_orientation_end = orth_end[0] * (point_y - end_val) -
                         orth_end[1] * (point_x - end_node);

  /* compute closest point on curve and set orientation of point w.r.t. curve */
  closest_pt = RHEA_ALLOC (double, 2);
  if (line_orientation_start > 0) { /* if pt. is left of line through start pt.
                                     * with direction normal to start deriv. */
    proj = - (start_node - point_x) * tangent_start[0]
           - (start_val - point_y) * tangent_start[1];
    closest_pt[0] = start_node + proj * tangent_start[0];
    closest_pt[1] = start_val + proj * tangent_start[1];

    *orientation_wrt_curve = SLABS_CURVE_ORIENT_BOTTOM_BACK;
  }
  else if (line_orientation_end < 0) { /* if pt. is right of line through end pt
                                        * with dir. normal to end deriv. */
    proj = - (end_node - point_x) * tangent_end[0]
           - (end_val - point_y) * tangent_end[1];
    closest_pt[0] = end_node + proj * tangent_end[0];
    closest_pt[1] = end_val + proj * tangent_end[1];

    if (point_y <= slabs_poly1 (point_x, poly1_right_coeff)) { /* if pt below */
      *orientation_wrt_curve = SLABS_CURVE_ORIENT_BOTTOM_LEFT;
    }
    else { /* if point above linear poly. */
      *orientation_wrt_curve = SLABS_CURVE_ORIENT_TOP_RIGHT;
    }
  }
  else { /* if pt. is between the lines */
    /* set initial guess for Newton's method */
    closest_pt[0] = (start_node + end_node) / 2.0;
    closest_pt[1] = slabs_poly2 (closest_pt[0], poly2_coeff);
    dist_sq_mid = (point_x - closest_pt[0]) * (point_x - closest_pt[0])
                  + (point_y - closest_pt[1]) * (point_y - closest_pt[1]);
    dist_sq_start = (point_x - start_node) * (point_x - start_node)
                    + (point_y - start_val) * (point_y - start_val);
    dist_sq_end = (point_x - end_node) * (point_x - end_node)
                  + (point_y - end_val) * (point_y - end_val);
    if (dist_sq_mid > dist_sq_start || dist_sq_mid > dist_sq_end) {
      if (dist_sq_start < dist_sq_end) {
        closest_pt[0] = start_node;
      }
      else {
        closest_pt[0] = end_node;
      }
    }

    /* find x-value of closest point via Newton's method in 1D */
    /*newton_num_iter = */
    slabs_closest_pt_newton (closest_pt, point_x, point_y,
                             SLABS_CLOSEST_PT_NEWTON_RTOL,
                             SLABS_CLOSEST_PT_NEWTON_MAXITER,
                             poly1_left_coeff, start_node,
                             poly2_coeff, end_node,
                             poly1_right_coeff);

    /* compute corresponding y-value on curve */
    if (closest_pt[0] <= start_node) { /* if at left lin. poly. */
      closest_pt[1] = slabs_poly1 (closest_pt[0], poly1_left_coeff);
    }
    else if (end_node <= closest_pt[0]) { /* if at right lin. poly. */
      closest_pt[1] = slabs_poly1 (closest_pt[0], poly1_right_coeff);
    }
    else { /* if in middle at quadratic poly. */
      closest_pt[1] = slabs_poly2 (closest_pt[0], poly2_coeff);
    }

    /* set orientation */
    if (point_y <= slabs_poly2 (point_x, poly2_coeff)) { /* if pt. below curve*/
      *orientation_wrt_curve = SLABS_CURVE_ORIENT_BOTTOM;
    }
    else { /* if pt. above curve */
      *orientation_wrt_curve = SLABS_CURVE_ORIENT_TOP;
    }
  }

  /* return distance */
  return closest_pt;
}

/* Computes distance and orientation of a point to a quadratic polynomial in 2D.*/
static double
slabs_compute_rot_on_poly2 (double *poly2_coeff, double *closest_pt)
{
  double              tangent[2];
  double              horiz[2] = {1.0, .0};
  double              rot;

  /* compute vectors of curve at start and end points */
  tangent[1] = slabs_poly2_deriv (closest_pt[0], poly2_coeff);
  tangent[0] = 1.0 / sqrt (1.0 + tangent[1] * tangent[1]);
  tangent[1] *= tangent[0];

  rot = acos(horiz[0]*tangent[0] + horiz[1]*tangent[1]);
  /* return distance */
  return rot;
}


