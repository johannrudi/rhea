#include <rhea_math.h>
#include <rhea_base.h>
#include <fenv.h>

double
rhea_math_exp (const double x)
{
  double              result;

  RHEA_ASSERT (isfinite (x));

  feclearexcept (FE_ALL_EXCEPT);
  result = exp (x);
  if (fetestexcept (FE_OVERFLOW)) { /* if argument of exp is too large */
    result = DBL_MAX;
  }

  return result;
}

double
rhea_math_log (const double x)
{
  double              result;

  RHEA_ASSERT (isfinite (x));
  RHEA_ASSERT (0.0 <= x);

  result = log (x);
  if (!isfinite (result)) {
    result = -DBL_MAX;
  }

  return result;
}

double
rhea_math_smin_logexp_nondim (const double x, const double y, const double p,
                              const double d)
{
  double              minimum, maximum;

  if (x < y) {
    minimum = x;
    maximum = y;
  }
  else if (y < x) {
    minimum = y;
    maximum = x;
  }
  else { /* otherwise x=y */
    return x;
  }

  RHEA_ASSERT (p < 0.0);
  RHEA_ASSERT (0.0 < d);
  RHEA_ASSERT ( isfinite (exp (p/d*(maximum - minimum))) );

  /* return smooth minimum */
  return minimum + log (1.0 + exp (p/d*(maximum - minimum))) * d/p;
}

double
rhea_math_smin_logexp (const double x, const double y, const double p)
{
  return rhea_math_smin_logexp_nondim (x, y, p, SC_MAX (fabs (x), fabs (y)));
}

double
rhea_math_smax_logexp_nondim (const double x, const double y, const double p,
                              const double d)
{
  double              minimum, maximum;

  if (x < y) {
    minimum = x;
    maximum = y;
  }
  else if (y < x) {
    minimum = y;
    maximum = x;
  }
  else { /* otherwise x=y */
    return x;
  }

  RHEA_ASSERT (0.0 < p);
  RHEA_ASSERT (0.0 < d);
  RHEA_ASSERT ( isfinite (exp (p/d*(minimum - maximum))) );

  /* return smooth maximum */
  return maximum + log (1.0 + exp (p/d*(minimum - maximum))) * d/p;
}

double
rhea_math_smax_logexp (const double x, const double y, const double p)
{
  return rhea_math_smax_logexp_nondim (x, y, p, SC_MAX (fabs (x), fabs (y)));
}

double
rhea_math_smin_gpm_nondim (const double x, const double y, const double p,
                           const double d)
{
  RHEA_ASSERT (p < 0.0);
  RHEA_ASSERT (0.0 < d);
  return d * pow (pow (x/d, p) + pow (y/d, p), 1.0/p);
}

double
rhea_math_smin_gpm_dx_nondim (const double x, const double y, const double p,
                              const double d)
{
  RHEA_ASSERT (p < 0.0);
  RHEA_ASSERT (0.0 < d);
  return pow (pow (x/d, p) + pow (y/d, p), 1.0/p - 1.0) *
         pow (x/d, p - 1.0);
}

double
rhea_math_smin_gpm_dx_impl_nondim (const double smin, const double y,
                                   const double p, const double d)
{
  RHEA_ASSERT (p < 0.0);
  RHEA_ASSERT (0.0 < d);
  return pow (pow (smin/d, p) - pow (y/d, p), 1.0 - 1.0/p) *
         pow (smin/d, 1.0-p);
}

double
rhea_math_smin_gpm (const double x, const double y, const double p)
{
  return rhea_math_smin_gpm_nondim (x, y, p, SC_MAX (fabs (x), fabs (y)));
}

double
rhea_math_smax_gpm_nondim (const double x, const double y, const double p,
                           const double d)
{
  RHEA_ASSERT (0.0 < p);
  RHEA_ASSERT (0.0 < d);
  return d * pow (pow (x/d, p) + pow (y/d, p), 1.0/p);
}

double
rhea_math_smax_gpm (const double x, const double y, const double p)
{
  return rhea_math_smax_gpm_nondim (x, y, p, SC_MAX (fabs (x), fabs (y)));
}
