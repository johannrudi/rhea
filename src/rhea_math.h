/** RHEA_MATH
 *
 * Basic mathematical operations.
 */

#ifndef RHEA_MATH_H
#define RHEA_MATH_H

/**
 * Computes exponential or (natural) logarithm and catches infinite output.
 */
double              rhea_math_exp (const double x);
double              rhea_math_log (const double x);

/**
 * Computes a smooth version of the minimum avoiding overflow and underflow:
 *
 *   smin (x, y) = log (exp (p*x) + exp(p*y)) / p
 *                       = min + log (1 + exp(p*(max - min))) / p
 *
 * where `min=min(x,y)`, `max=max(x,y)`, and the parameter `p<0` controls
 * "smoothness" such that
 *
 *   smin (x, y) -> min (x, y)  as  p<0 and p -> -infinity
 */
double              rhea_math_smin_logexp (const double x, const double y,
                                           const double p);

double              rhea_math_smin_logexp_nondim (const double x,
                                                  const double y,
                                                  const double p,
                                                  const double d);

/**
 * Computes a smooth version of the maximum avoiding overflow and underflow:
 *
 *   smax (x, y) = log (exp (p*x) + exp(p*y)) / p
 *                       = max + log (1 + exp(p*(min - max))) / p
 *
 * where `min=min(x,y)`, `max=max(x,y)`, and the parameter `p>0` controls
 * "smoothness" such that
 *
 *   smax (x, y) -> max (x, y)  as  p>0 and p -> +infinity
 */
double              rhea_math_smax_logexp (const double x, const double y,
                                           const double p);

double              rhea_math_smax_logexp_nondim (const double x,
                                                  const double y,
                                                  const double p,
                                                  const double d);

/**
 * Computes a smooth version of the minimum based on the generalized/power mean.
 */
double              rhea_math_smin_gpm (const double x, const double y,
                                        const double p);

double              rhea_math_smin_gpm_nondim (const double x, const double y,
                                               const double p, const double d);

double              rhea_math_smin_gpm_dx_nondim (const double x,
                                                  const double y,
                                                  const double p,
                                                  const double d);

double              rhea_math_smin_gpm_dx_impl_nondim (const double smin,
                                                       const double y,
                                                       const double p,
                                                       const double d);

/**
 * Computes a smooth version of the mximum based on the generalized/power mean.
 */
double              rhea_math_smax_gpm (const double x, const double y,
                                        const double p);

double              rhea_math_smax_gpm_nondim (const double x, const double y,
                                               const double p, const double d);

#endif /* RHEA_MATH_H */
