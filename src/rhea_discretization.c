/*
 */

#include <rhea_discretization.h>
#include <rhea_base.h>

/* default options */
#define RHEA_DISCRETIZATION_DEFAULT_LEVEL_MIN (1)
#define RHEA_DISCRETIZATION_DEFAULT_LEVEL_MAX (20)

int                 rhea_discretization_level_min =
  RHEA_DISCRETIZATION_DEFAULT_LEVEL_MIN;
int                 rhea_discretization_level_max =
  RHEA_DISCRETIZATION_DEFAULT_LEVEL_MAX;

void
rhea_discretization_add_options (ymir_options_t * opt_sup)
{
  const char         *opt_prefix = "Discretization";
  ymir_options_t     *opt = ymir_options_new ();

  /* *INDENT-OFF* */
  ymir_options_addv (opt,

  YMIR_OPTIONS_I, "level-min", '\0',
    &rhea_discretization_level_min, RHEA_DISCRETIZATION_DEFAULT_LEVEL_MIN,
    "Minumum level of mesh refinement",
  YMIR_OPTIONS_I, "level-max", '\0',
    &rhea_discretization_level_max, RHEA_DISCRETIZATION_DEFAULT_LEVEL_MAX,
    "Maximum level of mesh refinement",

  YMIR_OPTIONS_END_OF_LIST);
  /* *INDENT-ON* */

  /* add these options as sub-options */
  ymir_options_add_suboptions (opt_sup, opt, opt_prefix);
  ymir_options_destroy (opt);
}
