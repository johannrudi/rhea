/*
 */

#include <rhea.h>

void
rhea_add_suboptions_all (ymir_options_t * options)
{
  rhea_domain_add_suboptions (options);
  rhea_viscosity_add_suboptions (options);
}
