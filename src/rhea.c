/*
 */

#include <rhea.h>
#include <ymir.h>

void
rhea_add_options_all (ymir_options_t *options)
{
  rhea_domain_add_options (options);
  rhea_temperature_add_options (options);
  rhea_temperature_add_options_sinker (options);
  rhea_viscosity_add_options (options);
  rhea_discretization_add_options (options);
}

void
rhea_process_options_all (rhea_domain_options_t *domain_options,
                          rhea_temperature_options_t *temperature_options,
                          rhea_viscosity_options_t *viscosity_options,
                          rhea_discretization_options_t *discr_options)
{
  rhea_domain_process_options (domain_options);
  rhea_temperature_process_options (temperature_options, domain_options);
  rhea_viscosity_process_options (viscosity_options, domain_options);
  rhea_discretization_process_options (discr_options, domain_options);
}

int
rhea_get_production_run ()
{
  return ymir_get_production_run ();
}

void
rhea_set_production_run (const int is_production_run)
{
  ymir_set_production_run (is_production_run);;
}

