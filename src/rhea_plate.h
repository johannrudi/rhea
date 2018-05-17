/* RHEA_PLATE TODO write comment */

#ifndef RHEA_PLATE_H
#define RHEA_PLATE_H

#include <rhea_domain.h>

/* generic label for plate "none" */
#define RHEA_PLATE_NONE (-1)

/******************************************************************************
 * Options
 *****************************************************************************/

/* options for plates */
typedef struct rhea_plate_options
{
  /* binary/text files with coordinates of plate polygons */
  char               *polygons_file_path_txt;

  /* options & properties of the computational domain */
  rhea_domain_options_t *domain_options;
}
rhea_plate_options_t;

/**
 * Defines options and adds them as sub-options.
 */
void                rhea_plate_add_options (ymir_options_t * opt_sup);

/**
 * Processes options and stores them.
 */
void                rhea_plate_process_options (
                                        rhea_plate_options_t *opt,
                                        rhea_domain_options_t *domain_options);

/******************************************************************************
 * Plates of Cube Domain
 *****************************************************************************/

typedef enum
{
  RHEA_PLATE_CUBE_NONE = RHEA_PLATE_NONE, /* "none" (must come first) */
  RHEA_PLATE_CUBE_NW = 0,                 /* North-West (must start at zero) */
  RHEA_PLATE_CUBE_NE,                     /* North-East */
  RHEA_PLATE_CUBE_SE,                     /* South-East */
  RHEA_PLATE_CUBE_SW,                     /* South-West */
  RHEA_PLATE_CUBE_N                       /* (number of plates) */
}
rhea_plate_cube_label_t;

rhea_plate_cube_label_t   rhea_plate_cube_get_label (const double test_x,
                                                     const double test_y);

/******************************************************************************
 * Plates of Earth Domain
 *****************************************************************************/

typedef enum
{
  RHEA_PLATE_EARTH_NONE = RHEA_PLATE_NONE,  /* "none" (must come first) */

  /* MORVEL (25 plates) */
  RHEA_PLATE_EARTH_AM = 0,  /* Amur (must start at zero) */
  RHEA_PLATE_EARTH_AN,      /* Antarctic */
  RHEA_PLATE_EARTH_AR,      /* Arabia */
  RHEA_PLATE_EARTH_AU,      /* Australia */
  RHEA_PLATE_EARTH_CP,      /* Capricorn */
  RHEA_PLATE_EARTH_CA,      /* Caribbean */
  RHEA_PLATE_EARTH_CO,      /* Cocos */
  RHEA_PLATE_EARTH_EU,      /* Eurasia */
  RHEA_PLATE_EARTH_IN,      /* India */
  RHEA_PLATE_EARTH_JF,      /* Juan de Fuca */
  RHEA_PLATE_EARTH_LW,      /* Lwandle */
  RHEA_PLATE_EARTH_MQ,      /* Macquarie */
  RHEA_PLATE_EARTH_NZ,      /* Nazca */
  RHEA_PLATE_EARTH_NA,      /* North America */
  RHEA_PLATE_EARTH_NU,      /* Nubia */
  RHEA_PLATE_EARTH_PA,      /* Pacific */
  RHEA_PLATE_EARTH_PS,      /* Philippine Sea */
  RHEA_PLATE_EARTH_RI,      /* Rivera */
  RHEA_PLATE_EARTH_SW,      /* Sandwich */
  RHEA_PLATE_EARTH_SC,      /* Scotia */
  RHEA_PLATE_EARTH_SM,      /* Somalia */
  RHEA_PLATE_EARTH_SA,      /* South America */
  RHEA_PLATE_EARTH_SU,      /* Sundaland */
  RHEA_PLATE_EARTH_SR,      /* Sur */
  RHEA_PLATE_EARTH_YZ,      /* Yangtze */

  /* Bird, 2003 (31 plates) */
  RHEA_PLATE_EARTH_AS,      /* Aegean Sea */
  RHEA_PLATE_EARTH_AP,      /* Altiplano */
  RHEA_PLATE_EARTH_AT,      /* Anatolia */
  RHEA_PLATE_EARTH_BR,      /* Bahnoral Reef */
  RHEA_PLATE_EARTH_BS,      /* Banda Sea */
  RHEA_PLATE_EARTH_BH,      /* Birds Head */
  RHEA_PLATE_EARTH_BU,      /* Burma */
  RHEA_PLATE_EARTH_CL,      /* Caroline */
  RHEA_PLATE_EARTH_CR,      /* Conway Reef */
  RHEA_PLATE_EARTH_EA,      /* Easter */
  RHEA_PLATE_EARTH_FT,      /* Futuna */
  RHEA_PLATE_EARTH_GP,      /* Galapagos */
  RHEA_PLATE_EARTH_JZ,      /* Juan Fernandez */
  RHEA_PLATE_EARTH_KE,      /* Kermadec */
  RHEA_PLATE_EARTH_MN,      /* Manus */
  RHEA_PLATE_EARTH_MO,      /* Maoke */
  RHEA_PLATE_EARTH_MA,      /* Mariana */
  RHEA_PLATE_EARTH_MS,      /* Molucca Sea */
  RHEA_PLATE_EARTH_NH,      /* New Hebrides */
  RHEA_PLATE_EARTH_NI,      /* Niuafo'ou */
  RHEA_PLATE_EARTH_ND,      /* North Andes */
  RHEA_PLATE_EARTH_NB,      /* North Bismarck */
  RHEA_PLATE_EARTH_OK,      /* Okhotsk */
  RHEA_PLATE_EARTH_ON,      /* Okinawa */
  RHEA_PLATE_EARTH_PM,      /* Panama */
  RHEA_PLATE_EARTH_SL,      /* Shetland */
  RHEA_PLATE_EARTH_SS,      /* Solomon Sea */
  RHEA_PLATE_EARTH_SB,      /* South Bismarck */
  RHEA_PLATE_EARTH_TI,      /* Timor */
  RHEA_PLATE_EARTH_TO,      /* Tonga */
  RHEA_PLATE_EARTH_WL,      /* Woodlark */

  RHEA_PLATE_EARTH_N        /* (number of plates) */
}
rhea_plate_earth_label_t;

rhea_plate_earth_label_t  rhea_plate_earth_get_label (const double test_x,
                                                      const double test_y);

#endif /* RHEA_PLATE_H */
