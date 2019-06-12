/** RHEA_WEAKZONE_LABEL
 *
 * Weak zone labels that correspond to specific plate boudaries.
 */

#ifndef RHEA_WEAKZONE_LABEL_H
#define RHEA_WEAKZONE_LABEL_H

/* labels for distinguishing between weak zones */
typedef enum
{
  RHEA_WEAKZONE_LABEL_UNKNOWN        = -1,

  /* the "none" class is used when the specific type (ridge/fault) is unknown */
  RHEA_WEAKZONE_LABEL_CLASS_NONE     = 0,

  /* slabs are described by general surfaces/manifolds */
  RHEA_WEAKZONE_LABEL_CLASS_SLAB     = 1,

  /* ridges and fractures are described by lines on the domain's top surface
   * that are extended radially to a max depth */
  RHEA_WEAKZONE_LABEL_CLASS_RIDGE    = 2,
  RHEA_WEAKZONE_LABEL_CLASS_FRACTURE = 3,

  /* slabs of earth */
  RHEA_WEAKZONE_LABEL_EARTH_SL_ALU = 1001,
  RHEA_WEAKZONE_LABEL_EARTH_SL_CAL = 1002,
  RHEA_WEAKZONE_LABEL_EARTH_SL_CAM = 1003,
  RHEA_WEAKZONE_LABEL_EARTH_SL_CAR = 1004,
  RHEA_WEAKZONE_LABEL_EARTH_SL_CAS = 1005,
  RHEA_WEAKZONE_LABEL_EARTH_SL_COT = 1006,
  RHEA_WEAKZONE_LABEL_EARTH_SL_HAL = 1007,
  RHEA_WEAKZONE_LABEL_EARTH_SL_HEL = 1008,
  RHEA_WEAKZONE_LABEL_EARTH_SL_HIM = 1009,
  RHEA_WEAKZONE_LABEL_EARTH_SL_HIN = 1010,
  RHEA_WEAKZONE_LABEL_EARTH_SL_IZU = 1011,
  RHEA_WEAKZONE_LABEL_EARTH_SL_KER = 1012,
  RHEA_WEAKZONE_LABEL_EARTH_SL_KUR = 1013,
  RHEA_WEAKZONE_LABEL_EARTH_SL_MAK = 1014,
  RHEA_WEAKZONE_LABEL_EARTH_SL_MAN = 1015,
  RHEA_WEAKZONE_LABEL_EARTH_SL_MUE = 1016,
  RHEA_WEAKZONE_LABEL_EARTH_SL_PAM = 1017,
  RHEA_WEAKZONE_LABEL_EARTH_SL_PHI = 1018,
  RHEA_WEAKZONE_LABEL_EARTH_SL_PNG = 1019,
  RHEA_WEAKZONE_LABEL_EARTH_SL_PUY = 1020,
  RHEA_WEAKZONE_LABEL_EARTH_SL_RYU = 1021,
  RHEA_WEAKZONE_LABEL_EARTH_SL_SAM = 1022,
  RHEA_WEAKZONE_LABEL_EARTH_SL_SCO = 1023,
  RHEA_WEAKZONE_LABEL_EARTH_SL_SOL = 1024,
  RHEA_WEAKZONE_LABEL_EARTH_SL_SUL = 1025,
  RHEA_WEAKZONE_LABEL_EARTH_SL_SUM = 1026,
  RHEA_WEAKZONE_LABEL_EARTH_SL_VAN = 1027,

  /* ridges of earth */
  RHEA_WEAKZONE_LABEL_EARTH_RI_KE_AU   = 2001,
  RHEA_WEAKZONE_LABEL_EARTH_RI_TO_NI_1 = 2002,
  RHEA_WEAKZONE_LABEL_EARTH_RI_TO_NI_2 = 2003,
  RHEA_WEAKZONE_LABEL_EARTH_RI_TO_AU   = 2004,
  RHEA_WEAKZONE_LABEL_EARTH_RI_BR_NH   = 2005,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AU_CR   = 2006,
  RHEA_WEAKZONE_LABEL_EARTH_RI_NH_PA_1 = 2007,
  RHEA_WEAKZONE_LABEL_EARTH_RI_NH_PA_2 = 2008,
  RHEA_WEAKZONE_LABEL_EARTH_RI_CR_NH   = 2009,
  RHEA_WEAKZONE_LABEL_EARTH_RI_PA_BR   = 2010,
  RHEA_WEAKZONE_LABEL_EARTH_RI_ANDAMAN = 2011,
  RHEA_WEAKZONE_LABEL_EARTH_RI_PS_MA_1 = 2012,
  RHEA_WEAKZONE_LABEL_EARTH_RI_PS_MA_2 = 2013,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AYU     = 2014,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AF_AN_1 = 2015,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AF_AN_2 = 2016,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AF_AN_3 = 2017,
  RHEA_WEAKZONE_LABEL_EARTH_RI_SO_AN_1 = 2018,
  RHEA_WEAKZONE_LABEL_EARTH_RI_SO_AN_2 = 2019,
  RHEA_WEAKZONE_LABEL_EARTH_RI_SO_IN   = 2020,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AU_SO   = 2021,
  RHEA_WEAKZONE_LABEL_EARTH_RI_NA_AF   = 2022,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AF_SA_1 = 2023,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AF_SA_2 = 2024,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AF_SA_3 = 2025,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AU_AN_1 = 2026,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AU_AN_2 = 2027,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AU_AN_3 = 2028,
  RHEA_WEAKZONE_LABEL_EARTH_RI_PA_AN   = 2029,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AN_SA   = 2030,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AN_NZ   = 2031,
  RHEA_WEAKZONE_LABEL_EARTH_RI_EA_PA   = 2032,
  RHEA_WEAKZONE_LABEL_EARTH_RI_EA_NZ   = 2033,
  RHEA_WEAKZONE_LABEL_EARTH_RI_JZ_PA   = 2034,
  RHEA_WEAKZONE_LABEL_EARTH_RI_JZ_AN   = 2035,
  RHEA_WEAKZONE_LABEL_EARTH_RI_JZ_NZ   = 2036,
  RHEA_WEAKZONE_LABEL_EARTH_RI_NZ_PA   = 2037,
  RHEA_WEAKZONE_LABEL_EARTH_RI_CO_NZ_1 = 2038,
  RHEA_WEAKZONE_LABEL_EARTH_RI_NZ_PA_1 = 2039,
  RHEA_WEAKZONE_LABEL_EARTH_RI_NZ_PA_2 = 2040,
  RHEA_WEAKZONE_LABEL_EARTH_RI_CO_NZ_2 = 2041,
  RHEA_WEAKZONE_LABEL_EARTH_RI_GP_PA   = 2042,
  RHEA_WEAKZONE_LABEL_EARTH_RI_GP_NZ   = 2043,
  RHEA_WEAKZONE_LABEL_EARTH_RI_GP_CO   = 2044,
  RHEA_WEAKZONE_LABEL_EARTH_RI_CO_PA   = 2045,
  RHEA_WEAKZONE_LABEL_EARTH_RI_PA_CO   = 2046,
  RHEA_WEAKZONE_LABEL_EARTH_RI_RI_PA   = 2047,
  RHEA_WEAKZONE_LABEL_EARTH_RI_NA_PA   = 2048,
  RHEA_WEAKZONE_LABEL_EARTH_RI_PA_JF   = 2049,
  RHEA_WEAKZONE_LABEL_EARTH_RI_SW_SC   = 2050,
  RHEA_WEAKZONE_LABEL_EARTH_RI_NA_EU_1 = 2051,
  RHEA_WEAKZONE_LABEL_EARTH_RI_NA_EU_2 = 2052,
  RHEA_WEAKZONE_LABEL_EARTH_RI_EU_NA   = 2053,
  RHEA_WEAKZONE_LABEL_EARTH_RI_IN_SO   = 2054,
  RHEA_WEAKZONE_LABEL_EARTH_RI_IN_AR   = 2055,
  RHEA_WEAKZONE_LABEL_EARTH_RI_SO_AR   = 2056,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AF_AR   = 2057,
  RHEA_WEAKZONE_LABEL_EARTH_RI_AR_AF   = 2058,

  /* fractures of earth */
  RHEA_WEAKZONE_LABEL_EARTH_FZ_SB_WL   = 3001,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_PS_MA_1 = 3002,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_PS_MA_2 = 3003,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_AU_CR   = 3004,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_BR_CR   = 3005,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_BR_AU_1 = 3006,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_BR_AU_2 = 3007,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_PA_BR   = 3008,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_FT_PA   = 3009,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_KE_TO   = 3010,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_KE_AU   = 3011,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_PA_AU_1 = 3012,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_AU_PA   = 3013,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_PA_AU_2 = 3014,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_NZ_PM   = 3015,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_NA_RI   = 3016,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_NA_PA_1 = 3017,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_NA_PA_2 = 3018,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_OK_PS   = 3019,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_SC_AN   = 3020,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_SC_SA_1 = 3021,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_SC_SA_2 = 3022,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_AN_SA   = 3023,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_CA_SA   = 3024,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_CA_ND   = 3025,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_CA_NA   = 3026,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_EU_IN   = 3027,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_IN_EU   = 3028,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_AR_IN   = 3029,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_AR_AF   = 3030,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_AT_AR   = 3031,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_AR_EU   = 3032,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_EU_AR   = 3033,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_AT_EU   = 3034,
  RHEA_WEAKZONE_LABEL_EARTH_FZ_AS_EU   = 3035
}
rhea_weakzone_label_t;

/* number of label classes */
#define RHEA_WEAKZONE_LABEL_CLASS_N 3

/* label counts for earth */
#define RHEA_WEAKZONE_LABEL_EARTH_N_SL 27
#define RHEA_WEAKZONE_LABEL_EARTH_N_RI 58
#define RHEA_WEAKZONE_LABEL_EARTH_N_FZ 35
#define RHEA_WEAKZONE_LABEL_EARTH_N (RHEA_WEAKZONE_LABEL_EARTH_N_SL + \
                                     RHEA_WEAKZONE_LABEL_EARTH_N_RI + \
                                     RHEA_WEAKZONE_LABEL_EARTH_N_FZ)

static inline int
rhea_weakzone_label_is_class (const rhea_weakzone_label_t label)
{
  return (RHEA_WEAKZONE_LABEL_CLASS_NONE <= label &&
          label <= RHEA_WEAKZONE_LABEL_CLASS_N);
}

static inline rhea_weakzone_label_t
rhea_weakzone_label_get_class (const rhea_weakzone_label_t label)
{
  rhea_weakzone_label_t class_id;

  /* return itself if label corresponds to a class */
  if (RHEA_WEAKZONE_LABEL_CLASS_NONE <= label &&
      label <= RHEA_WEAKZONE_LABEL_CLASS_N) {
    return label;
  }

  /* calculate class from label */
  class_id = (rhea_weakzone_label_t) ((int) label / 1000);
  if (RHEA_WEAKZONE_LABEL_CLASS_NONE <= class_id &&
      class_id <= RHEA_WEAKZONE_LABEL_CLASS_N) {
    return class_id;
  }

  /* otherwise the class is unknown */
  return RHEA_WEAKZONE_LABEL_UNKNOWN;
}

static inline int
rhea_weakzone_label_is_valid (const rhea_weakzone_label_t label)
{
  return (RHEA_WEAKZONE_LABEL_UNKNOWN !=
          rhea_weakzone_label_get_class (label));
}

static inline int
rhea_weakzone_label_is_valid_int (const int label)
{
  return (RHEA_WEAKZONE_LABEL_UNKNOWN !=
          rhea_weakzone_label_get_class ((rhea_weakzone_label_t) label));
}

static inline int
rhea_weakzone_label_is_slab (const rhea_weakzone_label_t label)
{
  return (RHEA_WEAKZONE_LABEL_CLASS_SLAB ==
          rhea_weakzone_label_get_class (label));
}

static inline int
rhea_weakzone_label_is_ridge (const rhea_weakzone_label_t label)
{
  return (RHEA_WEAKZONE_LABEL_CLASS_RIDGE ==
          rhea_weakzone_label_get_class (label));
}

static inline int
rhea_weakzone_label_is_fracture (const rhea_weakzone_label_t label)
{
  return (RHEA_WEAKZONE_LABEL_CLASS_FRACTURE ==
          rhea_weakzone_label_get_class (label));
}

/**
 * Gets the array index in (0:RHEA_WEAKZONE_LABEL_EARTH_N) corresponding to a
 * label for earth.
 */
static inline int
rhea_weakzone_label_earth_get_idx (const rhea_weakzone_label_t label)
{
  const rhea_weakzone_label_t class_id = (rhea_weakzone_label_t)
                                         ((int) label / 1000);

  switch (class_id) {
  case RHEA_WEAKZONE_LABEL_CLASS_SLAB:
    return (label % (1000*RHEA_WEAKZONE_LABEL_CLASS_SLAB)) - 1;
  case RHEA_WEAKZONE_LABEL_CLASS_RIDGE:
    return RHEA_WEAKZONE_LABEL_EARTH_N_SL +
           (label % (1000*RHEA_WEAKZONE_LABEL_CLASS_RIDGE)) - 1;
  case RHEA_WEAKZONE_LABEL_CLASS_FRACTURE:
    return RHEA_WEAKZONE_LABEL_EARTH_N_SL +
           RHEA_WEAKZONE_LABEL_EARTH_N_RI +
           (label % (1000*RHEA_WEAKZONE_LABEL_CLASS_FRACTURE)) - 1;
  default: /* unknown label */
    return -1;
  }
}

/**
 * Gets the label corresponding to an array index in
 * (0:RHEA_WEAKZONE_LABEL_EARTH_N).
 */
static inline rhea_weakzone_label_t
rhea_weakzone_label_earth_get_label (const int idx)
{
  if (0 <= idx &&
      idx < RHEA_WEAKZONE_LABEL_EARTH_N_SL) { /* if slab */
    return (rhea_weakzone_label_t)
           (1000*RHEA_WEAKZONE_LABEL_CLASS_SLAB + idx + 1);
  }
  if (RHEA_WEAKZONE_LABEL_EARTH_N_SL <= idx &&
      idx < RHEA_WEAKZONE_LABEL_EARTH_N_RI) { /* if ridge */
    return (rhea_weakzone_label_t)
           (1000*RHEA_WEAKZONE_LABEL_CLASS_RIDGE + idx + 1);
  }
  if (RHEA_WEAKZONE_LABEL_EARTH_N_RI <= idx &&
      idx < RHEA_WEAKZONE_LABEL_EARTH_N_FZ) { /* if fracture */
    return (rhea_weakzone_label_t)
           (1000*RHEA_WEAKZONE_LABEL_CLASS_FRACTURE + idx + 1);
  }

  /* otherwise index cannot be assigned to a label */
  return RHEA_WEAKZONE_LABEL_UNKNOWN;
}

#endif /* RHEA_WEAKZONE_LABEL_H */
