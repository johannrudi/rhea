/*
 */

#include <rhea_base.h>
#include <ymir.h>

int                 rhea_package_id = -1;

void
rhea_init (sc_log_handler_t log_handler, int log_threshold)
{
  const int           w = 24;

  rhea_package_id = sc_package_register (log_handler, log_threshold,
                                         "rhea", "Global mantle convection");

  RHEA_GLOBAL_ESSENTIALF ("This is %s\n", RHEA_PACKAGE_STRING);
  RHEA_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "CC", RHEA_CC);
//RHEA_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "C_VERSION", RHEA_C_VERSION);
  RHEA_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "CFLAGS", RHEA_CFLAGS);
  RHEA_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "CPP", RHEA_CPP);
  RHEA_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "CPPFLAGS", RHEA_CPPFLAGS);
  RHEA_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "LDFLAGS", RHEA_LDFLAGS);
  RHEA_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "LIBS", RHEA_LIBS);
}

void
rhea_initialize (int argc, char **argv, MPI_Comm mpicomm)
{
  sc_log_handler_t    log_handler = NULL;
  int                 rhea_log_threshold = SC_LP_DEFAULT;
  ymir_log_threshold_set_t  ymir_log_threshold;

  /* initialize sub-packages */
  ymir_log_threshold.sc = SC_LP_INFO;
  ymir_log_threshold.p4est = SC_LP_INFO;
  ymir_log_threshold.mangll = SC_LP_INFO;
  ymir_log_threshold.ymir = SC_LP_DEFAULT;
  ymir_initialize_log (argc, argv, mpicomm, log_handler, &ymir_log_threshold);

  /* initialize this package */
  if (rhea_package_id < 0) {
    rhea_init (log_handler, rhea_log_threshold);
  }
}

void
rhea_finalize ()
{
  /* finalize this package */
  sc_package_unregister (rhea_package_id);
  rhea_package_id = -1;

  /* finalize sub-packages */
  ymir_finalize ();
}

#ifndef __cplusplus
#undef RHEA_GLOBAL_LOGF
#undef RHEA_LOGF
#undef RHEA_GLOBAL_TRACEF
#undef RHEA_GLOBAL_LDEBUGF
#undef RHEA_GLOBAL_VERBOSEF
#undef RHEA_GLOBAL_INFOF
#undef RHEA_GLOBAL_STATISTICSF
#undef RHEA_GLOBAL_PRODUCTIONF
#undef RHEA_GLOBAL_ESSENTIALF
#undef RHEA_GLOBAL_LERRORF
#undef RHEA_TRACEF
#undef RHEA_LDEBUGF
#undef RHEA_VERBOSEF
#undef RHEA_INFOF
#undef RHEA_STATISTICSF
#undef RHEA_PRODUCTIONF
#undef RHEA_ESSENTIALF
#undef RHEA_LERRORF
#endif

#ifndef SC_SPLINT

void
RHEA_GLOBAL_LOGF (int priority, const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  sc_logv ("<unknown>", 0, p4est_package_id, SC_LC_GLOBAL, priority, fmt, ap);
  va_end (ap);
}

void
RHEA_LOGF (int priority, const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  sc_logv ("<unknown>", 0, p4est_package_id, SC_LC_NORMAL, priority, fmt, ap);
  va_end (ap);
}

#define RHEA_LOG_IMP(n,p)                               \
  void                                                  \
  RHEA_GLOBAL_ ## n ## F (const char *fmt, ...)         \
  {                                                     \
    va_list             ap;                             \
    va_start (ap, fmt);                                 \
    sc_logv ("<unknown>", 0, p4est_package_id,          \
             SC_LC_GLOBAL, SC_LP_ ## p, fmt, ap);       \
    va_end (ap);                                        \
  }                                                     \
  void                                                  \
  RHEA_ ## n ## F (const char *fmt, ...)                \
  {                                                     \
    va_list             ap;                             \
    va_start (ap, fmt);                                 \
    sc_logv ("<unknown>", 0, p4est_package_id,          \
             SC_LC_NORMAL, SC_LP_ ## p, fmt, ap);       \
    va_end (ap);                                        \
  }

/* *INDENT-OFF* */
RHEA_LOG_IMP (TRACE, TRACE)
RHEA_LOG_IMP (LDEBUG, DEBUG)
RHEA_LOG_IMP (VERBOSE, VERBOSE)
RHEA_LOG_IMP (INFO, INFO)
RHEA_LOG_IMP (STATISTICS, STATISTICS)
RHEA_LOG_IMP (PRODUCTION, PRODUCTION)
RHEA_LOG_IMP (ESSENTIAL, ESSENTIAL)
RHEA_LOG_IMP (LERROR, ERROR)
/* *INDENT-ON* */

#endif /* !SC_SPLINT */
