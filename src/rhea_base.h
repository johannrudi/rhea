/*
 */

#ifndef RHEA_BASE_H
#define RHEA_BASE_H

/* include rhea config header */
#include <rhea_config.h>

/* include sc header */
#include <sc.h>

/* check MPI configuration */
#if \
  (defined (RHEA_ENABLE_MPI) && !defined (SC_ENABLE_MPI)) || \
  (!defined (RHEA_ENABLE_MPI) && defined (SC_ENABLE_MPI))
#error "MPI configured differently in rhea and libsc"
#endif
#if \
  (defined (RHEA_ENABLE_MPIIO) && !defined (SC_ENABLE_MPIIO)) || \
  (!defined (RHEA_ENABLE_MPIIO) && defined (SC_ENABLE_MPIIO))
#error "MPI I/O configured differently in rhea and libsc"
#endif

SC_EXTERN_C_BEGIN;

/* basic macros */
#define RHEA_NOOP() SC_NOOP ()
#define RHEA_ABORT(s) SC_ABORT (s)
#define RHEA_ABORT_NOT_REACHED() SC_ABORT_NOT_REACHED ()

/* checks & assertions */
#define RHEA_CHECK_ABORT(q,s) SC_CHECK_ABORT (q,s)
#ifdef RHEA_ENABLE_DEBUG
#define RHEA_ASSERT(c) SC_CHECK_ABORT ((c), "Assertion '" #c "'")
#else
#define RHEA_ASSERT(c) SC_NOOP ()
#endif

/* memory allocation (will abort if out of memory) */
#define RHEA_ALLOC(t,n)          (t *) sc_malloc (rhea_package_id,    \
                                                  (n) * sizeof(t))
#define RHEA_ALLOC_ZERO(t,n)     (t *) sc_calloc (rhea_package_id,    \
                                                  (size_t) (n), sizeof(t))
#define RHEA_REALLOC(p,t,n)      (t *) sc_realloc (rhea_package_id,   \
                                                   (p), (n) * sizeof(t))
#define RHEA_STRDUP(s)                 sc_strdup (rhea_package_id, (s))
#define RHEA_FREE(p)                   sc_free (rhea_package_id, (p))

/* log helper macros */
#define RHEA_GLOBAL_LOG(p,s)                           \
  SC_GEN_LOG (rhea_package_id, SC_LC_GLOBAL, (p), (s))
#define RHEA_LOG(p,s)                                  \
  SC_GEN_LOG (rhea_package_id, SC_LC_NORMAL, (p), (s))
void                RHEA_GLOBAL_LOGF (int priority, const char *fmt, ...)
  __attribute__ ((format (printf, 2, 3)));
void                RHEA_LOGF (int priority, const char *fmt, ...)
  __attribute__ ((format (printf, 2, 3)));
#ifndef __cplusplus
#define RHEA_GLOBAL_LOGF(p,f,...)                                      \
  SC_GEN_LOGF (rhea_package_id, SC_LC_GLOBAL, (p), (f), __VA_ARGS__)
#define RHEA_LOGF(p,f,...)                                             \
  SC_GEN_LOGF (rhea_package_id, SC_LC_NORMAL, (p), (f), __VA_ARGS__)
#endif

/* global log macros will only print if identifier <= 0 */
#define RHEA_GLOBAL_TRACE(s)      RHEA_GLOBAL_LOG (SC_LP_TRACE, (s))
#define RHEA_GLOBAL_LDEBUG(s)     RHEA_GLOBAL_LOG (SC_LP_DEBUG, (s))
#define RHEA_GLOBAL_VERBOSE(s)    RHEA_GLOBAL_LOG (SC_LP_VERBOSE, (s))
#define RHEA_GLOBAL_INFO(s)       RHEA_GLOBAL_LOG (SC_LP_INFO, (s))
#define RHEA_GLOBAL_STATISTICS(s) RHEA_GLOBAL_LOG (SC_LP_STATISTICS, (s))
#define RHEA_GLOBAL_PRODUCTION(s) RHEA_GLOBAL_LOG (SC_LP_PRODUCTION, (s))
#define RHEA_GLOBAL_ESSENTIAL(s)  RHEA_GLOBAL_LOG (SC_LP_ESSENTIAL, (s))
#define RHEA_GLOBAL_LERROR(s)     RHEA_GLOBAL_LOG (SC_LP_ERROR, (s))
void                RHEA_GLOBAL_TRACEF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_GLOBAL_LDEBUGF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_GLOBAL_VERBOSEF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_GLOBAL_INFOF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_GLOBAL_STATISTICSF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_GLOBAL_PRODUCTIONF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_GLOBAL_ESSENTIALF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_GLOBAL_LERRORF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
#ifndef __cplusplus
#define RHEA_GLOBAL_TRACEF(f,...)                             \
        RHEA_GLOBAL_LOGF (SC_LP_TRACE, (f), __VA_ARGS__)
#define RHEA_GLOBAL_LDEBUGF(f,...)                            \
        RHEA_GLOBAL_LOGF (SC_LP_DEBUG, (f), __VA_ARGS__)
#define RHEA_GLOBAL_VERBOSEF(f,...)                           \
        RHEA_GLOBAL_LOGF (SC_LP_VERBOSE, (f), __VA_ARGS__)
#define RHEA_GLOBAL_INFOF(f,...)                              \
        RHEA_GLOBAL_LOGF (SC_LP_INFO, (f), __VA_ARGS__)
#define RHEA_GLOBAL_STATISTICSF(f,...)                        \
        RHEA_GLOBAL_LOGF (SC_LP_STATISTICS, (f), __VA_ARGS__)
#define RHEA_GLOBAL_PRODUCTIONF(f,...)                        \
        RHEA_GLOBAL_LOGF (SC_LP_PRODUCTION, (f), __VA_ARGS__)
#define RHEA_GLOBAL_ESSENTIALF(f,...)                         \
        RHEA_GLOBAL_LOGF (SC_LP_ESSENTIAL, (f), __VA_ARGS__)
#define RHEA_GLOBAL_LERRORF(f,...)                            \
        RHEA_GLOBAL_LOGF (SC_LP_ERROR, (f), __VA_ARGS__)
#endif
#define RHEA_GLOBAL_NOTICE        RHEA_GLOBAL_STATISTICS
#define RHEA_GLOBAL_NOTICEF       RHEA_GLOBAL_STATISTICSF

/* global log macros for beginning and end of a function (XML type output) */
#define RHEA_GLOBAL_VERBOSE_FN_BEGIN(fn)  RHEA_GLOBAL_VERBOSEF ("<%s>\n", fn)
#define RHEA_GLOBAL_VERBOSE_FN_END(fn)    RHEA_GLOBAL_VERBOSEF ("</%s>\n", fn)
#define RHEA_GLOBAL_VERBOSE_FN_TAG(fn)    RHEA_GLOBAL_VERBOSEF ("<%s />\n", fn)
#define RHEA_GLOBAL_INFO_FN_BEGIN(fn)     RHEA_GLOBAL_INFOF ("<%s>\n", fn)
#define RHEA_GLOBAL_INFO_FN_END(fn)       RHEA_GLOBAL_INFOF ("</%s>\n", fn)
#define RHEA_GLOBAL_INFO_FN_TAG(fn)       RHEA_GLOBAL_INFOF ("<%s />\n", fn)
#define RHEA_GLOBAL_PRODUCTION_FN_BEGIN(fn) \
        RHEA_GLOBAL_PRODUCTIONF ("<%s>\n", fn)
#define RHEA_GLOBAL_PRODUCTION_FN_END(fn) \
        RHEA_GLOBAL_PRODUCTIONF ("</%s>\n", fn)
#define RHEA_GLOBAL_PRODUCTION_FN_TAG(fn) \
        RHEA_GLOBAL_PRODUCTIONF ("<%s />\n", fn)

#define RHEA_GLOBAL_VERBOSEF_FN_BEGIN(fn,a,...) \
        RHEA_GLOBAL_VERBOSEF ("<%s " a ">\n", fn, __VA_ARGS__)
#define RHEA_GLOBAL_VERBOSEF_FN_END(fn,a,...) \
        RHEA_GLOBAL_VERBOSEF ("</%s " a ">\n", fn, __VA_ARGS__)
#define RHEA_GLOBAL_VERBOSEF_FN_TAG(fn,a,...) \
        RHEA_GLOBAL_VERBOSEF ("<%s " a " />\n", fn, __VA_ARGS__)
#define RHEA_GLOBAL_INFOF_FN_BEGIN(fn,a,...) \
        RHEA_GLOBAL_INFOF ("<%s " a ">\n", fn, __VA_ARGS__)
#define RHEA_GLOBAL_INFOF_FN_END(fn,a,...) \
        RHEA_GLOBAL_INFOF ("</%s " a ">\n", fn, __VA_ARGS__)
#define RHEA_GLOBAL_INFOF_FN_TAG(fn,a,...) \
        RHEA_GLOBAL_INFOF ("<%s " a " />\n", fn, __VA_ARGS__)
#define RHEA_GLOBAL_PRODUCTIONF_FN_BEGIN(fn,a,...) \
        RHEA_GLOBAL_PRODUCTIONF ("<%s " a ">\n", fn, __VA_ARGS__)
#define RHEA_GLOBAL_PRODUCTIONF_FN_END(fn,a,...) \
        RHEA_GLOBAL_PRODUCTIONF ("</%s " a ">\n", fn, __VA_ARGS__)
#define RHEA_GLOBAL_PRODUCTIONF_FN_TAG(fn,a,...) \
        RHEA_GLOBAL_PRODUCTIONF ("<%s " a " />\n", fn, __VA_ARGS__)

/* log macros that are active on every processor */
#define RHEA_TRACE(s)      RHEA_LOG (SC_LP_TRACE, (s))
#define RHEA_LDEBUG(s)     RHEA_LOG (SC_LP_DEBUG, (s))
#define RHEA_VERBOSE(s)    RHEA_LOG (SC_LP_VERBOSE, (s))
#define RHEA_INFO(s)       RHEA_LOG (SC_LP_INFO, (s))
#define RHEA_STATISTICS(s) RHEA_LOG (SC_LP_STATISTICS, (s))
#define RHEA_PRODUCTION(s) RHEA_LOG (SC_LP_PRODUCTION, (s))
#define RHEA_ESSENTIAL(s)  RHEA_LOG (SC_LP_ESSENTIAL, (s))
#define RHEA_LERROR(s)     RHEA_LOG (SC_LP_ERROR, (s))
void                RHEA_TRACEF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_LDEBUGF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_VERBOSEF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_INFOF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_STATISTICSF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_PRODUCTIONF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_ESSENTIALF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                RHEA_LERRORF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
#ifndef __cplusplus
#define RHEA_TRACEF(f,...)                             \
        RHEA_LOGF (SC_LP_TRACE, (f), __VA_ARGS__)
#define RHEA_LDEBUGF(f,...)                            \
        RHEA_LOGF (SC_LP_DEBUG, (f), __VA_ARGS__)
#define RHEA_VERBOSEF(f,...)                           \
        RHEA_LOGF (SC_LP_VERBOSE, (f), __VA_ARGS__)
#define RHEA_INFOF(f,...)                              \
        RHEA_LOGF (SC_LP_INFO, (f), __VA_ARGS__)
#define RHEA_STATISTICSF(f,...)                        \
        RHEA_LOGF (SC_LP_STATISTICS, (f), __VA_ARGS__)
#define RHEA_PRODUCTIONF(f,...)                        \
        RHEA_LOGF (SC_LP_PRODUCTION, (f), __VA_ARGS__)
#define RHEA_ESSENTIALF(f,...)                         \
        RHEA_LOGF (SC_LP_ESSENTIAL, (f), __VA_ARGS__)
#define RHEA_LERRORF(f,...)                            \
        RHEA_LOGF (SC_LP_ERROR, (f), __VA_ARGS__)
#endif
#define RHEA_NOTICE            RHEA_STATISTICS
#define RHEA_NOTICEF           RHEA_STATISTICSF

/* log macros for beginning and end of a function (XML type output) */
#define RHEA_VERBOSE_FN_BEGIN(fn)     RHEA_VERBOSEF ("<%s>\n", fn)
#define RHEA_VERBOSE_FN_END(fn)       RHEA_VERBOSEF ("</%s>\n", fn)
#define RHEA_INFO_FN_BEGIN(fn)        RHEA_INFOF ("<%s>\n", fn)
#define RHEA_INFO_FN_END(fn)          RHEA_INFOF ("</%s>\n", fn)
#define RHEA_PRODUCTION_FN_BEGIN(fn)  RHEA_PRODUCTIONF ("<%s>\n", fn)
#define RHEA_PRODUCTION_FN_END(fn)    RHEA_PRODUCTIONF ("</%s>\n", fn)

#define RHEA_VERBOSEF_FN_BEGIN(fn,a,...) \
        RHEA_VERBOSEF ("<%s " a ">\n", fn, __VA_ARGS__)
#define RHEA_VERBOSEF_FN_END(fn,a,...) \
        RHEA_VERBOSEF ("</%s " a ">\n", fn, __VA_ARGS__)
#define RHEA_INFOF_FN_BEGIN(fn,a,...) \
        RHEA_INFOF ("<%s " a ">\n", fn, __VA_ARGS__)
#define RHEA_INFOF_FN_END(fn,a,...) \
        RHEA_INFOF ("</%s " a ">\n", fn, __VA_ARGS__)
#define RHEA_PRODUCTIONF_FN_BEGIN(fn,a,...) \
        RHEA_PRODUCTIONF ("<%s " a ">\n", fn, __VA_ARGS__)
#define RHEA_PRODUCTIONF_FN_END(fn,a,...) \
        RHEA_PRODUCTIONF ("</%s " a ">\n", fn, __VA_ARGS__)

/* extern declarations */
extern int          rhea_package_id;

/** Registers rhea with the SC Library and sets the logging behavior.
 * This function is optional.
 * If this function is not called or called with log_handler == NULL,
 * the default SC log handler will be used.
 * If this function is not called or called with log_threshold == SC_LP_DEFAULT,
 * the default SC log threshold will be used.
 * The default SC log settings can be changed with sc_set_log_defaults ().
 */
void                rhea_init (sc_log_handler_t log_handler,
                               int log_threshold);

/**
 * Initializes rhea and sub-packages.
 */
void                rhea_initialize (int argc, char **argv,
                                     sc_MPI_Comm mpicomm);

/**
 * Finalizes rhea and sub-packages.
 */
void                rhea_finalize ();

SC_EXTERN_C_END;

#endif /* RHEA_BASE_H */
