
# This file is part of the rhea library
# Makefile.am in src
# included non-recursively from toplevel directory

librhea_generated_headers = src/rhea_config.h
librhea_installed_headers = \
        src/rhea.h
librhea_internal_headers =
librhea_compiled_sources =

RHEA_CPPFLAGS =
RHEA_LDADD = @top_builddir@/src/librhea.la @RHEA_YMIR_LDADD@ \
             @RHEA_MANGLL_LDADD@ @RHEA_P4EST_LDADD@ @RHEA_SC_LDADD@
if RHEA_WITH_PETSC
RHEA_CPPFLAGS += @RHEA_PETSC_INCLUDES@
RHEA_LDADD += @RHEA_PETSC_LINK_LIBS@
endif

lib_LTLIBRARIES += src/librhea.la
src_librhea_la_SOURCES = \
        $(librhea_internal_headers) \
        $(librhea_compiled_sources)
src_librhea_la_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS)
LDADD += $(RHEA_LDADD)
EXTRA_src_librhea_la_DEPENDENCIES = @RHEA_YMIR_LDADD@ @RHEA_MANGLL_LDADD@ \
                                    @RHEA_P4EST_LDADD@ @RHEA_SC_LDADD@

nodist_include_HEADERS += $(librhea_generated_headers)
include_HEADERS += $(librhea_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src @RHEA_YMIR_CPPFLAGS@ @RHEA_MANGLL_CPPFLAGS@ \
               @RHEA_P4EST_CPPFLAGS@ @RHEA_SC_CPPFLAGS@

LINT_CSOURCES += $(librhea_compiled_sources)
