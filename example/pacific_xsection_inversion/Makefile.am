
# This file is part of rhea.
# Makefile.am in example/pacific_xsection_inversion
# Included non-recursively from toplevel directory.

rhea_pacific_xsection_inversion_headers =
rhea_pacific_xsection_inversion_sources =
LINT_CSOURCES += $(rhea_pacific_xsection_inversion_headers) \
                 $(rhea_pacific_xsection_inversion_sources)

## the pacific cross section inversion program
##
## program name:
##   example/pacific_xsection_inversion/rhea_pacific_xsection_inversion
## normalized name:
##   example_pacific_xsection_inversion_rhea_pacific_xsection_inversion
rhea_pacific_xsection_inversion_program_sources = \
        example/pacific_xsection_inversion/pacific_xsection_inversion.c
LINT_CSOURCES += $(rhea_pacific_xsection_inversion_program_sources)

bin_PROGRAMS += example/pacific_xsection_inversion/rhea_pacific_xsection_inversion
example_pacific_xsection_inversion_rhea_pacific_xsection_inversion_SOURCES = \
        $(rhea_pacific_xsection_inversion_headers) \
        $(rhea_pacific_xsection_inversion_sources) \
        $(rhea_pacific_xsection_inversion_program_sources)
example_pacific_xsection_inversion_rhea_pacific_xsection_inversion_CPPFLAGS = \
        $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
        -I@top_srcdir@/example/pacific_xsection_inversion
example_pacific_xsection_inversion_rhea_pacific_xsection_inversion_LDADD = \
        $(LDADD)
