
# This file is part of rhea.
# Makefile.am in example/earth_inversion
# Included non-recursively from toplevel directory.

rhea_earth_inversion_headers =
rhea_earth_inversion_sources =
LINT_CSOURCES += $(rhea_earth_inversion_headers) $(rhea_earth_inversion_sources)

## the earth inversion program
##
## program name:    example/earth_inversion/rhea_earth_inversion
## normalized name: example_earth_inversion_rhea_earth_inversion
rhea_earth_inversion_program_sources = \
        example/earth_inversion/earth_inversion.c
LINT_CSOURCES += $(rhea_earth_inversion_program_sources)

bin_PROGRAMS += example/earth_inversion/rhea_earth_inversion
example_earth_inversion_rhea_earth_inversion_SOURCES = \
        $(rhea_earth_inversion_headers) \
        $(rhea_earth_inversion_sources) \
        $(rhea_earth_inversion_program_sources)
example_earth_inversion_rhea_earth_inversion_CPPFLAGS = \
        $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
        -I@top_srcdir@/example/earth_inversion
example_earth_inversion_rhea_earth_inversion_LDADD = $(LDADD)
