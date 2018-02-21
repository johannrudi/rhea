
# This file is part of rhea.
# Makefile.am in example/earth
# Included non-recursively from toplevel directory.

rhea_earth_headers =
rhea_earth_sources =
LINT_CSOURCES += $(rhea_earth_headers) $(rhea_earth_sources)

## the earth program
##
## program name:    example/earth/rhea_earth
## normalized name: example_earth_rhea_earth
rhea_earth_program_sources = example/earth/earth.c
LINT_CSOURCES += $(rhea_earth_program_sources)

bin_PROGRAMS += example/earth/rhea_earth
example_earth_rhea_earth_SOURCES = $(rhea_earth_headers) \
														       $(rhea_earth_sources) \
															     $(rhea_earth_program_sources)
example_earth_rhea_earth_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
                                    -I@top_srcdir@/example/earth
example_earth_rhea_earth_LDADD = $(LDADD)
