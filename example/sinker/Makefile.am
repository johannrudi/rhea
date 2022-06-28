
# This file is part of rhea.
# Makefile.am in example/sinker
# Included non-recursively from toplevel directory.

rhea_sinker_headers =
rhea_sinker_sources =
LINT_CSOURCES += $(rhea_sinker_headers) $(rhea_sinker_sources)

## the sinker program
##
## program name:    example/sinker/rhea_sinker
## normalized name: example_sinker_rhea_sinker
rhea_sinker_program_sources = example/sinker/sinker.c
LINT_CSOURCES += $(rhea_sinker_program_sources)

bin_PROGRAMS += example/sinker/rhea_sinker
example_sinker_rhea_sinker_SOURCES = $(rhea_sinker_headers) \
                                     $(rhea_sinker_sources) \
                                     $(rhea_sinker_program_sources)
example_sinker_rhea_sinker_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
                                      -I@top_srcdir@/example/sinker
example_sinker_rhea_sinker_LDADD = $(LDADD)
