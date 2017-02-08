
# This file is part of rhea.
# Makefile.am in example/yielding
# Included non-recursively from toplevel directory.

rhea_yielding_headers =
rhea_yielding_sources =
LINT_CSOURCES += $(rhea_yielding_headers) $(rhea_yielding_sources)

## the yielding program
##
## program name:    example/yielding/rhea_yielding
## normalized name: example_yielding_rhea_yielding
rhea_yielding_program_sources = example/yielding/yielding.c
LINT_CSOURCES += $(rhea_yielding_program_sources)

bin_PROGRAMS += example/yielding/rhea_yielding
example_yielding_rhea_yielding_SOURCES = $(rhea_yielding_headers) \
																				 $(rhea_yielding_sources) \
																				 $(rhea_yielding_program_sources)
example_yielding_rhea_yielding_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
                                    			-I@top_srcdir@/example/yielding
example_yielding_rhea_yielding_LDADD = $(LDADD)
