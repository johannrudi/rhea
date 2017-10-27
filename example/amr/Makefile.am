
# This file is part of rhea.
# Makefile.am in example/amr
# Included non-recursively from toplevel directory.

rhea_amr_headers =
rhea_amr_sources =
LINT_CSOURCES += $(rhea_amr_headers) $(rhea_amr_sources)

## the amr program
##
## program name:    example/amr/rhea_amr
## normalized name: example_amr_rhea_amr
rhea_amr_program_sources = example/amr/amr.c
LINT_CSOURCES += $(rhea_amr_program_sources)

bin_PROGRAMS += example/amr/rhea_amr
example_amr_rhea_amr_SOURCES = $(rhea_amr_headers) \
														   $(rhea_amr_sources) \
															 $(rhea_amr_program_sources)
example_amr_rhea_amr_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
                                -I@top_srcdir@/example/amr
example_amr_rhea_amr_LDADD = $(LDADD)
