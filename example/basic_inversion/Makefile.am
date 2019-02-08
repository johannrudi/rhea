
# This file is part of rhea.
# Makefile.am in example/basic_inversion
# Included non-recursively from toplevel directory.

rhea_basic_inversion_headers =
rhea_basic_inversion_sources =
LINT_CSOURCES += $(rhea_basic_inversion_headers) $(rhea_basic_inversion_sources)

## the basic inversion program
##
## program name:    example/basic_inversion/rhea_basic_inversion
## normalized name: example_basic_inversion_rhea_basic_inversion
rhea_basic_inversion_program_sources = example/basic_inversion/basic_inversion.c
LINT_CSOURCES += $(rhea_basic_inversion_program_sources)

bin_PROGRAMS += example/basic_inversion/rhea_basic_inversion
example_basic_inversion_rhea_basic_inversion_SOURCES = \
        $(rhea_basic_inversion_headers) \
        $(rhea_basic_inversion_sources) \
        $(rhea_basic_inversion_program_sources)
example_basic_inversion_rhea_basic_inversion_CPPFLAGS = \
        $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
        -I@top_srcdir@/example/basic_inversion
example_basic_inversion_rhea_basic_inversion_LDADD = $(LDADD)
