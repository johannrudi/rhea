
# This file is part of rhea.
# Makefile.am in example/basic
# Included non-recursively from toplevel directory.

rhea_basic_headers =
rhea_basic_sources =
LINT_CSOURCES += $(rhea_basic_headers) $(rhea_basic_sources)

## the basic program
##
## program name:    example/basic/rhea_basic
## normalized name: example_basic_rhea_basic
rhea_basic_program_sources = example/basic/basic.c
LINT_CSOURCES += $(rhea_basic_program_sources)

bin_PROGRAMS += example/basic/rhea_basic
example_basic_rhea_basic_SOURCES = $(rhea_basic_headers) \
                                   $(rhea_basic_sources) \
                                   $(rhea_basic_program_sources)
example_basic_rhea_basic_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
                                    -I@top_srcdir@/example/basic
example_basic_rhea_basic_LDADD = $(LDADD)
