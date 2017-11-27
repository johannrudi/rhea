
# This file is part of rhea.
# Makefile.am in example/io
# Included non-recursively from toplevel directory.

rhea_io_headers =
rhea_io_sources =
LINT_CSOURCES += $(rhea_io_headers) $(rhea_io_sources)

## the io program
##
## program name:    example/io/rhea_io
## normalized name: example_io_rhea_io
rhea_io_program_sources = example/io/io.c
LINT_CSOURCES += $(rhea_io_program_sources)

bin_PROGRAMS += example/io/rhea_io
example_io_rhea_io_SOURCES = $(rhea_io_headers) \
                             $(rhea_io_sources) \
                             $(rhea_io_program_sources)
example_io_rhea_io_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
                              -I@top_srcdir@/example/io
example_io_rhea_io_LDADD = $(LDADD)
