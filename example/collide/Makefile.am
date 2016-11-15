
# This file is part of rhea.
# Makefile.am in example/collide
# Included non-recursively from toplevel directory.

rhea_collide_headers =
rhea_collide_sources =
LINT_CSOURCES += $(rhea_collide_headers) $(rhea_collide_sources)

## the collide program
##
## program name:    example/collide/rhea_collide
## normalized name: example_collide_rhea_collide
rhea_collide_program_sources = example/collide/collide.c
LINT_CSOURCES += $(rhea_collide_program_sources)

bin_PROGRAMS += example/collide/rhea_collide
example_collide_rhea_collide_SOURCES = $(rhea_collide_headers) \
                                       $(rhea_collide_sources) \
                                       $(rhea_collide_program_sources)
example_collide_rhea_collide_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
                                    -I@top_srcdir@/example/collide
example_collide_rhea_collide_LDADD = $(LDADD)
