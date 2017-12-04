
# This file is part of rhea.
# Makefile.am in example/slabweakzone
# Included non-recursively from toplevel directory.

rhea_slabweakzone_headers =
rhea_slabweakzone_sources =
LINT_CSOURCES += $(rhea_slabweakzone_headers) $(rhea_slabweakzone_sources)

## the slabweakzone program
##
## program name:    example/slabweakzone/rhea_slabweakzone
## normalized name: example_slabweakzone_rhea_slabweakzone
rhea_slabweakzone_program_sources = example/slabweakzone/slabweakzone.c
LINT_CSOURCES += $(rhea_slabweakzone_program_sources)

bin_PROGRAMS += example/slabweakzone/rhea_slabweakzone
example_slabweakzone_rhea_slabweakzone_SOURCES = $(rhea_slabweakzone_headers) \
                                       $(rhea_slabweakzone_sources) \
                                       $(rhea_slabweakzone_program_sources)
example_slabweakzone_rhea_slabweakzone_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
                                    -I@top_srcdir@/example/slabweakzone
example_slabweakzone_rhea_slabweakzone_LDADD = $(LDADD)
