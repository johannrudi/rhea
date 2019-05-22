
# This file is part of rhea.
# Makefile.am in example/sphere_inversion
# Included non-recursively from toplevel directory.

rhea_sphere_inversion_headers =
rhea_sphere_inversion_sources =
LINT_CSOURCES += $(rhea_sphere_inversion_headers) $(rhea_sphere_inversion_sources)

## the sphere inversion program
##
## program name:    example/sphere_inversion/rhea_sphere_inversion
## normalized name: example_sphere_inversion_rhea_sphere_inversion
rhea_sphere_inversion_program_sources = \
        example/sphere_inversion/sphere_inversion.c
LINT_CSOURCES += $(rhea_sphere_inversion_program_sources)

bin_PROGRAMS += example/sphere_inversion/rhea_sphere_inversion
example_sphere_inversion_rhea_sphere_inversion_SOURCES = \
        $(rhea_sphere_inversion_headers) \
        $(rhea_sphere_inversion_sources) \
        $(rhea_sphere_inversion_program_sources)
example_sphere_inversion_rhea_sphere_inversion_CPPFLAGS = \
        $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
        -I@top_srcdir@/example/sphere_inversion
example_sphere_inversion_rhea_sphere_inversion_LDADD = $(LDADD)
