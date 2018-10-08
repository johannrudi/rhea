
# This file is part of rhea.
# Makefile.am in example/subduction
# Included non-recursively from toplevel directory.

rhea_subduction_headers = \
      example/subduction/subduction_options.h \
      example/subduction/subduction_physics.h \
      example/subduction/subduction_test.h \
      example/subduction/subduction_TI.h \
      example/subduction/subduction_postp.h \
      example/subduction/subduction_adjoint.h \
      example/subduction/subduction_vtk.h \
      example/subduction/subduction_io.h \
      example/subduction/subduction_geometry.h \
      example/subduction/subduction_essential.h \
      example/subduction/subduction.h
rhea_subduction_sources = \
      example/subduction/subduction_options.c \
      example/subduction/subduction_physics.c \
      example/subduction/subduction_test.c \
      example/subduction/subduction_TI.c \
      example/subduction/subduction_postp.c \
      example/subduction/subduction_adjoint.c \
      example/subduction/subduction_vtk.c \
      example/subduction/subduction_io.c \
      example/subduction/subduction_geometry.c \
      example/subduction/subduction_essential.c
LINT_CSOURCES += $(rhea_subduction_headers) $(rhea_subduction_sources)

## the subduction program
##
## program name:    example/subduction/rhea_subduction
## normalized name: example_subduction_rhea_subduction
rhea_subduction_program_sources = example/subduction/subduction.c
LINT_CSOURCES += $(rhea_subduction_program_sources)

bin_PROGRAMS += example/subduction/rhea_subduction
example_subduction_rhea_subduction_SOURCES = $(rhea_subduction_headers) \
                                       $(rhea_subduction_sources) \
                                       $(rhea_subduction_program_sources)
example_subduction_rhea_subduction_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
                                    -I@top_srcdir@/example/subduction
example_subduction_rhea_subduction_LDADD = $(LDADD)
