
# This file is part of rhea.
# Makefile.am in example/adjoint
# Included non-recursively from toplevel directory.

rhea_adjoint_headers = \
      example/adjoint/adjoint_math.h \
      example/adjoint/adjoint_options.h \
      example/adjoint/adjoint_physics.h \
      example/adjoint/adjoint_test.h \
      example/adjoint/adjoint_TI.h \
      example/adjoint/adjoint_postp.h \
      example/adjoint/adjoint_adjoint.h \
      example/adjoint/adjoint_vtk.h \
      example/adjoint/adjoint_io.h \
      example/adjoint/adjoint_geometry.h \
      example/adjoint/adjoint_essential.h \
      example/adjoint/adjoint.h
rhea_adjoint_sources = \
      example/adjoint/adjoint_math.c \
      example/adjoint/adjoint_options.c \
      example/adjoint/adjoint_physics.c \
      example/adjoint/adjoint_test.c \
      example/adjoint/adjoint_TI.c \
      example/adjoint/adjoint_postp.c \
      example/adjoint/adjoint_adjoint.c \
      example/adjoint/adjoint_vtk.c \
      example/adjoint/adjoint_io.c \
      example/adjoint/adjoint_geometry.c \
      example/adjoint/adjoint_essential.c
LINT_CSOURCES += $(rhea_adjoint_headers) $(rhea_adjoint_sources)

## the adjoint program
##
## program name:    example/adjoint/rhea_adjoint
## normalized name: example_adjoint_rhea_adjoint
rhea_adjoint_program_sources = example/adjoint/adjoint.c
LINT_CSOURCES += $(rhea_adjoint_program_sources)

bin_PROGRAMS += example/adjoint/rhea_adjoint
example_adjoint_rhea_adjoint_SOURCES = $(rhea_adjoint_headers) \
                                       $(rhea_adjoint_sources) \
                                       $(rhea_adjoint_program_sources)
example_adjoint_rhea_adjoint_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
                                    -I@top_srcdir@/example/adjoint
example_adjoint_rhea_adjoint_LDADD = $(LDADD)
