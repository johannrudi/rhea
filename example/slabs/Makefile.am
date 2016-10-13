
# This file is part of rhea.
# Makefile.am in example/slabs
# Included non-recursively from toplevel directory.

rhea_slabs_headers = \
      example/slabs/nanoflann.hh \
      example/slabs/weakzone.hh \
      example/slabs/weakzone_wrapper.h \
      example/slabs/slabs_base.h \
      example/slabs/slabs_norm.h \
      example/slabs/slabs_stokes_state.h \
      example/slabs/slabs_physics.h \
      example/slabs/slabs_physics_extended.h \
      example/slabs/slabs_discretization.h \
      example/slabs/slabs_discretization_extended.h \
      example/slabs/slabs_linear_stokes_problem.h \
      example/slabs/slabs_nonlinear_stokes_problem.h \
      example/slabs/slabs_io.h \
      example/slabs/slabs_vtk.h \
      example/slabs/slabs_setup.h
rhea_slabs_sources = \
      example/slabs/weakzone.cc \
      example/slabs/weakzone_wrapper.cc \
      example/slabs/slabs_base.c \
      example/slabs/slabs_norm.c \
      example/slabs/slabs_stokes_state.c \
      example/slabs/slabs_physics.c \
      example/slabs/slabs_discretization.c \
      example/slabs/slabs_linear_stokes_problem.c \
      example/slabs/slabs_nonlinear_stokes_problem.c \
      example/slabs/slabs_io.c \
      example/slabs/slabs_vtk.c \
      example/slabs/slabs_setup.c
LINT_CSOURCES += $(rhea_slabs_headers) $(rhea_slabs_sources)

## slabs program
##
## program name:    example/slabs/rhea_slabs
## normalized name: example_slabs_rhea_slabs
rhea_slabs_program_sources = example/slabs/slabs.c
LINT_CSOURCES += $(rhea_slabs_program_sources)

bin_PROGRAMS += example/slabs/rhea_slabs
example_slabs_rhea_slabs_SOURCES = $(rhea_slabs_headers) \
																	 $(rhea_slabs_sources) \
																	 $(rhea_slabs_program_sources)
example_slabs_rhea_slabs_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
                                    -I@top_srcdir@/example/slabs
example_slabs_rhea_slabs_LDADD = $(LDADD)
