
# This file is part of the rhea library
# Makefile.am in src
# included non-recursively from toplevel directory

librhea_generated_headers = src/rhea_config.h
librhea_installed_headers = \
        src/rhea_base.h \
        src/rhea_domain.h \
        src/rhea_discretization.h \
        src/rhea_amr.h \
        src/rhea_pointcloud.hpp \
        src/rhea_pointcloud_adaptor.h \
        src/rhea_temperature.h \
        src/rhea_weakzone.h \
        src/rhea_viscosity.h \
        src/rhea_velocity.h \
        src/rhea_pressure.h \
        src/rhea_velocity_pressure.h \
        src/rhea_strainrate.h \
        src/rhea_stress.h \
        src/rhea_stokes_norm.h \
        src/rhea_stokes_problem.h \
        src/rhea_stokes_problem_amr.h \
        src/rhea_newton.h \
        src/rhea_newton_check.h \
        src/rhea_io_std.h \
        src/rhea_io_mpi.h \
        src/rhea_vtk.h \
        src/rhea_vis.h \
        src/rhea_vis_adaptor.h \
        src/rhea.h
librhea_internal_headers =
librhea_compiled_sources = \
        src/rhea_base.c \
        src/rhea_domain.c \
        src/rhea_discretization.c \
        src/rhea_amr.c \
        src/rhea_pointcloud.cpp \
        src/rhea_pointcloud_adaptor.cxx \
        src/rhea_temperature.c \
        src/rhea_weakzone.c \
        src/rhea_viscosity.c \
        src/rhea_velocity.c \
        src/rhea_pressure.c \
        src/rhea_velocity_pressure.c \
        src/rhea_strainrate.c \
        src/rhea_stress.c \
        src/rhea_stokes_norm.c \
        src/rhea_stokes_problem.c \
        src/rhea_stokes_problem_amr.c \
        src/rhea_newton.c \
        src/rhea_newton_check.c \
        src/rhea_io_std.c \
        src/rhea_io_mpi.c \
        src/rhea_vtk.c \
        src/rhea_vis.c \
        src/rhea_vis_adaptor.cxx \
        src/rhea.c
librhea_original_headers = \
        nanoflann/include/nanoflann.hpp

RHEA_CPPFLAGS =
RHEA_LDADD = @top_builddir@/src/librhea.la @RHEA_YMIR_LDADD@ \
             @RHEA_MANGLL_LDADD@ @RHEA_P4EST_LDADD@ @RHEA_SC_LDADD@
if RHEA_WITH_PETSC
RHEA_CPPFLAGS += @RHEA_PETSC_INCLUDES@
RHEA_LDADD += @RHEA_PETSC_LINK_LIBS@
endif

lib_LTLIBRARIES += src/librhea.la
src_librhea_la_SOURCES = \
        $(librhea_internal_headers) \
        $(librhea_compiled_sources)
src_librhea_la_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS)
LDADD += $(RHEA_LDADD)
EXTRA_src_librhea_la_DEPENDENCIES = @RHEA_YMIR_LDADD@ @RHEA_MANGLL_LDADD@ \
                                    @RHEA_P4EST_LDADD@ @RHEA_SC_LDADD@

nodist_include_HEADERS += $(librhea_generated_headers) \
                          $(librhea_original_headers)
include_HEADERS += $(librhea_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src @RHEA_YMIR_CPPFLAGS@ @RHEA_MANGLL_CPPFLAGS@ \
               @RHEA_P4EST_CPPFLAGS@ @RHEA_SC_CPPFLAGS@

LINT_CSOURCES += $(librhea_compiled_sources)
