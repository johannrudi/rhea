
# This file is part of rhea.
# Makefile.am in example/share
# Included non-recursively from toplevel directory.

rhea_example_share_internal_headers = \
      example/share/example_share_mesh.h \
      example/share/example_share_stokes.h \
      example/share/example_share_vtk.h
rhea_example_share_compiled_sources = \
      example/share/example_share_mesh.c \
      example/share/example_share_stokes.c \
      example/share/example_share_vtk.c

librhea_internal_headers += $(rhea_example_share_internal_headers)
librhea_compiled_sources += $(rhea_example_share_compiled_sources)
RHEA_CPPFLAGS += -I@top_srcdir@/example/share
