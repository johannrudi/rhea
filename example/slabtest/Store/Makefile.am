
# This file is part of rhea.
# Makefile.am in example/slabtest
# Included non-recursively from toplevel directory.

rhea_slabtest_headers = \
     example/slabtest/slabtest_options.h \
		 example/slabtest/slabtest_stress.h
rhea_slabtest_sources = \
     example/slabtest/slabtest_options.c \
		 example/slabtest/slabtest_stress.c
LINT_CSOURCES += $(rhea_slabtest_headers) $(rhea_slabtest_sources)

## slabtest program
##
## program name:    example/slabtest/rhea_slabtest
## normalized name: example_slabtest_rhea_slabtest
rhea_slabtest_program_sources = example/slabtest/slabtest.c
LINT_CSOURCES += $(rhea_slabtest_program_sources)

bin_PROGRAMS += example/slabtest/rhea_slabtest
example_slabtest_rhea_slabtest_SOURCES = $(rhea_slabtest_headers) \
																	 $(rhea_slabtest_sources) \
																	 $(rhea_slabtest_program_sources)
example_slabtest_rhea_slabtest_CPPFLAGS = $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
                                    -I@top_srcdir@/example/slabtest
example_slabtest_rhea_slabtest_LDADD = $(LDADD)
