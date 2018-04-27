
# This file is part of rhea.
# Makefile.am in example/newton_polynomial
# Included non-recursively from toplevel directory.

rhea_newton_polynomial_headers = \
        example/newton_polynomial/newton_polynomial_base.h
rhea_newton_polynomial_sources = \
        example/newton_polynomial/newton_polynomial_base.c
LINT_CSOURCES += \
        $(rhea_newton_polynomial_headers) \
				$(rhea_newton_polynomial_sources)

## the newton_polynomial program
##
## program name:    example/newton_polynomial/rhea_newton_polynomial
## normalized name: example_newton_polynomial_rhea_newton_polynomial
rhea_newton_polynomial_program_sources = \
        example/newton_polynomial/newton_polynomial.c
LINT_CSOURCES += $(rhea_newton_polynomial_program_sources)

bin_PROGRAMS += example/newton_polynomial/rhea_newton_polynomial
example_newton_polynomial_rhea_newton_polynomial_SOURCES = \
        $(rhea_newton_polynomial_headers) \
        $(rhea_newton_polynomial_sources) \
        $(rhea_newton_polynomial_program_sources)
example_newton_polynomial_rhea_newton_polynomial_CPPFLAGS = \
        $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
        -I@top_srcdir@/example/newton_polynomial
example_newton_polynomial_rhea_newton_polynomial_LDADD = $(LDADD)

## the newton_polynomial meshfree program
##
## program name:    example/newton_polynomial/rhea_newton_polynomial_meshfree
## normalized name: example_newton_polynomial_rhea_newton_polynomial_meshfree
rhea_newton_polynomial_meshfree_program_sources = \
        example/newton_polynomial/newton_polynomial_meshfree.c
LINT_CSOURCES += $(rhea_newton_polynomial_meshfree_program_sources)

bin_PROGRAMS += example/newton_polynomial/rhea_newton_polynomial_meshfree
example_newton_polynomial_rhea_newton_polynomial_meshfree_SOURCES = \
        $(rhea_newton_polynomial_headers) \
        $(rhea_newton_polynomial_sources) \
        $(rhea_newton_polynomial_meshfree_program_sources)
example_newton_polynomial_rhea_newton_polynomial_meshfree_CPPFLAGS = \
        $(AM_CPPFLAGS) $(RHEA_CPPFLAGS) \
        -I@top_srcdir@/example/newton_polynomial
example_newton_polynomial_rhea_newton_polynomial_meshfree_LDADD = $(LDADD)
