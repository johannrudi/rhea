
# This file is part of the ymir library
# Makefile.am in toplevel directory

ACLOCAL_AMFLAGS = -I config @YMIR_SC_AMFLAGS@
if YMIR_SC_MK_USE
@YMIR_SC_MK_INCLUDE@
endif
if YMIR_P4EST_MK_USE
@YMIR_P4EST_MK_INCLUDE@
endif
if YMIR_MANGLL_MK_USE
@YMIR_MANGLL_MK_INCLUDE@
endif

# initialize empty variables
AM_CPPFLAGS =
CLEANFILES =
DISTCLEANFILES =
EXTRA_DIST =
LDADD =
LINT_CSOURCES =
TESTS =
bin_SCRIPTS =
bin_PROGRAMS =
check_PROGRAMS =
include_HEADERS =
lib_LTLIBRARIES =
nodist_include_HEADERS =
noinst_HEADERS =
noinst_PROGRAMS =
sysconf_DATA =
#TODO variables below may need update (see mangll/p4est Makefile.am)
optfilesdir = $(datadir)/optfiles
dist_optfiles_DATA =

# use this if you want to link in ymir without autotools
sysconf_DATA += Makefile.ymir.mk
CLEANFILES += Makefile.ymir.mk
Makefile.ymir.mk : Makefile.ymir.pre
	cat $< | \
        sed -e 's,{\(\(.*prefix\|sysconfdir\)\)},{ymir_\1},g' \
            -e 's,^\(\(.*prefix\|sysconfdir\) *=\),ymir_\1,g' > $@

# install ymir m4 macros in the correct directory
ymiraclocaldir = $(datadir)/aclocal
dist_ymiraclocal_DATA = \
				config/ymir_include.m4 config/ymir_netcdf.m4 config/ymir_petsc.m4

# install ymir data files in the correct directory
ymirdatadir = $(datadir)/data

# setup test environment
if YMIR_MPIRUN
LOG_COMPILER = @YMIR_MPIRUN@
AM_LOG_FLAGS = @YMIR_MPI_TEST_FLAGS@
endif

# recurse only into subpackages
SUBDIRS = @YMIR_SC_SUBDIR@ @YMIR_P4EST_SUBDIR@ @YMIR_MANGLL_SUBDIR@
DIST_SUBDIRS = $(SUBDIRS)

# handle toplevel directory
EXTRA_DIST += \
        bootstrap ymirindent build-aux/git-version-gen build-aux/git2cl \
        ymir_albmap_jagjan12.p8c ymir_albmap_noholes.p8c
DISTCLEANFILES += \
        _configs.sed src/ymir_config.h @YMIR_DISTCLEAN@
.PHONY: lint ChangeLog

# non-recursive build
include src/Makefile.am
include test/Makefile.am
if YMIR_ENABLE_EXAMPLES
#include example/antarctica/Makefile.am
#include example/figscripts/Makefile.am
include example/gmg/Makefile.am
#include example/ismip/Makefile.am
#include example/meshsheet/Makefile.am
#include example/performance/Makefile.am
include example/sinker/Makefile.am
include example/slabs/Makefile.am
#include example/squaresheet/Makefile.am
#include example/splinesheet/Makefile.am
#include example/stiffness/Makefile.am
#include example/stress/Makefile.am
#include example/tv/Makefile.am
endif

# lint static syntax checker
ALL_LINT_FLAGS = $(LINT_FLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
                 $(MPI_INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) \
                 $(src_libymir_a_CPPFLAGS)

lint:
if LINT
if YMIR_ENABLE_COMPONENT_CHECKS
	@for subdir in $(SUBDIRS) ; do \
		echo "Making $@ in $$subdir"; \
		(cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) lint) ; \
	done
endif
	for f in $(LINT_CSOURCES) ; do \
		$(LINT) $(ALL_LINT_FLAGS) $(top_srcdir)/$$f || \
		echo "Lint check failed for $$f" ; \
	done
else
	@echo "Static source code check disabled by configure"
endif

if !YMIR_ENABLE_COMPONENT_CHECKS
check-recursive:
	$(MAKE) $(AM_MAKEFLAGS) $(check_PROGRAMS)
	$(MAKE) $(AM_MAKEFLAGS) check-TESTS
endif

# revision control and ChangeLog
ChangeLog:
	(GIT_DIR=@top_srcdir@/.git git log > .ChangeLog.tmp && \
         cat .ChangeLog.tmp | @top_srcdir@/build-aux/git2cl > ChangeLog) ; \
        rm -f .ChangeLog.tmp

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
	test "x$(VERSION)" = "x`@top_srcdir@/build-aux/git-version-gen\
              @top_srcdir@/.tarball-version`" || \
        ((echo "Stale version;" ; echo "Please run:" ; \
          echo "     cd @top_srcdir@ && ./bootstrap" ; \
          echo "before make dist") 1>&2 ; rm -r $(distdir) ; exit 1)
if YMIR_DIST_DENY
	@echo "-----------------------------------------------------"
	@echo "make dist does not work with external libsc"
	@echo "-----------------------------------------------------"
	@exit 1
endif

clean-local:
	rm -f ChangeLog
