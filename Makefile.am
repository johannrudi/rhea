
# This file is part of the rhea library
# Makefile.am in toplevel directory

ACLOCAL_AMFLAGS = -I config @RHEA_SC_AMFLAGS@
if RHEA_SC_MK_USE
@RHEA_SC_MK_INCLUDE@
endif
if RHEA_P4EST_MK_USE
@RHEA_P4EST_MK_INCLUDE@
endif
if RHEA_MANGLL_MK_USE
@RHEA_MANGLL_MK_INCLUDE@
endif
if RHEA_YMIR_MK_USE
@RHEA_YMIR_MK_INCLUDE@
endif

# initialize empty variables
AM_CPPFLAGS =
CLEANFILES =
DISTCLEANFILES =
EXTRA_DIST =
LDADD =
LINT_CSOURCES =
TESTS =
bin_PROGRAMS =
check_PROGRAMS =
include_HEADERS =
lib_LTLIBRARIES =
nodist_include_HEADERS =
noinst_HEADERS =
noinst_PROGRAMS =
sysconf_DATA =
dist_rheadata_DATA =

# use this if you want to link in rhea without autotools
sysconf_DATA += Makefile.rhea.mk
CLEANFILES += Makefile.rhea.mk
Makefile.rhea.mk : Makefile.rhea.pre
	cat $< | \
        sed -e 's,{\(\(.*prefix\|sysconfdir\)\)},{rhea_\1},g' \
            -e 's,^\(\(.*prefix\|sysconfdir\) *=\),rhea_\1,g' > $@

# install rhea m4 macros in the correct directory
rheaaclocaldir = $(datadir)/aclocal
dist_rheaaclocal_DATA = \
				config/rhea_include.m4

# install rhea data files in the correct directory
rheadatadir = $(datadir)/data

# setup test environment
if RHEA_MPIRUN
LOG_COMPILER = @RHEA_MPIRUN@
AM_LOG_FLAGS = @RHEA_MPI_TEST_FLAGS@
endif

# recurse only into subpackages
SUBDIRS = @RHEA_SC_SUBDIR@ @RHEA_P4EST_SUBDIR@ @RHEA_MANGLL_SUBDIR@ \
          @RHEA_YMIR_SUBDIR@
DIST_SUBDIRS = $(SUBDIRS)

# handle toplevel directory
EXTRA_DIST += \
        bootstrap build-aux/git-version-gen build-aux/git2cl
DISTCLEANFILES += \
        _configs.sed src/rhea_config.h @RHEA_DISTCLEAN@
.PHONY: lint ChangeLog

# non-recursive build
include src/Makefile.am
##include test/Makefile.am
##if RHEA_ENABLE_EXAMPLES
include example/share/Makefile.am
include example/basic/Makefile.am
include example/collide/Makefile.am
include example/slabweakzone/Makefile.am
include example/newton_polynomial/Makefile.am
include example/yielding/Makefile.am
include example/slabs/Makefile.am
##endif

# lint static syntax checker
ALL_LINT_FLAGS = $(LINT_FLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
                 $(MPI_INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) \
                 $(src_librhea_a_CPPFLAGS)

lint:
if LINT
##if RHEA_ENABLE_COMPONENT_CHECKS
	@for subdir in $(SUBDIRS) ; do \
		echo "Making $@ in $$subdir"; \
		(cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) lint) ; \
	done
##endif
	for f in $(LINT_CSOURCES) ; do \
		$(LINT) $(ALL_LINT_FLAGS) $(top_srcdir)/$$f || \
		echo "Lint check failed for $$f" ; \
	done
else
	@echo "Static source code check disabled by configure"
endif

##if !RHEA_ENABLE_COMPONENT_CHECKS
##check-recursive:
##	$(MAKE) $(AM_MAKEFLAGS) $(check_PROGRAMS)
##	$(MAKE) $(AM_MAKEFLAGS) check-TESTS
##endif

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
if RHEA_DIST_DENY
	@echo "-----------------------------------------------------"
	@echo "make dist does not work with external libsc"
	@echo "-----------------------------------------------------"
	@exit 1
endif

clean-local:
	rm -f ChangeLog
