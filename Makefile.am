SUBDIRS = libltdl src USERLIBS/src USERFUNCS/src
dist_doc_DATA = README USERLIBS/ReadME.USERLIBS USERFUNCS/ReadME.USERFUNCS
include examplefiles.mk
nobase_dist_pkgdata_DATA = $(example_files)

# Include the Python package in 'make dist'.
# We use dist-hook rather than EXTRA_DIST so that generated artifacts
# (__pycache__, *.egg-info, *.pyc) are excluded from the tarball.
dist-hook:
	cp -r $(srcdir)/python $(distdir)/python
	find $(distdir)/python -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
	find $(distdir)/python -name '*.egg-info' -type d -exec rm -rf {} + 2>/dev/null || true
	find $(distdir)/python -name '.pytest_cache' -type d -exec rm -rf {} + 2>/dev/null || true
	find $(distdir)/python -name '.hypothesis' -type d -exec rm -rf {} + 2>/dev/null || true
	find $(distdir)/python -name '*.pyc' -delete 2>/dev/null || true

.PHONY: dist-hook

# Python package installation
# Run 'make install-python' after 'make install' to install pyvartools.
# This records the vartools install path into the Python package so that
# pyvartools can find the binary without requiring it to be on PATH.
install-python:
	echo "VARTOOLS_INSTALL_PATH = '$(DESTDIR)$(bindir)/vartools'" \
	    > $(srcdir)/python/pyvartools/_install_path.py
	cd $(srcdir)/python && pip install .

uninstall-python:
	pip uninstall -y pyvartools

.PHONY: install-python uninstall-python

# Sanitizer-instrumented test pass.  This is intentionally NOT part of the
# regular `make check` flow because an ASan/UBSan build is 5-10x slower than
# a normal one — a full unittest.sh run takes ~5 minutes here vs. ~30 seconds
# without instrumentation.  Run it before tagging a release, or after editing
# memory-sensitive code (parsers, allocators, numerical loops).  See CLAUDE.md
# for usage notes and how to interpret the results.
#
# AM_CFLAGS / AM_LDFLAGS overrides preserve the configure-time CFLAGS while
# layering sanitizer flags on top.  halt_on_error=0 lets unittest.sh keep
# running past the first ASan/UBSan hit so a single run surfaces all issues.
ASAN_FLAGS = -fsanitize=address,undefined -fno-omit-frame-pointer \
             -fno-sanitize-recover=undefined -g -O1
ASAN_REPORT_PREFIX = $(top_builddir)/asan-report

check-asan:
	@echo "==> [check-asan] Rebuilding with -fsanitize=address,undefined ..."
	$(MAKE) clean
	$(MAKE) -k AM_CFLAGS="$(ASAN_FLAGS)" AM_LDFLAGS="$(ASAN_FLAGS)"
	@rm -f $(ASAN_REPORT_PREFIX).*
	@echo "==> [check-asan] Running unittest.sh under ASan/UBSan ..."
	@echo "    sanitizer reports (if any) will be written to $(ASAN_REPORT_PREFIX).<pid>"
	ASAN_OPTIONS="halt_on_error=0:abort_on_error=0:detect_leaks=0:print_stacktrace=1:log_path=$(ASAN_REPORT_PREFIX)" \
	UBSAN_OPTIONS="halt_on_error=0:print_stacktrace=1:log_path=$(ASAN_REPORT_PREFIX)" \
	LD_LIBRARY_PATH="$(top_builddir)/src/.libs:$$LD_LIBRARY_PATH" \
	bash $(top_srcdir)/unittest.sh $(top_builddir)/src/.libs/vartools \
	    || echo "    (unittest.sh exited non-zero — sanitizer reports below take precedence)"
	@hits=`ls $(ASAN_REPORT_PREFIX).* 2>/dev/null | wc -l`; \
	if test $$hits -gt 0; then \
	    echo ""; \
	    echo "==> [check-asan] FAIL: $$hits sanitizer report(s) written:"; \
	    ls $(ASAN_REPORT_PREFIX).*; \
	    echo "    Inspect each file for the offending stack trace."; \
	    exit 1; \
	else \
	    echo ""; \
	    echo "==> [check-asan] PASS: no ASan/UBSan reports."; \
	fi

.PHONY: check-asan
