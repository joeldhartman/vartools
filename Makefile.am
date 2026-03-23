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
