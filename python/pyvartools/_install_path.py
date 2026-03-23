# This file is overwritten by 'make install' (or 'make install-python') to
# record the path where the vartools binary was installed.  It is used as the
# third entry in the binary-discovery chain:
#   1. config file
#   2. VARTOOLS_BINARY environment variable
#   3. this file  (set at install time)
#   4. PATH
#
# If you are running from the source tree without installing, leave this
# unchanged and rely on VARTOOLS_BINARY or PATH instead.
VARTOOLS_INSTALL_PATH = ""
