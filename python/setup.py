# Minimal setup.py shim for compatibility with older pip versions (<21.3)
# that do not support editable installs from pyproject.toml alone.
# All package metadata is in pyproject.toml.
from setuptools import setup
setup()
