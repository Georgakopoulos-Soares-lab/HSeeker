"""
setup.py — declares the _hdna CPython C extension.

This file is intentionally minimal.  All project metadata lives in
pyproject.toml.  This file exists only to register the C extension so
that setuptools knows to compile it.

Build commands:
    pip install -e .                  # editable install (development)
    pip install .                     # regular install
    python -m build                   # build sdist + wheel
    cibuildwheel                      # build platform wheels for PyPI
"""

from setuptools import Extension, setup

hdna_ext = Extension(
    # The module will be importable as  hdna_hunter._hdna
    name="hdna_hunter._hdna",
    sources=["src/hdna_hunter/_hdna.c"],
    extra_compile_args=[
        "-O2",    # optimise — same flag used by the original Makefile
        "-Wall",  # keep warnings on to catch regressions
    ],
    libraries=["m"],  # link libm for math.h (no-op on Windows/MSVC)
)

setup(ext_modules=[hdna_ext])
