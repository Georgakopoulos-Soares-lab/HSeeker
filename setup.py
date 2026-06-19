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

import os
import platform

from setuptools import Extension, setup


def _is_cross_compiling():
    """Return True when cross-compiling (e.g. x86_64 wheel on an Apple M-series host).

    cibuildwheel sets ARCHFLAGS on macOS to declare the target architecture(s),
    e.g. ``-arch x86_64`` or ``-arch arm64 -arch x86_64`` (universal2).
    If any requested target differs from the native machine we are cross-compiling
    and must not pass -march=native to the compiler.
    """
    archflags = os.environ.get("ARCHFLAGS", "")
    if not archflags:
        return False
    parts = archflags.split()
    target_archs = [
        parts[i + 1]
        for i, p in enumerate(parts)
        if p == "-arch" and i + 1 < len(parts)
    ]
    native = platform.machine()  # e.g. "arm64" or "x86_64"
    return any(t != native for t in target_archs)


extra_compile_args = ["-O3", "-Wall"]
if not _is_cross_compiling():
    extra_compile_args.append("-march=native")

hdna_ext = Extension(
    # The module will be importable as  hseeker._hdna
    name="hseeker._hdna",
    sources=["src/hseeker/_hdna.c"],
    extra_compile_args=extra_compile_args,
    libraries=["m"],  # link libm for math.h (no-op on Windows/MSVC)
)

setup(ext_modules=[hdna_ext])
