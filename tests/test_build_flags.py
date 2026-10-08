"""Source distributions must not bake the builder's CPU features into wheels."""

import platform
import runpy
from pathlib import Path

import pytest
import setuptools


SETUP_PY = Path(__file__).resolve().parents[1] / "setup.py"


def extension_for(monkeypatch, system, *, native=False, archflags=None):
    captured = []
    monkeypatch.setattr(setuptools, "setup", lambda **kwargs: captured.append(kwargs))
    monkeypatch.setattr(platform, "system", lambda: system)
    monkeypatch.setattr(platform, "machine", lambda: "arm64")
    monkeypatch.delenv("HSEEKER_NATIVE_BUILD", raising=False)
    monkeypatch.delenv("ARCHFLAGS", raising=False)
    if native:
        monkeypatch.setenv("HSEEKER_NATIVE_BUILD", "1")
    if archflags is not None:
        monkeypatch.setenv("ARCHFLAGS", archflags)
    runpy.run_path(str(SETUP_PY))
    return captured[0]["ext_modules"][0]


@pytest.mark.parametrize("system", ["Linux", "Darwin"])
def test_default_extension_does_not_require_build_cpu_features(monkeypatch, system):
    ext = extension_for(monkeypatch, system)
    assert "-march=native" not in ext.extra_compile_args
    assert ext.extra_compile_args == ["-O3", "-Wall"]


def test_native_optimization_requires_explicit_opt_in(monkeypatch):
    ext = extension_for(monkeypatch, "Linux", native=True)
    assert "-march=native" in ext.extra_compile_args


def test_native_optimization_rejects_other_macos_target(monkeypatch):
    with pytest.raises(RuntimeError, match="different architecture"):
        extension_for(monkeypatch, "Darwin", native=True,
                      archflags="-arch x86_64 -arch arm64")
