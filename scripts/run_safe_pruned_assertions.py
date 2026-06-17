#!/usr/bin/env python3
"""Run SafePrunedHSeeker synthetic assertions without requiring pytest."""

from __future__ import annotations

import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


if __name__ == "__main__":
    runpy.run_path(str(ROOT / "tests" / "test_safe_pruned.py"), run_name="__main__")
