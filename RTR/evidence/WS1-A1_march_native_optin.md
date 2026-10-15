# WS1-A1 — make `-march=native` opt-in (reviewer R2.8)

**Status (2026-10-02):** resolved — RTR implementation superseded by Nikol's fix on main (merged)
**Date:** 2026-09-23 (RTR fix); updated 2026-10-02 after merging main
**Reviewer item R2.8 (verbatim):** "On most of the linux and mac distributions the
package will be compiled with -march=native. Actually ARCHFLAGS is never set on
Linux, so -march=native will always be set on Linux. Was that intentional?"

## Reviewer's reading confirmed exactly

`setup.py:42-49` (original code, before either fix):

```python
else:
    extra_compile_args = ["-O3", "-Wall"]
    if not _is_cross_compiling():
        extra_compile_args.append("-march=native")
```

`_is_cross_compiling()` decides solely from `ARCHFLAGS`, which **cibuildwheel sets
only on macOS**. On Linux it is unset, so the function returns `False` for every
ordinary build and `-march=native` was applied unconditionally — precisely what the
reviewer described.

**Why this matters beyond tidiness:** `-march=native` tunes the binary to the
building machine's instruction set. A binary so built is not portable; executing it
on a host without the same extensions can abort with SIGILL. This affects anyone
who builds from source (`pip install .`, the documented route) and any Linux wheel
built outside cibuildwheel's macOS-only ARCHFLAGS path.

## Resolution (current tree)

`-march=native` is **opt-in**, via `HSEEKER_NATIVE_BUILD=1` (`setup.py:42-55`):

- Default (non-Windows) build is `-O3 -Wall` only — portable.
- `-march=native` is appended only when `HSEEKER_NATIVE_BUILD` is exactly `"1"`.
- If opted in while `ARCHFLAGS` targets a different macOS architecture
  (`_is_cross_compiling()`, `setup.py:21-39`), the build **fails** with
  `RuntimeError("HSEEKER_NATIVE_BUILD=1 cannot target a different architecture")`
  rather than silently dropping the flag.
- README "Installation → Optional: CPU-native optimization" documents
  `HSEEKER_NATIVE_BUILD=1 pip install .`, states the result is not portable and must
  not be distributed, and notes the cross-arch error.

**Wheels unaffected:** the cibuildwheel config (`pyproject.toml:63-82`) and
`.github/workflows/` set no `HSEEKER_NATIVE_BUILD`, so CI wheels are built portably.

## Verification

`tests/test_build_flags.py` — 3 tests (4 cases; runs `setup.py` with
`setuptools.setup` captured, no compiler needed):

| Test | Checks |
|---|---|
| `test_default_extension_does_not_require_build_cpu_features[Linux/Darwin]` | default args == `["-O3", "-Wall"]`, no `-march=native` |
| `test_native_optimization_requires_explicit_opt_in` | `HSEEKER_NATIVE_BUILD=1` adds `-march=native` |
| `test_native_optimization_rejects_other_macos_target` | opt-in + `ARCHFLAGS=-arch x86_64 -arch arm64` on arm64 raises `RuntimeError` |

All 4 pass on the merged tree (2026-10-02). Full suite on the merged tree:
**211 passed, 3 skipped**.

## History (superseded RTR implementation)

RTR (2026-09-23) made the flag opt-in via `HSEEKER_NATIVE=1` (`_native_requested()`,
accepting `1/true/yes/on`), applied only when not cross-compiling — i.e. a
cross-compile was silently built without the flag. It was verified with a real
`build_ext` matrix (default 0, `HSEEKER_NATIVE=1` 1, `=0` 0, `=1`+foreign `ARCHFLAGS` 0
occurrences of `-march=native`; suite then 149 passed, 3 skipped). Nikol fixed the same
defect independently on main; on merging, the authors chose Nikol's version
(2026-10-02). It differs only in the variable name, the stricter `"1"` match, and
raising instead of silently ignoring a cross-arch opt-in. `_native_requested()` and
`HSEEKER_NATIVE` no longer exist in the tree.

## Reproduce

```
python -m pytest -q tests/test_build_flags.py                                   # 4 passed
rm -rf build/temp.linux* && touch src/hseeker/_hdna.c
python3 setup.py build_ext --force 2>&1 | grep -c march=native                  # -> 0
HSEEKER_NATIVE_BUILD=1 python3 setup.py build_ext --force 2>&1 | grep -c march=native  # -> 1
```
