# WS1-A1 — make `-march=native` opt-in (reviewer R2.8)

**Status:** done, verified
**Date:** 2026-09-23
**Reviewer item R2.8 (verbatim):** "On most of the linux and mac distributions the
package will be compiled with -march=native. Actually ARCHFLAGS is never set on
Linux, so -march=native will always be set on Linux. Was that intentional?"

## Reviewer's reading confirmed exactly

`setup.py:42-49` (before the fix):

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

## Fix

`-march=native` is now **opt-in**, via `HSEEKER_NATIVE=1`:

- Added `_native_requested()` reading `HSEEKER_NATIVE` (accepts `1/true/yes/on`,
  case-insensitive).
- Applied only when opted in **and** not cross-compiling — the cross-compile guard
  is retained as a second safety net rather than replaced.
- `_is_cross_compiling()` docstring now records that it cannot serve as the opt-in,
  because ARCHFLAGS is unset on Linux — so the original bug cannot be reintroduced
  by someone reading only that function.
- README "Installation → From source" documents the flag and states plainly that
  the result is not portable and must not be distributed.

## Verification — build matrix

Each row is a real `setup.py build_ext --force` with the C source touched, counting
`-march=native` in the emitted compiler command:

| Build environment | occurrences | expected |
|---|---|---|
| *(default — no env)* | **0** | 0 ✓ |
| `HSEEKER_NATIVE=1` | **1** | 1 ✓ |
| `HSEEKER_NATIVE=0` | **0** | 0 ✓ |
| `HSEEKER_NATIVE=1` + `ARCHFLAGS=-arch ppc64` | **0** | 0 ✓ (cross-compile guard wins) |

Default compile line, confirming a portable build:

```
gcc ... -DNDEBUG -g -O3 -Wall -fPIC -I<python> -c src/hseeker/_hdna.c -o ... -O3 -Wall
```

Opted-in compile line:

```
gcc ... -c src/hseeker/_hdna.c -o ... -O3 -Wall -march=native
```

Note the flag is appended **after** the source filename — an earlier grep of mine
truncated at `_hdna.c` and appeared to show the opt-in failing. It was the grep that
was wrong, not the build; re-verified with the full line above.

**Wheels unaffected:** cibuildwheel config (`pyproject.toml:64-74`) sets no
`HSEEKER_NATIVE`, so CI wheels are now built portably — which is the correct and
previously violated behaviour on Linux.

**Test suite after restoring a default (portable) build: 149 passed, 3 skipped.**

## Files changed

- `setup.py` (+`_native_requested()`, opt-in condition, docstring note)
- `README.md` (Installation → "Optional: CPU-native optimization")

## Reproduce

```
rm -rf build/temp.linux* && touch src/hseeker/_hdna.c
python3 setup.py build_ext --force 2>&1 | grep -c march=native        # -> 0
HSEEKER_NATIVE=1 python3 setup.py build_ext --force 2>&1 | grep -c march=native  # -> 1
```
