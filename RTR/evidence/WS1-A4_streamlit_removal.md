# WS1-A4 — remove the Streamlit front-end entirely (reviewer R2.11)

**Status:** done, verified (re-checked 2026-10-02 on the merged tree)
**Date:** 2026-09-23

**Merge note (2026-10-02):** Nikol's later edit to `src/hseeker/app.py` on main (an
overlap-strategy dropdown) was dropped along with the file during the merge of main
into RTR; Streamlit stays removed. Re-checked on the merged tree: no `streamlit`
reference outside `RTR/` (README, pyproject, src, webapp all clean), no `app.py` /
`streamlit_app.py` in `src/hseeker/`, `[app]` extra = `pandas`, `plotly`, Dockerfile
lines 15/42 unchanged, `webapp/main.py` parses; suite 211 passed, 3 skipped.
**Reviewer item R2.11 (verbatim):** "It seems there is some Streamlit demo, while
the app stands on FastAPI, giving rise to two app.py files. This might be
misleading. Apart from this, I really acknowledge the IT part of the project — the
care for the speed of the code as well as putting the package on PyPi and as a
standalone on-line tool. This might really help the user!"

## Resolution changed by author decision

The plan (WS1-A4) originally specified *renaming* `src/hseeker/app.py` to
`streamlit_app.py`. On 2026-09-23 the authors directed that **Streamlit be removed
from the project entirely**. That is the stronger answer to R2.11: rather than
making two front-ends distinguishable, there is now only one.

The rename was performed first and then superseded within the same session, so the
git history shows `app.py` → `streamlit_app.py` → deleted.

`RTR/HSeeker_Reviewer_Answer_Plan_v2.md` was updated in both places that described
the rename (the WS1-A4 row and the R2.11 draft response), because leaving them would
have put a false claim — "demo renamed" — into the response letter.

## What was removed

| Item | Action |
|---|---|
| `src/hseeker/app.py` (→ `streamlit_app.py`) | deleted |
| `streamlit>=1.30` in the `[app]` extra, `pyproject.toml:32` | removed; extra now `pandas`, `plotly` |
| README "Web interfaces" two-front-end table | replaced by a single **"Web interface"** section documenting only the FastAPI service |
| `.github/copilot-instructions.md:614` | `[app]` extra description corrected to "`pandas` and `plotly`" |

The `[app]` extra itself is retained: `pandas` and `plotly` are genuine
dependencies of the FastAPI webapp and the Dockerfile's first stage installs
`".[app]"`.

## Verification

**1. No Streamlit reference survives anywhere in the shipped tree** (excluding the
`RTR/` revision archive, which quotes the reviewer verbatim and must not be
altered):

```
grep -rn -i "streamlit" . --exclude-dir=.git --exclude-dir=RTR --exclude-dir=*egg-info
  NONE — clean
```

**2. The package builds and imports with Streamlit uninstalled:**

```
hseeker imports OK
hseeker.streamlit_app gone: True
hseeker.app gone         : True
[app] extra deps present : pandas, plotly OK
streamlit installed?     : False
```

`pip install -e ".[app]"` resolves cleanly with the dependency dropped.

**3. The remaining front-end is intact:**

```
webapp/main.py parses OK
Dockerfile:15  RUN pip install --no-cache-dir ".[app]"
Dockerfile:42  CMD uvicorn webapp.main:app --host 0.0.0.0 --port "${PORT:-8000}"
```

The deployment path is unchanged: the Docker image served the FastAPI app before
and still does.

**4. Full suite: 150 passed, 3 skipped.**

Note the web front-end is not exercised by the test suite (it was not before this
change either), so verification here is import resolution, dependency resolution,
syntax and the absence of stale references — not a live HTTP check.

## Files changed

- `src/hseeker/app.py` → `streamlit_app.py` → **deleted**
- `pyproject.toml` (`streamlit` dropped from the `[app]` extra)
- `README.md` (§4 "Web interface", single FastAPI entry)
- `.github/copilot-instructions.md` (`[app]` extra description)
- `RTR/HSeeker_Reviewer_Answer_Plan_v2.md` (WS1-A4 row + R2.11 draft response
  corrected from "renamed" to "removed")
