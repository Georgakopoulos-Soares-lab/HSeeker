"""
webapp/main.py — HSeeker Web Application
=========================================
Deploy:  uvicorn main:app --host 0.0.0.0 --port $PORT
"""
from __future__ import annotations

import csv
import json
import os
import tempfile
import threading
import time
import uuid
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

import hseeker

# ── App setup ─────────────────────────────────────────────────────────────────

BASE_DIR = Path(__file__).parent
app = FastAPI(title="HSeeker Web", docs_url=None, redoc_url=None)
app.mount("/static", StaticFiles(directory=BASE_DIR / "static"), name="static")
templates = Jinja2Templates(directory=BASE_DIR / "templates")
templates.env.globals["hseeker_version"] = hseeker.__version__

MAX_UPLOAD_MB = int(os.getenv("MAX_UPLOAD_MB", "200"))
MAX_UPLOAD_BYTES = MAX_UPLOAD_MB * 1024 * 1024
JOB_TTL_SECONDS = int(os.getenv("JOB_TTL_SECONDS", "7200"))  # 2 h

# ── Job registry ──────────────────────────────────────────────────────────────

_jobs: dict[str, dict[str, Any]] = {}
_jobs_lock = threading.Lock()


def _new_job(filename: str, params: dict) -> str:
    job_id = str(uuid.uuid4())
    with _jobs_lock:
        _jobs[job_id] = {
            "id": job_id,
            "status": "pending",
            "created": time.time(),
            "filename": filename,
            "params": params,
            "processed_records": 0,
            "hits": [],
            "stats": {},
            "charts": {},
            "seq_lengths": {},
            "tsv_path": None,
            "error": None,
        }
    return job_id


def _get_job(job_id: str) -> dict:
    with _jobs_lock:
        job = _jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if time.time() - job["created"] > JOB_TTL_SECONDS:
        with _jobs_lock:
            _jobs.pop(job_id, None)
        raise HTTPException(status_code=410, detail="Results expired (older than 2 h)")
    return job


# ── Plotly palette & helpers ──────────────────────────────────────────────────

_BG      = "#0d1117"
_CARD    = "#161b22"
_GRID    = "#21262d"
_TEAL    = "#39d5ab"
_PURPLE  = "#9b8dff"
_BLUE    = "#79c0ff"
_ORANGE  = "#f0883e"
_GREEN   = "#3fb950"
_TEXT    = "#c9d1d9"
_SUBTEXT = "#7d8590"


def _base_layout(**extra) -> dict:
    axis_defaults = dict(
        gridcolor=_GRID,
        zerolinecolor=_GRID,
        tickcolor=_SUBTEXT,
        linecolor=_GRID,
    )
    return dict(
        paper_bgcolor=_CARD,
        plot_bgcolor=_CARD,
        font=dict(
            family="Inter, system-ui, -apple-system, sans-serif",
            color=_TEXT,
            size=12,
        ),
        margin=dict(l=55, r=20, t=46, b=46),
        xaxis=axis_defaults,
        yaxis=axis_defaults,
        colorway=[_TEAL, _PURPLE, _BLUE, _ORANGE, _GREEN],
        **extra,
    )


def _title(text: str) -> dict:
    return dict(
        text=text,
        font=dict(size=13, color="#e6edf3"),
        x=0.02,
        xanchor="left",
        pad=dict(b=6),
    )


def _to_dict(fig: go.Figure) -> dict:
    """Convert Figure to a plain dict via to_json (handles NaN/Inf safely)."""
    return json.loads(fig.to_json())


# ── Chart generation ──────────────────────────────────────────────────────────

def generate_charts(
    hits: list[dict],
    seq_lengths: dict[str, int],
    params: dict,
) -> dict[str, dict]:
    """Return a dict of Plotly figure dicts, keyed by chart name.

    All histograms are pre-aggregated server-side so the JSON embedded in the
    page stays small regardless of hit count.
    """
    if not hits:
        return {}

    df = pd.DataFrame(hits)
    charts: dict[str, dict] = {}

    # ── 1. Hits per sequence ──────────────────────────────────────────────────
    try:
        counts = df["seq_id"].value_counts().head(25).reset_index()
        counts.columns = ["seq_id", "n"]
        fig = go.Figure(go.Bar(
            x=counts["n"].tolist(),
            y=counts["seq_id"].tolist(),
            orientation="h",
            marker=dict(
                color=counts["n"].tolist(),
                colorscale=[[0, _PURPLE], [1, _TEAL]],
                showscale=False,
                line_width=0,
            ),
            text=[f"{v:,}" for v in counts["n"]],
            textposition="outside",
            textfont=dict(size=10, color=_TEXT),
            hovertemplate="%{y}: <b>%{x:,}</b> sites<extra></extra>",
        ))
        fig.update_layout(
            **_base_layout(height=max(280, len(counts) * 30 + 80)),
            title=_title("H-DNA Sites per Sequence"),
            xaxis_title="Number of H-DNA Sites",
            bargap=0.22,
        )
        fig.update_yaxes(autorange="reversed")
        charts["hits_per_seq"] = _to_dict(fig)
    except Exception:
        pass

    # ── 2. Perfect vs Imperfect donut ─────────────────────────────────────────
    try:
        n_perfect = int(df["is_perfect"].sum())
        n_imp = len(df) - n_perfect
        fig = go.Figure(go.Pie(
            labels=["Perfect", "Imperfect"],
            values=[n_perfect, n_imp],
            hole=0.65,
            marker=dict(colors=[_TEAL, _PURPLE], line=dict(color=_BG, width=2)),
            textinfo="label+percent",
            textfont=dict(size=12, color="#000000"),
            hovertemplate="%{label}: %{value:,} (%{percent})<extra></extra>",
        ))
        fig.add_annotation(
            text=f"<b>{n_perfect:,}</b><br><span style='font-size:10px'>perfect</span>",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=15, color=_TEAL), align="center",
        )
        fig.update_layout(
            **_base_layout(height=310, showlegend=True),
            title=_title("Perfect vs Imperfect H-DNA"),
            legend=dict(orientation="h", y=-0.06, x=0.5, xanchor="center",
                        font=dict(size=11)),
        )
        fig.update_layout(margin=dict(l=20, r=20, t=46, b=46))
        fig.update_layout(plot_bgcolor="#000000")
        charts["perfect_donut"] = _to_dict(fig)
    except Exception:
        pass

    # ── 3. Arm length distribution ────────────────────────────────────────────
    try:
        arm_c = df["arm_length"].value_counts().sort_index()
        fig = go.Figure(go.Bar(
            x=arm_c.index.tolist(),
            y=arm_c.values.tolist(),
            marker=dict(
                color=arm_c.values.tolist(),
                colorscale=[[0, _PURPLE], [1, _TEAL]],
                showscale=False,
                line_width=0,
            ),
            hovertemplate="arm <b>%{x} bp</b>: %{y:,} sites<extra></extra>",
        ))
        fig.update_layout(
            **_base_layout(height=310),
            title=_title("Arm Length Distribution"),
            xaxis_title="Arm Length (bp)",
            yaxis_title="Count",
            bargap=0.1,
        )
        charts["arm_length_dist"] = _to_dict(fig)
    except Exception:
        pass

    # ── 4. Spacer length distribution ─────────────────────────────────────────
    try:
        sp_c = df["spacer_length"].value_counts().sort_index()
        fig = go.Figure(go.Bar(
            x=[str(v) for v in sp_c.index.tolist()],
            y=sp_c.values.tolist(),
            marker=dict(color=_PURPLE, line_width=0),
            hovertemplate="spacer <b>%{x} bp</b>: %{y:,} sites<extra></extra>",
        ))
        fig.update_layout(
            **_base_layout(height=310),
            title=_title("Spacer (Loop) Length Distribution"),
            xaxis_title="Spacer Length (bp)",
            yaxis_title="Count",
            bargap=0.2,
        )
        charts["spacer_dist"] = _to_dict(fig)
    except Exception:
        pass

    # ── 5. Mirror identity distribution ──────────────────────────────────────
    try:
        mir_hist, edges = np.histogram(
            df["mirror_identity"].values, bins=50, range=(0.0, 100.0)
        )
        centers = ((edges[:-1] + edges[1:]) / 2).tolist()
        min_mirror = (1.0 - params.get("mismatch", 0.10)) * 100.0
        fig = go.Figure(go.Bar(
            x=centers,
            y=mir_hist.tolist(),
            marker=dict(color=_BLUE, line_width=0),
            hovertemplate="mirror <b>%{x:.1f}%</b>: %{y:,} sites<extra></extra>",
        ))
        fig.add_vline(
            x=min_mirror,
            line_dash="dash", line_color=_ORANGE, line_width=1.5,
            annotation_text=f"threshold ({min_mirror:.0f}%)",
            annotation_font=dict(color=_ORANGE, size=10),
            annotation_position="top right",
        )
        fig.update_layout(
            **_base_layout(height=310),
            title=_title("Mirror Identity Distribution"),
            xaxis_title="Mirror Identity (%)",
            yaxis_title="Count",
            bargap=0.05,
        )
        charts["mirror_dist"] = _to_dict(fig)
    except Exception:
        pass

    # ── 6. GA% vs CT% composition scatter ────────────────────────────────────
    try:
        n_sample = min(len(df), 4000)
        sample = df.sample(n_sample, random_state=42) if len(df) > n_sample else df
        fig = go.Figure()
        for is_perf, color, name, opacity in [
            (True,  _TEAL,   "Perfect",   0.80),
            (False, _PURPLE, "Imperfect", 0.50),
        ]:
            mask = sample["is_perfect"] == is_perf
            if mask.any():
                s = sample[mask]
                fig.add_trace(go.Scatter(
                    x=s["ga_pct"].tolist(),
                    y=s["ct_pct"].tolist(),
                    mode="markers",
                    name=name,
                    marker=dict(color=color, size=4, opacity=opacity, line_width=0),
                    hovertemplate=(
                        f"<b>{name}</b><br>GA=%{{x:.1f}}%  CT=%{{y:.1f}}%"
                        "<extra></extra>"
                    ),
                ))
        # guide line GA + CT = 100
        fig.add_trace(go.Scatter(
            x=[0, 100], y=[100, 0],
            mode="lines",
            line=dict(color="#30363d", width=1, dash="dot"),
            showlegend=False, hoverinfo="skip",
        ))
        fig.update_layout(
            **_base_layout(height=375),
            title=_title(
                f"Sequence Composition: GA% vs CT% — right arm"
                + (f" (n = {n_sample:,} sampled)" if len(df) > n_sample else "")
            ),
            legend=dict(x=0.01, y=0.99, bgcolor="rgba(0,0,0,0)"),
        )
        fig.update_xaxes(title_text="GA% (Right Arm)", range=[-2, 105])
        fig.update_yaxes(title_text="CT% (Right Arm)", range=[-2, 105])
        charts["composition_scatter"] = _to_dict(fig)
    except Exception:
        pass

    # ── 7. Arm × Spacer 2-D density ──────────────────────────────────────────
    try:
        arm_v = df["arm_length"].values
        sp_v  = df["spacer_length"].values
        arm_bins = np.arange(int(arm_v.min()), int(arm_v.max()) + 2)
        sp_bins  = np.arange(int(sp_v.min()),  int(sp_v.max())  + 2)
        h2d, xedges, yedges = np.histogram2d(arm_v, sp_v, bins=[arm_bins, sp_bins])
        fig = go.Figure(go.Heatmap(
            z=h2d.T.tolist(),
            x=arm_bins[:-1].tolist(),
            y=sp_bins[:-1].tolist(),
            colorscale=[[0, _CARD], [0.4, _PURPLE], [1, _TEAL]],
            hovertemplate=(
                "arm=%{x} bp  spacer=%{y} bp<br>count=%{z:,}<extra></extra>"
            ),
            colorbar=dict(title="Count", tickfont=dict(size=10, color=_TEXT)),
        ))
        fig.update_layout(
            **_base_layout(height=360),
            title=_title("Arm Length × Spacer Length Density"),
            xaxis_title="Arm Length (bp)",
            yaxis_title="Spacer Length (bp)",
        )
        charts["arm_spacer_heatmap"] = _to_dict(fig)
    except Exception:
        pass

    # ── 8. Score distribution ──────────────────────────────────────────────────
    try:
        scores = df["total_score"].dropna().values
        if len(scores) > 0:
            score_hist, sedges = np.histogram(scores, bins=50)
            scenters = ((sedges[:-1] + sedges[1:]) / 2).tolist()
            fig = go.Figure(go.Bar(
                x=scenters,
                y=score_hist.tolist(),
                marker=dict(color=_TEAL, line_width=0),
                hovertemplate="score <b>%{{x:.1f}}</b>: %{{y:,}} sites<extra></extra>",
            ))
            fig.add_vline(
                x=60.0,
                line_dash="dash", line_color=_ORANGE, line_width=1.5,
                annotation_text="stable ≥ 60",
                annotation_font=dict(color=_ORANGE, size=10),
                annotation_position="top right",
            )
            fig.update_layout(
                **_base_layout(height=310),
                title=_title("Total Score Distribution"),
                xaxis_title="Total Score",
                yaxis_title="Count",
                bargap=0.05,
            )
            charts["score_dist"] = _to_dict(fig)
    except Exception:
        pass

    # ── 9. Score vs Mirror Identity scatter ────────────────────────────────────
    try:
        df_sc = df.dropna(subset=["total_score"])
        if len(df_sc) > 0:
            n_sample = min(len(df_sc), 4000)
            sample_sc = df_sc.sample(n_sample, random_state=42) if len(df_sc) > n_sample else df_sc
            fig = go.Figure(go.Scatter(
                x=sample_sc["mirror_identity"].tolist(),
                y=sample_sc["total_score"].tolist(),
                mode="markers",
                marker=dict(
                    color=_PURPLE,
                    size=3,
                    opacity=0.55,
                    line_width=0,
                ),
                hovertemplate=(
                    "mirror=%{{x:.1f}}%  score=%{{y:.1f}}<br>"
                    "<extra></extra>"
                ),
            ))
            fig.update_layout(
                **_base_layout(height=375),
                title=_title(
                    f"Score vs Mirror Identity"
                    + (f" (n = {n_sample:,} sampled)" if len(df_sc) > n_sample else "")
                ),
                xaxis_title="Mirror Identity (%)",
                yaxis_title="Total Score",
            )
            charts["score_vs_mirror"] = _to_dict(fig)
    except Exception:
        pass

    # ── 10. Locus map (top sequences) ────────────────────────────────────────
    try:
        top_seqs = df["seq_id"].value_counts().head(15).index.tolist()
        arm_min = int(df["arm_length"].min())
        arm_max = int(df["arm_length"].max())
        fig = go.Figure()

        for seq_id in top_seqs:
            sub = df[df["seq_id"] == seq_id].copy()
            seq_len = seq_lengths.get(seq_id, int(sub["end"].max()))

            # sequence backbone
            fig.add_trace(go.Scatter(
                x=[1, seq_len], y=[seq_id, seq_id],
                mode="lines",
                line=dict(color=_GRID, width=7),
                showlegend=False, hoverinfo="skip",
            ))

            # sample if dense
            display = sub.sample(600, random_state=42) if len(sub) > 600 else sub
            midpoints = ((display["start"] + display["end"]) / 2).tolist()
            _perf_span = "<span style='color:#39d5ab'>✓ perfect</span>"
            hover_texts = [
                f"<b>{r['seq_id']}</b><br>"
                f"pos {r['start']:,} \u2013 {r['end']:,}<br>"
                f"arm = {r['arm_length']} bp \u00a0 spacer = {r['spacer_length']} bp<br>"
                f"GA = {r['ga_pct']:.1f}% \u00a0 CT = {r['ct_pct']:.1f}%<br>"
                f"mirror = {r['mirror_identity']:.1f}%<br>"
                + (f"score = {r['total_score']:.1f}<br>" if r.get('total_score') is not None else "")
                + (_perf_span if r['is_perfect'] else 'imperfect')
                for _, r in display.iterrows()
            ]
            marker_sizes = (display["arm_length"].clip(upper=30) * 0.65).tolist()

            fig.add_trace(go.Scatter(
                x=midpoints,
                y=[seq_id] * len(display),
                mode="markers",
                name=seq_id,
                showlegend=False,
                marker=dict(
                    size=marker_sizes,
                    color=display["arm_length"].tolist(),
                    colorscale=[[0, _PURPLE], [1, _TEAL]],
                    cmin=arm_min, cmax=arm_max,
                    opacity=0.75,
                    line_width=0,
                    colorbar=dict(
                        title="Arm (bp)",
                        thickness=12,
                        len=0.6,
                        tickfont=dict(size=10, color=_TEXT),
                    ) if seq_id == top_seqs[0] else None,
                ),
                text=hover_texts,
                hovertemplate="%{text}<extra></extra>",
            ))

        dense_note = ""
        if any(len(df[df["seq_id"] == s]) > 600 for s in top_seqs):
            dense_note = " (600 sites sampled per sequence)"

        fig.update_layout(
            **_base_layout(
                height=max(320, min(700, len(top_seqs) * 50 + 90)),
                hovermode="closest",
            ),
            title=_title(
                f"H-DNA Locus Map — top {len(top_seqs)} sequences{dense_note}"
            ),
            xaxis_title="Genomic Position (bp)",
            yaxis_title="",
        )
        charts["locus_map"] = _to_dict(fig)
    except Exception:
        pass

    return charts


# ── Stats ─────────────────────────────────────────────────────────────────────

def compute_stats(hits: list[dict]) -> dict:
    empty = dict(
        total_hits=0, n_perfect=0, n_imperfect=0, pct_perfect=0.0,
        n_sequences=0, mean_arm_length=0.0, mean_spacer_length=0.0,
        mean_mirror_identity=0.0, max_arm_length=0, min_arm_length=0,
        mean_total_score=0.0, max_total_score=0.0, n_scored=0,
    )
    if not hits:
        return empty
    df = pd.DataFrame(hits)
    n_perfect = int(df["is_perfect"].sum())
    # Scoring stats — handle hits where scoring is None (disabled or failed)
    scores = [h.get("total_score") for h in hits if h.get("total_score") is not None]
    scored_count = len(scores)
    mean_score = round(float(np.mean(scores)), 1) if scores else 0.0
    max_score = round(float(np.max(scores)), 1) if scores else 0.0
    n_stable = sum(1 for s in scores if s >= 60)
    return dict(
        total_hits=len(hits),
        n_perfect=n_perfect,
        n_imperfect=len(hits) - n_perfect,
        pct_perfect=round(n_perfect / len(hits) * 100, 1),
        n_sequences=int(df["seq_id"].nunique()),
        mean_arm_length=round(float(df["arm_length"].mean()), 1),
        mean_spacer_length=round(float(df["spacer_length"].mean()), 1),
        mean_mirror_identity=round(float(df["mirror_identity"].mean()), 1),
        max_arm_length=int(df["arm_length"].max()),
        min_arm_length=int(df["arm_length"].min()),
        mean_total_score=mean_score,
        max_total_score=max_score,
        n_scored=scored_count,
        n_stable=n_stable,
    )


# ── TSV writer ────────────────────────────────────────────────────────────────

_TSV_FIELDS = [
    "seq_id", "source", "start", "end",
    "arm_length", "spacer_length", "total_length",
    "ga_pct", "ct_pct", "mirror_identity", "is_perfect",
    "left_arm", "spacer", "right_arm", "full_sequence",
    "stacking_score", "pairing_score", "total_score",
    "putative_triplex",
]


def _write_tsv(hits: list[dict], path: str) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=_TSV_FIELDS, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(hits)


# ── Background job processor ──────────────────────────────────────────────────

def _run_job(job_id: str, fasta_path: str) -> None:
    job = _jobs[job_id]
    try:
        job["status"] = "running"
        params = job["params"]
        all_hits: list[dict] = []
        seq_lengths: dict[str, int] = {}
        processed = 0

        for seq_id, seq, offset in hseeker.parse_fasta(fasta_path):
            processed += 1
            job["processed_records"] = processed
            seq_lengths[seq_id] = len(seq)

            hits = hseeker.scan_sequence(
                seq,
                minrep=params["minrep"],
                maxrep=params["maxrep"],
                maxspacer=params["maxspacer"],
                purity=params["purity"],
                mismatch=params["mismatch"],
                remove_overlaps=params["remove_overlaps"],
                seq_offset=offset,
                score=params.get("score", True),
            )
            for h in hits:
                h["seq_id"] = seq_id
                h["source"] = "findHDNA"
            all_hits.extend(hits)

        job["seq_lengths"] = seq_lengths
        job["hits"] = all_hits
        job["stats"] = compute_stats(all_hits)
        job["charts"] = generate_charts(all_hits, seq_lengths, params)

        fd, tsv_path = tempfile.mkstemp(
            suffix=".tsv", prefix=f"hseeker_{job_id[:8]}_"
        )
        os.close(fd)
        _write_tsv(all_hits, tsv_path)
        job["tsv_path"] = tsv_path
        job["status"] = "complete"

    except Exception as exc:
        job["status"] = "failed"
        job["error"] = str(exc)

    finally:
        Path(fasta_path).unlink(missing_ok=True)


# ── Routes ────────────────────────────────────────────────────────────────────

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse(request, "index.html")


@app.get("/about", response_class=HTMLResponse)
async def about(request: Request):
    return templates.TemplateResponse(request, "about.html")


@app.post("/submit")
async def submit(
    request: Request,
    file: Optional[UploadFile] = File(None),
    seq_text: str = Form(""),
    minrep: int = Form(10),
    maxrep: int = Form(1000),
    maxspacer: int = Form(10),
    purity: float = Form(0.90),
    mismatch: float = Form(0.10),
    remove_overlaps: int = Form(1),
    score: int = Form(1),
):
    has_file = file is not None and bool(file.filename)
    has_text = bool(seq_text.strip())
    if not has_file and not has_text:
        raise HTTPException(400, "Provide a FASTA file or paste FASTA text")

    fd, tmp_path = tempfile.mkstemp(suffix=".fa", prefix="hseeker_upload_")

    if has_text:
        # Write pasted text directly to temp file
        try:
            with os.fdopen(fd, "w") as fh:
                fh.write(seq_text)
        except Exception as exc:
            Path(tmp_path).unlink(missing_ok=True)
            raise HTTPException(500, f"Write error: {exc}") from exc
        filename = "pasted_sequence.fasta"
    else:
        suffix = Path(file.filename).suffix.lower()
        if suffix not in {".fa", ".fasta", ".fna", ".fas", ".txt"}:
            os.close(fd)
            Path(tmp_path).unlink(missing_ok=True)
            raise HTTPException(
                400, f"Unsupported file type '{suffix}'. Use .fa / .fasta / .fna"
            )
        # Stream to temp file with size guard
        total_bytes = 0
        try:
            with os.fdopen(fd, "wb") as fh:
                while True:
                    chunk = await file.read(1 << 20)  # 1 MB chunks
                    if not chunk:
                        break
                    total_bytes += len(chunk)
                    if total_bytes > MAX_UPLOAD_BYTES:
                        Path(tmp_path).unlink(missing_ok=True)
                        raise HTTPException(
                            413, f"File exceeds {MAX_UPLOAD_MB} MB limit"
                        )
                    fh.write(chunk)
        except HTTPException:
            raise
        except Exception as exc:
            Path(tmp_path).unlink(missing_ok=True)
            raise HTTPException(500, f"Upload error: {exc}") from exc
        filename = file.filename

    params = dict(
        minrep=minrep,
        maxrep=maxrep,
        maxspacer=maxspacer,
        purity=purity,
        mismatch=mismatch,
        remove_overlaps=bool(remove_overlaps),
        score=bool(score),
    )
    job_id = _new_job(filename, params)
    threading.Thread(target=_run_job, args=(job_id, tmp_path), daemon=True).start()
    return RedirectResponse(f"/job/{job_id}", status_code=303)


@app.get("/job/{job_id}", response_class=HTMLResponse)
async def job_page(request: Request, job_id: str):
    job = _get_job(job_id)
    return templates.TemplateResponse(request, "job.html", {
        "job": job,
        "charts_json": json.dumps(job.get("charts", {})),
        "stats": job.get("stats", {}),
    })


@app.get("/api/job/{job_id}/status")
async def job_status(job_id: str):
    job = _get_job(job_id)
    return {
        "status": job["status"],
        "processed_records": job["processed_records"],
        "total_hits": len(job["hits"]),
        "error": job["error"],
    }


@app.get("/api/job/{job_id}/hits")
async def job_hits(
    job_id: str,
    page: int = 1,
    per_page: int = 50,
    q: str = "",
    min_arm: int = 0,
    max_arm: int = 9999,
    min_mirror: float = 0.0,
    perfect_only: bool = False,
    seq_type: str = "",
    min_score: float = 0.0,
):
    job = _get_job(job_id)
    if job["status"] != "complete":
        raise HTTPException(400, "Job not complete yet")
    hits = job["hits"]
    if q:
        ql = q.lower()
        hits = [h for h in hits if ql in h.get("seq_id", "").lower()]
    if min_arm > 0:
        hits = [h for h in hits if h["arm_length"] >= min_arm]
    if max_arm < 9999:
        hits = [h for h in hits if h["arm_length"] <= max_arm]
    if min_mirror > 0:
        hits = [h for h in hits if h["mirror_identity"] >= min_mirror]
    if perfect_only:
        hits = [h for h in hits if h["is_perfect"]]
    if seq_type == "ga":
        hits = [h for h in hits if h["ga_pct"] >= 80]
    elif seq_type == "ct":
        hits = [h for h in hits if h["ct_pct"] >= 80]
    if min_score > 0:
        hits = [h for h in hits if h.get("total_score") is not None and h["total_score"] >= min_score]
    total = len(hits)
    start = (page - 1) * per_page
    return {
        "total": total,
        "page": page,
        "per_page": per_page,
        "pages": max(1, (total + per_page - 1) // per_page),
        "hits": hits[start: start + per_page],
    }


@app.get("/job/{job_id}/download")
async def download_tsv(job_id: str):
    job = _get_job(job_id)
    if job["status"] != "complete":
        raise HTTPException(400, "Job not complete yet")
    tsv_path = job.get("tsv_path")
    if not tsv_path or not Path(tsv_path).exists():
        raise HTTPException(404, "TSV not available")
    stem = Path(job["filename"]).stem
    return FileResponse(
        tsv_path,
        filename=f"{stem}_HDNA.tsv",
        media_type="text/tab-separated-values",
    )


@app.get("/health")
async def health():
    return {"status": "ok", "version": hseeker.__version__}


if __name__ == "__main__":
    import uvicorn
    port = int(os.getenv("PORT", "8000"))
    uvicorn.run("main:app", host="0.0.0.0", port=port, reload=False)
