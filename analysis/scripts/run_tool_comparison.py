#!/usr/bin/env python3
"""HSeeker vs triplex, Triplexator, NeSSie and nBMST on the frozen benchmark (R2.1, R2.7, R3.1, R3.6).

Implements RTR/TOOL_COMPARISON_PLAN_2026-10-11.md (§2-§5, with the changes in §11), which was committed
before any comparator was run. Every tool reads the same FASTA; a record is called forming when the tool
reports at least one qualifying hit in it. Settings per tool:
  A  documented defaults (primary comparison)
  B  geometry matched to HSeeker
  C  small grid, selected out of fold (5-fold x 20, leave-one-family-out, leave-one-study-out)
  manuscript / exploratory  reference runs, not used for claims
Wording rule (§2.4): "better" only if the family-cluster CI of the paired MCC difference under A excludes 0
and the difference has the same sign under B and C (leave-one-family-out); otherwise "comparable".

Tools are expected in $TOOLS_DIR (default $WORK/tools); R triplex in the conda env given by
$R_TRIPLEX_RSCRIPT. Raw outputs are cached per command.

Usage:  python3 analysis/scripts/run_tool_comparison.py [--tools triplex,triplexator,nessie,nbmst]
Output: analysis/results/tool_comparison/
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hseeker  # noqa: E402
import run_grouped_validation as gv  # noqa: E402
import run_parameter_robustness as rb  # noqa: E402

OUT = ROOT / "analysis" / "results" / "tool_comparison"
RAW = OUT / "raw"
TOOLS = Path(os.environ.get("TOOLS_DIR", Path(os.environ.get("WORK", "")) / "tools"))
RSCRIPT = os.environ.get("R_TRIPLEX_RSCRIPT",
                         "/home1/10899/kimopro/WORK/miniconda3/envs/r_triplex/bin/Rscript")
TRIPLEXATOR = TOOLS / "triplexator" / "bin" / "triplexator"
NESSIE = TOOLS / "nessie" / "nessie"
GFA = TOOLS / "non-B_gfa" / "gfa"
SEED = 20261011
N_BOOT = 2000
HOMOPOLYMERS = [("HDNA0053", "A" * 20, "non_forming"), ("HDNA0054", "C" * 30, "forming"),
                ("SYN0053", "A" * 15, "non_forming"), ("SYN0054", "A" * 30, "non_forming"),
                ("SYN0055", "T" * 25, "non_forming"), ("SYN0056", "T" * 35, "non_forming")]
GFA_MR_ONLY = ["-skipWGET", "-skipAPR", "-skipSTR", "-skipDR", "-skipIR", "-skipGQ", "-skipZ"]


# ---------------------------------------------------------------- settings

def settings() -> list[dict]:
    s = [dict(tool="HSeeker", id="defaults", tier="A", hit="any hit", args={})]
    # triplex: min_score, p_value, min_len, max_len, seq_type
    tp = lambda i, tier, ms, pv, mn, mx, st="eukaryotic": dict(  # noqa: E731
        tool="triplex", id=i, tier=tier, hit="any triplex.search hit", args=dict(ms=ms, pv=pv, mn=mn, mx=mx, st=st))
    s += [tp("defaults", "A", 15, 0.05, 6, 25), tp("matched", "B", 15, 0.05, 8, 50),
          tp("manuscript", "manuscript", 15, 0.05, 8, 50, "prokaryotic")]
    s += [tp(f"grid_s{ms}_p{pv}_l{mn}", "C", ms, pv, mn, 25) for ms in (10, 15, 20, 25) for pv in (0.01, 0.05)
          for mn in (6, 8)]
    # Triplexator (duplex target sites)
    tx = lambda i, tier, a: dict(tool="Triplexator", id=i, tier=tier, hit="triplex target site (-ds)", args=a)  # noqa: E731
    s += [tx("defaults", "A", []), tx("matched", "B", ["-l", "16", "-e", "10"]),
          tx("manuscript", "manuscript", ["-l", "10", "-L", "-1", "-e", "15", "-fr", "off"])]
    s += [tx(f"grid_l{l}_e{e}_fr{fr}", "C", ["-l", str(l), "-e", str(e), "-fr", fr])
          for l in (10, 12, 16) for e in (5, 10, 15) for fr in ("on", "off")]
    s += [dict(tool="Triplexator", id="tfo_tts_same_molecule", tier="exploratory",
               hit="TFO-TTS pair within the same record (-ss and -ds)", args=[])]
    # NeSSie: union of fixed lengths k..k+4
    ns = lambda i, tier, mode, k0, a, hit: dict(tool="NeSSie", id=i, tier=tier, hit=hit,  # noqa: E731
                                                args=dict(mode=mode, ks=list(range(k0, k0 + 5)), extra=a))
    s += [ns("T_defaults", "A", "-T", 16, [], "triplex-potential motif (-T)"),
          ns("T_matched", "B", "-T", 16, ["-p", "10", "-m", "10"], "triplex-potential motif (-T)"),
          ns("M_defaults", "A2", "-M", 16, [], "mirror motif (-M)"),
          ns("M_matched", "B2", "-M", 16, ["-m", "10", "-g", "10"], "mirror motif (-M)")]
    s += [ns(f"grid_T_k{k}_p{p}", "C", "-T", k, ["-p", str(p)], "triplex-potential motif (-T)")
          for k in (12, 16, 20) for p in (0, 10, 20)]
    s += [ns(f"grid_M_k{k}_m{m}_g{g}", "C2", "-M", k, ["-m", str(m), "-g", str(g)], "mirror motif (-M)")
          for k in (12, 16, 20) for m in (0, 10) for g in (0, 10)]
    # nBMST: mirror-repeat module; hit = triplex subset (Subset == 1) or any mirror repeat
    nb = lambda i, tier, a, sub: dict(tool="nBMST", id=i, tier=tier,  # noqa: E731
                                      hit="mirror repeat, triplex subset" if sub else "any mirror repeat",
                                      args=dict(extra=a, subset=sub))
    s += [nb("triplex_defaults", "A", [], True), nb("triplex_matched", "B",
          ["-minMRrep", "8", "-maxMRspacer", "10", "-maxTriplexSpacer", "11"], True),
          nb("MR_defaults", "A2", [], False),
          nb("MR_matched", "B2", ["-minMRrep", "8", "-maxMRspacer", "10", "-maxTriplexSpacer", "11"], False)]
    s += [nb(f"grid_rep{r}_sp{sp}_{'tri' if sub else 'MR'}", "C",
             ["-minMRrep", str(r), "-maxMRspacer", str(sp), "-maxTriplexSpacer", str(sp + 1)], sub)
          for r in (8, 10, 12) for sp in (5, 10) for sub in (True, False)]
    return s


# ---------------------------------------------------------------- runners

def _run(cmd: list[str], cwd: Path, timeout: int = 3600) -> None:
    (cwd / "cmd.txt").write_text(" ".join(map(str, cmd)) + "\n")
    t = time.perf_counter()
    p = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True, timeout=timeout)
    (cwd / "stderr.txt").write_text(p.stderr[-20000:])
    (cwd / "runtime.txt").write_text(f"{time.perf_counter() - t:.2f} s, exit {p.returncode}\n")
    if p.returncode != 0:
        raise RuntimeError(f"{cmd[0]} exited {p.returncode}: {p.stderr[-500:]}")


def setting_dir(st: dict, fasta: Path) -> Path:
    key = json.dumps([st["tool"], st["id"], st["args"], fasta.name], sort_keys=True)
    d = RAW / st["tool"] / f"{fasta.stem}__{st['id']}__{hashlib.sha1(key.encode()).hexdigest()[:8]}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def run_setting(st: dict, fasta: Path, seqs: dict[str, str]) -> dict[str, list[dict]]:
    """Return record_id -> list of hits (dicts with start, end, score)."""
    hits: dict[str, list[dict]] = {r: [] for r in seqs}
    tool = st["tool"]
    if tool == "HSeeker":
        for rid, s in seqs.items():
            for h in hseeker.scan_sequence(s):
                hits[rid].append(dict(start=h["start"], end=h["end"], score=h.get("total_score")))
        return hits
    d = setting_dir(st, fasta)
    done = d / "DONE"
    if tool == "triplex":
        out = d / "hits.csv"
        if not done.exists():
            a = st["args"]
            _run([RSCRIPT, ROOT / "analysis" / "scripts" / "triplex_search_params.R", fasta, out,
                  a["ms"], a["pv"], a["mn"], a["mx"], a["st"]], d)
            done.touch()
        for r in csv.DictReader(open(out)):
            hits[r["record_id"]].append(dict(start=int(r["start"]), end=int(r["end"]), score=float(r["score"])))
    elif tool == "Triplexator":
        out = d / "out"
        if not done.exists():
            if st["id"] == "tfo_tts_same_molecule":
                cmd = [TRIPLEXATOR, "-ss", fasta, "-ds", fasta, "-of", "0", "-od", d, "-o", "out"]
            else:
                cmd = [TRIPLEXATOR, "-ds", fasta, *st["args"], "-of", "0", "-od", d, "-o", "out"]
            _run(cmd, d)
            done.touch()
        for line in open(out):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if st["id"] == "tfo_tts_same_molecule":
                if f[0] != f[3]:
                    continue  # TFO and TTS from different records
                hits[f[3]].append(dict(start=int(f[4]) + 1, end=int(f[5]), score=float(f[6])))
            else:
                hits[f[0]].append(dict(start=int(f[1]) + 1, end=int(f[2]), score=float(f[3])))
    elif tool == "NeSSie":
        a = st["args"]
        for k in a["ks"]:
            out = d / f"out_k{k}"
            if not (d / f"DONE_k{k}").exists():
                _run([NESSIE, "-I", fasta, "-O", out, a["mode"], "-k", str(k), *a["extra"]], d)
                (d / f"DONE_k{k}").touch()
            rid = None
            for line in open(out):
                line = line.rstrip("\n")
                if line.startswith(">"):
                    rid = line[1:].split()[0]
                elif line.startswith("$|") and rid is not None:
                    length = int(line.split("|")[1])
                    hits[rid].append(dict(start=None, end=None, score=None, length=length, k=k))
                elif line.startswith("@indexes:") and rid is not None and hits[rid]:
                    idx = [int(x) for x in line.split(":", 1)[1].split("|") if x.strip()]
                    if idx:
                        hits[rid][-1]["start"] = idx[0] + 1
                        hits[rid][-1]["end"] = idx[0] + hits[rid][-1]["length"]
    elif tool == "nBMST":
        a = st["args"]
        out = d / "o_MR.tsv"
        if not done.exists():
            shutil.copy(fasta, d / "s.fa")
            _run([GFA, *GFA_MR_ONLY, "-seq", "s.fa", "-out", "o", *a["extra"]], d)  # short relative paths
            done.touch()
        if out.exists():
            for r in csv.DictReader(open(out), delimiter="\t"):
                if a["subset"] and r["Subset"] != "1":
                    continue
                hits[r["Sequence_name"]].append(dict(start=int(r["Start"]), end=int(r["Stop"]), score=None))
    return hits


# ---------------------------------------------------------------- statistics

def metrics(y, p, fam, rng) -> dict:
    c = rb.confusion(y, p)
    tp, fn, tn, fp = c["TP"], c["FN"], c["TN"], c["FP"]
    w = lambda k, n: [round(x, 4) for x in rb.wilson(k, n)]  # noqa: E731
    lo, hi = rb.bootstrap_mcc_ci(y, p, rng, N_BOOT)
    _, _, cl = gv.cluster_bootstrap(y, p, fam, rng, N_BOOT)
    r = lambda k, n: round(k / n, 4) if n else None  # noqa: E731
    return dict(n=len(y), TP=tp, FN=fn, TN=tn, FP=fp,
                sensitivity=r(tp, tp + fn), sensitivity_ci95=w(tp, tp + fn),
                specificity=r(tn, tn + fp), specificity_ci95=w(tn, tn + fp),
                PPV=r(tp, tp + fp), PPV_ci95=w(tp, tp + fp), NPV=r(tn, tn + fn), NPV_ci95=w(tn, tn + fn),
                MCC=round(rb.mcc_of(c), 4), MCC_ci95_record=[round(lo, 4), round(hi, 4)],
                MCC_ci95_family=[round(cl[0], 4), round(cl[1], 4)])


def _stats(y, p):
    c = rb.confusion(y, p)
    sens = c["TP"] / (c["TP"] + c["FN"]) if c["TP"] + c["FN"] else float("nan")
    spec = c["TN"] / (c["TN"] + c["FP"]) if c["TN"] + c["FP"] else float("nan")
    return rb.mcc_of(c), sens, spec


def paired(y, ph, pt, fam, rng) -> dict:
    """HSeeker minus tool: point differences, record and family bootstrap CIs, exact McNemar."""
    point = [a - b for a, b in zip(_stats(y, ph), _stats(y, pt))]
    by_f: dict[str, list[int]] = {}
    for i, g in enumerate(fam):
        by_f.setdefault(g, []).append(i)
    keys = list(by_f)
    out = {}
    for unit in ("record", "family"):
        ds = []
        for _ in range(N_BOOT):
            idx = ([rng.randrange(len(y)) for _ in y] if unit == "record"
                   else [i for _ in keys for i in by_f[rng.choice(keys)]])
            yy = [y[i] for i in idx]
            if len(set(yy)) < 2:
                continue
            a = _stats(yy, [ph[i] for i in idx])
            b = _stats(yy, [pt[i] for i in idx])
            ds.append([u - v for u, v in zip(a, b)])
        for j, name in enumerate(("MCC", "sensitivity", "specificity")):
            v = sorted(d[j] for d in ds if not math.isnan(d[j]))
            out[f"d{name}_ci95_{unit}"] = [round(v[int(0.025 * len(v))], 4),
                                           round(v[min(len(v) - 1, int(0.975 * len(v)))], 4)] if v else None
    b = sum(1 for a, h, t in zip(y, ph, pt) if h == a and t != a)  # HSeeker right, tool wrong
    c = sum(1 for a, h, t in zip(y, ph, pt) if h != a and t == a)  # tool right, HSeeker wrong
    n = b + c
    pval = min(1.0, 2 * sum(math.comb(n, i) for i in range(0, min(b, c) + 1)) / 2 ** n) if n else 1.0
    return dict(dMCC=round(point[0], 4), dsensitivity=round(point[1], 4), dspecificity=round(point[2], 4),
                **out, discordant_hseeker_right=b, discordant_tool_right=c, mcnemar_exact_p=round(pval, 6))


def select(train, grid_preds, y) -> int:
    best, bi = -2.0, 0
    for j, p in enumerate(grid_preds):
        m = rb.mcc_of(rb.confusion([y[i] for i in train], [p[i] for i in train]))
        if m > best + 1e-12:
            best, bi = m, j
    return bi


def out_of_fold(grid_preds, y, folds) -> tuple[list[int], list[int]]:
    pred, chosen = [None] * len(y), []
    for test in folds:
        ts = set(test)
        j = select([i for i in range(len(y)) if i not in ts], grid_preds, y)
        chosen.append(j)
        for i in test:
            pred[i] = grid_preds[j][i]
    return pred, chosen


# ---------------------------------------------------------------- main

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tools", default="triplex,Triplexator,NeSSie,nBMST")
    args = ap.parse_args()
    wanted = {"HSeeker", *args.tools.split(",")}
    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    csv_path = rb.DEFAULT_CSV
    recs = rb.load_balanced(csv_path)
    exp = {r["record_id"]: r for r in csv.DictReader(open(ROOT / "hdna_benchmark_experimental_v3.csv"))}
    bal = {r["record_id"]: r for r in csv.DictReader(open(csv_path))}
    ids = [r["record_id"] for r in recs]
    y = [r["y"] for r in recs]
    fam = gv.assign_groups(recs, "family")
    subsets = {"balanced": list(range(len(recs))),
               "experimental": [i for i, r in enumerate(recs) if r["record_type"] == "experimental"],
               "primary": [i for i, r in enumerate(recs) if r["record_type"] == "experimental"
                           and exp[r["record_id"]]["disposition"] == "primary"]}
    fasta = RAW / "benchmark130.fa"
    fasta.write_text("".join(f">{r['record_id']}\n{r['sequence']}\n" for r in recs))
    seqs = {r["record_id"]: r["sequence"] for r in recs}
    hp_fasta = RAW / "homopolymers.fa"
    hp_fasta.write_text("".join(f">{i}\n{s}\n" for i, s, _ in HOMOPOLYMERS))
    hp_seqs = {i: s for i, s, _ in HOMOPOLYMERS}

    sets = [s for s in settings() if s["tool"] in wanted]
    calls, records_rows, failures = {}, [], []
    for st in sets:
        key = (st["tool"], st["id"])
        try:
            h = run_setting(st, fasta, seqs)
        except Exception as e:  # recorded, never silently dropped
            failures.append(dict(tool=st["tool"], setting=st["id"], error=str(e)[:500]))
            print("FAILED", key, e, flush=True)
            continue
        calls[key] = [int(bool(h[i])) for i in ids]
        for i, rid in enumerate(ids):
            sc = [x["score"] for x in h[rid] if x.get("score") is not None]
            records_rows.append(dict(record_id=rid, label=recs[i]["label"], tool=st["tool"], setting=st["id"],
                                     tier=st["tier"], call=calls[key][i], n_hits=len(h[rid]),
                                     best_score=max(sc) if sc else "",
                                     hits=";".join(f"{x.get('start')}-{x.get('end')}" for x in h[rid][:20])))
        print(f"{st['tool']:12s} {st['id']:28s} {st['tier']:10s} calls={sum(calls[key])}", flush=True)

    # HSeeker setting C: its own 576-configuration grid
    hs_grid = rb.grid_configs()
    hs_pred = rb.predict_all(recs, hs_grid)
    grids = {"HSeeker": [hs_pred[rb.config_key(c)] for c in hs_grid]}
    for tool in ("triplex", "Triplexator", "NeSSie", "nBMST"):
        g = [calls[(s["tool"], s["id"])] for s in sets if s["tool"] == tool and s["tier"] in ("C", "C2")
             and (s["tool"], s["id"]) in calls]
        if g:
            grids[tool] = g

    # ---- summary metrics
    summary = []
    for (tool, sid), p in calls.items():
        tier = next(s["tier"] for s in sets if (s["tool"], s["id"]) == (tool, sid))
        for name, idx in subsets.items():
            summary.append(dict(tool=tool, setting=sid, tier=tier, set=name,
                                **metrics([y[i] for i in idx], [p[i] for i in idx], [fam[i] for i in idx], rng)))
    # ---- out-of-fold tuned (setting C)
    tuned, oof_lofo = [], {}
    schemes = {"5-fold x 20": [f for rep in range(20) for f in rb.stratified_kfold(recs, 5, SEED + rep)],
               "leave-one-family-out": [[i for i, g in enumerate(fam) if g == u] for u in sorted(set(fam))]}
    stud = gv.assign_groups(recs, "study")
    schemes["leave-one-study-out"] = [[i for i, g in enumerate(stud) if g == u] for u in sorted(set(stud))]
    for tool, g in grids.items():
        for sname, folds in schemes.items():
            if sname == "5-fold x 20":
                yy, pp = [], []
                for rep in range(20):
                    p, _ = out_of_fold(g, y, folds[rep * 5:(rep + 1) * 5])
                    yy += y
                    pp += p
                m = rb.metrics_of(yy, pp)
                tuned.append(dict(tool=tool, scheme=sname, grid_size=len(g), TP=m["TP"], FN=m["FN"], TN=m["TN"],
                                  FP=m["FP"], MCC=round(m["MCC"], 4)))
            else:
                p, chosen = out_of_fold(g, y, folds)
                if sname == "leave-one-family-out":
                    oof_lofo[tool] = p
                tuned.append(dict(tool=tool, scheme=sname, grid_size=len(g),
                                  **{k: v for k, v in metrics(y, p, fam, rng).items()
                                     if k in ("TP", "FN", "TN", "FP", "MCC", "MCC_ci95_family")}))
    # ---- paired differences vs HSeeker
    pairs_rows = []
    h_a = calls[("HSeeker", "defaults")]
    for (tool, sid), p in calls.items():
        if tool == "HSeeker":
            continue
        tier = next(s["tier"] for s in sets if (s["tool"], s["id"]) == (tool, sid))
        if tier in ("C", "C2"):
            continue
        for name, idx in subsets.items():
            pairs_rows.append(dict(tool=tool, setting=sid, tier=tier, set=name,
                                   **paired([y[i] for i in idx], [h_a[i] for i in idx], [p[i] for i in idx],
                                            [fam[i] for i in idx], rng)))
    for tool, p in oof_lofo.items():
        if tool == "HSeeker":
            continue
        pairs_rows.append(dict(tool=tool, setting="tuned out of fold (LOFO)", tier="C", set="balanced",
                               **paired(y, oof_lofo["HSeeker"], p, fam, rng)))
    # ---- wording verdict (§2.4)
    verdict = []
    for tool in ("triplex", "Triplexator", "NeSSie", "nBMST"):
        get = lambda tier: next((r for r in pairs_rows if r["tool"] == tool and r["tier"] == tier  # noqa: E731
                                 and r["set"] == "balanced"), None)
        a, b, c = get("A"), get("B"), get("C")
        if not a:
            continue
        lo, hi = a["dMCC_ci95_family"]
        signs = {math.copysign(1, r["dMCC"]) for r in (a, b, c) if r and r["dMCC"] != 0}
        if lo > 0 and len(signs) == 1 and signs == {1.0}:
            v = "HSeeker better"
        elif hi < 0 and len(signs) == 1 and signs == {-1.0}:
            v = f"{tool} better"
        else:
            v = "comparable"
        verdict.append(dict(tool=tool, dMCC_A=a["dMCC"], dMCC_A_ci95_family=a["dMCC_ci95_family"],
                            dMCC_B=b["dMCC"] if b else None, dMCC_C=c["dMCC"] if c else None, verdict=v))
    # ---- negatives by class, false negatives
    by_class = []
    for (tool, sid), p in calls.items():
        tier = next(s["tier"] for s in sets if (s["tool"], s["id"]) == (tool, sid))
        if tier not in ("A", "B", "A2", "B2", "manuscript"):
            continue
        cls: dict[str, list[int]] = {}
        for i, r in enumerate(recs):
            if r["y"] == 0:
                k = (bal[r["record_id"]].get("generation_method") or "experimental negative").split(":")[-1]
                cls.setdefault(k, []).append(p[i])
        by_class.append(dict(tool=tool, setting=sid, tier=tier,
                             **{f"FP {k}": f"{sum(v)}/{len(v)}" for k, v in sorted(cls.items())},
                             false_negatives=";".join(ids[i] for i in range(len(ids)) if y[i] == 1 and not p[i])))
    # ---- homopolymer panel (settings A and B)
    hp_rows = []
    for st in sets:
        if st["tier"] not in ("A", "B", "A2", "B2"):
            continue
        try:
            h = run_setting(st, hp_fasta, hp_seqs)
        except Exception as e:
            failures.append(dict(tool=st["tool"], setting=st["id"] + " (homopolymers)", error=str(e)[:500]))
            continue
        for rid, s, lab in HOMOPOLYMERS:
            hp_rows.append(dict(tool=st["tool"], setting=st["id"], record_id=rid, sequence=f"{s[0]}x{len(s)}",
                                label=lab, call="hit" if h[rid] else "no hit"))

    def write(name, rows):
        if not rows:
            return
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(OUT / name, "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow({k: (json.dumps(v) if isinstance(v, list) else v) for k, v in r.items()})

    write("per_record_calls.csv", records_rows)
    write("summary_metrics_v2.csv", summary)
    write("tuned_out_of_fold.csv", tuned)
    write("paired_differences.csv", pairs_rows)
    write("by_negative_class.csv", by_class)
    write("homopolymer_panel_tools.csv", hp_rows)
    write("verdict.csv", verdict)
    write("failures.csv", failures)
    meta = dict(date=time.strftime("%Y-%m-%d %H:%M:%S"), plan="RTR/TOOL_COMPARISON_PLAN_2026-10-11.md",
                input_csv=csv_path.name, input_csv_sha256=hashlib.sha256(csv_path.read_bytes()).hexdigest(),
                n_records=len(recs), seed=SEED, n_boot=N_BOOT, settings=sets,
                tools=dict(triplexator=str(TRIPLEXATOR), nessie=str(NESSIE), gfa=str(GFA), rscript=RSCRIPT,
                           hseeker=hseeker.__version__),
                git_commit=subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True,
                                          text=True).stdout.strip(),
                failures=failures)
    (OUT / "metadata.json").write_text(json.dumps(meta, indent=1, default=str))
    for r in summary:
        if r["set"] == "balanced" and r["tier"] in ("A", "B", "A2", "B2", "manuscript"):
            print(f"{r['tool']:12s} {r['setting']:24s} {r['tier']:10s} TP{r['TP']:3d} FN{r['FN']:3d} "
                  f"TN{r['TN']:3d} FP{r['FP']:3d} MCC {r['MCC']} fam{r['MCC_ci95_family']}")
    for v in verdict:
        print(v)


if __name__ == "__main__":
    main()
