# -*- coding: utf-8 -*-
"""Streaming 10x Flex reader for GSE334588, the CLE negative control.

NOTHING HERE OPENS REAL DATA AT IMPORT. Fetching is injected, so the whole reduction path is
exercised offline by selftest_flex_harness.py on a fixture whose answer is known by construction.

WHY A SECOND HARNESS. replication_harness reads SMART-seq2 plate files, a dense genes-by-cells
table per plate. GSE334588 is 10x Flex: barcodes, features and a MatrixMarket sparse matrix per
sample, up to 55 MB compressed. IT CANNOT BE READ DENSELY AND IT DOES NOT NEED TO BE.

THE ARITHMETIC IS DELIBERATELY THE SAME as the plate reader wherever the platform allows, so the
control differs from the disease runs in PLATFORM AND COHORT ONLY. Donor collapse is not
reimplemented here at all: it is imported from replication_harness, so the log1p-after-the-donor-
mean ordering stays pinned by the test that already guards it.
"""
import gzip, io
import numpy as np
from replication_harness import HarnessError, donor_matrix_and_labels  # ordering stays pinned there

def parse_features(raw):
    """features.tsv.gz bytes -> list of gene symbols, one per matrix ROW, in file order.
    10x writes id<TAB>symbol<TAB>type. The symbol is column 2. Rows are 1-based in the matrix."""
    txt = gzip.decompress(raw).decode("utf8", "replace")
    out = []
    for ln in txt.split("\n"):
        if not ln.strip(): continue
        p = ln.rstrip("\r").split("\t")
        out.append(p[1] if len(p) > 1 else p[0])
    if not out:
        raise HarnessError("features file parsed to zero rows")
    return out

def reduce_flex(features, mtx_raw, genes):
    """Stream a MatrixMarket matrix ONE TRIPLET AT A TIME and return the same four things the
    plate reader returns: (cpm, nz, n_cells, med_lib).

    Memory is one float per cell plus a handful of counters. The matrix is never materialized.

    A SYMBOL ON SEVERAL FEATURE ROWS IS SUMMED ACROSS THEM, matching the plate reader, where a
    gene on several index rows is summed rather than first-wins. A symbol absent from the
    feature list yields NaN rather than zero, because absent from the panel and measured as zero
    are different facts and the detection guard treats them so. ABCC6 is absent from this panel
    and that is the case this rule exists for.
    """
    want = {}
    for i, sym in enumerate(features, start=1):
        if sym in genes: want.setdefault(sym, []).append(i)
    rowmap = {}
    for sym, rows in want.items():
        for r in rows: rowmap[r] = sym

    lib = None; gsum = {s: 0.0 for s in want}; gcell = None
    seen_header = False
    with gzip.open(io.BytesIO(mtx_raw), "rt", encoding="utf8", errors="replace") as fh:
        for ln in fh:
            if not ln or ln[0] == "%": continue
            p = ln.split()
            if len(p) < 3: continue
            if not seen_header:
                nrow, ncol = int(p[0]), int(p[1])
                if nrow != len(features):
                    raise HarnessError("matrix has %d rows but the feature list has %d; the two "
                                       "files do not belong together" % (nrow, len(features)))
                lib = np.zeros(ncol, dtype=float)
                # ⛔ DETECTION IS COUNTED PER CELL, NEVER PER TRIPLET. A symbol on several
                # feature rows contributes several entries for the SAME cell, and incrementing
                # a counter per entry counts that cell more than once. That inflates DET, which
                # feeds guard_detected, so a gene could clear the 1 per cent floor it should not.
                # The plate reader avoids this by summing a gene's rows into a per-cell vector
                # before counting; a boolean per cell is the sparse equivalent. Caught by the
                # duplicated-symbol fixture in selftest_flex_harness.py, which is why it is there.
                gcell = {s: np.zeros(ncol, dtype=bool) for s in want}
                seen_header = True; continue
            r = int(p[0]); c = int(p[1]); v = float(p[2])
            lib[c-1] += v
            s = rowmap.get(r)
            if s is not None:
                gsum[s] += v
                if v > 0: gcell[s][c-1] = True
    if not seen_header:
        raise HarnessError("matrix header never seen; refusing to proceed")
    keep = lib > 0
    n_cells = int(keep.sum())
    if n_cells == 0:
        raise HarnessError("sample has no cell with a positive library size")
    total = float(lib[keep].sum())
    if not np.isfinite(total) or total <= 0:
        raise HarnessError("sample total count is %r" % total)
    cpm = {g: (gsum[g] / total * 1e6 if g in want else np.nan) for g in genes}
    nz = {g: (int(gcell[g].sum()) if g in want else 0) for g in genes}
    return cpm, nz, n_cells, float(np.median(lib[keep]))

def build_sample_table(records, fetch, genes, on_skip=None):
    """records: dicts with gsm, prefix, donor, arm. fetch(gsm, prefix, kind) -> raw bytes,
    kind being 'features' or 'matrix'. Returns (SL, DET, skipped), the same shape the plate
    reader returns, so everything downstream is unchanged."""
    import pandas as pd
    rows, nz_tot, n_tot, skipped = [], {g: 0 for g in genes}, 0, []
    for rec in records:
        try:
            feats = parse_features(fetch(rec["gsm"], rec["prefix"], "features"))
            mtx = fetch(rec["gsm"], rec["prefix"], "matrix")
        except Exception as e:
            skipped.append((rec["gsm"], rec["prefix"], type(e).__name__))
            if on_skip: on_skip(rec, e)
            continue
        cpm, nz, n_cells, med = reduce_flex(feats, mtx, genes)
        for g in genes: nz_tot[g] += nz[g]
        n_tot += n_cells
        row = dict(cpm); row.update({"donor": rec["donor"], "arm": rec["arm"],
                                     "ctrl": bool(rec.get("ctrl", False)), "cells": n_cells,
                                     "lib": med, "gsm": rec["gsm"], "prefix": rec["prefix"]})
        rows.append(row)
        del mtx, feats
    if not rows:
        raise HarnessError("no sample was read; refusing to return an empty table")
    DET = {g: (nz_tot[g] / n_tot if n_tot else 0.0) for g in genes}
    return pd.DataFrame(rows), DET, skipped
