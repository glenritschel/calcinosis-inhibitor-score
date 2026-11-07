# -*- coding: utf-8 -*-
"""Held-out scoring harness for GPL24676, under PREREG_REPLICATION v2 plus the section 3
addendum, the v4 section 6 addendum and the v5 section 3 addendum.

NOTHING HERE OPENS REAL DATA AT IMPORT. Fetching is injected, so the whole reduction path is
exercised by selftest_replication_harness.py on a synthetic fixture whose answer is known by
construction, WITHOUT consuming the one look.

THE PROCEDURE IS NOT INVENTED HERE. It is transcribed from the discovery run that produced
RESULT_03 on GPL18573, recovered from the notebook cell that ran it, so that the replication
differs from the discovery in PLATFORM AND DONORS ONLY and not in arithmetic. The two places
that order matters are called out below and are pinned by tests.
"""
import numpy as np, pandas as pd

# The same sixteen genes the discovery run carried. Detection is computed over all of them so
# the GPL24676-native survivor set (the pre-registered SECONDARY) can be read off afterwards.
GENES = ["ENPP1","ANKH","MGP","ABCC6","TNFRSF11B","GGCX","VKORC1","NT5E",
         "RUNX2","SP7","BMP2","MSX2","ALPL","IBSP","BGLAP","SPP1"]
FROZEN_SIX = ["ENPP1","ANKH","MGP","GGCX","VKORC1","NT5E"]   # the PRIMARY set, fixed
NONPERSON  = {"-","ctrl","ctrl_pool","hip","dssc"}
NONSSC     = ("gvhd","morphea")

class HarnessError(Exception):
    """A refusal by the harness. Never a result."""

# ---- 1. per-plate reduction, pure and testable -----------------------------------------
def reduce_plate(df, genes=GENES):
    """df: genes (index) x cells (columns) of raw counts, as the GEO per-sample file supplies.

    Returns (cpm, nz, n_cells, med_lib):
      cpm[g]   gene's share of the PLATE's total counts, times 1e6, over kept cells
      nz[g]    number of KEPT cells in which the gene is nonzero
      n_cells  kept cells, meaning cells whose library size is strictly positive
      med_lib  median library size over kept cells

    EMPTY WELLS ARE DROPPED BEFORE ANYTHING ELSE and are excluded from the plate total, from
    the per-gene sums and from the detection denominator. A 384-well plate is mostly a
    measuring device and its empty wells are not observations.

    A GENE APPEARING ON SEVERAL INDEX ROWS IS SUMMED ACROSS THEM, not silently reduced to the
    first. A gene absent from the index yields NaN rather than zero, because absent from the
    platform and measured as zero are different facts and the detection guard treats them so.
    """
    L = np.asarray(df.sum(axis=0), dtype=float)
    keep = L > 0
    n_cells = int(keep.sum())
    if n_cells == 0:
        raise HarnessError("plate has no cell with a positive library size")
    total = float(L[keep].sum())
    if not np.isfinite(total) or total <= 0:
        raise HarnessError("plate total count is %r" % total)
    cpm, nz = {}, {}
    for g in genes:
        if g not in df.index:
            cpm[g] = np.nan; nz[g] = 0; continue
        a = np.asarray(df.loc[g], dtype=float)
        v = (a.sum(axis=0) if a.ndim > 1 else a)[keep]
        nz[g] = int((v > 0).sum())
        cpm[g] = float(v.sum()) / total * 1e6
    return cpm, nz, n_cells, float(np.median(L[keep]))

# ---- 2. the plate table, with fetching injected ----------------------------------------
def build_plate_table(records, fetch, genes=GENES, on_skip=None):
    """records: iterable of dicts with gsm, title, donor, ctrl.
    fetch(gsm, title) -> DataFrame of raw counts. Injected so this is testable offline.

    Returns (PL, DET, skipped). DET[g] is the CELL-LEVEL detection rate, nonzero cells over
    ALL kept cells across every plate, which is the denominator guard_detected expects.
    """
    rows, nz_tot, n_cells_tot, skipped = [], {g: 0 for g in genes}, 0, []
    for rec in records:
        try:
            df = fetch(rec["gsm"], rec["title"])
        except Exception as e:
            skipped.append((rec["gsm"], rec["title"], type(e).__name__))
            if on_skip: on_skip(rec, e)
            continue
        cpm, nz, n_cells, med_lib = reduce_plate(df, genes)
        for g in genes: nz_tot[g] += nz[g]
        n_cells_tot += n_cells
        row = dict(cpm); row.update({"donor": rec["donor"], "ctrl": bool(rec["ctrl"]),
                                     "cells": n_cells, "lib": med_lib,
                                     "gsm": rec["gsm"], "title": rec["title"]})
        rows.append(row)
        del df
    if not rows:
        raise HarnessError("no plate was read; refusing to return an empty table")
    PL = pd.DataFrame(rows)
    DET = {g: (nz_tot[g] / n_cells_tot if n_cells_tot else 0.0) for g in genes}
    return PL, DET, skipped

# ---- 3. donor collapse. THE ORDER HERE IS THE PROCEDURE AND IS PINNED BY TEST ----------
def donor_matrix_and_labels(PL, genes=GENES):
    """Returns (DM, LAB). DM is genes x donors.

    ⛔ log1p IS APPLIED AFTER THE DONOR MEAN, NOT BEFORE. The discovery run computed
    PL.groupby("donor")[G].mean().T and then np.log1p on that. Taking log1p per plate and then
    averaging is a DIFFERENT STATISTIC, because the log of a mean is not the mean of logs, and
    donors contribute between one and five plates here so the gap does not cancel. Swapping the
    two would silently change the held-out statistic relative to the discovery one, which is the
    single thing a replication may not do. selftest_replication_harness.py asserts the two
    orders disagree on the fixture, so this line cannot be reordered without a test failing.
    """
    D = PL.groupby("donor")[list(genes)].mean().T
    LAB = PL.groupby("donor")["ctrl"].first()
    DM = np.log1p(D)
    return DM, LAB

# ---- 4. the held-out donor set, rebuilt from the v5 addendum alone ---------------------
def option_b_primary(ids_24676, ids_18573):
    """Option B, fixed in PREREG_REPLICATION_v5_ADDENDUM_SECTION3 at 23:00 on 2026-10-21.

    Takes the raw donor-id sets of the two platforms' CD90+ skin plates, selected by EXACT
    string equality on tissue and selection marker, and applies the three exclusions the
    addendum names: the non-person ids, the non-SSc donors, and the donors present on both
    platforms. Returns the sorted primary ids.
    """
    keep = lambda s: s.lower() not in NONPERSON and not s.lower().startswith(NONSSC)
    a = {d for d in ids_24676 if keep(d)}
    b = {d for d in ids_18573 if keep(d)}
    return sorted(a - b)

def check_primary(primary, is_control):
    """Refuse to proceed unless the set the addendum fixes is the set in hand."""
    n = len(primary)
    nc = sum(1 for d in primary if is_control(d))
    if (n, nc, n - nc) != (50, 25, 25):
        raise HarnessError("the primary set is %d donors, %d control and %d patient; the "
                           "addendum fixes 50 at 25 and 25, so the set in hand is NOT the "
                           "pre-registered one and nothing may be scored" % (n, nc, n - nc))
    return True
