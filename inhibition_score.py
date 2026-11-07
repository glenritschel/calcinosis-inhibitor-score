# -*- coding: utf-8 -*-
"""Scoring and guards for PREREG_INHIBITION_v1_2026-10-20 (sha256 be4d5f0c).

Pure numpy + pandas so it runs on the device shell, which has no scipy.
NOTHING HERE OPENS REAL DATA. The gene sets are the committed ones and are not editable
at call time on purpose: passing your own set is a deviation and must be declared.

THE GUARDS ARE THE POINT. This program has already pre-registered a primary statistic that
was mathematically incapable of varying, while its power guard certified it because the guard
relabeled outputs instead of re-scoring inputs. Every guard here is therefore exercised on a
corrupted copy in the selftest and must be seen to FAIL there before any real run is trusted.
"""
import numpy as np, pandas as pd
from math import comb

# ---- committed gene sets, prereg section 3. Do not edit. -------------------------------
INHIBITOR = ["ENPP1","ANKH","MGP","ABCC6","TNFRSF11B","GGCX","VKORC1","NT5E"]
OSTEOGENIC = ["RUNX2","SP7","BMP2","MSX2","ALPL","IBSP","BGLAP"]
AMBIGUOUS  = ["SPP1"]          # reported, never scored into a primary
PREREG_SHA = "be4d5f0c"

# Sentinel so set_score can tell "caller said nothing" from "caller explicitly chose a floor",
# while the floor constant itself stays defined with the other detection thresholds below.
MIN_SET_SURVIVORS_SENTINEL = object()

class Void(Exception):
    """Raised when a guard says the DESIGN is void. Not a result. Not a null."""

# ---- per-donor collapse and scoring ----------------------------------------------------
def donor_matrix(expr, donors):
    """expr: genes x samples. donors: sample -> donor id. Returns genes x donors."""
    return expr.groupby(donors.reindex(expr.columns).values, axis=1).mean()

def set_score(dmat, genes, min_survivors=MIN_SET_SURVIVORS_SENTINEL):
    """z-score each gene across donors, then mean over the genes that actually survive.
    Returns (scores per donor, coverage fraction, genes used).

    ⛔ THE SURVIVOR FLOOR IS RE-CHECKED HERE, AFTER THE ZERO-VARIANCE DROP, AND AN EARLIER
    VERSION DID NOT CHECK IT AT ALL. guard_detected enforces the floor on the genes that clear
    the DETECTION gates; item 14 of the score spec then drops any gene whose donor-level SD is
    exactly zero. Nothing re-counted in between, so a set could clear a floor of five, lose two
    genes to zero variance, and return a score computed on THREE genes with no void raised.
    Worse, coverage was len(present)/len(genes), computed from the gate list, so it reported
    1.0000 while a third of the set had silently left the mean: the one number that would have
    shown the shrinkage was blind to it. Found 2026-10-26 by Enock Mararo, reading the
    specification rather than the code, as point 5 of his review. Demonstrated before fixing.

    min_survivors: the floor. Defaults to MIN_SET_SURVIVORS. A DIAGNOSTIC that deliberately
    scores a subset below the floor -- the leave-one-out, which scores single genes -- must
    pass its own value EXPLICITLY, so that scoring below the floor is always a declared act and
    never a silent one.
    """
    if min_survivors is MIN_SET_SURVIVORS_SENTINEL:
        min_survivors = MIN_SET_SURVIVORS
    present = [g for g in genes if g in dmat.index]
    if not present:
        raise Void("no genes from the set are present; the set cannot be scored here")
    sub = dmat.loc[present]
    sd = sub.std(axis=1, ddof=1).replace(0, np.nan)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0)
    z = z.dropna(how="all")
    if z.empty:
        raise Void("every gene in the set is constant across donors")
    used = list(z.index)
    if len(used) < min_survivors:
        raise Void("%d of %d genes reach the mean after the zero-variance drop (%s), below the "
                   "floor of %d; the SET IS VOID. The detection gates passed %d, so the loss is "
                   "to zero variance across donors, not to dropout."
                   % (len(used), len(genes), ", ".join(used), min_survivors, len(present)))
    return z.mean(axis=0), len(used)/len(genes), used

# ---- GUARD 1: can the statistic vary at all? -------------------------------------------
def guard_can_vary(scores, min_sd=1e-6):
    sd = float(np.std(np.asarray(scores, float), ddof=1))
    if not np.isfinite(sd) or sd < min_sd:
        raise Void("between-donor SD of the set score is %.3g; the statistic cannot vary, "
                   "so the DESIGN is void rather than the hypothesis" % sd)
    return sd

# ---- GUARD 2: does an injected effect actually move the statistic? ---------------------
def guard_responds(dmat, genes, disease_donors, delta=1.0, min_frac=0.25):
    """Inject delta (in units of each gene's own SD) into the set's genes in disease donors
    and confirm the statistic moves IN THAT DIRECTION by a substantial amount.

    ⛔ IT DOES NOT CHECK THAT THE MOVEMENT EQUALS delta, AND AN EARLIER VERSION DID.
    set_score z-scores every gene across donors, so it is scale-normalizing by construction:
    the injection inflates the very SD the score then divides by, and exact recovery is
    arithmetically impossible. The first version demanded moved ~= delta and VOIDED A VALID
    DESIGN on its second positive control. A guard that fails for a reason of its own making
    is worse than no guard, because it looks like evidence. What can be asserted is direction
    and substantiality, so that is what is asserted. Found and fixed 2026-10-20.
    """
    base, _, used = set_score(dmat, genes)
    spiked = dmat.copy()
    for g in used:
        s = spiked.loc[g].std(ddof=1)
        spiked.loc[g, disease_donors] = spiked.loc[g, disease_donors] + delta*s
    after, _, _ = set_score(spiked, genes)
    d_before = base[disease_donors].mean() - base.drop(disease_donors).mean()
    d_after  = after[disease_donors].mean() - after.drop(disease_donors).mean()
    moved = d_after - d_before
    if not np.isfinite(moved) or moved < min_frac*delta:
        raise Void("injected %.2f SD and the statistic moved only %.3f (floor %.3f); the "
                   "scoring path does not respond to its own inputs"
                   % (delta, moved, min_frac*delta))
    return moved

# ---- GUARD 3: is the criterion arithmetically reachable? -------------------------------
def guard_reachable(n_disease, n_control, alpha=0.05):
    """A permutation p below alpha is impossible when there are too few distinct donor-label
    assignments. This is the check that was missing when a top-tertile criterion turned out
    to be unreachable about 40 percent of the time."""
    n = n_disease + n_control
    distinct = comb(n, n_disease)
    min_p = 1.0/distinct
    if min_p > alpha:
        raise Void("only %d distinct donor labelings exist, so the smallest attainable "
                   "permutation p is %.3f and the %.2f criterion CANNOT BE MET. The design is "
                   "void at this n, which is not evidence about the hypothesis."
                   % (distinct, min_p, alpha))
    return distinct, min_p

# ---- GUARD 4: is the gene actually DETECTED, or is a zero a dropout? -------------------
# Thresholds FIXED 2026-10-20 before any real expression value was read. Declared in
# DEVIATION_03. Changing one after seeing a result is fitting the guard to the answer.
MIN_CELL_DETECT   = 0.01   # gene must be nonzero in >=1% of cells in the dataset
MIN_DONOR_NONZERO = 0.50   # and nonzero in the pseudobulk of >=50% of donors
MIN_SET_SURVIVORS = 5      # a set with fewer surviving genes is VOID for that dataset

def guard_detected(dmat, genes, cell_detect_rate=None):
    """Separate 'the inhibitor is low' from 'we could not see the inhibitor'.

    THE THREE EARLIER GUARDS DO NOT CATCH DROPOUT. In sparse data a gene that is never
    detected produces a confident, stable, entirely meaningless set score, and a failure of
    inhibition and a failure of measurement are then the SAME NUMBER. That is the tertile
    error wearing different clothes.

    cell_detect_rate: {gene: fraction of cells nonzero}, computed when the pseudobulk is
    built. Omit it only when the input was never cell-level.
    Returns (usable_genes, dropped) where dropped is {gene: reason}.
    """
    usable, dropped = [], {}
    for g in genes:
        if g not in dmat.index:
            dropped[g] = "absent from the platform's feature list"; continue
        if cell_detect_rate is not None:
            r = cell_detect_rate.get(g)
            if r is None:
                dropped[g] = "no detection rate supplied"; continue
            if r < MIN_CELL_DETECT:
                dropped[g] = "detected in %.3f%% of cells, below the %.1f%% floor" % (
                    100*r, 100*MIN_CELL_DETECT); continue
        row = dmat.loc[g]
        # ⛔ (row > 0), NEVER (row != 0). NaN != 0 EVALUATES TRUE IN PANDAS, so the earlier
        # form counted a donor whose gene is ABSENT FROM THE PLATFORM as a donor in which the
        # gene was DETECTED, and a gene absent for 60 percent of donors scored frac_nz = 1.00
        # and sailed through the one gate written to catch exactly that. replication_harness
        # computes the same gate as (DM.loc[g] > 0).mean() and therefore DISAGREED with this
        # module: on a gene absent for 6 of 10 donors the harness said 0.400 and this said
        # 1.000. Every published figure came through the harness path, so no result moves.
        # Found 2026-10-26 via Enock Mararo's point 2.
        frac_nz = float((row > 0).sum()) / max(len(row), 1)
        if frac_nz < MIN_DONOR_NONZERO:
            dropped[g] = "nonzero in only %.0f%% of donors, below the %.0f%% floor" % (
                100*frac_nz, 100*MIN_DONOR_NONZERO); continue
        usable.append(g)
    if len(usable) < MIN_SET_SURVIVORS:
        raise Void("only %d of %d genes survive detection (%s); the SET IS VOID for this "
                   "dataset. This is not evidence that inhibition is intact."
                   % (len(usable), len(genes), "; ".join("%s: %s" % kv for kv in dropped.items())))
    return usable, dropped

# ---- donor-level label null ------------------------------------------------------------
def label_null(scores, is_disease, n_perm=10000, direction="lower", seed=0):
    """Permute the disease label ACROSS DONORS. A gene-set null is not substitutable:
    the claim is about disease, so the label is what must be permuted."""
    rng = np.random.default_rng(seed)
    s = np.asarray(scores, float); lab = np.asarray(is_disease, bool)
    obs = s[lab].mean() - s[~lab].mean()
    perm = np.empty(n_perm)
    for i in range(n_perm):
        p = rng.permutation(lab)
        perm[i] = s[p].mean() - s[~p].mean()
    if direction == "lower":
        p_one = (np.sum(perm <= obs) + 1) / (n_perm + 1)
    else:
        p_one = (np.sum(perm >= obs) + 1) / (n_perm + 1)
    p_two = (np.sum(np.abs(perm) >= abs(obs)) + 1) / (n_perm + 1)
    return dict(observed=obs, p_one_sided=p_one, p_two_sided=p_two,
                direction=direction, n_perm=n_perm)

# ---- verdict, prereg section 6 ---------------------------------------------------------
def verdict(inh, ost, alpha=0.05):
    a = inh["p_one_sided"] < alpha                 # inhibitor LOWER in disease
    b = ost["p_one_sided"] < alpha                 # osteogenic HIGHER in disease
    if a and not b: return "MODEL A SUPPORTED (failure of inhibition)"
    if b and not a: return "MODEL B SUPPORTED (excess of osteogenic stimulus)"
    if a and b:     return "MIXED (both move)"
    return "NEITHER SUPPORTED IN THIS TISSUE (declared in advance as the most likely outcome)"
