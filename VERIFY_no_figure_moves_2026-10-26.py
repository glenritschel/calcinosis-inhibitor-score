# -*- coding: utf-8 -*-
"""Addendum v2 asserts that tonight's three fixes to inhibition_score.py move no published
figure. THAT WAS AN ARGUMENT, NOT A MEASUREMENT, because the plate table was still inside the
Colab session. It is on disk now, so the argument is replaced by a run.

The test: re-derive RESULT_11's observed difference from the plate table under the PATCHED
module and under the PRE-PATCH BACKUP, and require both to equal the recorded -0.2939.
A fix that silently moved the statistic fails here.
"""
import importlib.util, sys, hashlib
import numpy as np, pandas as pd

PLATES   = "PLATE_TABLE_GSE195452_optionB_2026-10-26.csv"
EXPECTED = -0.2939
TOL      = 0.00005          # the figure is recorded to four decimals
SIX      = ["ENPP1","ANKH","MGP","GGCX","VKORC1","NT5E"]

print("plate table sha256:", hashlib.sha256(open(PLATES,"rb").read()).hexdigest()[:16])
PL = pd.read_csv(PLATES)
print("plates: %d   donors: %d" % (len(PL), PL["donor"].nunique()))
assert len(PL) == 108 and PL["donor"].nunique() == 50

def load(path, name):
    """importlib needs a .py suffix to infer a loader, and the pre-patch module is a .bak file.
    The copy goes to the SYSTEM temp dir, not this folder: deleting inside a connected folder is
    off by default, so a scratch file written here could not be cleaned up afterward."""
    import shutil, tempfile, os
    if not path.endswith(".py"):
        tmp = os.path.join(tempfile.mkdtemp(), name + ".py")
        shutil.copyfile(path, tmp); path = tmp
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); sys.modules[name] = m
    spec.loader.exec_module(m); return m

def observed(mod, genes):
    # the committed donor collapse: mean of plate values per donor, THEN log1p (never before)
    D   = PL.groupby("donor")[list(genes)].mean().T
    dm  = np.log1p(D)
    lab = PL.groupby("donor")["ctrl"].first().reindex(dm.columns)
    assert not lab.isna().any(), "a donor has no control label"
    ctrl = np.asarray(lab.values, dtype=bool)
    sc, cov, used = mod.set_score(dm, list(genes))
    v = np.asarray(sc, float)
    return float(v[~ctrl].mean() - v[ctrl].mean()), cov, used

rows = []
for label, path in (("PATCHED  (tonight)", "inhibition_score.py"),
                    ("PRE-PATCH (backup)", "inhibition_score.py.bak_pre_mararo_2026-10-26")):
    mod = load(path, "m_" + label.split()[0].lower())
    diff, cov, used = observed(mod, SIX)
    rows.append((label, diff, cov, len(used)))
    print("%-20s  observed %+.6f   coverage %.4f   genes %d" % (label, diff, cov, len(used)))

print("\nRESULT_11 recorded   %+.4f" % EXPECTED)
a, b = rows[0][1], rows[1][1]
print("patched minus pre-patch: %.3e" % abs(a - b))

fail = []
if abs(a - EXPECTED) > TOL: fail.append("patched does not reproduce RESULT_11")
if abs(b - EXPECTED) > TOL: fail.append("pre-patch does not reproduce RESULT_11")
if abs(a - b) > 1e-12:      fail.append("THE FIXES MOVED THE STATISTIC")
if fail:
    print("\n⛔ FAIL: " + "; ".join(fail)); sys.exit(1)
print("\n✅ BOTH MODULES REPRODUCE RESULT_11 AND AGREE TO MACHINE PRECISION.")
print("   The claim in addendum v2 is now measured rather than argued.")
print("\n⚠ WHAT THIS DOES NOT COVER: the per-gene DETECTION figures. The plate table carries")
print("   cpm, donor, ctrl, cells, lib, gsm and title, but NOT the per-gene nonzero cell counts,")
print("   so DET cannot be rebuilt from it and the 0.0137 / 0.5291 figures are not re-derived")
print("   here. Those come through reduce_plate, which the (row > 0) fix does not touch.")
