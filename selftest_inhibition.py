# -*- coding: utf-8 -*-
"""SELFTEST. Every guard must be SEEN TO FAIL on a corrupted copy before any real run.
A check is worth nothing until it has returned the answer you do not want, on a concrete
input, along the code path the real run uses."""
import numpy as np, pandas as pd, sys
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from inhibition_score import (INHIBITOR, OSTEOGENIC, set_score, donor_matrix, guard_can_vary,
                              guard_responds, guard_reachable, guard_detected, label_null,
                              verdict, Void, MIN_CELL_DETECT)

rng = np.random.default_rng(11)
GENES = INHIBITOR + OSTEOGENIC + ["SPP1"] + ["FILLER%d" % i for i in range(60)]

def synth(n_dis=20, n_con=20, inh_shift=0.0, ost_shift=0.0, seed=1):
    r = np.random.default_rng(seed)
    donors = ["D%02d" % i for i in range(n_dis+n_con)]
    # ⛔ CENTERED AT 5, NOT 0, AND THAT IS NOT CosMETIC. The real dmat is np.log1p of a gene's
    # share of plate counts, so IT CANNOT BE NEGATIVE. A fixture centered at zero puts half its
    # values below zero, and the MIN_DONOR_NONZERO gate then passed VACUOUSLY for years: under
    # the old (row != 0) form every negative read as "detected", so the clean fixture asserting
    # "all-detected set passes" was passing on data the gate will never see. Centering at 5
    # leaves every z-scored statistic identical, because z-scoring is location invariant, while
    # putting the gate's input in the domain the gate was written for. Changed 2026-10-26.
    X = pd.DataFrame(r.normal(5,1,(len(GENES), len(donors))), index=GENES, columns=donors)
    assert (X.values >= 0).all(), "fixture produced a negative value; dmat is log1p and cannot"
    dis = donors[:n_dis]
    X.loc[INHIBITOR, dis] += inh_shift
    X.loc[OSTEOGENIC, dis] += ost_shift
    lab = pd.Series([d in dis for d in donors], index=donors)
    return X, lab, dis

def run(X, lab, dis, alpha=0.05, n_perm=4000):
    guard_reachable(int(lab.sum()), int((~lab).sum()), alpha)
    si,_,_ = set_score(X, INHIBITOR); guard_can_vary(si); guard_responds(X, INHIBITOR, dis)
    so,_,_ = set_score(X, OSTEOGENIC); guard_can_vary(so); guard_responds(X, OSTEOGENIC, dis)
    inh = label_null(si, lab, n_perm=n_perm, direction="lower")
    ost = label_null(so, lab, n_perm=n_perm, direction="higher")
    return verdict(inh, ost, alpha), inh, ost

ok = True
def check(name, cond, extra=""):
    global ok
    print(("  ok   " if cond else "  FAIL ") + name + ("  " + extra if extra else ""))
    ok = ok and cond

print("=== POSITIVE CONTROLS: the code must find what was planted ===")
v,i,o = run(*synth(inh_shift=-1.2)); check("Model A planted -> Model A found", v.startswith("MODEL A"), v)
v,i,o = run(*synth(ost_shift=+1.2)); check("Model B planted -> Model B found", v.startswith("MODEL B"), v)
v,i,o = run(*synth(inh_shift=-1.2, ost_shift=+1.2)); check("both planted -> MIXED", v.startswith("MIXED"), v)

print("=== NEGATIVE CONTROL: nothing planted must return nothing ===")
v,i,o = run(*synth()); check("no effect -> NEITHER", v.startswith("NEITHER"), v)

print("=== GUARD 1 must FAIL on a statistic that cannot vary ===")
X,lab,dis = synth(); X.loc[INHIBITOR,:] = 5.0          # constant across donors
try:
    s,_,_ = set_score(X, INHIBITOR); guard_can_vary(s); check("guard_can_vary fired", False, "it did not")
except Void as e: check("guard_can_vary fired", True, str(e)[:60])

print("=== GUARD 2 must FAIL when scoring ignores its inputs ===")
import inhibition_score as m
X,lab,dis = synth(); real = m.set_score
m.set_score = lambda d,g: (pd.Series(0.0, index=d.columns), 1.0, list(g))   # deaf scorer
try:
    m.guard_responds(X, INHIBITOR, dis); check("guard_responds fired", False, "it did not")
except Void as e: check("guard_responds fired", True, str(e)[:60])
finally: m.set_score = real

print("=== GUARD 3 must FAIL when the criterion is unreachable ===")
try:
    guard_reachable(2, 2, 0.05); check("guard_reachable fired at n=2v2", False, "it did not")
except Void as e: check("guard_reachable fired at n=2v2", True, str(e)[:70])
d,p = guard_reachable(20,20,0.05); check("guard_reachable passes at n=20v20", p < 0.05, "min p=%.2g" % p)

print("=== GUARD 4: detection ===")
X,lab,dis = synth()
usable,dropped = guard_detected(X, INHIBITOR, {g:0.9 for g in INHIBITOR})
check("all-detected set passes", len(usable)==len(INHIBITOR) and not dropped)

rates = {g:0.9 for g in INHIBITOR}; rates["ENPP1"]=0.001
usable,dropped = guard_detected(X, INHIBITOR, rates)
check("one dropout gene is dropped, set survives", "ENPP1" in dropped and len(usable)==len(INHIBITOR)-1,
      dropped.get("ENPP1",""))

print("=== GUARD 4 must FAIL when the set is mostly undetected ===")
try:
    guard_detected(X, INHIBITOR, {g:0.0001 for g in INHIBITOR}); check("guard_detected fired", False, "it did not")
except Void as e: check("guard_detected fired", True, str(e)[:70])

print("=== GUARD 4 must catch structural zeros even with no cell rates ===")
X2,_,_ = synth(); X2.loc[INHIBITOR[:6]] = 0.0
try:
    guard_detected(X2, INHIBITOR); check("guard_detected fired on donor zeros", False, "it did not")
except Void as e: check("guard_detected fired on donor zeros", True, str(e)[:70])

print("=== the ambiguous gene is never scored into a primary ===")
check("SPP1 absent from both sets", "SPP1" not in INHIBITOR and "SPP1" not in OSTEOGENIC)
check("no collagen in either set", not any(g.startswith("COL") for g in INHIBITOR+OSTEOGENIC))

# ======================================================================================
# FIXTURES ADDED 2026-10-26 after Enock Mararo's review of the score specification.
# Each defect is exercised in BOTH directions: it must fire where it should and stay
# silent where it should not. A rule with only a fire fixture passes when it is broken.
# ======================================================================================
import numpy as _np, pandas as _pd

def _fix_dmat(n_donors=10, seed=0, genes=("ENPP1","ANKH","MGP","GGCX","VKORC1","NT5E")):
    r = _np.random.default_rng(seed)
    d = [f"D{i}" for i in range(n_donors)]
    return _pd.DataFrame({g: r.uniform(1, 5, n_donors) for g in genes}, index=d).T

# ---- POINT 5: the survivor floor must be re-checked after the zero-variance drop ------
_X = _fix_dmat(); _X.loc["MGP"] = 3.0; _X.loc["NT5E"] = 2.0      # two constant genes -> 4 left
try:
    set_score(_X, list(_X.index))
    check("set_score voids below the floor after zero-variance drop", False,
          "it returned a score on 4 genes under a floor of 5")
except Void as e:
    check("set_score voids below the floor after zero-variance drop", True, str(e)[:70])

_X1 = _fix_dmat(); _X1.loc["MGP"] = 3.0                          # one constant gene -> 5 left
try:
    sc, cov, used = set_score(_X1, list(_X1.index))
    check("set_score does NOT void at exactly the floor", len(used) == 5, "used=%d" % len(used))
    check("coverage counts the genes USED, not the gate list", abs(cov - 5/6) < 1e-12,
          "coverage reported %.4f and 5 of 6 genes entered the mean" % cov)
except Void as e:
    check("set_score does NOT void at exactly the floor", False, "it voided: %s" % str(e)[:50])

# a declared diagnostic may score below the floor, but only by saying so
_X2 = _fix_dmat()
try:
    sc, cov, used = set_score(_X2, ["VKORC1"], min_survivors=1)
    check("a DECLARED diagnostic may score one gene", len(used) == 1, "used=%d" % len(used))
except Void as e:
    check("a DECLARED diagnostic may score one gene", False, str(e)[:60])
try:
    set_score(_X2, ["VKORC1"])
    check("an UNDECLARED single-gene score is refused", False, "it was allowed silently")
except Void as e:
    check("an UNDECLARED single-gene score is refused", True, str(e)[:60])

# ---- POINT 2: the donor gate must treat an ABSENT donor value as not detected --------
_X3 = _fix_dmat()
_X3.loc["ENPP1"] = [_np.nan]*6 + [4.0, 4.1, 3.9, 4.2]            # absent for 6 of 10 donors
_rates = {g: 0.9 for g in _X3.index}
_dropped = {}
try:
    _u, _dropped = guard_detected(_X3, list(_X3.index), _rates)
except Void:
    _dropped = {"ENPP1": "set voided"}
check("NaN donor values do not count as detected", "ENPP1" in _dropped,
      "ENPP1 absent for 6 of 10 donors was %s" % ("dropped" if "ENPP1" in _dropped else "KEPT"))

_X4 = _fix_dmat()                                                 # nothing absent
_u4, _d4 = guard_detected(_X4, list(_X4.index), {g: 0.9 for g in _X4.index})
check("a fully present set is not dropped by the NaN fix", len(_u4) == 6, "usable=%d" % len(_u4))


print("\n" + ("SELFTEST PASS: every guard fires when it should and stays silent when it should not."
               if ok else "SELFTEST FAIL"))
sys.exit(0 if ok else 1)
