# -*- coding: utf-8 -*-
"""Selftest for replication_harness.py. NO NETWORK. NO REAL DATA. CONSUMES NO LOOK.

Every check is on a fixture whose answer is known by construction, and EVERY GUARD IS SHOWN
FAILING ON A CASE IT SHOULD FAIL before it is trusted when it passes. Two guards shipped this
month passed their selftests for weeks and still let a defect through, in both cases because
the test only ever asked them to say yes.
"""
import numpy as np, pandas as pd, replication_harness as H

FAILS = []
def chk(label, got, want):
    ok = (got == want) if not isinstance(want, float) else abs(got-want) < 1e-9
    if not ok: FAILS.append(label)
    print("  %-4s %-54s got %-22r want %r" % ("ok" if ok else "FAIL", label, got, want))
def raises(label, fn, exc=H.HarnessError):
    try:
        fn(); FAILS.append(label); print("  %-4s %-54s DID NOT RAISE" % ("FAIL", label))
    except exc as e:
        print("  %-4s %-54s raised %s" % ("ok", label, type(e).__name__))
    except Exception as e:
        FAILS.append(label); print("  %-4s %-54s raised WRONG type %s" % ("FAIL", label, type(e).__name__))

# ---- fixture: 4 wells, one empty, one duplicated gene row, one absent gene -------------
PLATE = pd.DataFrame(
    [[10,0,5,0],[2,0,0,0],[3,0,1,0],[5,0,4,0],[0,0,0,0]],
    index=["ENPP1","MGP","MGP","GGCX","VKORC1"], columns=["c1","c2","c3","c4"])
# column sums 20, 0, 10, 0 -> keep c1 and c3, n_cells 2, plate total 30, median library 15
# ENPP1 = (10+5)/30*1e6 = 500000 ; MGP summed over both rows = (5+1)/30*1e6 = 200000

print("1. PER-PLATE REDUCTION")
cpm, nz, n_cells, med = H.reduce_plate(PLATE, ["ENPP1","MGP","GGCX","VKORC1","ABCC6"])
chk("empty wells dropped from the cell count", n_cells, 2)
chk("median library over KEPT cells only", med, 15.0)
chk("CPM is the gene's share of the plate total", cpm["ENPP1"], 500000.0)
chk("a gene on two index rows is SUMMED, not first-wins", cpm["MGP"], 200000.0)
chk("  and its detection counts both rows", nz["MGP"], 2)
chk("a gene that is all zero has CPM 0", cpm["VKORC1"], 0.0)
chk("  and is detected in no cell", nz["VKORC1"], 0)
chk("a gene ABSENT from the platform is NaN, not 0", np.isnan(cpm["ABCC6"]), True)
chk("  and contributes no detection", nz["ABCC6"], 0)
chk("empty wells excluded from detection too", nz["ENPP1"], 2)
raises("an all-empty plate is REFUSED, not scored as zero",
       lambda: H.reduce_plate(pd.DataFrame([[0,0]], index=["ENPP1"], columns=["a","b"]), ["ENPP1"]))

print("\n2. PLATE TABLE AND THE DETECTION DENOMINATOR")
recs = [dict(gsm="G1", title="P1", donor="Ctrl01", ctrl=True),
        dict(gsm="G2", title="P2", donor="Ctrl01", ctrl=True),
        dict(gsm="G3", title="P3", donor="pt01",  ctrl=False),
        dict(gsm="G4", title="BAD", donor="pt01", ctrl=False)]
def fetch(gsm, title):
    if title == "BAD": raise IOError("simulated fetch failure")
    return PLATE
PL, DET, skipped = H.build_plate_table(recs, fetch, ["ENPP1","MGP","ABCC6"])
chk("a failing fetch is SKIPPED and recorded, not fatal", len(skipped), 1)
chk("  and the skip names the plate", skipped[0][1], "BAD")
chk("the surviving plates are in the table", len(PL), 3)
chk("DET denominator is kept cells across ALL plates", DET["ENPP1"], 1.0)
chk("  an absent gene has detection 0", DET["ABCC6"], 0.0)
raises("a table with no readable plate is REFUSED",
       lambda: H.build_plate_table([recs[3]], fetch, ["ENPP1"]))

print("\n3. THE log1p ORDER, WHICH IS THE PROCEDURE")
PL2 = pd.DataFrame([{"donor":"d1","ctrl":True,"MGP":0.0},{"donor":"d1","ctrl":True,"MGP":100.0}])
DM, LAB = H.donor_matrix_and_labels(PL2, ["MGP"])
after  = float(DM.loc["MGP","d1"])
before = float(np.mean([np.log1p(0.0), np.log1p(100.0)]))
chk("log1p is applied AFTER the donor mean", round(after, 6), round(float(np.log1p(50.0)), 6))
chk("  which is NOT log1p then mean", round(after, 6) == round(before, 6), False)
chk("  and the two genuinely differ on this fixture", abs(after-before) > 1.0, True)
chk("labels survive the collapse", bool(LAB["d1"]), True)

print("\n4. THE OPTION B DONOR SET")
ids24 = {"Ctrl01","Ctrl02","pt01","pt02","GVHD_pt001","Morphea_pt001","Ctrl","Hip","dSSc","Ctrl_pool","-"}
ids18 = {"Ctrl02","pt02","Ctrl_pool","-"}
chk("non-person, non-SSc and overlap ids all removed",
    H.option_b_primary(ids24, ids18), ["Ctrl01","pt01"])
chk("  a platform with nothing left yields nothing",
    H.option_b_primary({"Hip","dSSc"}, set()), [])
isc = lambda d: d.lower().startswith("ctrl")
fifty = ["Ctrl%02d" % i for i in range(1,26)] + ["pt%02d" % i for i in range(1,26)]
chk("check_primary PASSES on exactly 50 at 25 and 25", H.check_primary(fifty, isc), True)
raises("  and REFUSES 49 donors", lambda: H.check_primary(fifty[:49], isc))
raises("  and REFUSES a 26 to 24 split",
       lambda: H.check_primary(["Ctrl%02d" % i for i in range(1,27)] + ["pt%02d" % i for i in range(1,25)], isc))

print("\n5. END TO END ON A COHORT WHOSE ANSWER IS KNOWN")
# two donors, one plate each, one gene; the donor matrix must reproduce by hand
def mk(v): return pd.DataFrame([[v,v],[1,1]], index=["ENPP1","FILLER"], columns=["a","b"])
recs2 = [dict(gsm="A", title="A", donor="Ctrl01", ctrl=True),
         dict(gsm="B", title="B", donor="pt01",  ctrl=False)]
PL3, DET3, _ = H.build_plate_table(recs2, lambda g,t: mk(9 if g=="A" else 1), ["ENPP1","FILLER"])
DM3, LAB3 = H.donor_matrix_and_labels(PL3, ["ENPP1","FILLER"])
# plate A: ENPP1 18 of total 20 -> 900000 ; plate B: 2 of total 4 -> 500000
chk("control donor CPM then log1p", round(float(DM3.loc["ENPP1","Ctrl01"]),6), round(float(np.log1p(900000.0)),6))
chk("patient donor CPM then log1p", round(float(DM3.loc["ENPP1","pt01"]),6), round(float(np.log1p(500000.0)),6))
chk("every gene detected in every cell here", DET3["ENPP1"], 1.0)

print("\n%s  (%d failures)" % ("SELFTEST PASS" if not FAILS else "SELFTEST FAILED: " + ", ".join(FAILS), len(FAILS)))
raise SystemExit(1 if FAILS else 0)
