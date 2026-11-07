# -*- coding: utf-8 -*-
"""Selftest for flex_harness.py. NO NETWORK. NO REAL DATA. CONSUMES NO LOOK.

Every figure is known by construction and EVERY GUARD IS SHOWN FAILING ON A CASE IT SHOULD FAIL.
The seam is tested too, not only the module: the first attempt at the GPL24676 run died on an
untested seam between a tested harness and the notebook, and that is recorded at decisions 976.
"""
import gzip, io
import numpy as np, flex_harness as F

FAILS = []
def chk(label, got, want):
    ok = (got == want) if not isinstance(want, float) else abs(got - want) < 1e-6
    if not ok: FAILS.append(label)
    print("  %-4s %-54s got %-20r want %r" % ("ok" if ok else "FAIL", label, got, want))
def raises(label, fn, exc=F.HarnessError):
    try:
        fn(); FAILS.append(label); print("  %-4s %-54s DID NOT RAISE" % ("FAIL", label))
    except exc: print("  %-4s %-54s raised %s" % ("ok", label, exc.__name__))
    except Exception as e:
        FAILS.append(label); print("  %-4s %-54s raised WRONG type %s" % ("FAIL", label, type(e).__name__))

gz = lambda s: gzip.compress(s.encode())
# features: MGP appears TWICE, ABCC6 is absent from the panel exactly as in the real Flex data
FEAT = gz("ENSG1\tENPP1\tGene Expression\nENSG2\tMGP\tGene Expression\n"
          "ENSG3\tMGP\tGene Expression\nENSG4\tFILLER\tGene Expression\n")
# 4 features x 3 cells. cell3 is EMPTY (no entries at all).
#   ENPP1 row1: c1=10, c2=5          MGP rows 2 and 3: c1=2+3=5, c2=0+1=1
#   FILLER row4: c1=8, c2=4
#   library: c1 = 10+2+3+8 = 23, c2 = 5+0+1+4 = 10, c3 = 0
#   kept cells 2, total 33, median library 16.5
#   ENPP1 cpm = 15/33*1e6 ; MGP cpm = 6/33*1e6
MTX = gz("%%MatrixMarket matrix coordinate integer general\n%\n4 3 7\n"
         "1 1 10\n1 2 5\n2 1 2\n3 1 3\n3 2 1\n4 1 8\n4 2 4\n")
GENES = ["ENPP1", "MGP", "ABCC6", "VKORC1"]

print("1. FEATURE PARSING")
feats = F.parse_features(FEAT)
chk("symbol is taken from column 2, not the ensembl id", feats[0], "ENPP1")
chk("one entry per matrix row, in file order", len(feats), 4)
chk("a repeated symbol keeps both rows", feats.count("MGP"), 2)
raises("an empty features file is REFUSED", lambda: F.parse_features(gz("")))

print("\n2. STREAMING REDUCTION")
cpm, nz, n_cells, med = F.reduce_flex(feats, MTX, GENES)
chk("empty cells excluded from the cell count", n_cells, 2)
chk("median library over KEPT cells only", med, 16.5)
chk("CPM is the gene's share of the sample total", cpm["ENPP1"], 15/33*1e6)
chk("a symbol on two feature rows is SUMMED", cpm["MGP"], 6/33*1e6)
chk("  and its detection counts cells across both rows", nz["MGP"], 2)
chk("a gene ABSENT from the panel is NaN, not 0", bool(np.isnan(cpm["ABCC6"])), True)
chk("  ABCC6 is the real case this rule exists for", nz["ABCC6"], 0)
chk("a panel gene with no counts is also NaN here", bool(np.isnan(cpm["VKORC1"])), True)
chk("detection counts only nonzero cells", nz["ENPP1"], 2)

print("\n3. GUARDS SHOWN FAILING ON CASES THEY SHOULD FAIL")
raises("a matrix whose row count disagrees with features is REFUSED",
       lambda: F.reduce_flex(feats, gz("%%MatrixMarket\n9 3 1\n1 1 5\n"), GENES))
raises("a matrix with no header is REFUSED",
       lambda: F.reduce_flex(feats, gz("%%MatrixMarket matrix coordinate integer general\n%\n"), GENES))
raises("an all-empty sample is REFUSED, not scored as zero",
       lambda: F.reduce_flex(feats, gz("%%MatrixMarket\n4 2 0\n"), GENES))

print("\n4. THE SAMPLE TABLE AND THE DETECTION DENOMINATOR")
recs = [dict(gsm="G1", prefix="CLE1", donor="CLE1", arm="CLE", ctrl=False),
        dict(gsm="G2", prefix="H1", donor="H1", arm="Healthy", ctrl=True),
        dict(gsm="G3", prefix="BAD", donor="H2", arm="Healthy", ctrl=True)]
def fetch(gsm, prefix, kind):
    if prefix == "BAD": raise IOError("simulated fetch failure")
    return FEAT if kind == "features" else MTX
SL, DET, skipped = F.build_sample_table(recs, fetch, GENES)
chk("a failing fetch is SKIPPED and recorded, not fatal", len(skipped), 1)
chk("  and the skip names the sample", skipped[0][1], "BAD")
chk("surviving samples are in the table", len(SL), 2)
chk("DET denominator is kept cells across ALL samples", DET["ENPP1"], 1.0)
chk("  an off-panel gene has detection 0", DET["ABCC6"], 0.0)
raises("a table with no readable sample is REFUSED",
       lambda: F.build_sample_table([recs[2]], fetch, GENES))

print("\n5. THE SEAM: DONOR COLLAPSE IS IMPORTED, NOT REIMPLEMENTED")
import replication_harness as H
chk("flex_harness uses the SAME donor collapse", F.donor_matrix_and_labels is H.donor_matrix_and_labels, True)
DM, LAB = F.donor_matrix_and_labels(SL, ["ENPP1"])
chk("  so log1p still follows the donor mean", round(float(DM.loc["ENPP1", "CLE1"]), 6),
    round(float(np.log1p(15/33*1e6)), 6))
chk("  and labels survive", bool(LAB["H1"]), True)

print("\n%s  (%d failures)" % ("SELFTEST PASS" if not FAILS else "SELFTEST FAILED: " + ", ".join(FAILS), len(FAILS)))
raise SystemExit(1 if FAILS else 0)
