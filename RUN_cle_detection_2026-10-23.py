# CLE CONTROL: DETECTION GATE. STOPS AT THE DETECTION VERDICT. NOTHING IS SCORED.
# Required by PREREG_CLE_CONTROL_v1_2026-10-23 section 4: DEVIATION_03 thresholds unchanged,
# MIN_CELL_DETECT 0.01, MIN_DONOR_NONZERO 0.50, MIN_SET_SURVIVORS 5, and the VOID rule applies.
#
# ⭐ PREDICTION, FIXED HERE BEFORE THE RUN: I PREDICT DETECTION CLEARS.
# RESULT_05 already read 13 samples and 61,646 cells from THIS dataset and found 6 of 8 inhibitors
# above the floor, the same frozen six as GSE195452. This reads a SUBSET of that dataset, so the
# only way it voids is if the CLE and Healthy arms are much shallower than the pooled set was.
# That is a weak prediction and it is recorded anyway, because a surprise would be informative.
#
# Needs replication_harness, donor_group, flex_harness in the session (payload v2).
import gzip, re, io, collections, hashlib
import numpy as np, pandas as pd, requests
import replication_harness as H, donor_group as DG, flex_harness as F
EXPECT = {"replication_harness.py":"9dd1fa94c7807225","donor_group.py":"f61d761a8a7f72a1","flex_harness.py":"512fe8ab13b7460f"}
hv = lambda f: hashlib.sha256(open(f,"rb").read()).hexdigest()[:16]
got = {f: hv(f) for f in EXPECT}
assert got == EXPECT, "MODULE HASHES DIFFER FROM THE LEDGER: %r" % got
print("module integrity OK")
MIN_CELL_DETECT, MIN_DONOR_NONZERO, MIN_SET_SURVIVORS = 0.01, 0.50, 5
print("DEVIATION_03: cell %.2f | donor %.2f | survivors %d" % (MIN_CELL_DETECT, MIN_DONOR_NONZERO, MIN_SET_SURVIVORS))
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE334nnn/GSE334588/matrix/"
UA = {"User-Agent": "Mozilla/5.0 (research; ORCID 0009-0006-5341-8840)"}
txt = gzip.decompress(requests.get(BASE+"GSE334588_series_matrix.txt.gz", headers=UA, timeout=300).content).decode("utf8","replace")
assert "!series_matrix_table_begin" in txt, "table marker missing"
head = txt[:txt.index("!series_matrix_table_begin")]
L = lambda k: [l for l in head.split("\n") if l.startswith(k)]
V = lambda k: [x.strip('"') for x in L(k)[0].split("\t")[1:]] if L(k) else []
ch = [[x.strip('"') for x in l.split("\t")[1:]] for l in head.split("\n") if l.startswith("!Sample_characteristics")]
gsm, tit = V("!Sample_geo_accession"), V("!Sample_title")
rec = [dict([("gsm",g),("title",t)]+[(c[i].split(": ")[0].strip().lower(), c[i].split(": ",1)[-1]) for c in ch if ": " in c[i]]) for i,(g,t) in enumerate(zip(gsm,tit))]
did = lambda t: re.sub(r"\s*-\s*batch\s*\d+\s*$","",t.split(",")[-1],flags=re.I).strip()
armf = lambda d: "Healthy" if re.fullmatch(r"H\d+",d,re.I) else ("CLE" if re.fullmatch(r"CLE\d+",d,re.I) else ("NL_CLE" if re.fullmatch(r"NL_CLE\d+",d,re.I) else ("DM" if re.fullmatch(r"DM\d+",d,re.I) else ("NL_DM" if re.fullmatch(r"NL_DM\d+",d,re.I) else "?"))))
for r in rec: r["donor"]=did(r["title"]); r["arm"]=armf(r["donor"]); r["run"]=(r.get("batch") or "").strip()
assert all(r["arm"]!="?" for r in rec) and all(r["run"] for r in rec), "cohort parse failed"
# ---- supplementary file PREFIXES, read from the series matrix. NOT GUESSED. ----
sup = {}
for l in L("!Sample_supplementary_file"):
    for i, x in enumerate([y.strip('"') for y in l.split("\t")[1:]]):
        if x and x != "NONE": sup.setdefault(i, []).append(x.rsplit("/",1)[-1])
pref = {}
for i, names in sup.items():
    base = [re.sub(r"_(barcodes|features|matrix)\.(tsv|mtx)\.gz$","",n) for n in names if re.search(r"_(barcodes|features|matrix)\.",n)]
    base = sorted(set(base))
    if len(base) == 1: pref[rec[i]["gsm"]] = re.sub(r"^%s_" % rec[i]["gsm"], "", base[0])
print("\nprefixes resolved for %d of %d samples" % (len(pref), len(rec)))
assert len(pref) == len(rec), "some sample has no resolvable 10x prefix; STOP rather than guess"
D = collections.defaultdict(set)
for r in rec: D[(r["arm"],r["donor"])].add(r["run"])
nC = len([1 for a,d in D if a=="CLE"]); nH = len([1 for a,d in D if a=="Healthy"])
assert (nC, nH) == (9, 5), "cohort is %r not the 9 CLE and 5 Healthy RESULT_14 verified" % ((nC,nH),)
print("cohort matches RESULT_14: CLE %d donors, Healthy %d donors" % (nC, nH))
# ---- Aditi's ruling: H1 drops. Then H2 lands in run2 which holds no CLE, so it leaves too. ----
TEST = {"CLE1","CLE2","CLE3","CLE4","CLE5","CLE6","CLE7","CLE8","CLE9","H3","H4","H5"}
WIDER = TEST | {"H2"}
sel = lambda keep: [dict(gsm=r["gsm"], prefix=pref[r["gsm"]], donor=r["donor"], arm=r["arm"], ctrl=(r["arm"]=="Healthy")) for r in rec if r["donor"] in keep]
print("\nTEST SET  (H1 dropped by ruling, H2 leaves with single-arm run2): %d donors, %d samples" % (len(TEST), len(sel(TEST))))
print("WIDER SET (H1 dropped only, H2 retained):                         %d donors, %d samples" % (len(WIDER), len(sel(WIDER))))
sbase = lambda g: "https://ftp.ncbi.nlm.nih.gov/geo/samples/" + g[:-3] + "nnn/" + g + "/suppl/"
def fetch(g, p, kind):
    fn = "%s_%s_%s.%s.gz" % (g, p, kind, "tsv" if kind=="features" else "mtx")
    return requests.get(sbase(g)+fn, headers=UA, timeout=600).content
def detect(records, label):
    print("\n--- %s: reading %d samples ---" % (label, len(records)))
    SL, DET, skipped = F.build_sample_table(records, fetch, H.GENES)
    print("read %d | skipped %d | cells %d | median lib %.1f" % (len(SL), len(skipped), int(SL["cells"].sum()), float(SL["lib"].median())))
    if skipped: print("SKIPPED:", skipped)
    DM_, LAB = H.donor_matrix_and_labels(SL)
    nz = {g: float((DM_.loc[g] > 0).mean()) if g in DM_.index else 0.0 for g in H.GENES}
    keep = lambda g: DET[g] >= MIN_CELL_DETECT and nz[g] >= MIN_DONOR_NONZERO
    print("  %-12s %9s %9s  %s" % ("gene","cell_det","donor_nz","verdict"))
    for g in H.GENES[:8]: print("  %-12s %9.4f %9.4f  %s" % (g, DET[g], nz[g], "keep" if keep(g) else "drop"))
    inh = [g for g in H.GENES[:8] if keep(g)]
    ost = [g for g in H.GENES[8:] if keep(g)]
    six = [g for g in H.FROZEN_SIX if keep(g)]
    print("  INHIBITOR %d of 8 -> %s" % (len(inh), inh))
    print("  FROZEN SIX %d of 6 -> %s" % (len(six), six))
    print("  OSTEOGENIC %d of 8 -> %s" % (len(ost), ost))
    print("  VERDICT: %s" % ("VOID, below MIN_SET_SURVIVORS" if len(inh) < MIN_SET_SURVIVORS else "NOT VOID, the inhibitor set is testable here"))
    return len(inh), len(ost)
iT, oT = detect(sel(TEST), "TEST SET, the 12 donors that enter the stratified test")
iW, oW = detect(sel(WIDER), "WIDER SET, 13 donors including H2, reported so the choice is visible")
print("\n*** PREDICTION ***")
print("I predicted detection CLEARS. Test set inhibitor survivors %d of 8. PREDICTION %s" % (iT, "CORRECT" if iT >= MIN_SET_SURVIVORS else "WRONG"))
print("\n⭐ THE OSTEOGENIC ARM IS THE ONE TO LOOK AT. It has NEVER been evaluable: void at 2 of 7 on")
print("   GPL18573 and not evaluable on GPL24676, so the composite criterion in section 6 of the")
print("   original pre-registration has never once been testable. Here it reads %d of 8." % oT)
print("\nSTOP. NOTHING SCORED. NO STATISTIC. NO PERMUTATION. A void is not a null (DEVIATION_03).")
