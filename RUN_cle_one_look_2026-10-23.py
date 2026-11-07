"""⛔ THIS SCRIPT CONSUMES THE ONE LOOK ON THE CLE CONTROL AND CANNOT BE UNDONE.

Section 5 of PREREG_CLE_CONTROL_v1_2026-10-23 binds on whatever this returns:
  IF p < 0.05 IN THE LOWER DIRECTION, THE LINE IS DEAD. No reanalysis, no subgroup, no second
  cohort may rescue it, because permitting one is how a stopping rule becomes decoration.
  IF p IS NOT BELOW 0.05, specificity is NOT DEMONSTRATED. Only that a deficit was not detected.

Governing documents, all fixed BEFORE this file could produce a number:
  PREREG_CLE_CONTROL_v1_2026-10-23   48295fa514dbbced  criterion, direction, rulings binding
  RESULT_14_cle_reachability          dbdee15a5d2bf775  60 labelings, min p 0.016667
  RESULT_15_cle_detection             dac5cac8066d4567  6 of 8 inhibitors, frozen six 6 of 6
  selftest_cle_arms.py                1251a0e3d5686416  arm selection and strata, 22 checks
  DEVIATION_03                                          detection floors, VOID rule
"""
import gzip, re, collections, hashlib
import numpy as np, pandas as pd, requests
import replication_harness as H, donor_group as DG, flex_harness as F
import inhibition_score as INS, stratified_null as SN
EXPECT = {"replication_harness.py":"9dd1fa94c7807225","donor_group.py":"f61d761a8a7f72a1",
          "flex_harness.py":"512fe8ab13b7460f","inhibition_score.py":None,"stratified_null.py":None}
hv = lambda f: hashlib.sha256(open(f,"rb").read()).hexdigest()[:16]
for f, want in EXPECT.items():
    got = hv(f)
    assert want is None or got == want, "%s is %s not %s" % (f, got, want)
print("modules present and the three pinned ones match the ledger")
UA = {"User-Agent": "Mozilla/5.0 (research; ORCID 0009-0006-5341-8840)"}
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE334nnn/GSE334588/matrix/"
txt = gzip.decompress(requests.get(BASE+"GSE334588_series_matrix.txt.gz", headers=UA, timeout=300).content).decode("utf8","replace")
assert "!series_matrix_table_begin" in txt
head = txt[:txt.index("!series_matrix_table_begin")]
L = lambda k: [l for l in head.split("\n") if l.startswith(k)]
V = lambda k: [x.strip('"') for x in L(k)[0].split("\t")[1:]] if L(k) else []
ch = [[x.strip('"') for x in l.split("\t")[1:]] for l in head.split("\n") if l.startswith("!Sample_characteristics")]
gsm, tit = V("!Sample_geo_accession"), V("!Sample_title")
rec = [dict([("gsm",g),("title",t)]+[(c[i].split(": ")[0].strip().lower(), c[i].split(": ",1)[-1]) for c in ch if ": " in c[i]]) for i,(g,t) in enumerate(zip(gsm,tit))]
did = lambda t: re.sub(r"\s*-\s*batch\s*\d+\s*$","",t.split(",")[-1],flags=re.I).strip()
af = lambda d: "Healthy" if re.fullmatch(r"H\d+",d,re.I) else ("CLE" if re.fullmatch(r"CLE\d+",d,re.I) else ("NL_CLE" if re.fullmatch(r"NL_CLE\d+",d,re.I) else ("DM" if re.fullmatch(r"DM\d+",d,re.I) else ("NL_DM" if re.fullmatch(r"NL_DM\d+",d,re.I) else "?"))))
for r in rec: r["donor"]=did(r["title"]); r["arm"]=af(r["donor"]); r["run"]=(r.get("batch") or "").strip()
assert all(r["arm"]!="?" and r["run"] for r in rec), "cohort parse failed"
sup = {}
for l in L("!Sample_supplementary_file"):
    for i,x in enumerate([y.strip('"') for y in l.split("\t")[1:]]):
        if x and x!="NONE": sup.setdefault(i,[]).append(x.rsplit("/",1)[-1])
pref = {}
for i,names in sup.items():
    b = sorted({re.sub(r"_(barcodes|features|matrix)\.(tsv|mtx)\.gz$","",n) for n in names if re.search(r"_(barcodes|features|matrix)\.",n)})
    if len(b)==1: pref[rec[i]["gsm"]] = re.sub(r"^%s_"%rec[i]["gsm"],"",b[0])
assert len(pref)==len(rec), "a prefix did not resolve; STOP rather than guess"
D = collections.defaultdict(set)
for r in rec: D[(r["arm"],r["donor"])].add(r["run"])
assert (len([1 for a,d in D if a=="CLE"]), len([1 for a,d in D if a=="Healthy"]))==(9,5), "cohort is not RESULT_14's"
TEST = {"CLE1","CLE2","CLE3","CLE4","CLE5","CLE6","CLE7","CLE8","CLE9","H3","H4","H5"}
recs = [dict(gsm=r["gsm"], prefix=pref[r["gsm"]], donor=r["donor"], arm=r["arm"], ctrl=(r["arm"]=="Healthy")) for r in rec if r["donor"] in TEST]
assert len(recs)==12, "test set is %d samples not 12" % len(recs)
print("test set: 12 donors, 12 samples, H1 dropped by the ruling and H2 with single-arm run2")
sbase = lambda g: "https://ftp.ncbi.nlm.nih.gov/geo/samples/"+g[:-3]+"nnn/"+g+"/suppl/"
fetch = lambda g,p,kind: requests.get(sbase(g)+"%s_%s_%s.%s.gz"%(g,p,kind,"tsv" if kind=="features" else "mtx"), headers=UA, timeout=600).content
SL, DET, skipped = F.build_sample_table(recs, fetch, H.GENES)
assert not skipped, "a sample was skipped: %r" % (skipped,)
DM, LAB = H.donor_matrix_and_labels(SL)
nz = {g: float((DM.loc[g]>0).mean()) if g in DM.index else 0.0 for g in H.GENES}
keep = lambda g: DET[g] >= INS.MIN_CELL_DETECT and nz[g] >= INS.MIN_DONOR_NONZERO
six = [g for g in H.FROZEN_SIX if keep(g)]
print("GUARD 4 detection: frozen six %d of 6 -> %s" % (len(six), six))
assert len(six)==6, "detection does not reproduce RESULT_15's 6 of 6; STOP, do not score"
scores, cov, used = INS.set_score(DM, H.FROZEN_SIX)
sd = INS.guard_can_vary(scores); print("GUARD 1 can_vary: between-donor SD %.4f" % sd)
cle = [d for d in scores.index if af(d)=="CLE"]
moved = INS.guard_responds(DM, H.FROZEN_SIX, cle); print("GUARD 2 responds: injected 1.0 SD moved %.4f, floor 0.25" % moved)
runof = {d: sorted(D[(af(d),d)])[0] for d in scores.index}
arm = np.array([af(d) for d in scores.index], dtype=object)
strata = np.array([runof[d] for d in scores.index], dtype=object)
nlab, minp = SN.guard_reachable_stratified(arm, strata, "CLE", alpha=0.05)
print("GUARD 3 reachable: %d labelings, min attainable p %.6f" % (nlab, minp))
assert nlab==60, "labelings %d not the 60 RESULT_14 fixed" % nlab
print("\nALL FOUR GUARDS PASSED. SCORING NOW. THIS IS THE ONE LOOK.\n")
res = SN.stratified_label_null(scores, arm, strata, "CLE", direction="lower", n_perm=20000, seed=20261020)
print("observed, CLE minus Healthy  %.4f" % res["observed"])
print("one-sided p, LOWER           %.6f" % res["p_one_sided"])
print("two-sided p                  %.6f" % res["p_two_sided"])
print("null size %d | exact enumeration %s | labelings %d | donors used %d" % (res["n_null"], res["exact"], res["n_labelings"], res["n_used"]))
print("dropped strata %r" % (res["dropped_strata"],))
print("\nCLE mean %.4f | Healthy mean %.4f" % (float(scores[[d for d in scores.index if af(d)=='CLE']].mean()), float(scores[[d for d in scores.index if af(d)=='Healthy']].mean())))
print("\nFOR COMPARISON, THE DISEASE FIGURES (effect size is the reported quantity, section 1):")
print("  GPL18573 discovery  -0.3637   GPL24676 held out  -0.2939")
print("\n*** SECTION 5 STOPPING RULE ***")
print(("⛔ p < 0.05 IN THE LOWER DIRECTION. THE LINE IS DEAD. The inhibitor deficit is not specific\n"
       "   to a calcifying disease. No reanalysis, no subgroup, no second cohort may rescue it.")
      if res["p_one_sided"] < 0.05 else
      ("p IS NOT BELOW 0.05. The line continues. SPECIFICITY IS NOT DEMONSTRATED, only that a\n"
       "   deficit was not detected at this n. Section 1 bars reporting this as specificity shown."))
print("\nTHE ONE LOOK IS SPENT. The CLE arm has been scored and is not scored again.")
