# DETECTION ONLY, GPL24676 CD45+ immune plates, the Deviation-2 specificity check.
# RESULT_11 section 4 named this as not run on this platform. RESULT_12 fixed a PREDICTION
# before this cell was written: I predict it VOIDS, as RESULT_04 did on GPL18573, because
# RESULT_04 had 86 plates and 33,009 cells and still cleared only 3 of 8, and this has 36 plates.
# IT STOPS AFTER DETECTION. NOTHING IS SCORED. NO STATISTIC IS COMPUTED.
# Requires replication_harness.py and donor_group.py uploaded to this Colab session.
import io, zlib, hashlib, collections, gzip, re
import numpy as np, pandas as pd, requests
import replication_harness as H, donor_group as DG
# --- INTEGRITY: the modules must be the ones the ledger records, not some other upload ---
hv = lambda f: hashlib.sha256(open(f,"rb").read()).hexdigest()[:16]
EXPECT = {"replication_harness.py": "9dd1fa94c7807225", "donor_group.py": "f61d761a8a7f72a1"}
got = {f: hv(f) for f in EXPECT}
assert got == EXPECT, "UPLOADED MODULE HASHES DIFFER FROM THE LEDGER. got %r expected %r. STOP." % (got, EXPECT)
print("module integrity OK:", got)
print("FROZEN_SIX from the harness:", H.FROZEN_SIX)
print("full inhibitor set:", H.GENES[:8])
# --- thresholds, FIXED IN DEVIATION_03 BEFORE ANY DATA WAS READ, restated not reinvented ---
MIN_CELL_DETECT, MIN_DONOR_NONZERO, MIN_SET_SURVIVORS = 0.01, 0.50, 5
print("DEVIATION_03 thresholds: MIN_CELL_DETECT %.2f | MIN_DONOR_NONZERO %.2f | MIN_SET_SURVIVORS %d" % (MIN_CELL_DETECT, MIN_DONOR_NONZERO, MIN_SET_SURVIVORS))
# --- cohort, rebuilt exactly as the preflight did ---
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE195nnn/GSE195452/matrix/"
UA = {"User-Agent": "Mozilla/5.0 (research; ORCID 0009-0006-5341-8840)"}
get = lambda fn: gzip.decompress(requests.get(BASE+fn, headers=UA, timeout=300).content).decode("utf8","replace")
SM = {p: get("GSE195452-%s_series_matrix.txt.gz" % p) for p in ("GPL24676","GPL18573")}
assert all("!series_matrix_table_begin" in t for t in SM.values()), "table marker missing; refusing to parse"
cut = lambda t: t[:t.index("!series_matrix_table_begin")]
row = lambda t,k: [x.strip('"') for x in [l for l in cut(t).split("\n") if l.startswith(k)][0].split("\t")[1:]]
chs = lambda t: [[x.strip('"') for x in l.split("\t")[1:]] for l in cut(t).split("\n") if l.startswith("!Sample_characteristics")]
build = lambda t: [dict([("gsm",g),("title",ti)] + [(c[i].split(": ")[0].strip().lower(), c[i].split(": ",1)[-1]) for c in chs(t) if ": " in c[i]]) for i,(g,ti) in enumerate(zip(row(t,"!Sample_geo_accession"), row(t,"!Sample_title")))]
REC24 = build(SM["GPL24676"]); REC18 = build(SM["GPL18573"])
nm = lambda r: DG._norm(r.get("patient id",""))
ex = lambda R,mk: [r for r in R if r.get("tissue") == "Skin" and r.get("selection marker") == mk]
E24c = ex(REC24,"CD90+"); E18c = ex(REC18,"CD90+")
ids24 = {nm(r) for r in E24c}; ids18 = {nm(r) for r in E18c}
assert (len(E24c),len(ids24)) == (159,71) and (len(E18c),len(ids18)) == (184,101), "CD90+ does not reconstruct; STOP"
isc = lambda d: d.lower().startswith("ctrl")
PRIMARY = set(H.option_b_primary(ids24, ids18))
H.check_primary(sorted(PRIMARY), isc)
I24 = ex(REC24,"CD45+")
RECS = [dict(gsm=r["gsm"], title=r["title"], donor=nm(r), ctrl=isc(nm(r))) for r in I24 if nm(r) in PRIMARY]
assert len(RECS) == 36, "expected the 36 CD45+ plates RESULT_12 counted, got %d; STOP" % len(RECS)
print("\ncohort matches RESULT_12: %d CD45+ plates, %d donors, %d control %d patient" % (len(RECS), len({r["donor"] for r in RECS}), len({r["donor"] for r in RECS if r["ctrl"]}), len({r["donor"] for r in RECS if not r["ctrl"]})))
# --- fetch, identical to the runner's ---
sbase = lambda g: "https://ftp.ncbi.nlm.nih.gov/geo/samples/" + g[:-3] + "nnn/" + g + "/suppl/"
fetch = lambda gsm, title: pd.read_csv(io.BytesIO(zlib.decompress(requests.get(sbase(gsm) + gsm + "_" + title + ".txt.gz", headers=UA, timeout=240).content, 31)), sep="\t", index_col=0)
print("\nreading %d plates. detection only." % len(RECS))
PL, DET, skipped = H.build_plate_table(RECS, fetch)
print("plates read %d | skipped %d | total cells %d" % (len(PL), len(skipped), int(PL["cells"].sum())))
if skipped: print("SKIPPED:", skipped)
print("median library size across plates: %.1f" % float(PL["lib"].median()))
# --- the guard, per gene ---
DM, LAB = H.donor_matrix_and_labels(PL)
nzfrac = {g: float((DM.loc[g] > 0).mean()) if g in DM.index else 0.0 for g in H.GENES}
print("\nper-gene detection, CD45+ immune, GPL24676")
print("  %-12s %10s %10s  %s" % ("gene","cell_det","donor_nz","verdict"))
verdict = lambda g: "keep" if (DET[g] >= MIN_CELL_DETECT and nzfrac[g] >= MIN_DONOR_NONZERO) else "drop"
for g in H.GENES[:8]: print("  %-12s %10.4f %10.4f  %s" % (g, DET[g], nzfrac[g], verdict(g)))
surv = [g for g in H.GENES[:8] if verdict(g) == "keep"]
six = [g for g in H.FROZEN_SIX if verdict(g) == "keep"]
print("\nINHIBITOR SET survivors: %d of 8 -> %s" % (len(surv), surv))
print("FROZEN SIX survivors:    %d of 6 -> %s" % (len(six), six))
print("\n*** VERDICT ***")
print("VOID: %d survivors is below MIN_SET_SURVIVORS of %d. The specificity check CANNOT be run" % (len(surv), MIN_SET_SURVIVORS) if len(surv) < MIN_SET_SURVIVORS else "NOT VOID: %d survivors meets MIN_SET_SURVIVORS of %d. The check IS runnable, 9 control against 21 patient." % (len(surv), MIN_SET_SURVIVORS))
print("PREDICTION IN RESULT_12 WAS **%s**" % ("CORRECT" if len(surv) < MIN_SET_SURVIVORS else "WRONG"))
print("\nSTOP HERE. NOTHING HAS BEEN SCORED. NO STATISTIC COMPUTED. NO PERMUTATION RUN.")
print("A void is not a null and must never be reported as one (DEVIATION_03).")
