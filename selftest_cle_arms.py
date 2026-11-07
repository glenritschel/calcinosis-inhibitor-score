"""Fixture for the CLE control's ARM SELECTION and STRATIFICATION path.

flex_harness has been exercised on GSE334588 for DETECTION (RESULT_05) but never for a two-arm
contrast. This pins the part PREREG_CLE_CONTROL_v1_2026-10-23 actually relies on, with no network
and no expression data: which samples become which donors, which donors enter the test, how they
map to strata, and what the stratified null can reach.

The donor-to-run map below is METADATA VERIFIED IN RESULT_14 against the 'batch' characteristic,
cross-checked against the independent 'disease state' field for all 48 samples.
"""
import re, collections
import numpy as np
import stratified_null as SN

FAILS = []
def chk(label, got, want):
    ok = got == want
    print("  %-4s %-62s got %r" % ("ok" if ok else "FAIL", label, got) + ("" if ok else "  want %r" % (want,)))
    if not ok: FAILS.append(label)
def raises(label, fn, exc=Exception):
    try:
        fn(); print("  FAIL %-62s did not raise" % label); FAILS.append(label)
    except exc as e:
        print("  ok   %-62s refused: %s" % (label, type(e).__name__))

# ---- WHAT IS VERIFIED, AND WHAT IS DELIBERATELY ABSENT -------------------------------------
# ⛔ THE FIRST VERSION OF THIS FIXTURE FABRICATED THE DM RUN ASSIGNMENTS to pad the sheet to 48.
# RESULT_14 verified the run for the CLE and Healthy arms ONLY; the per-sample run of the DM and
# NL arms was never printed and I invented it. The fixture caught it on the run census.
# A FIXTURE MAY NOT CONTAIN INVENTED DATA, so the DM and NL arms appear here WITHOUT RUNS: they
# exist to test ARM EXCLUSION, which is all ruling 3 needs them for, and the test never reads
# their run. The run census assert is gone with them.
CLE_HEALTHY = (
    [("Cutaneous lupus erythematosus lesional skin, CLE%d" % i, "run3") for i in (1, 2, 3, 4)]
    + [("Cutaneous lupus erythematosus lesional skin, CLE%d" % i, "run4") for i in (7, 8, 9)]
    + [("Cutaneous lupus erythematosus lesional skin, CLE%d" % i, "run5") for i in (5, 6)]
    + [("Healthy control skin, H1 - batch1", "run1"), ("Healthy control skin, H1 - batch3", "run3"),
       ("Healthy control skin, H1 - batch4", "run4"), ("Healthy control skin, H1 - batch5", "run5"),
       ("Healthy control skin, H2 - batch2", "run2"), ("Healthy control skin, H2 - batch3", "run3"),
       ("Healthy control skin, H3 - batch3", "run3"), ("Healthy control skin, H4 - batch 4", "run4"),
       ("Healthy control skin, H5 - batch5", "run5")]
)
EXCLUDED_NO_RUN = (
    [("Dermatomyositis lesional skin, DM%d" % i, None) for i in range(1, 23)]
    + [("Cutaneous lupus erythematosus non-lesional skin, NL_CLE%d" % i, None) for i in range(1, 7)]
    + [("Dermatomyositis non-lesional skin, NL_DM%d" % i, None) for i in (1, 2)]
)
SAMPLES = CLE_HEALTHY + EXCLUDED_NO_RUN

donor_of = lambda t: re.sub(r"\s*-\s*batch\s*\d+\s*$", "", t.split(",")[-1], flags=re.I).strip()
arm_of = lambda d: ("Healthy" if re.fullmatch(r"H\d+", d, re.I) else
                    "CLE" if re.fullmatch(r"CLE\d+", d, re.I) else
                    "NL_CLE" if re.fullmatch(r"NL_CLE\d+", d, re.I) else
                    "DM" if re.fullmatch(r"DM\d+", d, re.I) else
                    "NL_DM" if re.fullmatch(r"NL_DM\d+", d, re.I) else "?")

print("1. THE SAMPLE SHEET ITSELF")
chk("48 samples", len(SAMPLES), 48)
chk("no donor id fails to classify", sum(1 for t, _ in SAMPLES if arm_of(donor_of(t)) == "?"), 0)
chk("CLE and Healthy samples carry a verified run", sum(1 for _, r in CLE_HEALTHY if r), 18)
chk("excluded arms carry NO run, because theirs was never verified",
    sum(1 for _, r in EXCLUDED_NO_RUN if r), 0)
chk("CLE and Healthy run census, verified in RESULT_14",
    dict(collections.Counter(r for _, r in CLE_HEALTHY)),
    {"run3": 7, "run4": 5, "run5": 4, "run1": 1, "run2": 1})
# run1 H1 | run2 H2 | run3 CLE1-4 + H1,H2,H3 | run4 CLE7-9 + H1,H4 | run5 CLE5,6 + H1,H5
# Consistent with the full 48 census in RESULT_14: run4 is 5 samples in total and all five are
# CLE or Healthy, so no excluded-arm sample is in run4 at all.

print("\n2. DONOR COLLAPSE, WHICH IS DEVIATION_01's DEFECT")
D = collections.defaultdict(set)
for t, r in SAMPLES:
    if r: D[donor_of(t)].add(r)
    else: D[donor_of(t)]
chk("H1 is ONE donor spanning four runs, not four donors", sorted(D["H1"]), ["run1", "run3", "run4", "run5"])
chk("H2 is ONE donor spanning two runs", sorted(D["H2"]), ["run2", "run3"])
chk("Healthy donors", sum(1 for d in D if arm_of(d) == "Healthy"), 5)
chk("CLE donors", sum(1 for d in D if arm_of(d) == "CLE"), 9)
chk("counting SAMPLES as donors would give 9 healthy, which is the defect",
    sum(1 for t, _ in SAMPLES if arm_of(donor_of(t)) == "Healthy"), 9)

print("\n3. ARM SELECTION, AND THE SUBSTRING HAZARD")
chk("NL_CLE does NOT classify as CLE", arm_of("NL_CLE1"), "NL_CLE")
chk("NL_DM does NOT classify as DM", arm_of("NL_DM1"), "NL_DM")
chk("DM excluded by ruling 3", sum(1 for d in D if arm_of(d) == "DM"), 22)
sel = sorted(d for d in D if arm_of(d) in ("CLE", "Healthy"))
chk("the test takes CLE and Healthy only, 14 donors", len(sel), 14)

print("\n4. STRATIFICATION AFTER ADITI'S RULING, H1 NOT FORCED")
keep = [d for d in sel if d != "H1"]
asg = {d: (sorted(D[d])[0] if len(D[d]) == 1 else sorted(collections.Counter(sorted(D[d])).items())[0][0]) for d in keep}
arm = np.array([arm_of(d) for d in keep], dtype=object)
strata = np.array([asg[d] for d in keep], dtype=object)
chk("H1 is out", "H1" in keep, False)
chk("H2 lands in run2", asg["H2"], "run2")
n, minp = SN.guard_reachable_stratified(arm, strata, "CLE", alpha=0.05)
chk("labelings", n, 60)
chk("minimum attainable p", round(minp, 6), 0.016667)
chk("REACHABLE at 0.05", minp <= 0.05, True)

print("\n5. THE SINGLE-ARM STRATUM IS DROPPED, NOT AN ERROR")
s, a, b, dropped, arms = SN._prepare(np.zeros(len(keep)), arm, strata, "CLE")
chk("run2 dropped as single-arm", "run2" in dropped, True)
chk("three strata survive", sorted(set(b.tolist())), ["run3", "run4", "run5"])
chk("12 donors survive: 14 minus H1 by the ruling, minus H2 with single-arm run2", len(s), 12)

print("\n6. CASES THAT MUST BE REFUSED")
raises("a contrast with only one arm present is REFUSED",
       lambda: SN._prepare(np.zeros(3), np.array(["CLE"] * 3, dtype=object), np.array(["run3"] * 3, dtype=object), "CLE"))
raises("a target that is not one of the arms is REFUSED",
       lambda: SN._prepare(np.zeros(2), np.array(["CLE", "Healthy"], dtype=object), np.array(["run3"] * 2, dtype=object), "DM"))
raises("mismatched lengths are REFUSED",
       lambda: SN._prepare(np.zeros(3), np.array(["CLE", "Healthy"], dtype=object), np.array(["run3"] * 2, dtype=object), "CLE"))
raises("every stratum single-arm is REFUSED, not silently empty",
       lambda: SN._prepare(np.zeros(2), np.array(["CLE", "Healthy"], dtype=object), np.array(["run3", "run4"], dtype=object), "CLE"))

print("\n7. FORCING H1 INTO run1 WOULD NOT CHANGE THE TEST, WHICH IS WHY THE RULING IS SAFE")
k2 = sel
a2 = np.array([arm_of(d) for d in k2], dtype=object)
s2 = np.array([("run1" if d == "H1" else asg.get(d, sorted(D[d])[0])) for d in k2], dtype=object)
n2, _ = SN.guard_reachable_stratified(a2, s2, "CLE", alpha=0.05)
chk("H1 forced into run1 still gives 60, because run1 holds no CLE donor", n2, 60)

print("\n%s  (%d failure%s)" % ("SELFTEST PASS" if not FAILS else "SELFTEST FAILED: " + ", ".join(FAILS),
                                len(FAILS), "" if len(FAILS) == 1 else "s"))
