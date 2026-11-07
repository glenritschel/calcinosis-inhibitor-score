"""STEP 5, GLEN'S SIDE. The code-to-score table, from the blinded matrix.

⛔ THIS SCRIPT REIMPLEMENTS NOTHING. It imports set_score from inhibition_score, which is the
committed scoring function that produced the discovery numbers, and that is the whole point:
Mararo's independent reimplementation is a test OF this code, so this side must run this code.
Prereg section 4 is titled "Procedure, reusing code that already exists" and governs this side.

⛔ NO DONOR LABEL IS READ AND NO ARM CONTRAST IS COMPUTED. Input is the blinded matrix with
D01-D15 codes. The look is not spent here.
"""
import sys, os, hashlib, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inhibition_score import set_score, INHIBITOR, MIN_SET_SURVIVORS, Void
from replication_harness import FROZEN_SIX

MATRIX = "GSE292979_donor_matrix_BLINDED_log1p_fetched_2026-10-29_2156.csv"

def main():
    raw = open(MATRIX,"rb").read()
    print("matrix   %s" % MATRIX)
    print("sha256   %s" % hashlib.sha256(raw).hexdigest())
    dm = pd.read_csv(MATRIX, index_col=0)
    dm.index.name = None

    # ---- guards on the input, before anything is scored -------------------------------
    assert dm.shape[1] == 15, "expected 15 donor columns, got %d" % dm.shape[1]
    cols = list(dm.columns)
    assert cols == ["D%02d" % i for i in range(1,16)], "columns are not D01..D15: %s" % cols
    # the blinding must still be intact: no real donor identifier anywhere in the header
    for bad in ("Pt","HC","GSM","ctrl","control","patient"):
        assert not any(bad.lower() in c.lower() for c in cols), "LABEL LEAK in columns: %s" % bad
    assert dm.notna().all().all(), "the matrix has missing values"
    print("columns  %s" % " ".join(cols))
    print("genes    %d rows: %s" % (dm.shape[0], " ".join(dm.index)))
    print()

    # ---- ⛔ THE DETECTION GATE IS NOT set_score's JOB AND MUST BE SUPPLIED ------------
    # set_score's own docstring says guard_detected enforces the DETECTION gates and that
    # set_score only drops zero-variance genes. An earlier version of THIS script passed the
    # full set and let set_score score every gene present in the matrix, which scored the
    # secondary on eight genes when step 4b had cleared only seven. Same shape as the defect
    # Mararo found as his point 5: a floor enforced in one place and assumed in another.
    # MIN_CELL_DETECT needs cell-level data, which this matrix does not carry, so the survivor
    # list CANNOT be recomputed here and must come from the 4b run. No default, no guess.
    DETECTION_SURVIVORS = {
        # set name -> list of genes that cleared BOTH detection gates in the 4b run,
        # or None when the per-gene outcome is not in hand.
        "PRIMARY frozen six": list(FROZEN_SIX),   # 4b: 6 of 6 CLEAR, recorded in protocol step 4b
        # ⭐ RESOLVED 2026-10-30 20:5x. ABCC6 IS THE GENE THAT FAILED, out at 0.2058 percent
        # against the 1 percent MIN_CELL_DETECT floor, recorded in section 5 of
        # PREREG_GSE292979_v1_ADDENDUM_SENSITIVITY_GSM8872105_2026-10-29.md.
        # ⛔ I TWICE WROTE THAT THIS WAS "NOT ON DISK ANYWHERE IN THE LEDGER". It was in a
        # file I had already read and quoted from the same day. All 15 genes pass the DONOR
        # gate (checked against the raw blinded matrix), so the failure had to be on
        # MIN_CELL_DETECT, and that is where the recorded figure sits.
        "SECONDARY inhibitor eight": [g for g in INHIBITOR if g != "ABCC6"],
    }

    out = {}
    for name, genes, floor in (("PRIMARY frozen six", FROZEN_SIX, 5),
                               ("SECONDARY inhibitor eight", INHIBITOR, 5)):
        print("=" * 78)
        print("%s   floor %d of %d" % (name, floor, len(genes)))
        missing = [g for g in genes if g not in dm.index]
        if missing:
            print("  NOT IN MATRIX: %s" % " ".join(missing))
        surv = DETECTION_SURVIVORS[name]
        if surv is None:
            print("  ⛔ NOT SCORED. The detection-gate outcome for this set is not on disk.")
            print("     Step 4b cleared 7 of 8 and WHICH gene failed is only in the Colab run")
            print("     output. Scoring all eight here would score a gene the gates excluded,")
            print("     and the disagreement with Mararo would look like his bug, not mine.")
            out[name] = {"not_scored": "detection-gate outcome for this set not in hand"}
            continue
        assert set(surv) <= set(genes), "survivor list is not a subset of the set"
        assert len(surv) >= floor, "survivor list is below the floor; the set is VOID at 4b"
        print("  detection survivors supplied from 4b: %d of %d (%s)"
              % (len(surv), len(genes), ", ".join(surv)))
        try:
            scores, coverage, used = set_score(dm, surv, min_survivors=floor)
        except Void as e:
            print("  VOID: %s" % e)
            print("  ⛔ VOID IS NOT A NULL AND IS NOT REPORTED AS ONE.")
            out[name] = {"void": str(e)}
            continue
        print("  genes reaching the mean: %d of %d  (%s)" % (len(used), len(genes), ", ".join(used)))
        print("  coverage: %.4f" % coverage)
        print("  between-donor SD of the score: %.12g" % float(np.std(scores.values, ddof=1)))
        print()
        print("  code   score")
        for c in cols:
            print("  %s    %.12g" % (c, scores[c]))
        out[name] = {"detection_survivors": surv, "used": used, "coverage": coverage,
                     "scores": {c: float(scores[c]) for c in cols}}
        print()

    # ---- what must NOT be in this output ----------------------------------------------
    blob = json.dumps(out)
    for bad in ("Pt1","HC1","GSM","is_control","ctrl"):
        assert bad not in blob, "LEAK in output: %s" % bad
    with open("STEP5_GLEN_code_to_score_2026-10-30.json","w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
    print("=" * 78)
    print("wrote STEP5_GLEN_code_to_score_2026-10-30.json")

if __name__ == "__main__":
    main()
