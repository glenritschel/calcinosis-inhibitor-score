# -*- coding: utf-8 -*-
"""THE ONE LOOK. Scores the held-out GPL24676 primary set under the frozen pre-registration.

⛔ THIS SCRIPT CONSUMES THE ONE LOOK AND CANNOT BE UNDONE. It refuses to run without the
explicit flag --i-am-taking-the-one-look, so it cannot be triggered by an import, a stray
exec or a rerun-all. Section 7 attaches a program-stopping decision to a failure here.

Governing documents, all recorded BEFORE this file could produce a number:
  PREREG_REPLICATION_v2_2026-10-21                 51f155ec   gene set, direction, VOID rule
  PREREG_REPLICATION_v3_SECTION3 + addendum        28e9d649 / 4e69cd91   counting basis
  PREREG_REPLICATION_v4_ADDENDUM_SECTION6          bee4f50a   spanning-donor rule
  PREREG_REPLICATION_v5_ADDENDUM_SECTION3          869b423d   Option B, the 50 donors
  RESULT_10_confound_test                          a51dc8e8   p=0.2945, UNRESTRICTED null
"""
import sys, io, zlib, time
import numpy as np, pandas as pd
import replication_harness as H
import inhibition_score as INS

N_PERM, SEED, DIRECTION = 20000, 20261020, "lower"
UA = {"User-Agent": "Mozilla/5.0 (research; ORCID 0009-0006-5341-8840)"}

def make_fetch(requests):
    base = lambda g: "https://ftp.ncbi.nlm.nih.gov/geo/samples/" + g[:-3] + "nnn/" + g + "/suppl/"
    def fetch(gsm, title):
        raw = zlib.decompress(requests.get(base(gsm) + gsm + "_" + title + ".txt.gz",
                                           headers=UA, timeout=240).content, 31)
        return pd.read_csv(io.BytesIO(raw), sep="\t", index_col=0)
    return fetch

def records_for(E24, primary, norm):
    """E24: the GPL24676 CD90+ skin sample records. primary: the 50 ids. norm: id normalizer."""
    keep = set(primary)
    isc = lambda d: d.lower().startswith("ctrl")
    return [dict(gsm=r["gsm"], title=r["title"], donor=norm(r.get("patient id","")),
                 ctrl=isc(norm(r.get("patient id",""))))
            for r in E24 if norm(r.get("patient id","")) in keep]

def run(E24, ids_24676, ids_18573, norm, requests):
    isc = lambda d: d.lower().startswith("ctrl")
    primary = H.option_b_primary(ids_24676, ids_18573)
    H.check_primary(primary, isc)
    print("PRIMARY SET RECONSTRUCTED FROM THE ADDENDUM: %d donors, %d control, %d patient"
          % (len(primary), sum(1 for d in primary if isc(d)), sum(1 for d in primary if not isc(d))))

    recs = records_for(E24, primary, norm)
    print("plates to read: %d" % len(recs))
    t0 = time.time()
    PL, DET, skipped = H.build_plate_table(recs, make_fetch(requests))
    print("read %d plates, %d cells, %.0fs; SKIPPED %d %r"
          % (len(PL), int(PL.cells.sum()), time.time()-t0, len(skipped), skipped[:5]))
    if skipped:
        raise H.HarnessError("%d plates failed to read; a partial held-out set is not the "
                             "pre-registered set and must not be scored" % len(skipped))

    DM, LAB = H.donor_matrix_and_labels(PL)
    got = (DM.shape[1], int(LAB.sum()), int((~LAB).sum()))
    print("donor matrix %r | control %d | patient %d" % (DM.shape, got[1], got[2]))
    if got != (50, 25, 25):
        raise H.HarnessError("donor matrix is %r, not 50 at 25 and 25" % (got,))

    print("\nCELL DETECTION RATES")
    for g in H.GENES: print("   %-10s %.4f" % (g, DET[g]))

    print("\nPRIMARY, THE FROZEN SIX")
    try:
        used, dropped = INS.guard_detected(DM, H.FROZEN_SIX, DET)
        print("   usable %r" % (used,))
        for k, v in dropped.items(): print("   dropped %s -> %s" % (k, v))
    except INS.Void as e:
        print("   *** PRIMARY IS VOID: %s" % e)
        print("   A VOID PRIMARY IS NOT A FAILURE TO REPLICATE. It states that this platform")
        print("   cannot measure the set. The native-survivor secondary is then the only result.")
        used = None

    if used:
        sc, cov, genes_used = INS.set_score(DM, used)
        dis = [d for d in DM.columns if not bool(LAB[d])]
        print("   GUARD 1 can_vary  SD %.4f | coverage %.2f | genes %r"
              % (INS.guard_can_vary(sc), cov, genes_used))
        print("   GUARD 2 responds  moved %.3f against a 0.25 floor"
              % INS.guard_responds(DM, used, dis))
        print("   GUARD 3 reachable %r" % (INS.guard_reachable(len(dis), DM.shape[1]-len(dis)),))
        # label_null returns a dict: observed, p_one_sided, p_two_sided, direction, n_perm.
        # is_disease is aligned to DM's column order explicitly; relying on the Series order
        # matching the matrix is the kind of silent mismatch that produces a confident number.
        is_dis = (~LAB.reindex(DM.columns).astype(bool)).values
        if int(is_dis.sum()) != 25 or int((~is_dis).sum()) != 25:
            raise H.HarnessError("label vector is %d disease and %d control after alignment"
                                 % (is_dis.sum(), (~is_dis).sum()))
        res = INS.label_null(sc.reindex(DM.columns).values, is_dis,
                             n_perm=N_PERM, direction=DIRECTION, seed=SEED)
        print("   LABEL NULL n_perm %d seed %d direction %s" % (N_PERM, SEED, DIRECTION))
        print("   observed difference, disease minus control : %+.4f" % res["observed"])
        print("   one-sided p, pre-registered direction      : %.5f" % res["p_one_sided"])
        print("   two-sided p                                : %.5f" % res["p_two_sided"])
        print("   mean score control %+.4f | disease %+.4f"
              % (sc[[d for d in DM.columns if bool(LAB[d])]].mean(), sc[dis].mean()))
        print("\n   ⛔ THIS IS THE INHIBITOR CLAUSE ONLY. Section 6 of the original")
        print("   pre-registration requires BOTH that the inhibitor set is lower AND that the")
        print("   osteogenic set is not elevated. On GPL18573 the osteogenic set was VOID at")
        print("   2 of 7 genes, so that clause could not be evaluated there either. Reporting")
        print("   this as a model being supported would be quoting half a criterion.")

    print("\nSECONDARY, GPL24676-NATIVE SURVIVORS OF THE FULL INHIBITOR SET")
    try:
        nat, natdrop = INS.guard_detected(DM, INS.INHIBITOR, DET)
        print("   native survivors %r" % (nat,))
        for k, v in natdrop.items(): print("   dropped %s -> %s" % (k, v))
    except INS.Void as e:
        print("   secondary VOID: %s" % e)

    print("\nTHE ONE LOOK IS NOW SPENT.")
    return PL, DET, DM, LAB

if __name__ == "__main__":
    if "--i-am-taking-the-one-look" not in sys.argv:
        sys.exit("REFUSED. This script consumes the one look. Pass --i-am-taking-the-one-look "
                 "to run it, and only on Glen's explicit word.")
    sys.exit("REFUSED here: run() must be called from the notebook, where E24, the donor id "
             "sets, the id normalizer and requests already exist. This guard exists so the "
             "file cannot score anything by being executed.")
