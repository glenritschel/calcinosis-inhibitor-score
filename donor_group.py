"""
Donor grouping for PREREG_INHIBITION_v1 and PREREG_REPLICATION_v1.

Separate from inhibition_score.py on purpose: PREREG_REPLICATION_v1 commits that module as
unchanged at sha256 16284d6d, so the RESULT_08 fix lives here instead of modifying it.

THE DEFECT THIS EXISTS TO PREVENT (RESULT_08, 2026-10-20):
  GSE195452 GPL18573 CD90+ skin contains two plates whose patient id is the literal string
  "-", meaning the donor was not recorded. Grouping by that string created a donor that
  does not exist, pooling two unrelated plates into one. The disease test was
  str.lower().startswith("ctrl"), which is False for "-", so the fabricated donor was
  silently assigned to the DISEASE arm. It became the 72nd of 72 "patients".

  The failure was not the placeholder. It was that an unmatched value had somewhere to go.
  A classification rule with no fallback does not fail loudly; it files the exception under
  whatever the default happens to be. This module raises instead.
"""
import re

PLACEHOLDERS = {"", "-", "--", "n/a", "na", "none", "null", "nan", "unknown", "not recorded", "?"}

class UnrecordedDonor(Exception):
    pass

class AmbiguousLabel(Exception):
    pass

def _norm(x):
    return re.sub(r"\s+", " ", str(x)).strip()

def check_donor_ids(ids):
    """Raise if any donor id is missing or a placeholder. Returns the cleaned ids."""
    out = []
    bad = {}
    for i, x in enumerate(ids):
        s = _norm(x)
        if s.lower() in PLACEHOLDERS:
            bad.setdefault(s, []).append(i)
        out.append(s)
    if bad:
        raise UnrecordedDonor(
            "unrecorded donor id(s) present; grouping would fabricate a donor: "
            + "; ".join("%r at rows %s" % (k, v[:8] + (["..."] if len(v) > 8 else []))
                        for k, v in sorted(bad.items()))
            + ". Drop these rows explicitly and record the drop, or supply a real id."
        )
    return out

def is_control(donor_id, control_prefix="ctrl", patient_prefixes=("pt", "gvhd", "morphea")):
    """
    Classify a donor. Raises on anything it does not positively recognize, rather than
    letting an unmatched value fall into the disease arm.
    """
    s = _norm(donor_id).lower()
    if s in PLACEHOLDERS:
        raise UnrecordedDonor("cannot classify an unrecorded donor id %r" % donor_id)
    if s.startswith(control_prefix):
        return True
    for p in patient_prefixes:
        if s.startswith(p):
            return False
    raise AmbiguousLabel(
        "donor id %r matches neither the control prefix %r nor any patient prefix %r. "
        "Add the prefix deliberately after checking what it is; do NOT let it default."
        % (donor_id, control_prefix, patient_prefixes)
    )

def group_plates(df, donor_col="patient id", drop_unrecorded=False):
    """
    Return (clean_df, dropped_df). With drop_unrecorded=False an unrecorded id raises.
    With True the rows are REMOVED and returned, so the drop is explicit and countable,
    never silent.
    """
    ids = [_norm(x) for x in df[donor_col]]
    mask = [s.lower() in PLACEHOLDERS for s in ids]
    if any(mask):
        if not drop_unrecorded:
            check_donor_ids(df[donor_col])
        keep = [not m for m in mask]
        return df[keep].copy(), df[mask].copy()
    return df.copy(), df[[False] * len(df)].copy()
