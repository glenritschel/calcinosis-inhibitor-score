"""
Batch-stratified label null for PREREG_INHIBITION_v1 under DEVIATION_04.

The unstratified label_null in inhibition_score.py permutes disease labels across all
donors. In GSE334588 the arms are not balanced across sequencing runs, so that null
attributes run-to-run variation to disease. This module permutes labels WITHIN stratum.

Strata containing only one arm carry no information about the contrast and are DROPPED
BY CONSTRUCTION, not by choice: a stratum with one arm admits exactly one labeling.

Statistic: scores are centered within stratum, then the difference of arm means is taken
over the retained samples. Centering is what removes the stratum offset; permuting within
stratum is what makes the null respect it. Both are required, and either alone is wrong.
"""
import numpy as np
from itertools import product, combinations

class Void(Exception):
    pass

def _prepare(scores, arm, strata, target):
    s = np.asarray(scores, dtype=float)
    a = np.asarray(arm, dtype=object)
    b = np.asarray(strata, dtype=object)
    if not (len(s) == len(a) == len(b)):
        raise Void("scores, arm and strata must be the same length")
    arms = sorted(set(a.tolist()))
    if len(arms) != 2:
        raise Void("exactly two arms required, got %r" % (arms,))
    if target not in arms:
        raise Void("target %r is not one of the arms %r" % (target, arms))
    keep = np.zeros(len(s), dtype=bool)
    dropped = {}
    for st in sorted(set(b.tolist())):
        m = (b == st)
        n_t = int((a[m] == target).sum())
        n_o = int(m.sum()) - n_t
        if n_t > 0 and n_o > 0:
            keep |= m
        else:
            dropped[st] = (n_t, n_o)
    if not keep.any():
        raise Void("no stratum contains both arms; the contrast is not estimable")
    return s[keep], a[keep], b[keep], dropped, arms

def _centered_stat(s, a, b, target):
    out = s.astype(float).copy()
    for st in set(b.tolist()):
        m = (b == st)
        out[m] = out[m] - out[m].mean()
    t = (a == target)
    return float(out[t].mean() - out[~t].mean())

def n_labelings(a, b, target):
    """Exact count of distinct within-stratum relabelings."""
    tot = 1
    for st in sorted(set(b.tolist())):
        m = (b == st)
        n = int(m.sum()); k = int((a[m] == target).sum())
        tot *= _comb(n, k)
    return tot

def _comb(n, k):
    from math import comb
    return comb(n, k)

def guard_reachable_stratified(a, b, target, alpha=0.05):
    """Can the stratified null attain p < alpha at all? Returns (n_labelings, min_p)."""
    n = n_labelings(a, b, target)
    return n, 1.0 / n

def stratified_label_null(scores, arm, strata, target, direction="lower",
                          n_perm=20000, seed=0, exact_limit=200000):
    """
    direction: 'lower' means the one-sided alternative is that `target` scores LOWER.
    Enumerates exactly when the number of distinct relabelings is small enough,
    otherwise samples. Which path was taken is reported in the result.
    """
    if direction not in ("lower", "higher"):
        raise Void("direction must be 'lower' or 'higher'")
    s, a, b, dropped, arms = _prepare(scores, arm, strata, target)
    obs = _centered_stat(s, a, b, target)
    n_lab = n_labelings(a, b, target)
    if n_lab < 2:
        raise Void("only one distinct relabeling exists; no null is possible")
    strata_ids = sorted(set(b.tolist()))
    idx = {st: np.flatnonzero(b == st) for st in strata_ids}
    ktab = {st: int((a[idx[st]] == target).sum()) for st in strata_ids}
    other = [x for x in arms if x != target][0]

    def stat_from(assign):
        lab = np.array([other] * len(s), dtype=object)
        for st, pos in assign:
            lab[idx[st][list(pos)]] = target
        return _centered_stat(s, lab, b, target)

    exact = n_lab <= exact_limit
    null = []
    if exact:
        per = [[(st, c) for c in combinations(range(len(idx[st])), ktab[st])] for st in strata_ids]
        for combo in product(*per):
            null.append(stat_from(combo))
    else:
        rng = np.random.default_rng(seed)
        for _ in range(n_perm):
            combo = [(st, rng.choice(len(idx[st]), ktab[st], replace=False)) for st in strata_ids]
            null.append(stat_from(combo))
    null = np.asarray(null, dtype=float)
    if direction == "lower":
        p1 = float((null <= obs).sum() + (0 if exact else 1)) / (len(null) + (0 if exact else 1))
    else:
        p1 = float((null >= obs).sum() + (0 if exact else 1)) / (len(null) + (0 if exact else 1))
    p2 = float((np.abs(null) >= abs(obs)).sum() + (0 if exact else 1)) / (len(null) + (0 if exact else 1))
    return dict(observed=obs, p_one_sided=p1, p_two_sided=p2, n_null=len(null),
                exact=exact, n_labelings=n_lab, dropped_strata=dropped,
                n_used=len(s), direction=direction, target=target, other=other)
