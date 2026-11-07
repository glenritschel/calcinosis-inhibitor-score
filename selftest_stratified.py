"""Selftest for stratified_null.py. Must pass before the module produces any real number."""
import numpy as np, sys
import stratified_null as S

ARM = np.array(['DM']*5+['CLE']*4+['DM']*2+['CLE']*2, dtype=object)
BAT = np.array(['run3']*9+['run5']*4, dtype=object)
F = []
def chk(name, cond, detail=""):
    print(("  PASS  " if cond else "  FAIL  ")+name+(("  "+detail) if detail else ""))
    if not cond: F.append(name)

print("TEST 1: planted ARM effect (DM lower by 2 SD) is recovered")
rng=np.random.default_rng(1)
s = rng.normal(0,1,13); s[ARM=='DM'] -= 2.0
r = S.stratified_label_null(s, ARM, BAT, 'DM', direction='lower')
chk("exact enumeration used", r['exact'], "n_null=%d"%r['n_null'])
chk("all 756 labelings enumerated", r['n_null']==756)
chk("p is small", r['p_one_sided'] < 0.05, "p=%.5f obs=%.3f"%(r['p_one_sided'],r['observed']))

print("\nTEST 2: PURE BATCH EFFECT, no arm effect. THE TEST THAT MATTERS.")
print("  run5 shifted +5 SD. run5 is 50%% DM, run3 is 56%% DM, so an unstratified")
print("  null sees a spurious arm difference. The stratified null must not.")
rng=np.random.default_rng(2)
ps=[]
for trial in range(200):
    s = rng.normal(0,1,13); s[BAT=='run5'] += 5.0
    r = S.stratified_label_null(s, ARM, BAT, 'DM', direction='lower')
    ps.append(r['p_one_sided'])
ps=np.array(ps)
chk("stratified p is not systematically small under a pure batch effect",
    (ps<0.05).mean() < 0.12, "share p<0.05 = %.3f over 200 trials"%(ps<0.05).mean())
# A design where the confound ALIGNS with the arm. The first formulation of this check
# used the real run3/run5 design and asserted the spurious effect pushed DM higher. It
# does not: run5 is 2 of 6 CLE against 2 of 7 DM, so it is slightly CLE enriched and the
# spurious effect pushes DM LOWER. The check was wrong, not the module. Rather than flip
# the direction to obtain a pass, the design below is built so the confound is unambiguous.
A3 = np.array(['DM']*5+['CLE']*1+['DM']*1+['CLE']*5, dtype=object)
B3 = np.array(['hi']*6+['lo']*6, dtype=object)
B1 = np.array(['one']*12, dtype=object)
rng=np.random.default_rng(3)
nun=[];nst=[]
for _ in range(200):
    s = rng.normal(0,0.3,12); s[B3=='hi'] += 5.0
    nun.append(S.stratified_label_null(s,A3,B1,'DM',direction='higher')['p_one_sided'])
    nst.append(S.stratified_label_null(s,A3,B3,'DM',direction='higher')['p_one_sided'])
nun=np.array(nun); nst=np.array(nst)
print("    stratum 'hi' is 5/6 DM and is shifted +5 SD; no real arm effect exists")
print("    unstratified: share p<0.05 = %.3f   stratified: share p<0.05 = %.3f"%((nun<0.05).mean(),(nst<0.05).mean()))
chk("unstratified null IS fooled by the aligned batch effect", (nun<0.05).mean() > 0.5)
chk("stratified null is NOT fooled", (nst<0.05).mean() < 0.12)

print("\nTEST 3: pure null gives no significance and p is roughly uniform")
rng=np.random.default_rng(4)
ps=[S.stratified_label_null(rng.normal(0,1,13),ARM,BAT,'DM',direction='lower')['p_one_sided'] for _ in range(300)]
ps=np.array(ps)
chk("false positive rate near alpha", 0.0 <= (ps<0.05).mean() <= 0.10, "share p<0.05 = %.3f"%(ps<0.05).mean())
chk("median p near 0.5", 0.35 < np.median(ps) < 0.65, "median=%.3f"%np.median(ps))

print("\nTEST 4: single-arm strata are dropped by construction")
A2 = np.concatenate([ARM, np.array(['DM']*15,dtype=object)])
B2 = np.concatenate([BAT, np.array(['run2']*15,dtype=object)])
s2 = np.concatenate([np.zeros(13), np.full(15, 99.0)])
r = S.stratified_label_null(s2, A2, B2, 'DM', direction='lower')
chk("run2 reported as dropped", 'run2' in r['dropped_strata'], str(r['dropped_strata']))
chk("only 13 samples used", r['n_used']==13)
chk("labeling count unchanged at 756", r['n_labelings']==756)
chk("the 99.0 outlier stratum did not move the statistic", abs(r['observed'])<1e-9,
    "obs=%.3g"%r['observed'])

print("\nTEST 5: guards refuse impossible designs")
try:
    S.stratified_label_null(np.zeros(4), np.array(['DM','DM','CLE','CLE'],dtype=object),
                            np.array(['a','b','c','d'],dtype=object), 'DM')
    chk("all-single-arm strata raise Void", False)
except S.Void as e:
    chk("all-single-arm strata raise Void", True, str(e)[:48])
try:
    S.stratified_label_null(np.zeros(3), np.array(['DM','CLE','X'],dtype=object),
                            np.array(['a','a','a'],dtype=object), 'DM')
    chk("three arms raise Void", False)
except S.Void as e:
    chk("three arms raise Void", True, str(e)[:48])
try:
    S.stratified_label_null(np.zeros(3), np.array(['DM','CLE','DM'],dtype=object),
                            np.array(['a','a'],dtype=object), 'DM')
    chk("length mismatch raises Void", False)
except S.Void as e:
    chk("length mismatch raises Void", True, str(e)[:48])

print("\n"+("ALL SELFTESTS PASS" if not F else "FAILURES: "+", ".join(F)))
sys.exit(1 if F else 0)
