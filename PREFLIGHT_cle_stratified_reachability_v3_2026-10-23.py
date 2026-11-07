# CLE STRATIFIED REACHABILITY, v3. THE RUN IS THE 'batch' CHARACTERISTIC, NOT THE TITLE SUFFIX.
# v1 never resolved donor ids and reported a VOID that was a parser artifact.
# v2 fixed donor ids and dumped the characteristics, WHICH PROVED batch EXISTS FOR ALL 48 SAMPLES
#    (run1 1, run2 16, run3 16, run4 5, run5 10), but then still built its donor table from the
#    TITLE, so the CLE arm printed runs=NONE and I wrongly concluded stratification might be
#    unexecutable. IT IS EXECUTABLE. The prereg's 'CLE 0,0,4,3,2' stands.
# v3 reads the run from the characteristics for EVERY arm and computes reachability.
# METADATA ONLY. NO EXPRESSION IS READ. NOTHING IS SCORED.
import gzip, re, requests, collections, math
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE334nnn/GSE334588/matrix/"
UA = {"User-Agent": "Mozilla/5.0 (research; ORCID 0009-0006-5341-8840)"}
txt = gzip.decompress(requests.get(BASE + "GSE334588_series_matrix.txt.gz", headers=UA, timeout=300).content).decode("utf8","replace")
assert "!series_matrix_table_begin" in txt, "table marker missing; refusing to parse"
head = txt[:txt.index("!series_matrix_table_begin")]
line = lambda k: [l for l in head.split("\n") if l.startswith(k)]
vals = lambda k: [x.strip('"') for x in line(k)[0].split("\t")[1:]] if line(k) else []
tit = vals("!Sample_title"); gsm = vals("!Sample_geo_accession")
ch = [[x.strip('"') for x in l.split("\t")[1:]] for l in head.split("\n") if l.startswith("!Sample_characteristics")]
rec = [dict([("gsm",g),("title",t)] + [(c[i].split(": ")[0].strip().lower(), c[i].split(": ",1)[-1]) for c in ch if ": " in c[i]]) for i,(g,t) in enumerate(zip(gsm,tit))]
did = lambda t: re.sub(r"\s*-\s*batch\s*\d+\s*$", "", t.split(",")[-1], flags=re.I).strip()
arm = lambda d: "Healthy" if re.fullmatch(r"H\d+", d, re.I) else ("CLE" if re.fullmatch(r"CLE\d+", d, re.I) else ("NL_CLE" if re.fullmatch(r"NL_CLE\d+", d, re.I) else ("DM" if re.fullmatch(r"DM\d+", d, re.I) else ("NL_DM" if re.fullmatch(r"NL_DM\d+", d, re.I) else "?"))))
for r in rec: r["donor"]=did(r["title"]); r["arm"]=arm(r["donor"]); r["run"]=(r.get("batch") or "").strip()
assert all(r["arm"]!="?" for r in rec), "a donor id did not classify; STOP"
assert all(r["run"] for r in rec), "a sample has no batch; STOP"
assert all(r["arm"]==r["disease state"] or (r["arm"].startswith("NL") and r["disease state"].startswith("NL")) for r in rec), "title arm and disease state disagree; STOP"
print("ALL 48 SAMPLES CARRY A RUN. title-derived arm agrees with the 'disease state' field for all 48.")
print("run census, all arms:", dict(collections.Counter(r["run"] for r in rec)))
D = collections.defaultdict(list)
for r in rec: D[(r["arm"], r["donor"])].append(r["run"])
print("\n*** DONOR TABLE, Healthy and CLE, RUNS FROM THE CHARACTERISTICS ***")
for (a,d) in sorted(D, key=lambda k:(k[0],k[1])):
    if a in ("Healthy","CLE"): print("   %-7s %-6s runs=%-26s n=%d" % (a, d, sorted(set(D[(a,d)])), len(D[(a,d)])))
cnt = lambda a: collections.Counter(r for (aa,d) in D if aa==a for r in set(D[(aa,d)]))
print("\nDONOR counts per run (a multi-run donor counted once per run it appears in)")
for a in ("CLE","Healthy"): print("   %-7s %r" % (a, dict(sorted(cnt(a).items()))))
print("\n   prereg section 2 asserts CLE 0,0,4,3,2 and Healthy 1,1,3,2,2 across run1..run5")
C = lambda n,k: math.comb(n,k) if 0<=k<=n else 0
HL = {d:sorted(set(v)) for (a,d),v in D.items() if a=="Healthy"}
CL = {d:sorted(set(v)) for (a,d),v in D.items() if a=="CLE"}
def reach(drop, label):
    H = {d:v for d,v in HL.items() if d not in drop}
    asg = {}
    for d,v in H.items():
        c = collections.Counter(x for (a,dd),vv in D.items() if a=="Healthy" and dd==d for x in vv)
        asg[d] = v[0] if len(v)==1 else sorted(c.items(), key=lambda kv:(-kv[1],kv[0]))[0][0]
    for d,v in CL.items(): asg[("CLE",d)] = v[0] if len(v)==1 else sorted(collections.Counter(v).items(), key=lambda kv:(-kv[1],kv[0]))[0][0]
    by = collections.defaultdict(lambda:[0,0])
    for d in H: by[asg[d]][0] += 1
    for d in CL: by[asg[("CLE",d)]][1] += 1
    tot, used, drp = 1, [], []
    for r in sorted(by):
        h,c = by[r]
        (used if (h and c) else drp).append((r,h,c))
        if h and c: tot *= C(h+c, h)
    p = 1.0/tot if tot else float("nan")
    print("\n%s" % label)
    print("   dropped: %r" % (sorted(drop) or "none"))
    print("   strata used (run, nHealthy, nCLE): %r" % used)
    print("   single-arm strata dropped: %r" % (drp or "none"))
    print("   labelings %d | min attainable p %.6f | %s" % (tot, p, "REACHABLE at 0.05" if p<=0.05 else "VOID ON REACHABILITY"))
    return tot
print("\n*** REACHABILITY UNDER THE THREE READINGS ***")
tA = reach({"H1"}, "A  drop H1 only, as she named it")
tB = reach({d for d,v in HL.items() if len(v)>2}, "B  drop donors spanning more than two runs")
noMaj = {d for d,v in HL.items() if len(v)>1 and max(collections.Counter(x for (a,dd),vv in D.items() if a=="Healthy" and dd==d for x in vv).values())*2 <= sum(1 for (a,dd),vv in D.items() if a=="Healthy" and dd==d for _ in vv)}
tC = reach(noMaj, "C  drop donors with no run majority")
print("\n*** DOES THE AMBIGUITY MATTER ***")
print("   A=%d  B=%d  C=%d  ->  %s" % (tA,tB,tC, "THEY AGREE. The H1/H2 reading does not change the test." if tA==tB==tC else "THEY DIFFER. Aditi must choose."))
print("\nNO EXPRESSION VALUE HAS BEEN READ. NOTHING SCORED. NO LOOK SPENT.")
