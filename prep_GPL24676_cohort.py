import gzip, requests
import donor_group as DG, replication_harness as H
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
ex = lambda R: [r for r in R if r.get("tissue") == "Skin" and r.get("selection marker") == "CD90+"]
E24 = ex(REC24); E18 = ex(REC18)
nm = lambda r: DG._norm(r.get("patient id",""))
ids24 = {nm(r) for r in E24}; ids18 = {nm(r) for r in E18}
print("CD90+ skin, exact match: GPL24676 %d plates %d ids | GPL18573 %d plates %d ids" % (len(E24), len(ids24), len(E18), len(ids18)))
assert (len(E18), len(ids18)) == (184, 101), "GPL18573 is %r, not the 184 plates and 101 ids RESULT_03 recorded" % ((len(E18), len(ids18)),)
assert (len(E24), len(ids24)) == (159, 71), "GPL24676 is %r, not the 159 plates and 71 ids the prereg recorded" % ((len(E24), len(ids24)),)
print("  both match the figures recorded before today")
keep = lambda s: s.lower() not in H.NONPERSON and not s.lower().startswith(H.NONSSC)
a = {d for d in ids24 if keep(d)}; b = {d for d in ids18 if keep(d)}
PRIMARY = H.option_b_primary(ids24, ids18)
isc = lambda d: d.lower().startswith("ctrl")
print("OPTION B: GPL18573 %d | GPL24676 %d | overlap %d | PRIMARY %d" % (len(b), len(a), len(a & b), len(PRIMARY)))
print("v5 addendum states: 97 | 56 | 6 | 50")
assert (len(b), len(a), len(a & b), len(PRIMARY)) == (97, 56, 6, 50), "the cohort does not reconstruct"
print("control %d, patient %d" % (sum(1 for d in PRIMARY if isc(d)), sum(1 for d in PRIMARY if not isc(d))))
H.check_primary(PRIMARY, isc)
print("\ncheck_primary PASSED. THE SET IN HAND IS THE PRE-REGISTERED SET.")
print("plates belonging to the primary 50: %d" % sum(1 for r in E24 if nm(r) in set(PRIMARY)))
print("\nNO EXPRESSION VALUE HAS BEEN READ. THE ONE LOOK IS NOT TAKEN.")
