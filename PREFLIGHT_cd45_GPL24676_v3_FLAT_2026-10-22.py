# PREFLIGHT v3: can the Deviation-2 specificity check be run on GPL24676 at all?
# SELF-CONTAINED AND ZERO-INDENT. v1 died on ModuleNotFoundError in a fresh Colab kernel.
# v2 fixed that but reintroduced indented def blocks, which is the Colab auto-indent defect that
# prep_GPL24676_cohort.py was deliberately written flat to avoid. v3 is flat.
# The two borrowed functions are TRANSLITERATED to lambdas, not rewritten from memory. Originals:
#   donor_group._norm:            def _norm(x): return re.sub(r"\s+", " ", str(x)).strip()
#   replication_harness:          NONPERSON = {"-","ctrl","ctrl_pool","hip","dssc"}
#                                 NONSSC    = ("gvhd","morphea")
#   replication_harness.option_b_primary:
#                                 keep = lambda s: s.lower() not in NONPERSON and not s.lower().startswith(NONSSC)
#                                 a = {d for d in ids_24676 if keep(d)}
#                                 b = {d for d in ids_18573 if keep(d)}
#                                 return sorted(a - b)
# THE ASSERTS BELOW ARE THE CHECK ON THAT TRANSLITERATION: 184/101, 159/71 and 50 at 25/25 all
# come from the pre-registration, so if a lambda diverges the cell stops instead of reporting.
# METADATA ONLY. NO EXPRESSION VALUE IS READ. THE ONE LOOK IS NOT TOUCHED.
import gzip, re, requests, collections
_norm = lambda x: re.sub(r"\s+", " ", str(x)).strip()
NONPERSON = {"-","ctrl","ctrl_pool","hip","dssc"}
NONSSC = ("gvhd","morphea")
keep = lambda s: s.lower() not in NONPERSON and not s.lower().startswith(NONSSC)
option_b_primary = lambda i24, i18: sorted({d for d in i24 if keep(d)} - {d for d in i18 if keep(d)})
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
nm = lambda r: _norm(r.get("patient id",""))
cen = lambda R: collections.Counter(r.get("selection marker") for r in R if r.get("tissue") == "Skin")
print("GPL24676 skin marker census:", dict(cen(REC24)))
print("GPL18573 skin marker census:", dict(cen(REC18)))
print("  exact-match discipline: 'CD45+' only, never a substring of 'CD45+CD90+' or 'CD90-CD45-'")
ex = lambda R,mk: [r for r in R if r.get("tissue") == "Skin" and r.get("selection marker") == mk]
E24c = ex(REC24,"CD90+"); E18c = ex(REC18,"CD90+")
ids24 = {nm(r) for r in E24c}; ids18 = {nm(r) for r in E18c}
assert (len(E18c), len(ids18)) == (184, 101), "GPL18573 CD90+ is %r not (184,101); stop" % ((len(E18c), len(ids18)),)
assert (len(E24c), len(ids24)) == (159, 71), "GPL24676 CD90+ is %r not (159,71); stop" % ((len(E24c), len(ids24)),)
PRIMARY = set(option_b_primary(ids24, ids18))
isc = lambda d: d.lower().startswith("ctrl")
assert (len(PRIMARY), sum(1 for d in PRIMARY if isc(d))) == (50, 25), "PRIMARY is %d donors %d control, not 50 at 25; the transliterated option_b_primary does NOT reproduce the pre-registered set, stop" % (len(PRIMARY), sum(1 for d in PRIMARY if isc(d)))
print("\nCD90+ reconstructs 184/101 and 159/71, and PRIMARY rebuilds as 50 at 25 and 25.")
print("THE TRANSLITERATION REPRODUCES THE PRE-REGISTERED SET. Proceeding to CD45+.")
I24 = ex(REC24,"CD45+"); I18 = ex(REC18,"CD45+")
iid24 = {nm(r) for r in I24}
print("\nCD45+ skin plates: GPL24676 %d plates, %d donors | GPL18573 %d plates (RESULT_04 used 86)" % (len(I24), len(iid24), len(I18)))
ov = iid24 & PRIMARY
print("\n*** THE QUESTION THIS CELL EXISTS TO ANSWER ***")
print("CD45+ donors on GPL24676 ALSO in the held-out primary 50: %d of 50" % len(ov))
print("  of those, control %d, patient %d" % (sum(1 for d in ov if isc(d)), sum(1 for d in ov if not isc(d))))
print("  CD45+ plates belonging to those donors: %d" % sum(1 for r in I24 if nm(r) in ov))
print("\nREAD BEFORE RUNNING ANYTHING ELSE:")
print("  Deviation 2 asks whether the inhibitor deficit ALSO appears in immune cells.")
print("  That needs BOTH arms. If either is small the check is VOID ON COHORT GROUNDS and no")
print("  gene needs counting to know it. RESULT_04 voided the GPL18573 version on DETECTION,")
print("  after reading expression. This voids or clears it one step earlier.")
print("\nNO EXPRESSION VALUE HAS BEEN READ. NOTHING SCORED. THE ONE LOOK IS NOT TOUCHED.")
