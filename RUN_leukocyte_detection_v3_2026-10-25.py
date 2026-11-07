"""Detection only. PREREG_LEUKOCYTE_DETECTION_v1_2026-10-23 sha256 786c3444165be313.
NO HYPOTHESIS, NO EFFECT, NO p, NO LOOK SPENT. Controls print BEFORE test genes.

v1 VOIDED CORRECTLY on a parsing fault: Monaco row ids are Ensembl, not symbols, so every
gene read NOT IN MATRIX. The positive control turned that into VOID instead of a false
'0 of 8 detected'. v2 matched symbol OR Ensembl and then FAILED on the HPA arm.

v3, 2026-10-25, fixes two defects found while repairing the v2 HPA failure. The second is
much worse than the crash that exposed it.

DEFECT A, THE REPORTED ONE. v2 hunted for immune nTPM columns inside the WIDE
proteinatlas.tsv with a regex, and matched 'RNA blood cell specific nTPM' and 'RNA blood
lineage specific nTPM'. ⛔ THOSE ARE SPECIFICITY ANNOTATION COLUMNS holding strings like
'Group enriched', not per-cell-type expression. v3 does not look for columns at all: HPA
publishes rna_immune_cell.tsv.zip in LONG format, one row per gene per immune cell, with an
explicit nTPM field. There is no column-finding step left to get wrong.

DEFECT B, WHICH THE CRASH CONCEALED. With every value a non-numeric string, v2's
mat[sym].get(c, 0.0) on a defaultdict CREATED a key for all fourteen genes, each an EMPTY
dict. The positive control then read
    sum(1 for x in mat['PTPRC'].values() if x >= FLOOR) >= 0.8 * len(mat['PTPRC'])
which on zero data is 0 >= 0.0, i.e. ⛔ PASS. The purity control, written with a median
instead, raised IndexError on the empty list two lines later and stopped the run. HAD IT
BEEN WRITTEN LIKE THE POSITIVE CONTROL, v2 WOULD HAVE PRINTED PASS / PASS ON AN EMPTY
MATRIX AND THEN REPORTED '0 of 8 DETECTED' -- exactly the false negative the positive
control exists to prevent. The crash was luck, not a safeguard. A void is not a null, and
an empty matrix is not a measurement of absence.

v3 therefore: (1) builds no phantom keys; (2) gates on DATA SUFFICIENCY before reading any
control, so a control with fewer than MIN_OBS observed values is VOID and never PASS;
(3) strips quotes before float; (4) states the observed count on every line so thin data is
visible rather than inferred."""
import io, zipfile, gzip, json, urllib.request, urllib.parse, collections, re

TEST = ["ENPP1","ANKH","MGP","GGCX","VKORC1","NT5E","ABCC6","TNFRSF11B"]
POS  = ["PTPRC","ACTB","B2M"]
PUR  = ["COL1A1","COL1A2","THY1"]
ALL  = TEST + POS + PUR
FLOOR, MIN_SUBSETS = 1.0, 3
MIN_OBS = 5          # a control carrying fewer observed values than this cannot be read

def get(url, tries=3, to=240):
    last=None
    for i in range(tries):
        try: return urllib.request.urlopen(url, timeout=to).read()
        except Exception as e:
            last=e; print("   attempt %d: %s" % (i+1, e))
    raise last

def ensembl_map():
    """symbol -> ENSG, from mygene.info. Printed so it is checkable."""
    q = urllib.parse.urlencode({"q":",".join(ALL),"scopes":"symbol","fields":"ensembl.gene",
                                "species":"human"}).encode()
    r = json.loads(urllib.request.urlopen(
        urllib.request.Request("https://mygene.info/v3/query", data=q), timeout=120).read())
    m={}
    for h in r:
        g=h.get("query"); e=h.get("ensembl")
        if not g or not e: continue
        ids=[e["gene"]] if isinstance(e,dict) else [x["gene"] for x in e]
        m[g]=ids
    print("   symbol -> Ensembl: " + ", ".join("%s=%s"%(k,"/".join(v)) for k,v in sorted(m.items())))
    missing=[g for g in ALL if g not in m]
    if missing: print("   ⛔ no Ensembl id resolved for: %s" % ", ".join(missing))
    return m

def _num(x):
    """float or None. Strips quotes and whitespace. Never raises, never returns 0.0 for junk."""
    try: return float(str(x).strip().strip('"').strip("'"))
    except (ValueError, TypeError): return None

def build(rows, cols, idmap):
    """rows: (rowid, [values]). Match rowid as symbol, or as ENSG (version stripped).
    ⛔ Creates NO key for a gene with no numeric value: a phantom empty entry is what let
    v2's positive control pass on zero data."""
    e2s = {e.split(".")[0]: s for s, es in idmap.items() for e in es}
    mat = {}; hits = collections.Counter(); junk = 0
    for rid, vals in rows:
        r = str(rid).strip().strip('"')
        sym = None
        if r in ALL: sym, how = r, "symbol"
        else:
            base = r.split(".")[0]
            if base in e2s: sym, how = e2s[base], "ensembl"
        if sym is None: continue
        hits[how] += 1
        for c, v in zip(cols, vals):
            f = _num(v)
            if f is None:
                junk += 1; continue
            d = mat.setdefault(sym, {})
            d[c] = max(d.get(c, 0.0), f)
    print("   matched rows: %s" % (dict(hits) or "NONE"))
    print("   genes carrying at least one numeric value: %d of %d" % (len(mat), len(ALL)))
    if junk: print("   ⛔ %d non-numeric cells skipped (not counted as zero)" % junk)
    return mat

def build_long(triples, idmap):
    """(gene_id, subset, value) -> matrix. For HPA's long-format immune cell file."""
    e2s = {e.split(".")[0]: s for s, es in idmap.items() for e in es}
    mat = {}; hits = collections.Counter(); subs = set(); junk = 0
    for gid, sub, val in triples:
        g = str(gid).strip().strip('"'); sub = str(sub).strip().strip('"')
        sym = None
        if g in ALL: sym, how = g, "symbol"
        else:
            base = g.split(".")[0]
            if base in e2s: sym, how = e2s[base], "ensembl"
        if sym is None: continue
        f = _num(val)
        if f is None:
            junk += 1; continue
        hits[how] += 1; subs.add(sub)
        d = mat.setdefault(sym, {})
        d[sub] = max(d.get(sub, 0.0), f)
    print("   matched rows: %s" % (dict(hits) or "NONE"))
    print("   distinct subsets seen: %d" % len(subs))
    print("   genes carrying at least one numeric value: %d of %d" % (len(mat), len(ALL)))
    if junk: print("   ⛔ %d non-numeric values skipped (not counted as zero)" % junk)
    return mat, sorted(subs)

def report(name, mat, subsets, unit):
    print("\n" + "="*76)
    print("SOURCE: %s   |   %d %s" % (name, len(subsets), unit))
    print("="*76)
    def vals(g): return sorted(mat.get(g, {}).values())
    def line(g):
        v = vals(g)
        if not v: return "%-11s  NOT IN MATRIX  (0 observed values)" % g
        n = sum(1 for x in v if x >= FLOOR)
        top = sorted(mat[g].items(), key=lambda kv: -kv[1])[:3]
        return "%-11s  median %8.2f  max %9.2f  >=%.1f in %3d/%3d (%4.0f%%)  top: %s" % (
            g, v[len(v)//2], v[-1], FLOOR, n, len(v), 100.0*n/len(v),
            ", ".join("%s %.0f" % (k[:20], x) for k, x in top))

    # DATA SUFFICIENCY GATE, v3. Read BEFORE any control, because a control computed on an
    # empty list is not a failed control, it is an unread one, and v2 read 0 >= 0 as PASS.
    thin = [(g, len(mat.get(g, {}))) for g in POS + PUR if len(mat.get(g, {})) < MIN_OBS]
    if thin:
        print("\n⛔ DATA SUFFICIENCY: %d of %d control genes carry fewer than %d observed "
              "values:" % (len(thin), len(POS+PUR), MIN_OBS))
        for g, n in thin: print("     %-11s %d observed" % (g, n))
        print("\n  *** VOID. The controls cannot be read, so they are NOT reported as "
              "passing or failing. Nothing is reported about the test genes. ***")
        print("  (This is the v2 defect: on zero data its positive control returned "
              "0 >= 0.8*0, i.e. PASS.)")
        return None

    print("\nPOSITIVE CONTROL  (must be HIGH or the data are not leukocytes)")
    for g in POS: print("  " + line(g))
    print("\nPURITY CONTROL    (must be LOW or the sort carries stroma and this is VOID)")
    for g in PUR: print("  " + line(g))

    pv = vals("PTPRC")
    pos_ok = len(pv) >= MIN_OBS and sum(1 for x in pv if x >= FLOOR) >= 0.8*len(pv)
    cols = [vals(g) for g in ("COL1A1","COL1A2") if len(mat.get(g, {})) >= MIN_OBS]
    pur_ok = bool(cols) and all(v[len(v)//2] < FLOOR for v in cols)
    print("\n  POSITIVE CONTROL: %s      PURITY CONTROL: %s" % (
        "PASS" if pos_ok else "FAIL", "PASS" if pur_ok else "FAIL"))
    if not (pos_ok and pur_ok):
        print("\n  *** VOID. Controls did not pass. Nothing reported about the test genes. ***")
        return None
    print("\nTEST SET")
    for g in TEST: print("  " + line(g))
    det = [g for g in TEST
           if sum(1 for x in vals(g) if x >= FLOOR) >= MIN_SUBSETS]
    print("\n  DETECTED (>=%.1f nTPM in >=%d %s): %d of %d -> %s" % (
        FLOOR, MIN_SUBSETS, unit, len(det), len(TEST), ", ".join(det) or "NONE"))
    print("  ⚠ the >=%d rule is applied to %s. Read the PERCENTAGE column as the comparable figure." % (MIN_SUBSETS, unit))
    print("  PRE-REGISTERED PREDICTION was NT5E and ENPP1 clear.   NT5E %s   ENPP1 %s" % (
        "CLEARS" if "NT5E" in det else "does NOT clear", "CLEARS" if "ENPP1" in det else "does NOT clear"))
    return set(det)

def monaco(idmap):
    base="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107011/suppl/"
    fn=re.findall(r'href="([^"]*[Tt][Pp][Mm][^"]*)"', get(base).decode("utf-8","replace"))[0].split("/")[-1]
    print("   file: %s" % fn)
    raw=get(base+fn); txt=(gzip.decompress(raw) if fn.endswith(".gz") else raw).decode("utf-8","replace").splitlines()
    sep="\t" if "\t" in txt[0] else ","
    hdr=[h.strip().strip('"') for h in txt[0].split(sep)]
    cols=hdr[1:]
    print("   first row ids: %s" % [l.split(sep)[0].strip('"') for l in txt[1:6]])
    print("   first columns: %s" % cols[:6])
    rows=[(l.split(sep)[0], l.split(sep)[1:]) for l in txt[1:] if len(l.split(sep))==len(hdr)]
    return build(rows, cols, idmap), cols, "sample columns"

def hpa(idmap):
    """HPA immune cell nTPM, LONG format: one row per gene per immune cell.
    Columns per proteinatlas.org/humanproteome/single+cell/immune+cell/data:
    Gene (Ensembl id), Immune cell, TPM, pTPM, nTPM. ⛔ NOT the wide proteinatlas.tsv,
    whose 'RNA blood cell specific nTPM' columns are specificity ANNOTATIONS, not values."""
    url="https://www.proteinatlas.org/download/tsv/rna_immune_cell.tsv.zip"
    z=zipfile.ZipFile(io.BytesIO(get(url)))
    nm=[n for n in z.namelist() if n.lower().endswith(".tsv")]
    if not nm: raise RuntimeError("no .tsv inside %s (members: %s)" % (url, z.namelist()))
    txt=z.read(nm[0]).decode("utf-8","replace").splitlines()
    hdr=[h.strip().strip('"') for h in txt[0].split("\t")]
    print("   file: %s   header: %s" % (nm[0], hdr))
    def need(*names):
        for n in names:
            for i,h in enumerate(hdr):
                if h.lower()==n.lower(): return i
        raise RuntimeError("no %s column in %s; header is %s" % (names[0], nm[0], hdr))
    gi, si, vi = need("Gene"), need("Immune cell","Blood cell","Cell type"), need("nTPM")
    print("   using columns: Gene=%r subset=%r value=%r" % (hdr[gi], hdr[si], hdr[vi]))
    tri=[]
    for l in txt[1:]:
        f=l.split("\t")
        if len(f)<=max(gi,si,vi): continue
        tri.append((f[gi], f[si], f[vi]))
    mat, subs = build_long(tri, idmap)
    return mat, subs, "immune cell types"

if __name__ == "__main__":
    print("PREREG_LEUKOCYTE_DETECTION_v1_2026-10-23   floor nTPM >= %.1f in >= %d subsets" % (FLOOR, MIN_SUBSETS))
    print("DETECTION ONLY. No hypothesis, no effect, no p, no look spent.\n")
    print("Resolving gene identifiers ...")
    try: IDM = ensembl_map()
    except Exception as e:
        print("   ⛔ mygene.info failed (%s); symbol matching only" % e); IDM = {}
    res={}
    for nm_, fn in (("Monaco GSE107011", monaco), ("Human Protein Atlas", hpa)):
        print("\nFetching %s ..." % nm_)
        try:
            mat, cols, unit = fn(IDM)
            res[nm_] = report(nm_, mat, cols, unit)
        except Exception as e:
            print("   FAILED: %s: %s" % (type(e).__name__, e)); res[nm_]="FAILED"
    print("\n" + "="*76)
    ok={k:v for k,v in res.items() if isinstance(v,set)}
    if len(ok)==2:
        a,b=list(ok.values())
        print("CROSS-SOURCE: both clear %s" % (", ".join(sorted(a&b)) or "NOTHING"))
        print(("⛔ DISAGREE on %s -- per section 5 that disagreement IS the result." % ", ".join(sorted(a^b))) if a^b else "No disagreement between sources.")
    elif len(ok)==1: print("ONE SOURCE ONLY. Section 5 wanted two; a single-source verdict is weaker and is reported as such.")
    else: print("NO SOURCE USABLE. Nothing is concluded.")
    print("="*76)
