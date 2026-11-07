"""⛔ THIS IS A DIAGNOSTIC, NOT A RE-RUN OF THE PRE-REGISTERED TEST. NO LOOK IS SPENT.

THE QUESTION. RESULT_20 found the same Monaco samples FAIL the purity control in the GEO
per-sample rendering (COL1A1 median 1.89 in 96/127) and PASS in the HPA donor-aggregated
rendering (median 0.00, 0/30). Processing and aggregation changed together, so the cause is
unidentified. If AGGREGATION is the cause, the aggregated passes are weak evidence and the
per-sample fail is the more sensitive test. If PROCESSING is the cause, the GEO author TPM is
the outlier.

THE DISCRIMINATING DATA. HPA publishes rna_immune_cell_sample.tsv.zip, ITS OWN pipeline at
SAMPLE level. Same processing as the aggregated files, different aggregation. So:

    HPA sample-level FAILS purity  -> the flip is AGGREGATION
    HPA sample-level PASSES purity -> the flip is PROCESSING

v6, 2026-10-25, after v5 REFUSED the sample-level file and was RIGHT TO. v5 printed:
    header: ['Sample ID','Donor','Immune cell','ENSG ID','Gene name','TPM','pTPM','nTPM']
    ⛔ REFUSED: need Gene, a subset column and pTPM.
⭐ THE REFUSAL CAUGHT TWO DOCUMENTATION ERRORS AT ONCE. The HPA download page describes this
file as "Gene, Immune cell, Donor, TPM, pTPM". Both halves are wrong: the gene column is called
ENSG ID, not Gene, AND THE FILE DOES CARRY nTPM. A script that guessed a column would have
matched nothing and reported a clean per-sample purity PASS from an EMPTY MATRIX, which is
precisely the v2 defect (decision 1189). Refusing by name and printing the header turned a
silent wrong answer into a one-line fix.

⛔ UNITS. v5 planned to compare pTPM because the documentation said the sample file had no
nTPM. It does. So v6 compares BOTH, with nTPM as the primary because that is the unit
RESULT_20's figures are in and the unit the pre-registration fixes its floor in. ⛔ THIS IS
STILL NOT THE PRE-REGISTERED TEST: purity genes only, no test gene is evaluated, and no
detection verdict is produced. It answers one question: does HPA's own per-sample view carry
stroma.
"""
import io, zipfile, urllib.request, collections, statistics

PUR = ["COL1A1","COL1A2","THY1"]
POS = ["PTPRC"]
WATCH = PUR + POS
IDS = {"COL1A1":"ENSG00000108821","COL1A2":"ENSG00000164692",
       "THY1":"ENSG00000154096","PTPRC":"ENSG00000081237"}
E2S = {v:k for k,v in IDS.items()}
BASE = "https://www.proteinatlas.org/download/tsv/"
FILES = [("HPA consensus, aggregated",        "rna_immune_cell.tsv.zip"),
         ("HPA/Monaco, aggregated",           "rna_immune_cell_monaco.tsv.zip"),
         ("HPA SAMPLE LEVEL, not aggregated", "rna_immune_cell_sample.tsv.zip")]

def get(url, tries=3, to=240):
    last=None
    for i in range(tries):
        try: return urllib.request.urlopen(url, timeout=to).read()
        except Exception as e:
            last=e; print("   attempt %d: %s" % (i+1, e))
    raise last

def num(x):
    try: return float(str(x).strip().strip('"').strip("'"))
    except (ValueError, TypeError): return None

def read(label, fn):
    print("\n" + "="*76); print("SOURCE: %s" % label); print("="*76)
    z=zipfile.ZipFile(io.BytesIO(get(BASE+fn)))
    nm=[n for n in z.namelist() if n.lower().endswith(".tsv")][0]
    txt=z.read(nm).decode("utf-8","replace").splitlines()
    hdr=[h.strip().strip('"') for h in txt[0].split("\t")]
    print("   file: %s" % nm); print("   header: %s" % hdr)
    def find(*names):
        for n in names:
            for i,h in enumerate(hdr):
                if h.lower()==n.lower(): return i
        return None
    gi=find("Gene","ENSG ID","Ensembl","Gene ID"); ni=find("Gene name")
    si=find("Immune cell","Blood cell","Cell type")
    di=find("Donor","Sample")
    vi=find("nTPM")                      # present in all three, despite the docs
    pi=find("pTPM")
    if gi is None or si is None or vi is None:
        print("   ⛔ REFUSED: need a gene column, a subset column and nTPM. Header above.")
        return None
    print("   using nTPM at column %d%s" % (vi, "" if di is None else
          ", donor/sample column %d present -> this rendering is PER SAMPLE" % di))
    vals=collections.defaultdict(list); pvals=collections.defaultdict(list); names={}
    for l in txt[1:]:
        f=l.split("\t")
        if len(f)<=max(x for x in (gi,si,vi,di or 0,ni or 0)): continue
        base=f[gi].strip().strip('"').split(".")[0]
        sym=E2S.get(base)
        if sym is None: continue
        if ni is not None: names.setdefault(sym, f[ni].strip().strip('"'))
        v=num(f[vi])
        if v is not None: vals[sym].append(v)
        if pi is not None and len(f)>pi:
            pv=num(f[pi])
            if pv is not None: pvals[sym].append(pv)
    bad=[(s,g) for s,g in names.items() if g.upper()!=s.upper()]
    print("   IDENTITY CONTROL: %d/%d agree%s" % (len(names)-len(bad), len(names),
          "" if not bad else "  ⛔⛔ DISAGREEMENT: %s" % bad))
    if not vals: print("   ⛔ no watched gene found"); return None
    n_obs = len(vals.get("COL1A1", []))
    print("   observations per gene: %d  (%s)" % (n_obs,
          "aggregated cell types" if di is None else "individual samples"))
    print()
    out={}
    for g in WATCH:
        v=sorted(vals.get(g,[]))
        if not v: print("   %-9s NOT FOUND" % g); continue
        over=sum(1 for x in v if x>=1.0)
        out[g]=(statistics.median(v), v[-1], over, len(v))
        pv=sorted(pvals.get(g,[]))
        extra=("   |  pTPM median %8.2f  >=1.0 in %4d/%4d" %
               (statistics.median(pv), sum(1 for x in pv if x>=1.0), len(pv))) if pv else ""
        print("   %-9s nTPM median %9.2f  max %9.2f  >=1.0 in %4d/%4d (%4.0f%%)%s" % (
              g, statistics.median(v), v[-1], over, len(v), 100.0*over/len(v), extra))
    return out

res={}
for label, fn in FILES:
    print("\nFetching %s ..." % fn)
    try: res[label]=read(label, fn)
    except Exception as e: print("   FAILED: %s: %s" % (type(e).__name__, e)); res[label]=None

print("\n" + "="*76); print("VERDICT"); print("="*76)
samp=[k for k in res if "SAMPLE LEVEL" in k and res[k]]
agg =[k for k in res if "SAMPLE LEVEL" not in k and res[k]]
if not samp:
    print("⛔ The sample-level rendering was not usable. The question is UNANSWERED and the")
    print("   RESULT_20 caveat stands in full: the flip may be aggregation or processing.")
else:
    k=samp[0]; med,mx,over,n = res[k]["COL1A1"]
    print("HPA per-sample COL1A1: nTPM median %.2f, >=1.0 in %d/%d (%.0f%%)" % (med,over,n,100.0*over/n))
    for a in agg:
        m2,x2,o2,n2 = res[a]["COL1A1"]
        print("   vs %-34s median %6.2f, %d/%d" % (a, m2, o2, n2))
    if over > 0.5*n:
        print("\n⭐ THE FLIP IS AGGREGATION. HPA's OWN per-sample view carries COL1A1 where its")
        print("   aggregated view does not. ⛔ Then every purity PASS computed on aggregated data")
        print("   in this program is weak evidence, RESULT_20 included, and the per-sample GEO")
        print("   FAIL in RESULT_18 is the more sensitive test rather than the outlier.")
    elif over == 0:
        print("\n⭐ THE FLIP IS PROCESSING. HPA's per-sample view is clean, so aggregation is not")
        print("   hiding stroma and the GEO author TPM rendering is the outlier. The aggregated")
        print("   purity passes stand.")
    else:
        print("\n⛔ INTERMEDIATE: %d/%d of samples carry COL1A1 at >=1.0 nTPM. Neither clean" % (over,n))
        print("   explanation holds. Report the figure, draw no verdict, and say which samples.")
print("\n⛔ NO LOOK SPENT. This touched purity genes only and did not apply the")
print("   pre-registered floor as a detection verdict or evaluate any test gene.")
print("="*76)
