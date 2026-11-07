"""Selftest for donor_group.py. Must pass before GPL24676 is touched."""
import sys, pandas as pd
import donor_group as DG
F=[]
def chk(n,c,d=""):
    print(("  PASS  " if c else "  FAIL  ")+n+(("  "+d) if d else "")); 
    if not c: F.append(n)

print("TEST 1: the exact RESULT_08 case raises instead of fabricating a donor")
df=pd.DataFrame({"patient id":["Ctrl060","pt01050","-","-","Ctrl003"],"v":[1,2,3,4,5]})
try:
    DG.check_donor_ids(df["patient id"]); chk("'-' raises UnrecordedDonor", False)
except DG.UnrecordedDonor as e:
    chk("'-' raises UnrecordedDonor", True, str(e)[:60])

print("\nTEST 2: '-' is never classified into an arm")
try:
    DG.is_control("-"); chk("is_control('-') raises", False)
except DG.UnrecordedDonor:
    chk("is_control('-') raises", True)

print("\nTEST 3: real ids still classify correctly")
chk("Ctrl060 is control", DG.is_control("Ctrl060") is True)
chk("ctrl003 is control", DG.is_control("ctrl003") is True)
chk("pt01050 is patient", DG.is_control("pt01050") is False)
chk("GVHD_pt008 is patient", DG.is_control("GVHD_pt008") is False)
chk("Morphea_pt001 is patient", DG.is_control("Morphea_pt001") is False)

print("\nTEST 4: an UNRECOGNIZED prefix raises, it does not default to disease")
print("  this is the whole point: the old rule made every unmatched value a patient")
try:
    DG.is_control("Subject_17"); chk("unknown prefix raises AmbiguousLabel", False)
except DG.AmbiguousLabel as e:
    chk("unknown prefix raises AmbiguousLabel", True, str(e)[:52])

print("\nTEST 5: other placeholder spellings are caught too")
for p in ["","NA","n/a","None","NaN","unknown","  -  "]:
    try:
        DG.is_control(p); chk("placeholder %r raises"%p, False)
    except DG.UnrecordedDonor:
        chk("placeholder %r raises"%p, True)

print("\nTEST 6: explicit drop is explicit and countable")
clean,drop=DG.group_plates(df,drop_unrecorded=True)
chk("clean has 3 rows", len(clean)==3, "ids=%s"%list(clean["patient id"]))
chk("dropped has 2 rows and is returned", len(drop)==2)
chk("no fabricated donor survives", "-" not in set(clean["patient id"]))

print("\nTEST 7: a clean frame passes through untouched")
ok=pd.DataFrame({"patient id":["Ctrl001","pt02002"],"v":[1,2]})
c2,d2=DG.group_plates(ok)
chk("clean frame unchanged", len(c2)==2 and len(d2)==0)

print("\n"+("ALL SELFTESTS PASS" if not F else "FAILURES: "+", ".join(F)))
sys.exit(1 if F else 0)
