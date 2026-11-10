import itertools
import numpy as np
import pandas as pd

log = pd.read_csv(
    "GSE292979_n14_donor_matrix_TRUEKEY_log1p.csv",
    index_col=0
)

primary_genes = [
    "ENPP1", "ANKH", "MGP", "GGCX", "VKORC1", "NT5E"
]

secondary_genes = [
    "ENPP1", "ANKH", "MGP", "TNFRSF11B",
    "GGCX", "VKORC1", "NT5E"
]

mapping = {
    "D01": ("GSM8872107", "HC", 2),
    "D02": ("GSM8872104", "SSc", 2),
    "D03": ("GSM8872097", "HC", 1),
    "D04": ("GSM8872100", "SSc", 2),
    "D05": ("GSM8872096", "SSc", 1),
    "D06": ("GSM8872108", "HC", 2),
    "D07": ("GSM8872103", "SSc", 2),
    "D08": ("GSM8872102", "SSc", 2),
    "D09": ("GSM8872095", "SSc", 1),
    "D10": ("GSM8872101", "SSc", 2),
    "D11": ("GSM8872106", "SSc", 2),
    "D12": ("GSM8872094", "SSc", 1),
    "D13": ("GSM8872099", "SSc", 2),
    "D14": ("GSM8872098", "SSc", 2),
}

def make_scores(genes):
    x = log.loc[genes].astype(float)
    z = x.sub(x.mean(axis=1), axis=0).div(
        x.std(axis=1, ddof=1), axis=0
    )
    return z.mean(axis=0)

def exact_test(genes):
    scores = make_scores(genes)

    rows = []
    for code, (gsm, arm, batch) in mapping.items():
        rows.append({
            "code": code,
            "gsm": gsm,
            "arm": arm,
            "batch": batch,
            "score": scores[gsm],
        })

    df = pd.DataFrame(rows)

    observed = (
        df.loc[df["arm"] == "SSc", "score"].mean()
        - df.loc[df["arm"] == "HC", "score"].mean()
    )

    batch1 = df.index[df["batch"] == 1].tolist()
    batch2 = df.index[df["batch"] == 2].tolist()

    null = []

    for hc1 in itertools.combinations(batch1, 1):
        for hc2 in itertools.combinations(batch2, 2):
            hc = set(hc1 + hc2)
            null.append(
                df.drop(index=list(hc))["score"].mean()
                - df.loc[list(hc), "score"].mean()
            )

    null = np.asarray(null)
    extreme = np.sum(null <= observed)
    p = extreme / len(null)

    return observed, extreme, len(null), p

for name, genes in [
    ("PRIMARY", primary_genes),
    ("SECONDARY", secondary_genes),
]:
    observed, extreme, total, p = exact_test(genes)

    print(name)
    print("SSc - HC =", observed)
    print(f"Exact p = {extreme}/{total} = {p}")
    print()
