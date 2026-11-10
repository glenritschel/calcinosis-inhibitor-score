import itertools
import pandas as pd

mapping = {
    "D01": "GSM8872107",
    "D02": "GSM8872104",
    "D03": "GSM8872097",
    "D04": "GSM8872100",
    "D05": "GSM8872096",
    "D06": "GSM8872108",
    "D07": "GSM8872103",
    "D08": "GSM8872102",
    "D09": "GSM8872095",
    "D10": "GSM8872101",
    "D11": "GSM8872106",
    "D12": "GSM8872094",
    "D13": "GSM8872099",
    "D14": "GSM8872098",
    "D15": "GSM8872105",
}

scores = {
    "D01": -0.494404638573,
    "D02": -0.561628879266,
    "D03": 0.456419221526,
    "D04": -0.132812637824,
    "D05": 1.372108008403,
    "D06": 0.136892343403,
    "D07": 0.084392785647,
    "D08": -0.546635530699,
    "D09": 0.368031640898,
    "D10": -0.936432533541,
    "D11": 0.644300378978,
    "D12": 0.370901491897,
    "D13": -0.138427324202,
    "D14": -0.557349529992,
    "D15": -0.065354796655,
}

metadata = {
    "GSM8872094": ("SSc", 1),
    "GSM8872095": ("SSc", 1),
    "GSM8872096": ("SSc", 1),
    "GSM8872097": ("HC", 1),
    "GSM8872098": ("SSc", 2),
    "GSM8872099": ("SSc", 2),
    "GSM8872100": ("SSc", 2),
    "GSM8872101": ("SSc", 2),
    "GSM8872102": ("SSc", 2),
    "GSM8872103": ("SSc", 2),
    "GSM8872104": ("SSc", 2),
    "GSM8872105": ("SSc", 2),
    "GSM8872106": ("SSc", 2),
    "GSM8872107": ("HC", 2),
    "GSM8872108": ("HC", 2),
}

rows = []
for code, gsm in mapping.items():
    arm, batch = metadata[gsm]
    rows.append({
        "code": code,
        "gsm": gsm,
        "arm": arm,
        "batch": batch,
        "score": scores[code],
    })

df = pd.DataFrame(rows).sort_values("code")

assert len(df) == 15
assert (df["arm"] == "SSc").sum() == 12
assert (df["arm"] == "HC").sum() == 3
assert len(df[(df["batch"] == 1) & (df["arm"] == "HC")]) == 1
assert len(df[(df["batch"] == 2) & (df["arm"] == "HC")]) == 2

ssc_mean = df.loc[df["arm"] == "SSc", "score"].mean()
hc_mean = df.loc[df["arm"] == "HC", "score"].mean()
observed = ssc_mean - hc_mean

batch1 = df.index[df["batch"] == 1].tolist()
batch2 = df.index[df["batch"] == 2].tolist()

null_effects = []

for hc_batch1 in itertools.combinations(batch1, 1):
    for hc_batch2 in itertools.combinations(batch2, 2):
        hc_indices = set(hc_batch1 + hc_batch2)
        ssc_indices = [i for i in df.index if i not in hc_indices]

        null_effects.append(
            df.loc[ssc_indices, "score"].mean()
            - df.loc[list(hc_indices), "score"].mean()
        )

assert len(null_effects) == 220

extreme = sum(x <= observed for x in null_effects)
p_exact = extreme / len(null_effects)

print("SSc mean:", ssc_mean)
print("HC mean:", hc_mean)
print("SSc - HC:", observed)
print("Admissible labelings:", len(null_effects))
print("Number <= observed:", extreme)
print("Exact one-sided p:", p_exact)
