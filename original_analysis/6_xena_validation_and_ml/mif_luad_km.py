import xenaPython as xena
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

host = 'https://toil.xenahubs.net'
surv_dataset = 'TCGA_survival_data'

with open('xena_phenotype_full.json') as f:
    metadata = json.load(f)
with open('all_gene_expression.json') as f:
    expr_data = json.load(f)

cancer_type = "Lung Adenocarcinoma"
tumor_ids = [s for s, m in metadata.items()
             if m['_study'] == 'TCGA' and m['_sample_type'] == 'Primary Tumor'
             and m['primary disease or tissue'] == cancer_type]
print(f"Tumor samples for {cancer_type}: {len(tumor_ids)}")

print("Fetching survival data (OS, OS.time)...")
surv_values = xena.dataset_probe_values(host, surv_dataset, tumor_ids, ['OS', 'OS.time'])
os_status_raw = surv_values[1][0]
os_time_raw = surv_values[1][1]

records = []
for i, s in enumerate(tumor_ids):
    try:
        os_stat = float(os_status_raw[i])
        os_days = float(os_time_raw[i])
    except (ValueError, TypeError):
        continue
    if np.isnan(os_stat) or np.isnan(os_days):
        continue
    if s not in expr_data["MIF"]:
        continue
    try:
        mif_val = float(expr_data["MIF"][s])
    except (ValueError, TypeError):
        continue
    if np.isnan(mif_val):
        continue
    records.append({"sample": s, "OS": os_stat, "OS_time_months": os_days / 30.44, "MIF_expr": mif_val})

df = pd.DataFrame(records)
print(f"Samples with complete survival + MIF expression data: {len(df)}")
print(df[["OS", "OS_time_months", "MIF_expr"]].describe())

median_mif = df["MIF_expr"].median()
df["group"] = np.where(df["MIF_expr"] >= median_mif, "High MIF", "Low MIF")
high = df[df["group"] == "High MIF"]
low = df[df["group"] == "Low MIF"]

lr = logrank_test(high["OS_time_months"], low["OS_time_months"], high["OS"], low["OS"])
print(f"\nLog-rank test: p = {lr.p_value:.4f}")

cph = CoxPHFitter()
cph.fit(df.rename(columns={"MIF_expr": "expr", "OS_time_months": "time", "OS": "event"})[["time", "event", "expr"]],
        duration_col="time", event_col="event")
hr = cph.hazard_ratios_["expr"]
ci = cph.confidence_intervals_.loc["expr"]
hr_lo, hr_hi = np.exp(ci.iloc[0]), np.exp(ci.iloc[1])
print(f"Cox HR (continuous MIF): {hr:.3f} (95% CI {hr_lo:.3f}-{hr_hi:.3f})")

fig, ax = plt.subplots(figsize=(7, 6))
kmf_h = KaplanMeierFitter()
kmf_h.fit(high["OS_time_months"], high["OS"], label=f"High MIF (n={len(high)})")
kmf_h.plot_survival_function(ax=ax, color="#C00000")
kmf_l = KaplanMeierFitter()
kmf_l.fit(low["OS_time_months"], low["OS"], label=f"Low MIF (n={len(low)})")
kmf_l.plot_survival_function(ax=ax, color="#4472C4")

ax.set_xlabel("Time (months)")
ax.set_ylabel("Overall survival probability")
ax.set_title(f"{cancer_type}: Overall survival by MIF expression\n(UCSC Xena TCGA TARGET GTEx cohort)")
ax.text(0.98, 0.05, f"Log-rank p = {lr.p_value:.4f}\nCox HR = {hr:.2f} (95% CI {hr_lo:.2f}-{hr_hi:.2f})",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray"))
ax.legend(loc="upper right")
plt.tight_layout()
plt.savefig("MIF_LUAD_KM_xena.png", dpi=300, bbox_inches="tight")
plt.savefig("MIF_LUAD_KM_xena.eps", format="eps", bbox_inches="tight")
print("\nSaved MIF_LUAD_KM_xena.png/.eps")
