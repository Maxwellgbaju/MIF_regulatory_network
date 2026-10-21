import xenaPython as xena
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
import pandas as pd

host = 'https://toil.xenahubs.net'
surv_dataset = 'TCGA_survival_data'

with open('xena_phenotype_full.json') as f:
    metadata = json.load(f)
with open('all_gene_expression.json') as f:
    expr_data = json.load(f)

# Get tumor samples for Lung Squamous Cell Carcinoma
cancer_type = "Lung Squamous Cell Carcinoma"
tumor_ids = [s for s, m in metadata.items()
             if m['_study'] == 'TCGA' and m['_sample_type'] == 'Primary Tumor'
             and m['primary disease or tissue'] == cancer_type]
print(f"Tumor samples for {cancer_type}: {len(tumor_ids)}")

# Fetch survival data for these samples
print("Fetching survival data (OS, OS.time)...")
surv_values = xena.dataset_probe_values(host, surv_dataset, tumor_ids, ['OS', 'OS.time'])
os_status = surv_values[1][0]
os_time = surv_values[1][1]

df = pd.DataFrame({
    "sample": tumor_ids,
    "OS": [float(v) if v not in (None, "NaN") else np.nan for v in os_status],
    "OS_time": [float(v) if v not in (None, "NaN") else np.nan for v in os_time],
    "MIF_expr": [float(expr_data["MIF"][s]) if s in expr_data["MIF"] else np.nan for s in tumor_ids]
})
df = df.dropna()
print(f"Samples with complete survival + expression data: {len(df)}")
print(df[["OS","OS_time","MIF_expr"]].describe())

# Median split
median_mif = df["MIF_expr"].median()
df["group"] = np.where(df["MIF_expr"] >= median_mif, "High MIF", "Low MIF")

high = df[df["group"] == "High MIF"]
low = df[df["group"] == "Low MIF"]

lr = logrank_test(high["OS_time"], low["OS_time"], high["OS"], low["OS"])
print(f"\nLog-rank test: p = {lr.p_value:.4f}")

cph = CoxPHFitter()
cph.fit(df.rename(columns={"MIF_expr":"expr","OS_time":"time","OS":"event"})[["time","event","expr"]],
        duration_col="time", event_col="event")
hr = cph.hazard_ratios_["expr"]
ci = cph.confidence_intervals_.loc["expr"]
print(f"Cox HR (continuous MIF): {hr:.3f} (95% CI {np.exp(ci.iloc[0]):.3f}-{np.exp(ci.iloc[1]):.3f})")

fig, ax = plt.subplots(figsize=(7,6))
kmf_h = KaplanMeierFitter()
kmf_h.fit(high["OS_time"], high["OS"], label=f"High MIF (n={len(high)})")
kmf_h.plot_survival_function(ax=ax, color="#C00000")
kmf_l = KaplanMeierFitter()
kmf_l.fit(low["OS_time"], low["OS"], label=f"Low MIF (n={len(low)})")
kmf_l.plot_survival_function(ax=ax, color="#4472C4")

ax.set_xlabel("Time (months)")
ax.set_ylabel("Overall survival probability")
ax.set_title(f"{cancer_type}: Overall survival by MIF expression\n(UCSC Xena TCGA TARGET GTEx cohort)")
ax.text(0.98, 0.05, f"Log-rank p = {lr.p_value:.4f}\nCox HR = {hr:.2f}", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=10, bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray"))
plt.tight_layout()
plt.savefig("MIF_LUSC_KM_xena.png", dpi=300, bbox_inches="tight")
plt.savefig("MIF_LUSC_KM_xena.eps", format="eps", bbox_inches="tight")
print("\nSaved MIF_LUSC_KM_xena.png/.eps")
